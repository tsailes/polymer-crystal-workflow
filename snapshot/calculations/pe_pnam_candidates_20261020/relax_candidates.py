"""Restartable, fixed-cell Pnam GFN2-xTB preoptimization of generated PE candidates.

Run with the isolated WSL Python environment. Every selected input reaches the
fine stage; coarse energies never eliminate candidates. Final sampling beyond
2x3x8 remains a separate convergence check, not an implied result of this script.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import time
import traceback

# Set before importing BLAS or tblite; orchestrator sets this for each worker.
THREADS = int(os.environ.get("PE_DFT_THREADS", "4"))
os.environ["OMP_NUM_THREADS"] = str(THREADS)
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OMP_MAX_ACTIVE_LEVELS"] = "1"

import numpy as np
import spglib
from ase.calculators.singlepoint import SinglePointCalculator
from ase.constraints import FixConstraint
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.optimize import BFGS
from ase.spacegroup.symmetrize import symmetrize_rank1
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parent
WORKFLOW = ROOT.parent / "pe_gfn2_first_run_20261020/convergence_relaxation/workflow.py"
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location("verified_pe_gfn2_workflow", WORKFLOW)
workflow = importlib.util.module_from_spec(spec)
spec.loader.exec_module(workflow)
GENERATOR = None
TOPOLOGY = None
SYMPREC = 1e-5


def stamp():
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, data):
    path = Path(path)
    temp = path.with_suffix(path.suffix + f".{os.getpid()}.tmp")
    temp.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                    encoding="utf-8")
    temp.replace(path)


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def prepare_geometry_context(input_path):
    global GENERATOR, TOPOLOGY
    generator_path = ROOT / "generate_candidates.py"
    if GENERATOR is None:
        generator_spec = importlib.util.spec_from_file_location("pe_pnam_generator", generator_path)
        GENERATOR = importlib.util.module_from_spec(generator_spec)
        generator_spec.loader.exec_module(GENERATOR)
    topology_path = Path(input_path).parent / "reference_topology.json"
    TOPOLOGY = load(topology_path)
    return {"generator_sha256": digest(generator_path),
            "reference_topology_path": str(topology_path),
            "reference_topology_sha256": digest(topology_path)}


def repeats(value):
    result = tuple(int(x) for x in value.replace("x", ",").split(","))
    if len(result) != 3 or min(result) < 1:
        raise argparse.ArgumentTypeError("Expected three positive integers, e.g. 2,2,6")
    return result


def target_operations(manifest):
    target = manifest.get("target_symmetry", {})
    if "rotations" in target and "translations" in target:
        rotations = np.asarray(target["rotations"], dtype=float)
        translations = np.asarray(target["translations"], dtype=float)
        source = "generation_manifest.target_symmetry"
    else:
        target = spglib.get_symmetry_from_database(297)
        rotations = np.asarray(target["rotations"], dtype=float)
        translations = np.asarray(target["translations"], dtype=float)
        source = "spglib Hall 297: Pnam, choice a-cb"
    if rotations.shape != (8, 3, 3) or translations.shape != (8, 3):
        raise ValueError("Exactly eight target Pnam operations are required")
    if not np.allclose(rotations, np.rint(rotations), atol=1e-12, rtol=0):
        raise ValueError("Noninteger fractional symmetry rotations")
    rotations = np.rint(rotations).astype(int)
    translations %= 1
    # Reject duplicate or non-closed supplied operation sets.
    def locate(rotation, translation):
        dt = translations - translation
        dt -= np.rint(dt)
        return np.flatnonzero(np.all(rotations == rotation, axis=(1, 2))
                              & (np.max(np.abs(dt), axis=1) < 1e-9))
    for r, t in zip(rotations, translations):
        if len(locate(r, t)) != 1:
            raise ValueError("Duplicate target operation")
        for r2, t2 in zip(rotations, translations):
            if len(locate(r @ r2, r @ t2 + t)) != 1:
                raise ValueError("Target operation set is not closed")
    return rotations, translations, source


def operation_maps(atoms, rotations, translations, tolerance=SYMPREC):
    """Bijective, element-preserving map i -> operation(i), in the fixed cell."""
    scaled = atoms.get_scaled_positions()
    numbers = atoms.numbers
    maps, largest = [], 0.0
    for rotation, translation in zip(rotations, translations):
        cart_rotation = atoms.cell.array.T @ rotation @ np.linalg.inv(atoms.cell.array.T)
        if not np.allclose(cart_rotation.T @ cart_rotation, np.eye(3), atol=1e-10):
            raise ValueError("Target operation is not an isometry of this cell")
        transformed = scaled @ rotation.T + translation
        displacement = scaled[None, :, :] - transformed[:, None, :]
        displacement -= np.rint(displacement)
        distance = np.linalg.norm(displacement @ atoms.cell.array, axis=2)
        distance[numbers[:, None] != numbers[None, :]] = 1e20
        rows, columns = linear_sum_assignment(distance)
        maximum = float(distance[rows, columns].max())
        if maximum > tolerance:
            raise ValueError(f"Target Pnam operation mismatch: {maximum:.6g} A")
        mapping = np.empty(len(atoms), dtype=int)
        mapping[rows] = columns
        maps.append(mapping)
        largest = max(largest, maximum)
    return maps, largest


class TargetPnam(FixConstraint):
    """Project only the supplied target group, not accidental extra symmetry."""

    def __init__(self, atoms, rotations, translations):
        self.rotations = rotations
        self.translations = translations
        self.maps, self.initial_error = operation_maps(atoms, rotations, translations)

    def project(self, atoms, vectors):
        return symmetrize_rank1(atoms.cell.array, atoms.cell.reciprocal().T,
                               np.asarray(vectors), self.rotations,
                               self.translations, self.maps)

    def adjust_positions(self, atoms, proposed):
        proposed[:] = atoms.positions + self.project(atoms, proposed - atoms.positions)

    def adjust_forces(self, atoms, forces):
        forces[:] = self.project(atoms, forces)


class ConstrainedReplicaCalculator(workflow.ReplicatedGFN2):
    """Keep constraints on primitive coordinates, never on the expanded cell."""

    def calculate(self, atoms=None, properties=("energy",), system_changes=workflow.all_changes):
        plain = (atoms if atoms is not None else self.atoms).copy()
        plain.set_constraint()
        super().calculate(plain, properties, system_changes)


def clean_copy(atoms):
    result = atoms.copy()
    result.set_constraint()
    result.calc = None
    return result


def structural_audit(atoms, initial_cell, rotations, translations):
    errors = []
    if len(atoms) != 12 or atoms.get_chemical_formula() != "C4H8":
        errors.append("Expected exactly 12 atoms, C4H8")
    if not all(atoms.pbc):
        errors.append("Three-dimensional periodicity required")
    if not np.array_equal(atoms.cell.array, initial_cell):
        errors.append("Fixed cell changed")
    if not np.allclose(atoms.cell.array, np.diag(np.diag(atoms.cell.array)), atol=1e-12):
        errors.append("This executor requires the orthorhombic a,b,c axis setting")
    if not np.isfinite(atoms.positions).all():
        errors.append("Nonfinite atom positions")
    if "chain_id" not in atoms.arrays or set(atoms.arrays["chain_id"]) != {0, 1}:
        errors.append("Two chain_id labels, 0 and 1, required")
    chemistry = None
    symmetry_error = None
    if not errors:
        if GENERATOR is None or TOPOLOGY is None:
            errors.append("Generator and immutable periodic bond topology are required")
        else:
            chemistry = GENERATOR.audit_geometry(atoms, TOPOLOGY)
            # Preserve the predefined periodic graph. A coarse stage can move
            # across a geometric screening boundary; it must still reach fine.
            # Stop only on a broken chemical structure or lost symmetry here.
            for bond in chemistry["bonds"]:
                bounds = (1.35, 1.75) if bond["kind"] == "CC" else (0.95, 1.22)
                if not bounds[0] < bond["distance_A"] < bounds[1]:
                    errors.append(f"Reference bond outside safety bounds: {bond}")
            if not chemistry["symmetry"]["contains_target_group"]:
                errors.append("Generator audit: target Pnam group lost")
        try:
            _, symmetry_error = operation_maps(atoms, rotations, translations)
        except ValueError as exc:
            errors.append(str(exc))
    return {"passed": not errors, "errors": errors, "chemistry": chemistry,
            "target_operation_max_error_A": symmetry_error}


def geometric_projection_test(atoms, constraint):
    """No electronic calculation: orthogonal-projector and tangent checks."""
    rng = np.random.default_rng(1739)
    u, v = rng.normal(size=(2, len(atoms), 3))
    pu, pv = constraint.project(atoms, u), constraint.project(atoms, v)
    idempotence = float(np.max(np.abs(constraint.project(atoms, pu) - pu)))
    adjoint = abs(float(np.vdot(u, pv) - np.vdot(pu, v)))
    displaced = clean_copy(atoms)
    displaced.positions += pu * 1e-4
    _, symmetry = operation_maps(displaced, constraint.rotations, constraint.translations)
    # Analytic quadratic potential supplies an independent directional derivative.
    h = 1e-5
    x = atoms.positions.copy()
    target = rng.normal(size=x.shape)
    eplus = 0.5 * np.sum((x + h * pv - target) ** 2)
    eminus = 0.5 * np.sum((x - h * pv - target) ** 2)
    numerical = (eplus - eminus) / (2 * h)
    analytic = float(np.vdot(constraint.project(atoms, x - target), pv))
    derivative_error = abs(float(numerical - analytic))
    passed = idempotence < 1e-10 and adjoint < 1e-10 and derivative_error < 1e-7
    return dict(passed=passed, idempotence_max=idempotence, adjoint_error=adjoint,
                symmetry_after_tangent_step_A=symmetry,
                quadratic_directional_derivative_error=derivative_error)


def package_versions():
    return {p: importlib.metadata.version(p)
            for p in ("ase", "tblite", "numpy", "scipy", "spglib", "gemmi")}


def resolve_input(entry):
    path = (ROOT / entry["input_path"]).resolve()
    if not path.is_relative_to(ROOT):
        raise ValueError("Candidate input must be inside the candidate directory")
    return path


def run_stage(atoms, name, repeat, fmax, steps, attempt, rotations, translations,
              original_cell, update):
    folder = attempt / name
    folder.mkdir()
    constraint = TargetPnam(atoms, rotations, translations)
    atoms.set_constraint(constraint)
    atoms.calc = ConstrainedReplicaCalculator(repeat, accuracy=0.01)
    trajectory = Trajectory(str(folder / "optimization.traj"), "w")
    optimizer = BFGS(atoms, logfile=str(folder / "optimizer.log"), maxstep=0.05,
                     restart=str(folder / "bfgs_restart.json"))
    started = time.monotonic()
    history = []

    def record():
        energy = float(atoms.get_potential_energy())
        raw = atoms.get_forces(apply_constraint=False)
        projected = atoms.get_forces()
        _, symmetry_error = operation_maps(atoms, rotations, translations)
        if not np.array_equal(atoms.cell.array, original_cell):
            raise RuntimeError("Cell changed")
        if not np.isfinite([energy]).all() or not np.isfinite(raw).all():
            raise RuntimeError("Nonfinite energy or forces")
        row = dict(step=int(optimizer.nsteps), energy_eV_per_cell=energy,
                   energy_eV_per_CH2=energy / 4,
                   raw_fmax_eV_A=float(np.linalg.norm(raw, axis=1).max()),
                   projected_fmax_eV_A=float(np.linalg.norm(projected, axis=1).max()),
                   projection_residual_eV_A=float(np.linalg.norm(raw - projected, axis=1).max()),
                   replica_force_deviation_eV_A=atoms.calc.last_folded[
                       "replica_force_max_deviation_eV_per_A"],
                   symmetry_error_A=symmetry_error, elapsed_seconds=time.monotonic() - started,
                   raw_forces_eV_A=raw.tolist(), projected_forces_eV_A=projected.tolist())
        history.append(row)
        dump(folder / "history.json", history)
        frame = clean_copy(atoms)
        frame.info.update(step=int(optimizer.nsteps), energy_eV_per_cell=energy,
                          energy_eV_per_CH2=energy / 4, stage=name)
        frame.set_array("raw_forces", raw.copy())
        frame.set_array("projected_forces", projected.copy())
        write(folder / "latest.extxyz", frame)
        write(folder / "optimization.extxyz", frame, append=True)
        frame.calc = SinglePointCalculator(frame, energy=energy, forces=projected)
        trajectory.write(frame)
        update(status="running", stage=name, step=optimizer.nsteps,
               checkpoint=str(folder / "latest.extxyz"), last_step=row)
        print("STEP", json.dumps({k: row[k] for k in ("step", "energy_eV_per_CH2",
                    "raw_fmax_eV_A", "projected_fmax_eV_A", "elapsed_seconds")}), flush=True)

    optimizer.attach(record)
    try:
        converged = bool(optimizer.run(fmax=fmax, steps=steps))
    finally:
        trajectory.close()
    audit = structural_audit(atoms, original_cell, rotations, translations)
    result = dict(stage=name, repeat=list(repeat), fmax_target_eV_A=fmax,
                  converged=converged, steps=optimizer.nsteps,
                  elapsed_seconds=time.monotonic() - started, final=history[-1],
                  structural_audit=audit)
    dump(folder / "result.json", result)
    write(folder / "final.extxyz", clean_copy(atoms))
    if not audit["passed"]:
        raise RuntimeError(f"{name}: geometry audit failed: {audit['errors']}")
    return clean_copy(atoms), result


def worker(args, manifest, entry):
    attempt = Path(args.attempt_dir)
    candidate = entry["candidate_id"]
    relaxation = attempt.parent.parent
    status_path = relaxation / "status.json"
    previous = load(status_path) if status_path.exists() else {}
    initial_path = resolve_input(entry)
    initial = read(initial_path)
    original_cell = initial.cell.array.copy()
    rotations, translations, op_source = target_operations(manifest)
    config = dict(candidate_id=candidate, parameters=entry.get("parameters"),
                  input_path=str(initial_path), input_sha256=digest(initial_path),
                  script_sha256=digest(__file__), workflow_path=str(WORKFLOW),
                  workflow_sha256=digest(WORKFLOW), manifest_sha256=digest(args.manifest),
                  method_settings=workflow.SETTINGS, packages=package_versions(),
                  fine_repeat=[2, 3, 8], fine_fmax_eV_A=args.fine_fmax,
                  continue_completed=args.continue_completed,
                  coarse_repeat=list(args.coarse_repeat) if args.coarse_repeat else None,
                  coarse_fmax_eV_A=0.03, max_steps_per_stage=args.steps,
                  maxstep_A=0.05, threads=THREADS, python=sys.executable,
                  symmetry_source=op_source, rotations=rotations.tolist(),
                  translations=translations.tolist(), fixed_cell_A=original_cell.tolist(),
                  sampling_scope="Gamma exact replicas; new-candidate dense-grid recheck pending",
                  resume_mode="geometry checkpoint, fresh BFGS Hessian",
                  no_energy_elimination=True)
    state = dict(candidate_id=candidate, attempt_directory=str(attempt), status="starting",
                 pid=os.getpid(), hostname=socket.gethostname(), started_utc=stamp(),
                 input_sha256=config["input_sha256"])

    def update(**kwargs):
        state.update(kwargs, updated_utc=stamp())
        dump(attempt / "status.json", state)
        dump(status_path, state)

    try:
        config.update(prepare_geometry_context(initial_path))
        snapshot = attempt / "source_snapshot"
        snapshot.mkdir(exist_ok=True)
        for source in (Path(__file__), WORKFLOW, ROOT / "generate_candidates.py",
                       args.manifest, Path(config["reference_topology_path"])):
            shutil.copy2(source, snapshot / source.name)
        if config["packages"]["tblite"] != "0.7.0":
            raise RuntimeError("This executor requires validated tblite 0.7.0")
        initial_audit = structural_audit(initial, original_cell, rotations, translations)
        dump(attempt / "input_audit.json", initial_audit)
        if not initial_audit["passed"]:
            raise ValueError(f"Invalid generated candidate: {initial_audit['errors']}")
        if not initial_audit["chemistry"]["passed"]:
            raise ValueError("Generated input no longer passes its physical survivor filters")
        state["reference_topology_sha256"] = config["reference_topology_sha256"]
        projection_test = geometric_projection_test(initial, TargetPnam(initial, rotations, translations))
        dump(attempt / "projection_validation.json", projection_test)
        if not projection_test["passed"]:
            raise RuntimeError("Target-group projection validation failed")
        atoms = clean_copy(initial)
        resumed_stage = None
        if not args.fresh and previous.get("checkpoint"):
            checkpoint = Path(previous["checkpoint"])
            if previous.get("input_sha256") != config["input_sha256"]:
                raise ValueError("Input changed; use --fresh to create an explicitly new attempt")
            if previous.get("reference_topology_sha256") != config["reference_topology_sha256"]:
                raise ValueError("Reference topology changed; use --fresh")
            if checkpoint.exists():
                resumed = read(checkpoint)
                audit = structural_audit(resumed, original_cell, rotations, translations)
                if not audit["passed"]:
                    raise ValueError("Previous checkpoint failed geometry validation")
                atoms = clean_copy(resumed)
                resumed_stage = previous.get("stage")
                config["resumed_from"] = str(checkpoint)
        dump(attempt / "configuration.json", config)
        write(attempt / "input.extxyz", atoms)
        update(status="running", stage="initialization")
        if args.coarse_repeat and resumed_stage not in ("fine", "final_validation", "finished"):
            atoms, coarse = run_stage(atoms, "coarse", args.coarse_repeat, 0.03, args.steps,
                                      attempt, rotations, translations, original_cell, update)
            # Hitting coarse maxsteps still proceeds to fine; no energy cutoff.
        atoms, fine = run_stage(atoms, "fine", (2, 3, 8), args.fine_fmax, args.steps,
                                attempt, rotations, translations, original_cell, update)
        update(stage="final_validation", checkpoint=str(attempt / "fine/final.extxyz"))
        # Independent electronic initialization and tighter accuracy at the same geometry.
        atoms.set_constraint(TargetPnam(atoms, rotations, translations))
        atoms.calc = ConstrainedReplicaCalculator((2, 3, 8), accuracy=0.001)
        energy = float(atoms.get_potential_energy())
        raw = atoms.get_forces(apply_constraint=False)
        projected = atoms.get_forces()
        audit = structural_audit(atoms, original_cell, rotations, translations)
        delta_e = energy / 4 - fine["final"]["energy_eV_per_CH2"]
        delta_f = float(np.linalg.norm(raw - np.asarray(
            fine["final"]["raw_forces_eV_A"]), axis=1).max())
        raw_fmax = float(np.linalg.norm(raw, axis=1).max())
        projected_fmax = float(np.linalg.norm(projected, axis=1).max())
        passed = (fine["converged"] and audit["passed"] and raw_fmax < args.fine_fmax
                  and projected_fmax < args.fine_fmax and abs(delta_e) < 0.001 and delta_f < 0.005)
        result = dict(passed=bool(passed), status="preoptimized" if passed else "validation_failed",
                      energy_eV_per_cell=energy, energy_eV_per_CH2=energy / 4,
                      raw_fmax_eV_A=raw_fmax, projected_fmax_eV_A=projected_fmax,
                      raw_forces_eV_A=raw.tolist(), projected_forces_eV_A=projected.tolist(),
                      accuracy_energy_change_eV_per_CH2=delta_e,
                      accuracy_force_change_eV_A=delta_f,
                      structural_audit=audit, folded_result=atoms.calc.last_folded,
                      final_original_geometric_filters_passed=audit["chemistry"]["passed"],
                      scientific_scope="Fixed-cell target-Pnam constrained preoptimization",
                      sampling_status="2x3x8 validated on pilot; candidate dense-grid recheck pending",
                      deduplication_status="Not deduplicated by this executor; retain every candidate")
        final = clean_copy(atoms)
        final.info.update(candidate_id=candidate, energy_eV_per_CH2=energy / 4,
                          optimization_status=result["status"])
        final.set_array("raw_forces", raw)
        final.set_array("projected_forces", projected)
        GENERATOR.export_structure(final, attempt, stem="final")
        dump(attempt / "final_result.json", result)
        update(status=result["status"], stage="finished", result=str(attempt / "final_result.json"),
               checkpoint=str(attempt / "final.extxyz"), finished_utc=stamp())
        return 0 if passed else 2
    except Exception as exc:
        update(status="failed", error=str(exc), traceback=traceback.format_exc())
        print(traceback.format_exc(), flush=True)
        return 2


def launch(args, entry):
    candidate = entry["candidate_id"]
    if not candidate or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in candidate):
        raise ValueError("Unsafe candidate_id")
    relaxation = ROOT / "candidates" / candidate / "relaxation"
    relaxation.mkdir(parents=True, exist_ok=True)
    status_path = relaxation / "status.json"
    if status_path.exists() and not args.fresh:
        previous = load(status_path)
        if previous.get("status") == "preoptimized" and not args.continue_completed:
            if previous.get("input_sha256") != digest(resolve_input(entry)):
                raise ValueError("Previously accepted input changed; use --fresh")
            return dict(candidate_id=candidate, status="preoptimized", skipped=True,
                        result=previous.get("result"))
    lock = relaxation / "worker.lock"
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError:
        existing = load(lock)
        alive = True
        try:
            os.kill(int(existing["pid"]), 0)
        except ProcessLookupError:
            alive = False
        if alive:
            return dict(candidate_id=candidate, status="busy", lock=existing)
        lock.unlink()
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    with os.fdopen(fd, "w") as handle:
        json.dump(dict(pid=os.getpid(), hostname=socket.gethostname(), time=stamp()), handle)
    try:
        attempt = relaxation / "attempts" / stamp()
        attempt.mkdir(parents=True, exist_ok=False)
        command = [sys.executable, str(Path(__file__).resolve()), "--worker",
                   "--candidate", candidate, "--manifest", str(args.manifest),
                   "--attempt-dir", str(attempt), "--steps", str(args.steps),
                   "--fine-fmax", str(args.fine_fmax)]
        if args.continue_completed:
            command += ["--continue-completed"]
        if args.coarse_repeat:
            command += ["--coarse-repeat", ",".join(map(str, args.coarse_repeat))]
        if args.fresh:
            command += ["--fresh"]
        threads = args.threads if args.threads else (2 if args.workers > 1 else 4)
        with (attempt / "native_output.log").open("w", encoding="utf-8") as output:
            process = subprocess.Popen(command, stdout=output, stderr=subprocess.STDOUT,
                                       env={**os.environ, "PE_DFT_THREADS": str(threads)})
            dump(lock, dict(pid=process.pid, launcher_pid=os.getpid(),
                            hostname=socket.gethostname(), time=stamp()))
            return_code = process.wait()
        state = load(status_path) if status_path.exists() else dict(status="worker_failed_before_status")
        valid_state = state.get("attempt_directory") == str(attempt)
        return dict(candidate_id=candidate, exit_code=return_code,
                    status=state["status"] if valid_state else "worker_failed_before_status",
                    attempt_directory=str(attempt),
                    result=state.get("result"))
    finally:
        lock.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choose = parser.add_mutually_exclusive_group(required=True)
    choose.add_argument("--candidate")
    choose.add_argument("--all", action="store_true")
    parser.add_argument("--manifest", type=Path, default=ROOT / "generation_manifest.json")
    parser.add_argument("--coarse-repeat", type=repeats, default=None,
                        help="Optional warm-up, e.g. 2,2,6; every candidate still reaches 2x3x8")
    parser.add_argument("--steps", type=int, default=100, help="Maximum steps in each stage")
    parser.add_argument("--workers", type=int, default=1, choices=(1, 2, 4))
    parser.add_argument("--threads", type=int, choices=(1, 2, 4))
    parser.add_argument("--fresh", action="store_true", help="New attempt from generated input")
    parser.add_argument("--fine-fmax", type=float, default=.01,
                        help="Fine optimization force target, at most 0.01 eV/A")
    parser.add_argument("--continue-completed", action="store_true",
                        help="Continue a completed geometry checkpoint in a new preserved attempt")
    parser.add_argument("--validate-only", action="store_true", help="Geometry/projection checks, no xTB")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--attempt-dir", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.steps < 1:
        parser.error("--steps must be positive")
    if not 0 < args.fine_fmax <= .01:
        parser.error("--fine-fmax must be positive and no looser than 0.01")
    if args.threads and args.threads * args.workers > 8:
        parser.error("At most eight computational threads in total")
    manifest = load(args.manifest)
    entries = manifest["unique_survivors"]
    if args.candidate:
        entries = [entry for entry in entries if entry["candidate_id"] == args.candidate]
        if len(entries) != 1:
            parser.error("Candidate not found uniquely in unique_survivors")
    if len({entry["candidate_id"] for entry in entries}) != len(entries):
        raise ValueError("Duplicate candidate IDs in manifest")
    if args.worker:
        return worker(args, manifest, entries[0])
    if args.validate_only:
        rotations, translations, source = target_operations(manifest)
        records = []
        for entry in entries:
            input_path = resolve_input(entry)
            prepare_geometry_context(input_path)
            atoms = read(input_path)
            audit = structural_audit(atoms, atoms.cell.array.copy(), rotations, translations)
            projection = geometric_projection_test(atoms, TargetPnam(atoms, rotations, translations)) if audit["passed"] else None
            records.append(dict(candidate_id=entry["candidate_id"], audit=audit, projection=projection))
        print(json.dumps(dict(symmetry_source=source, validation=records), indent=2))
        return 0 if all(r["audit"]["passed"] and r["projection"]["passed"] for r in records) else 2
    batch = ROOT / "relaxation_batches" / stamp()
    batch.mkdir(parents=True, exist_ok=False)
    dump(batch / "configuration.json", dict(manifest=str(args.manifest),
          manifest_sha256=digest(args.manifest), script_sha256=digest(__file__),
          workflow_sha256=digest(WORKFLOW), candidate_ids=[e["candidate_id"] for e in entries],
          workers=args.workers, threads_per_worker=args.threads or (2 if args.workers > 1 else 4),
          coarse_repeat=args.coarse_repeat, fine_repeat=[2, 3, 8], steps=args.steps,
          fine_fmax_eV_A=args.fine_fmax, continue_completed=args.continue_completed))
    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(launch, args, entry): entry for entry in entries}
        for future in as_completed(futures):
            entry = futures[future]
            try:
                result = future.result()
            except Exception as exc:
                result = dict(candidate_id=entry["candidate_id"],
                              status="orchestrator_failed", error=str(exc))
            results.append(result)
            dump(batch / "summary.json", results)
            print(json.dumps(result), flush=True)
    return 0 if all(r["status"] == "preoptimized" for r in results) else 2


if __name__ == "__main__":
    raise SystemExit(main())
