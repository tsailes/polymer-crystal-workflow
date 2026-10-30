"""Auditable resume of fixed-cell Pnam PE relaxation; identical on every MPI rank.

No calculation is started by importing this module or using --check-input.
An existing output directory is never overwritten.  A resume object pairs the
unevaluated pending geometry with the saved LBFGS history in a fresh attempt.
Original physical inputs stay unchanged; electronic states always start cold.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import sys
import time
import traceback

import numpy as np
import spglib
from ase import Atoms
from ase.calculators.mixing import SumCalculator
from ase.calculators.singlepoint import SinglePointCalculator
from ase.constraints import FixSymmetry
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.optimize import LBFGS
from ase.parallel import broadcast
from ase.units import Bohr, GPa, Hartree
from gpaw.mpi import world


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def convert(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(type(value).__name__)


def root_write(path, obj, append=False):
    """Only root writes JSON; propagate root I/O errors to all ranks."""
    error = None
    if world.rank == 0:
        try:
            with Path(path).open("a" if append else "w", encoding="utf-8") as handle:
                handle.write(json.dumps(obj, default=convert, allow_nan=False,
                                        indent=None if append else 2) + "\n")
                handle.flush()
        except Exception:
            error = traceback.format_exc()
    error = broadcast(error, root=0, comm=world)
    if error:
        raise OSError(error)


def project_root(path):
    for candidate in [path, *path.parents]:
        if (candidate / "scripts/run_point.py").is_file() and (candidate / "protocol.json").is_file():
            return candidate
    raise ValueError("Cannot locate workflow root containing scripts/run_point.py and protocol.json")


def validate_spec(spec):
    for key in ("method", "structure", "ecut_eV", "kmesh", "gpts", "xc_stencil"):
        if key not in spec:
            raise ValueError(f"Explicit input required: {key}")
    if spec["method"] != "PBE-D3(BJ)":
        raise ValueError("This fixed-cell runner is restricted to PBE-D3(BJ)")
    if not np.isfinite(float(spec["ecut_eV"])) or float(spec["ecut_eV"]) <= 0:
        raise ValueError("ecut_eV must be positive")
    for key in ("kmesh", "gpts"):
        value = np.asarray(spec[key])
        if value.shape != (3,) or np.any(value <= 0) or np.any(value != value.astype(int)):
            raise ValueError(f"{key} must contain three positive integers")
    if int(spec["xc_stencil"]) != spec["xc_stencil"] or spec["xc_stencil"] < 1:
        raise ValueError("xc_stencil must be a positive integer")
    if spec.get("h_A") is not None or spec.get("reference_cell_A") is not None:
        raise ValueError("Use explicit fixed gpts and an unscaled cutoff; h_A/reference_cell_A are unsupported here")
    fmax = float(spec.get("fmax_eV_A", spec.get("fmax", .002)))
    if not 0 < fmax <= .002:
        raise ValueError("Atomic-force target must be positive and no looser than 0.002 eV/A")
    if spec.get("relax_cell", False) is not False:
        raise ValueError("This runner never changes the cell")
    for key, maximum in (("scf_energy_eV_electron", 1e-9),
                         ("scf_density_electron_electron", 1e-9),
                         ("scf_eigenstates_eV2_electron", 1e-11)):
        if key not in spec or not 0 < float(spec[key]) <= maximum:
            raise ValueError(f"Explicit {key} must be positive and <= {maximum}")
    for key in ("source_geometry_sha256_expected", "protocol_sha256_expected", "point_runner_sha256"):
        if not isinstance(spec.get(key), str) or len(spec[key]) != 64:
            raise ValueError(f"Required SHA256 input: {key}")
    steps = int(spec.get("maxsteps", 200))
    if not 1 <= steps <= 200:
        raise ValueError("maxsteps must be in 1..200")
    if not 0 < float(spec.get("maxstep", .05)) <= .05:
        raise ValueError("Optimizer maxstep must be in (0, 0.05]")
    if not 0 < float(spec.get("symprec_A", 1e-5)) <= 1e-4:
        raise ValueError("symprec_A must be in (0, 1e-4]")
    return fmax, steps


def symmetry_info(atoms, symprec, hall_number=297):
    data = spglib.get_symmetry_dataset(
        (atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers),
        symprec=symprec, hall_number=hall_number)
    if data is None:
        return {"valid": False, "requested_hall_number": hall_number}
    return {
        "valid": True, "number": int(data.number), "international": data.international,
        "hall_number": int(data.hall_number), "hall": data.hall, "choice": data.choice,
        "transformation_matrix": data.transformation_matrix.tolist(),
        "origin_shift": data.origin_shift.tolist(), "operations": len(data.rotations),
    }


def pnam_closure(atoms):
    """Test Hall297 database operations in the original, unchanged input axes."""
    from scipy.optimize import linear_sum_assignment
    operations = spglib.get_symmetry_from_database(297)
    scaled = atoms.get_scaled_positions(wrap=False)
    maximum = 0.0
    for rotation, translation in zip(operations["rotations"], operations["translations"]):
        transformed = scaled @ rotation.T + translation
        for number in set(atoms.numbers):
            indices = np.flatnonzero(atoms.numbers == number)
            delta = transformed[indices, None, :] - scaled[None, indices, :]
            delta -= np.rint(delta)
            distances = np.linalg.norm(delta @ atoms.cell.array, axis=2)
            rows, columns = linear_sum_assignment(distances)
            maximum = max(maximum, float(distances[rows, columns].max()))
    return maximum


def terminal_symmetry(atoms, symprec):
    hall = symmetry_info(atoms, symprec)
    standard = symmetry_info(atoms, symprec, hall_number=0)
    passed = (
        hall.get("valid", False) and hall.get("number") == 62
        and hall.get("operations") == 8 and standard.get("number") == 62
        and standard.get("operations") == 8
        and np.allclose(hall["transformation_matrix"], np.eye(3), atol=1e-8, rtol=0)
        and pnam_closure(atoms) <= symprec)
    return {"passed": bool(passed), "Hall297_input_setting": hall,
            "default_standard_setting": standard,
            "policy": "Require #62 with eight operations; a higher-symmetry endpoint is recorded but not automatically accepted."}


def clean_copy(atoms):
    result = atoms.copy()
    result.set_constraint([])
    result.calc = None
    return result


def prepare_atoms(structure, symprec):
    # Atoms input is used only for validated pending resume coordinates.
    atoms = structure.copy() if isinstance(structure, Atoms) else read(structure)
    # File-based ASE I/O remains collective: every MPI rank enters this call.
    atoms.pbc = True
    atoms.set_constraint([])
    original = clean_copy(atoms)
    before = symmetry_info(atoms, symprec)
    default_before = symmetry_info(atoms, symprec, hall_number=0)
    if (not before["valid"] or before["number"] != 62 or before["operations"] != 8
            or not np.allclose(before["transformation_matrix"], np.eye(3), atol=1e-8, rtol=0)
            or default_before.get("number") != 62 or default_before.get("operations") != 8):
        raise ValueError("Initial structure must be #62 in Pnam/Hall297 input axes, with exactly 8 operations")
    if set(atoms.numbers) != {1, 6} or sum(atoms.numbers == 1) != 2 * sum(atoms.numbers == 6):
        raise ValueError("Expected periodic PE with C:H = 1:2")
    constraint = FixSymmetry(atoms, symprec=symprec, adjust_positions=True,
                             adjust_cell=False, verbose=False)
    cell_refinement_error = float(np.max(np.abs(atoms.cell.array - original.cell.array)))
    if cell_refinement_error > 1e-10:
        raise ValueError("Symmetry preparation would change the frozen input cell")
    # Remove floating-point-only lattice changes made inside FixSymmetry's
    # initialization; scale_atoms=False preserves its refined positions.
    atoms.set_cell(original.cell.array.copy(), scale_atoms=False, apply_constraint=False)
    atoms.set_constraint(constraint)
    after = symmetry_info(atoms, symprec)
    default_after = symmetry_info(atoms, symprec, hall_number=0)
    if (not after["valid"] or after["number"] != 62 or len(constraint.rotations) != 8
            or default_after.get("number") != 62 or default_after.get("operations") != 8):
        raise ValueError("FixSymmetry refinement changed or raised the initial symmetry")
    closure = pnam_closure(atoms)
    if closure > symprec:
        raise ValueError("Refined structure fails original-axis Pnam operation closure")
    movement = atoms.positions - original.positions
    refinement = {
        "before_Hall297": before, "after_Hall297": after,
        "before_default_standard_setting": default_before,
        "after_default_standard_setting": default_after,
        "max_initial_refinement_displacement_A": float(np.linalg.norm(movement, axis=1).max()),
        "initial_refinement_cell_delta_A": atoms.cell.array - original.cell.array,
        "Pnam_operation_max_residual_A": closure,
        "fixed_rotations": np.asarray(constraint.rotations),
        "fixed_translations": np.asarray(constraint.translations),
        "note": "Pnam is the input setting; Pnma may be reported as the equivalent standard setting. No axis permutation is applied.",
    }
    return atoms, original, refinement


def prepare_resume(spec, initial_atoms, original, symprec, maxsteps):
    """Validate a pending-geometry/LBFGS pair without any energy or force calls."""
    if "resume" not in spec:
        return initial_atoms, None
    resume = spec["resume"]
    required = {"geometry_json", "geometry_sha256", "optimizer_restart",
                "optimizer_restart_sha256", "previous_evaluation_json",
                "previous_evaluation_sha256", "optimizer_step_offset", "previous_attempt"}
    if not isinstance(resume, dict) or set(resume) != required:
        raise ValueError("Resume object must contain exactly the documented eight fields")
    offset = resume["optimizer_step_offset"]
    if type(offset) is not int or not 1 <= offset < maxsteps:
        raise ValueError("Resume step offset must be an integer inside the original total-step budget")
    pairs = (("geometry_json", "geometry_sha256"),
             ("optimizer_restart", "optimizer_restart_sha256"),
             ("previous_evaluation_json", "previous_evaluation_sha256"))
    for key, hashkey in pairs:
        if not isinstance(resume[key], str) or not Path(resume[key]).is_absolute():
            raise ValueError("Resume paths must be absolute: " + key)
        if sha(resume[key]) != resume[hashkey]:
            raise ValueError("Resume asset hash mismatch: " + key)
    previous_attempt = Path(resume["previous_attempt"])
    if not previous_attempt.is_absolute():
        raise ValueError("previous_attempt must be absolute")
    old_spec = json.loads((previous_attempt / "input.json").read_text())
    physics = {k: v for k, v in spec.items() if k != "resume"}
    old_physics = {k: v for k, v in old_spec.items() if k != "resume"}
    if physics != old_physics:
        raise ValueError("Resume must preserve every original input field outside resume")
    geometry = json.loads(Path(resume["geometry_json"]).read_text())
    previous = json.loads(Path(resume["previous_evaluation_json"]).read_text())
    if (geometry.get("optimizer_step") != offset or geometry.get("evaluation_id") != offset
            or previous.get("optimizer_step") != offset - 1
            or previous.get("evaluation_id") != offset - 1):
        raise ValueError("Resume evaluation/optimizer step numbering mismatch")
    if geometry.get("phase") != "initial" or previous.get("phase") != "initial":
        raise ValueError("Only interrupted initial optimization phases are supported for this resume")
    if previous.get("SCF_converged") is not True:
        raise ValueError("Previous completed evaluation must have converged SCF")
    for key, name in (("geometry_json", f"evaluation_input_{offset:04d}.json"),
                      ("optimizer_restart", "optimizer_restart.json"),
                      ("previous_evaluation_json", f"evaluation_{offset-1:04d}.json")):
        if sha(previous_attempt / name) != sha(resume[key]):
            raise ValueError("Frozen resume asset differs from the paused attempt: " + key)
    positions = np.asarray(geometry["positions_A"], dtype=float)
    cell = np.asarray(geometry["cell_A"], dtype=float)
    previous_positions = np.asarray(previous["positions_A"], dtype=float)
    previous_forces = np.asarray(previous["projected"]["forces_eV_A"], dtype=float)
    n = len(original)
    for label, value, shape in (("pending positions", positions, (n, 3)),
                                ("pending cell", cell, (3, 3)),
                                ("previous positions", previous_positions, (n, 3)),
                                ("previous projected forces", previous_forces, (n, 3))):
        if value.shape != shape or not np.all(np.isfinite(value)):
            raise ValueError("Invalid resume array: " + label)
    if (not np.array_equal(geometry["atomic_numbers"], original.numbers)
            or not np.array_equal(geometry["pbc"], original.pbc)
            or not np.array_equal(cell, original.cell.array)
            or not np.array_equal(previous["cell_A"], original.cell.array)):
        raise ValueError("Resume must preserve original atomic order/elements/PBC and exact cell")
    pending = clean_copy(original)
    pending.set_positions(positions, apply_constraint=False)
    atoms, _, refinement = prepare_atoms(pending, symprec)
    refinement_error = float(np.max(np.abs(atoms.positions - positions)))
    if refinement_error > 1e-10:
        raise ValueError("FixSymmetry would significantly change pending resume coordinates")
    # Preserve the exact pending coordinates, including unwrapped atom images.
    atoms.set_positions(positions.copy(), apply_constraint=False)
    if not terminal_symmetry(atoms, symprec)["passed"]:
        raise ValueError("Pending geometry does not satisfy fixed-cell Pnam contract")
    previous_atoms = initial_atoms.copy()
    previous_atoms.set_positions(previous_positions, apply_constraint=False)
    if not terminal_symmetry(previous_atoms, symprec)["passed"]:
        raise ValueError("Previous completed geometry does not satisfy Pnam")
    # Constructor reads restart only; no calculator is attached and no quantum
    # or optimizer step is executed.  logfile=None avoids any output file.
    probe = LBFGS(atoms.copy(), restart=resume["optimizer_restart"],
                  logfile=None, trajectory=None, maxstep=float(spec.get("maxstep", .05)), comm=world)
    if not probe.load_restart or probe.state.iteration != offset:
        raise ValueError("ASE LBFGS restart iteration differs from offset")
    r0, f0 = np.asarray(probe.r0), np.asarray(probe.f0)
    if r0.shape != (3*n,) or f0.shape != (3*n,) or not np.all(np.isfinite([r0, f0])):
        raise ValueError("Invalid LBFGS r0/f0 shape or values")
    r0_error = float(np.max(np.abs(r0 - previous_positions.ravel())))
    f0_error = float(np.max(np.abs(f0 - previous_forces.ravel())))
    if r0_error > 1e-12 or f0_error > 1e-12:
        raise ValueError("LBFGS r0/f0 do not match previous positions/projected forces")
    history_length = min(offset - 1, probe.state.memory)
    if not (len(probe.state.s) == len(probe.state.y) == len(probe.state.rho) == history_length):
        raise ValueError("LBFGS history length does not match iteration")
    rho_error = 0.0
    for s, y, rho in zip(probe.state.s, probe.state.y, probe.state.rho):
        s, y = np.asarray(s), np.asarray(y)
        if (s.shape != (3*n,) or y.shape != (3*n,)
                or not np.all(np.isfinite(s)) or not np.all(np.isfinite(y)) or not np.isfinite(rho)):
            raise ValueError("Invalid LBFGS history arrays")
        sy = float(np.dot(s, y))
        rho_error = max(rho_error, abs(float(rho) * sy - 1.0))
        if sy <= 0 or rho <= 0 or abs(float(rho) * sy - 1.0) > 1e-10:
            raise ValueError("LBFGS curvature history is invalid")
    # Reconstruct the already-proposed pending step.  Restart iteration is one
    # beyond the compute_step invocation that produced it; decrement only this
    # disposable in-memory probe.  The source restart is never written.
    probe.state.iteration = offset - 1
    proposed_step = probe.determine_step(probe.state.compute_step(-f0)) * probe.damping
    previous_atoms.set_positions(previous_positions + proposed_step.reshape(n, 3))
    replay_error = float(np.max(np.abs(previous_atoms.positions - positions)))
    if replay_error > 1e-10:
        raise ValueError("Pending geometry is not the LBFGS step implied by the saved history")
    probe.close()
    return atoms, {
        "validated": True, "optimizer_step_offset": offset,
        "first_evaluation_id": offset, "remaining_step_budget": maxsteps-offset,
        "electronic_restart": "cold_start_no_density_or_wavefunction_reuse",
        "geometry_sha256": resume["geometry_sha256"],
        "optimizer_restart_sha256": resume["optimizer_restart_sha256"],
        "previous_evaluation_sha256": resume["previous_evaluation_sha256"],
        "previous_attempt": str(previous_attempt), "history_length": history_length,
        "r0_previous_positions_max_delta_A": r0_error,
        "f0_previous_projected_forces_max_delta_eV_A": f0_error,
        "rho_sy_max_error": rho_error, "pending_step_reconstruction_max_delta_A": replay_error,
        "resume_preparation_max_delta_A": refinement_error,
        "pending_coordinates_preserved_exactly": True, "symmetry_refinement": refinement,
    }


class AuditedSumCalculator(SumCalculator):
    """Always compute and archive a complete raw property set at each geometry."""
    def __init__(self, dft, dispersion, context, callback):
        dispersion.implemented_properties = list(dispersion.implemented_properties)
        if "free_energy" not in dispersion.implemented_properties:
            dispersion.implemented_properties.append("free_energy")
        super().__init__([dft, dispersion])
        self.context = context
        self.callback = callback

    def calculate(self, atoms=None, properties=None, system_changes=None):
        from ase.calculators.calculator import all_changes
        fixed_cell_error = float(np.max(np.abs(atoms.cell.array - self.context["fixed_cell_A"])))
        if fixed_cell_error > 1e-10:
            raise RuntimeError(f"Frozen cell changed by {fixed_cell_error} A")
        self.context["evaluation_id"] += 1
        self.context["evaluating_cell_A"] = atoms.cell.array.copy()
        self.context["evaluating_positions_A"] = atoms.positions.copy()
        root_write(self.context["work"] / f"evaluation_input_{self.context['evaluation_id']:04d}.json",
                   {"evaluation_id": self.context["evaluation_id"], "phase": self.context["phase"],
                    "optimizer_step": getattr(self.context["optimizer"], "nsteps", 0),
                    "cell_A": atoms.cell.array, "positions_A": atoms.positions,
                    "atomic_numbers": atoms.numbers, "pbc": atoms.pbc})
        super().calculate(atoms, ["energy", "free_energy", "forces", "stress"],
                          all_changes if system_changes is None else system_changes)
        self.callback(atoms, self)


def load_point_module(path):
    module_spec = importlib.util.spec_from_file_location("frozen_pe_point_runner", path)
    module = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(module)
    return module


def check_dispersion(dispersion):
    params = dispersion.parameters
    cutoffs = params["realspace_cutoff"]
    required = {"disp2": 60, "disp3": 40, "cn": 60, "width2": .05, "width3": .05}
    if params.get("damping") != "d3bj" or params.get("params_tweaks") != {"method": "PBE", "atm": True}:
        raise ValueError("D3(BJ) PBE+ATM parameters do not match the frozen method")
    for name, expected in required.items():
        if not np.isclose(cutoffs[name], expected * Bohr, atol=1e-12, rtol=0):
            raise ValueError(f"Unexpected D3 cutoff {name}: {cutoffs[name]}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_json", type=Path)
    parser.add_argument("--check-input", action="store_true",
                        help="Validate structure/symmetry and construct interfaces; no SCF, optimization, or output files")
    args = parser.parse_args()
    input_path = args.input_json.resolve()
    spec = json.loads(input_path.read_text())
    force_target, maxsteps = validate_spec(spec)
    work = input_path.parent
    root = project_root(Path(__file__).resolve().parent)
    structure = Path(spec["structure"])
    if not structure.is_absolute():
        structure = (work / structure).resolve()
    symprec = float(spec.get("symprec_A", 1e-5))
    point_path = Path(spec.get("point_runner", root / "scripts/run_point.py")).resolve()
    if spec.get("point_runner_sha256") and sha(point_path) != spec["point_runner_sha256"]:
        raise ValueError("Point-runner SHA256 differs from the requested version")
    if sha(structure) != spec["source_geometry_sha256_expected"]:
        raise ValueError("Source geometry SHA256 differs from the frozen input")
    if sha(root / "protocol.json") != spec["protocol_sha256_expected"]:
        raise ValueError("Protocol SHA256 differs from the frozen input")
    if int(spec.get("mpi_ranks", world.size)) != world.size and not args.check_input:
        raise ValueError("Actual MPI size differs from requested mpi_ranks")
    atoms, original, refinement = prepare_atoms(structure, symprec)
    atoms, resume_validation = prepare_resume(spec, atoms, original, symprec, maxsteps)
    step_offset = 0 if resume_validation is None else resume_validation["optimizer_step_offset"]
    if args.check_input:
        point = load_point_module(point_path)
        dft = point.make_calculator(spec, os.devnull)
        dispersion = point.dispersion_calculator(spec["method"])
        check_dispersion(dispersion)
        if world.rank == 0:
            print(json.dumps({"status": "input_and_interfaces_validated_no_quantum_calculation",
                              "refinement": refinement, "resume_validation": resume_validation,
                              "dft_class": str(type(dft)),
                              "dispersion": dict(dispersion.parameters)}, default=convert, indent=2))
        return

    # Only the designated new attempt directory receives outputs.
    guarded = ["result.json", "status.json", "trajectory.traj", "gpaw.txt",
               "scf_iterations.jsonl", "evaluations.jsonl", "optimizer.log",
               "optimizer_restart.json", "gpaw_final_cold.txt", "final_cold_scf.jsonl",
               "resume_contract.json", "resume_geometry.json", "resume_previous_evaluation.json",
               "resume_optimizer_restart_seed.json"]
    if (any((work / name).exists() for name in guarded)
            or any(work.glob("evaluation_*.json")) or any(work.glob("evaluation_input_*.json"))):
        raise FileExistsError("Existing relaxation outputs found; create a new attempt directory")
    archive_error = None
    if world.rank == 0:
        try:
            for source, name in [(Path(__file__).resolve(), "run_relax_snapshot.py"),
                                 (point_path, "run_point_snapshot.py"),
                                 (root / "protocol.json", "protocol_snapshot.json"),
                                 (root / "environment/environment_verified.json", "environment_snapshot.json")]:
                target = work / name
                if source != target:
                    if target.exists() and sha(target) != sha(source):
                        raise FileExistsError(f"Different existing provenance snapshot: {target}")
                    if not target.exists():
                        shutil.copy2(source, target)
            if resume_validation is not None:
                resume = spec["resume"]
                for key, hashkey, name in (
                    ("geometry_json", "geometry_sha256", "resume_geometry.json"),
                    ("previous_evaluation_json", "previous_evaluation_sha256", "resume_previous_evaluation.json"),
                    ("optimizer_restart", "optimizer_restart_sha256", "resume_optimizer_restart_seed.json"),
                    ("optimizer_restart", "optimizer_restart_sha256", "optimizer_restart.json")):
                    if sha(resume[key]) != resume[hashkey]:
                        raise ValueError("Resume asset changed before copy: " + key)
                    shutil.copyfile(resume[key], work / name)
                    if sha(work / name) != resume[hashkey]:
                        raise ValueError("Resume copy hash mismatch: " + name)
        except Exception:
            archive_error = traceback.format_exc()
    archive_error = broadcast(archive_error, root=0, comm=world)
    if archive_error:
        raise OSError(archive_error)
    world.barrier()
    point = load_point_module(work / "run_point_snapshot.py")
    started = time.time()
    context = {"evaluation_id": step_offset - 1, "phase": "initial", "optimizer": None, "work": work,
               "fixed_cell_A": original.cell.array.copy()}
    state = {
        "status": "prepared", "input": spec, "input_sha256": sha(input_path),
        "runner_sha256": sha(__file__), "point_runner_sha256": sha(work / "run_point_snapshot.py"),
        "protocol_sha256": sha(work / "protocol_snapshot.json"),
        "environment_manifest_sha256": sha(work / "environment_snapshot.json"),
        "source_geometry_sha256": sha(structure), "source_geometry_path": str(structure),
        "software": {**point.versions(), "spglib": spglib.__version__},
        "mpi_world_size": world.size, "start_unix": started,
        "units": {"energy": "eV/cell", "forces": "eV/Angstrom",
                  "stress": "eV/Angstrom^3; ASE tensile-positive",
                  "stress_voigt_order": ["xx", "yy", "zz", "yz", "xz", "xy"]},
        "symmetry_refinement": refinement, "resume_validation": resume_validation,
        "optimizer": {"name": "LBFGS", "maxstep": float(spec.get("maxstep", .05)),
                      "maxsteps": maxsteps, "filter": None, "cell_fixed": True},
        "acceptance": {"raw_atomic_fmax_eV_A": force_target,
                       "stress_role": "Recorded only; fixed cell need not have zero stress",
                       "fixed_cell_max_delta_A": 1e-10, "requires_final_cold_check": True,
                       "requires_optimizer_convergence": True, "requires_Pnam_closure": True},
        "cell_comparison": "Frozen Tashiro cell; this run optimizes only atomic positions within Pnam.",
    }
    root_write(work / "status.json", state)
    root_write(work / "provenance.json", state)
    if resume_validation is not None:
        root_write(work / "resume_contract.json", resume_validation)
    write(work / "input_structure.extxyz", original)
    write(work / "symmetrized_initial.extxyz", clean_copy(atoms))
    dft = dispersion = trajectory = None
    latest = {}
    exit_code = 1
    try:
        dft = point.make_calculator(spec, work / "gpaw.txt")
        dispersion = point.dispersion_calculator(spec["method"])
        check_dispersion(dispersion)
        trajectory = Trajectory(work / "trajectory.traj", "w",
                                properties=["energy", "free_energy", "forces", "stress"],
                                comm=world, master=world.rank == 0)

        def record_iteration():
            # All ranks enter this observer.  Only root writes the local
            # Hamiltonian scalars; root_write broadcasts any I/O exception.
            item = {"evaluation_id": context["evaluation_id"], "phase": context["phase"],
                    "optimizer_step": getattr(context["optimizer"], "nsteps", 0),
                    "SCF_iteration": getattr(dft.scf, "niter", None),
                    "wall_seconds": time.time() - started}
            for name in ("e_total_free", "e_total_extrapolated", "e_kinetic",
                         "e_coulomb", "e_zero", "e_external", "e_xc", "e_entropy"):
                value = getattr(dft.hamiltonian, name, None)
                if value is not None:
                    item[name + "_eV"] = float(value * Hartree)
            scf_log = "final_cold_scf.jsonl" if context["phase"] == "final_cold" else "scf_iterations.jsonl"
            root_write(work / scf_log, item, append=True)

        dft.attach(record_iteration, 1)

        def capture(current, combined):
            raw = {key: np.array(combined.results[key]).copy() if key in ("forces", "stress")
                   else float(combined.results[key]) for key in ("energy", "free_energy", "forces", "stress")}
            if not bool(dft.scf.converged) or not all(np.isfinite(v).all() for v in raw.values()):
                raise RuntimeError("Nonfinite properties or unconverged SCF during relaxation")
            if not np.isclose(dispersion.results["free_energy"], dispersion.results["energy"],
                              rtol=0, atol=1e-12):
                raise RuntimeError("Dispersion free_energy must equal its geometry-only energy")
            projected_forces, projected_stress = raw["forces"].copy(), raw["stress"].copy()
            for constraint in current.constraints:
                constraint.adjust_forces(current, projected_forces)
                if hasattr(constraint, "adjust_stress"):
                    constraint.adjust_stress(current, projected_stress)
            closure = pnam_closure(current)
            item = {
                "evaluation_id": context["evaluation_id"], "phase": context["phase"],
                "optimizer_step": getattr(context["optimizer"], "nsteps", 0),
                "wall_seconds": time.time() - started,
                "SCF_converged": bool(dft.scf.converged), "SCF_iterations": int(dft.scf.niter),
                "cell_A": current.cell.array.copy(), "cell_lengths_A": current.cell.lengths(),
                "cell_angles_deg": current.cell.angles(), "positions_A": current.positions.copy(),
                "volume_A3": current.get_volume(), "natoms": len(current),
                "total": {"energy_extrapolated_eV": raw["energy"], "energy_free_eV": raw["free_energy"],
                          "forces_eV_A": raw["forces"], "stress_eV_A3": raw["stress"],
                          "stress_kbar": raw["stress"] / GPa * 10},
                "projected": {"forces_eV_A": projected_forces,
                              "stress_eV_A3": projected_stress, "stress_kbar": projected_stress / GPa * 10},
                "raw_max_force_norm_eV_A": float(np.linalg.norm(raw["forces"], axis=1).max()),
                "raw_max_stress_component_kbar": float(np.abs(raw["stress"] / GPa * 10).max()),
                "Pnam_operation_max_residual_A": closure,
                "fixed_cell_max_delta_A": float(np.max(np.abs(current.cell.array - context["fixed_cell_A"]))),
                "dft": {"energy_extrapolated_eV": float(dft.results["energy"]),
                        "energy_free_eV": float(dft.results["free_energy"]),
                        "forces_eV_A": np.array(dft.results["forces"]).copy(),
                        "stress_eV_A3": np.array(dft.results["stress"]).copy()},
                "dispersion": {"energy_eV": float(dispersion.results["energy"]),
                               "free_energy_eV": float(dispersion.results["free_energy"]),
                               "forces_eV_A": np.array(dispersion.results["forces"]).copy(),
                               "stress_eV_A3": np.array(dispersion.results["stress"]).copy(),
                               "calculator_parameters": dict(dispersion.parameters)},
            }
            local_basis = point.basis_metadata(dft)
            actual_xc = dft.hamiltonian.xc
            semilocal_xc = getattr(actual_xc, "semilocal_xc", actual_xc)
            local_basis["actual_XC_stencil"] = {
                "requested": spec["xc_stencil"],
                "xc_stencil_range": getattr(actual_xc, "stencil_range", None),
                "semilocal_stencil_range": getattr(semilocal_xc, "stencil_range", None)}
            item["basis_by_MPI_rank"] = [
                broadcast(local_basis if world.rank == rank else None, root=rank, comm=world)
                for rank in range(world.size)]
            if "datasets" not in state:
                state["datasets"] = [
                    {"symbol": setup.symbol, "type": setup.type, "filename": setup.filename,
                     "setupname": getattr(getattr(setup, "data", None), "setupname", None),
                     "fingerprint": getattr(setup, "fingerprint", None),
                     "file_sha256": sha(setup.filename)} for setup in dft.wfs.setups]
                root_write(work / "datasets.json", state["datasets"])
            root_write(work / "evaluations.jsonl", item, append=True)
            root_write(work / f"evaluation_{context['evaluation_id']:04d}.json", item)
            snapshot = clean_copy(current)
            snapshot.calc = SinglePointCalculator(snapshot, **raw)
            trajectory.write(snapshot)  # Every rank, with cached raw properties.
            write(work / "last_evaluated_structure.extxyz", snapshot)
            latest.clear()
            latest.update(item)
            state["last_evaluation_id"] = context["evaluation_id"]
            state["status"] = "final_cold_verification" if context["phase"] == "final_cold" else "optimizing"
            if world.rank == 0:
                print(json.dumps({"event": "evaluation_complete", "phase": context["phase"],
                                  "evaluation_id": context["evaluation_id"],
                                  "optimizer_step": item["optimizer_step"],
                                  "energy_eV_cell": raw["energy"],
                                  "raw_fmax_eV_A": item["raw_max_force_norm_eV_A"],
                                  "wall_seconds": item["wall_seconds"]}), flush=True)
            root_write(work / "status.json", {"status": state["status"],
                       "last_evaluation_id": context["evaluation_id"],
                       "optimizer_step": item["optimizer_step"],
                       "raw_fmax_eV_A": item["raw_max_force_norm_eV_A"],
                       "raw_stress_max_kbar": item["raw_max_stress_component_kbar"]})
            if closure > symprec:
                raise RuntimeError("Pnam operation closure lost during optimization")

        atoms.calc = AuditedSumCalculator(dft, dispersion, context, capture)
        optimizer = LBFGS(atoms, logfile=str(work / "optimizer.log"),
                          restart=str(work / "optimizer_restart.json"), trajectory=None,
                          maxstep=float(spec.get("maxstep", .05)), comm=world)
        if resume_validation is not None:
            if not optimizer.load_restart or optimizer.state.iteration != step_offset:
                raise ValueError("Copied optimizer restart failed to restore global iteration")
            # ASE restart restores LBFGS history, but not Optimizer.nsteps.
            # Set it before any SCF; run() receives only the remaining budget.
            optimizer.nsteps = step_offset
        context["optimizer"] = optimizer
        phases = []
        for phase in range(4):
            current_fmax = force_target / 2**phase
            context["phase"] = "initial" if phase == 0 else f"tightened_{phase}"
            remaining = maxsteps - optimizer.nsteps
            if remaining <= 0:
                break
            optimizer_converged = bool(optimizer.run(fmax=current_fmax, steps=remaining))
            # Every rank makes these calls, and apply_constraint=False bypasses
            # FixSymmetry.  Complete raw properties were already cached.
            raw_forces = atoms.get_forces(apply_constraint=False)
            raw_stress = atoms.get_stress(apply_constraint=False)
            raw_fmax = float(np.linalg.norm(raw_forces, axis=1).max())
            raw_smax = float(np.abs(raw_stress / GPa * 10).max())
            symmetry_check = terminal_symmetry(atoms, symprec)
            cell_error = float(np.max(np.abs(atoms.cell.array - context["fixed_cell_A"])))
            accepted = (optimizer_converged and raw_fmax <= force_target
                        and cell_error <= 1e-10 and symmetry_check["passed"])
            phase_result = {"phase": context["phase"], "optimizer_fmax_target": current_fmax,
                            "fixed_cell_max_delta_A": cell_error,
                            "optimizer_converged": optimizer_converged,
                            "total_optimizer_steps": optimizer.nsteps,
                            "raw_fmax_eV_A": raw_fmax, "raw_stress_max_kbar": raw_smax,
                            "terminal_symmetry": symmetry_check,
                            "physical_acceptance_pass": accepted}
            phases.append(phase_result)
            root_write(work / "optimization_phases.json", phases)
            if accepted:
                state["status"], exit_code = "optimized_pending_final_cold", 0
                break
            if not optimizer_converged:
                break
        if exit_code:
            state["status"] = ("maxsteps_reached" if optimizer.nsteps >= maxsteps
                               else "failed_physical_acceptance")
        state["optimization_phases"] = phases
        state["optimizer_steps"] = optimizer.nsteps
        state["optimizer_steps_in_this_attempt"] = optimizer.nsteps - step_offset
        state["optimized_evaluation"] = dict(latest)
        if exit_code == 0:
            # Start an entirely fresh SCF calculator: do not reuse orbitals,
            # density or a GPW checkpoint from the optimization.
            before_cold = dict(latest)
            dft.write(str(work / "optimized_warm.gpw"), mode="all")
            context["phase"] = "final_cold"
            dft = point.make_calculator(spec, work / "gpaw_final_cold.txt")
            dispersion = point.dispersion_calculator(spec["method"])
            check_dispersion(dispersion)
            dft.attach(record_iteration, 1)
            atoms.calc = AuditedSumCalculator(dft, dispersion, context, capture)
            atoms.get_potential_energy(force_consistent=True)
            cold_forces = atoms.get_forces(apply_constraint=False)
            cold_symmetry = terminal_symmetry(atoms, symprec)
            cold_fmax = float(np.linalg.norm(cold_forces, axis=1).max())
            cell_error = float(np.max(np.abs(atoms.cell.array - context["fixed_cell_A"])))
            cold_pass = bool(dft.scf.converged and cold_fmax <= force_target
                             and cell_error <= 1e-10 and cold_symmetry["passed"])
            state["final_cold_verification"] = {
                "passed": cold_pass, "independent_density_and_wavefunctions": True,
                "raw_fmax_eV_A": cold_fmax, "fixed_cell_max_delta_A": cell_error,
                "SCF_converged": bool(dft.scf.converged), "terminal_symmetry": cold_symmetry,
                "energy_change_eV_cell": latest["total"]["energy_free_eV"] - before_cold["total"]["energy_free_eV"],
                "max_force_vector_change_eV_A": float(np.linalg.norm(
                    cold_forces - np.asarray(before_cold["total"]["forces_eV_A"]), axis=1).max()),
                "stress_change_kbar": np.asarray(latest["total"]["stress_kbar"]) - np.asarray(before_cold["total"]["stress_kbar"]),
                "stress_acceptance_required": False,
            }
            root_write(work / "final_cold_verification.json", state["final_cold_verification"])
            state["status"] = "passed_relaxation" if cold_pass else "needs_attention_final_cold_acceptance"
            exit_code = 0 if cold_pass else 1
        state["final_evaluation"] = dict(latest)
        state["final_symmetry"] = terminal_symmetry(atoms, symprec)
        final = clean_copy(atoms)
        # All ranks enter ASE I/O and GPAW checkpoint operations in this order.
        write(work / "final_structure.extxyz", final)
        write(work / "final_POSCAR", final, format="vasp", direct=True, sort=False, vasp5=True)
        write(work / "final.cif", final, format="cif")
        state["CIF_representation"] = "Full-cell explicit atoms (ASE P1 header); actual Pnam closure is separately verified."
        dft.write(str(work / "final.gpw"), mode="all")
        state["final_geometry_sha256"] = sha(work / "final_structure.extxyz")
        if resume_validation is not None:
            for key, hashkey in (("geometry_json", "geometry_sha256"),
                                 ("optimizer_restart", "optimizer_restart_sha256"),
                                 ("previous_evaluation_json", "previous_evaluation_sha256")):
                if sha(spec["resume"][key]) != spec["resume"][hashkey]:
                    raise ValueError("Frozen resume asset changed during calculation: " + key)
        if sha(structure) != spec["source_geometry_sha256_expected"]:
            raise ValueError("Frozen source geometry changed during calculation")
        if sha(root / "protocol.json") != spec["protocol_sha256_expected"]:
            raise ValueError("Frozen protocol changed during calculation")
    except Exception as exc:
        state["status"] = "failed_calculation_or_io"
        state["error"] = repr(exc)
        state["traceback"] = traceback.format_exc()
        state["last_successful_evaluation"] = latest
        exit_code = 1
        if world.rank == 0:
            print(state["traceback"], file=sys.stderr, flush=True)
    finally:
        if trajectory is not None:
            trajectory.close()
        state["wall_seconds"] = time.time() - started
        state["exit_code"] = exit_code
        root_write(work / "result.json", state)
        root_write(work / "status.json", {"status": state["status"], "exit_code": exit_code,
                   "wall_seconds": state["wall_seconds"],
                   "optimizer_steps": state.get("optimizer_steps"),
                   "last_evaluation_id": state.get("last_evaluation_id")})
    if world.rank == 0:
        print(json.dumps({"status": state["status"], "exit_code": exit_code,
                          "optimizer_steps": state.get("optimizer_steps"),
                          "wall_seconds": state["wall_seconds"]}), flush=True)
    if exit_code:
        raise SystemExit(exit_code)


if __name__ == "__main__":
    main()
