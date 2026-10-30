"""Auditable Pnam PE cell relaxation; launch identically on every MPI rank.

No calculation is started by importing this module or using --check-input.
An existing output directory is never overwritten.  Use a fresh attempt with
the previous final_structure.extxyz as input to continue an incomplete run.
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
from ase.calculators.mixing import SumCalculator
from ase.calculators.singlepoint import SinglePointCalculator
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
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
        raise ValueError("This cell-comparison runner is restricted to PBE-D3(BJ)")
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
    stress = float(spec.get("stress_max_kbar", .05))
    if not 0 < stress <= .05:
        raise ValueError("Stress target must be positive and no looser than 0.05 kbar")
    steps = int(spec.get("maxsteps", 200))
    if not 1 <= steps <= 200:
        raise ValueError("maxsteps must be in 1..200")
    if not 0 < float(spec.get("maxstep", .05)) <= .05:
        raise ValueError("Optimizer maxstep must be in (0, 0.05]")
    if not 0 < float(spec.get("symprec_A", 1e-5)) <= 1e-4:
        raise ValueError("symprec_A must be in (0, 1e-4]")
    return fmax, stress, steps


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
    atoms = read(structure)  # Collective ASE I/O: every MPI rank calls this.
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
                             adjust_cell=True, verbose=False)
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
    force_target, stress_target, maxsteps = validate_spec(spec)
    work = input_path.parent
    root = project_root(Path(__file__).resolve().parent)
    structure = Path(spec["structure"])
    if not structure.is_absolute():
        structure = (work / structure).resolve()
    symprec = float(spec.get("symprec_A", 1e-5))
    point_path = Path(spec.get("point_runner", root / "scripts/run_point.py")).resolve()
    if spec.get("point_runner_sha256") and sha(point_path) != spec["point_runner_sha256"]:
        raise ValueError("Point-runner SHA256 differs from the requested version")
    atoms, original, refinement = prepare_atoms(structure, symprec)
    if args.check_input:
        point = load_point_module(point_path)
        dft = point.make_calculator(spec, os.devnull)
        dispersion = point.dispersion_calculator(spec["method"])
        check_dispersion(dispersion)
        if world.rank == 0:
            print(json.dumps({"status": "input_and_interfaces_validated_no_quantum_calculation",
                              "refinement": refinement, "dft_class": str(type(dft)),
                              "dispersion": dict(dispersion.parameters)}, default=convert, indent=2))
        return

    # Only the designated new attempt directory receives outputs.
    guarded = ["result.json", "status.json", "trajectory.traj", "gpaw.txt",
               "scf_iterations.jsonl", "evaluations.jsonl", "optimizer.log",
               "optimizer_restart.json"]
    if any((work / name).exists() for name in guarded):
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
        except Exception:
            archive_error = traceback.format_exc()
    archive_error = broadcast(archive_error, root=0, comm=world)
    if archive_error:
        raise OSError(archive_error)
    world.barrier()
    point = load_point_module(work / "run_point_snapshot.py")
    started = time.time()
    context = {"evaluation_id": -1, "phase": "initial", "optimizer": None, "work": work}
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
        "symmetry_refinement": refinement,
        "optimizer": {"name": "LBFGS", "maxstep": float(spec.get("maxstep", .05)),
                      "maxsteps": maxsteps, "filter": "FrechetCellFilter",
                      "mask": [1, 1, 1, 0, 0, 0], "exp_cell_factor": 1.0},
        "acceptance": {"raw_atomic_fmax_eV_A": force_target, "raw_stress_max_kbar": stress_target,
                       "requires_optimizer_convergence": True, "requires_Pnam_closure": True},
        "cell_comparison": "This run alone does not establish 0.5% agreement; compare independently accepted runs at two precisions.",
    }
    root_write(work / "status.json", state)
    root_write(work / "provenance.json", state)
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
            root_write(work / "scf_iterations.jsonl", item, append=True)

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
            latest.clear()
            latest.update(item)
            state["last_evaluation_id"] = context["evaluation_id"]
            state["status"] = "optimizing"
            root_write(work / "status.json", {"status": state["status"],
                       "last_evaluation_id": context["evaluation_id"],
                       "optimizer_step": item["optimizer_step"],
                       "raw_fmax_eV_A": item["raw_max_force_norm_eV_A"],
                       "raw_stress_max_kbar": item["raw_max_stress_component_kbar"]})
            if closure > symprec:
                raise RuntimeError("Pnam operation closure lost during optimization")

        atoms.calc = AuditedSumCalculator(dft, dispersion, context, capture)
        cell_filter = FrechetCellFilter(atoms, mask=[1, 1, 1, 0, 0, 0], exp_cell_factor=1.0)
        optimizer = LBFGS(cell_filter, logfile=str(work / "optimizer.log"),
                          restart=str(work / "optimizer_restart.json"), trajectory=None,
                          maxstep=float(spec.get("maxstep", .05)), comm=world)
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
            accepted = (optimizer_converged and raw_fmax <= force_target
                        and raw_smax <= stress_target and symmetry_check["passed"])
            phase_result = {"phase": context["phase"], "filter_fmax_target": current_fmax,
                            "optimizer_converged": optimizer_converged,
                            "total_optimizer_steps": optimizer.nsteps,
                            "raw_fmax_eV_A": raw_fmax, "raw_stress_max_kbar": raw_smax,
                            "terminal_symmetry": symmetry_check,
                            "physical_acceptance_pass": accepted}
            phases.append(phase_result)
            root_write(work / "optimization_phases.json", phases)
            if accepted:
                state["status"], exit_code = "passed_relaxation", 0
                break
            if not optimizer_converged:
                break
        if exit_code:
            state["status"] = ("maxsteps_reached" if optimizer.nsteps >= maxsteps
                               else "failed_physical_acceptance")
        state["optimization_phases"] = phases
        state["optimizer_steps"] = optimizer.nsteps
        state["final_evaluation"] = latest
        state["final_symmetry"] = terminal_symmetry(atoms, symprec)
        final = clean_copy(atoms)
        # All ranks enter ASE I/O and GPAW checkpoint operations in this order.
        write(work / "final_structure.extxyz", final)
        write(work / "final_POSCAR", final, format="vasp", direct=True, sort=False, vasp5=True)
        write(work / "final.cif", final, format="cif")
        state["CIF_representation"] = "Full-cell explicit atoms (ASE P1 header); actual Pnam closure is separately verified."
        dft.write(str(work / "final.gpw"), mode="all")
        state["final_geometry_sha256"] = sha(work / "final_structure.extxyz")
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
