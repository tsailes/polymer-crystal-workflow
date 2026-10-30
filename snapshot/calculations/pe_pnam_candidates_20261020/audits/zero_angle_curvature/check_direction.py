"""Two single-point checks along one Pnam-preserving rigid-chain rotation."""
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

os.environ["OMP_NUM_THREADS"] = "4"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OMP_MAX_ACTIVE_LEVELS"] = "1"
sys.dont_write_bytecode = True
import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parents[2]
AUDIT = Path(__file__).resolve().parent
WORKFLOW = ROOT.parent / "pe_gfn2_first_run_20261020/convergence_relaxation/workflow.py"
spec = importlib.util.spec_from_file_location("pnam_generation", ROOT / "generate_candidates.py")
generator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(generator)


def load(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, data):
    Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def main():
    folder = AUDIT / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    folder.mkdir()
    candidate = ROOT / "candidates/cand_0001"
    status = load(candidate / "relaxation/status.json")
    if status["status"] != "preoptimized":
        raise ValueError("Expected completed zero-angle preoptimization")
    attempt = Path(status["attempt_directory"])
    source = attempt / "final.extxyz"
    atoms = read(source)
    atoms.calc = None
    atoms.set_constraint()
    atoms.info.clear()
    for name in ("raw_forces", "projected_forces"):
        if name in atoms.arrays:
            del atoms.arrays[name]
    topology_path = candidate / "reference_topology.json"
    topology = load(topology_path)
    construction = load(candidate / "construction.json")
    centers = np.column_stack([np.asarray(construction["chain_centers_fractional"]),
                               np.zeros(2)]) @ atoms.cell.array
    wrap = np.asarray(atoms.arrays["wrap_image"])
    unwrapped = atoms.positions + wrap @ atoms.cell.array
    chain_ids = atoms.arrays["chain_id"]
    reference_audit = generator.audit_geometry(atoms, topology)
    reference_bonds = np.array([b["distance_A"] for b in reference_audit["bonds"]])
    fine_source = attempt / "fine/result.json"
    fine = load(fine_source)
    cold_source = attempt / "final_result.json"
    cold = load(cold_source)
    zero_energy = fine["final"]["energy_eV_per_cell"]  # accuracy=0.01, identical positions
    zero_cold_energy = cold["energy_eV_per_cell"]  # accuracy=0.001 independent initialization
    fine_geometry = read(attempt / "fine/final.extxyz")
    if not np.array_equal(fine_geometry.positions, atoms.positions):
        raise ValueError("Stored accuracy=0.01 baseline is a different geometry")
    if not reference_audit["passed"]:
        raise ValueError("Baseline failed structure audit")
    delta_deg = 1.0
    results = {}
    metadata = dict(
        candidate="cand_0001", original_status_preserved=True, source_path=str(source),
        source_sha256=sha(source), topology_path=str(topology_path), topology_sha256=sha(topology_path),
        baseline_energy_path=str(fine_source), baseline_energy_sha256=sha(fine_source),
        baseline_cold_path=str(cold_source), baseline_cold_sha256=sha(cold_source),
        baseline_energy_accuracy_0p01_eV_cell=zero_energy,
        baseline_energy_accuracy_0p001_eV_cell=zero_cold_energy,
        baseline_accuracy_difference_eV_cell=zero_cold_energy-zero_energy,
        script_sha256=sha(__file__), workflow_sha256=sha(WORKFLOW),
        generator_sha256=sha(ROOT / "generate_candidates.py"),
        chain_centers_A=centers.tolist(), delta_degrees=delta_deg,
        perturbation="chain 0 rotates +delta, chain 1 -delta, about fixed c axes; unwrap then reapply original image representation",
        fixed_cell_A=atoms.cell.array.tolist(), reference_geometry_audit=reference_audit,
        scope="One finite angular direction; no full Hessian, phonons, or further optimization")
    dump(folder / "protocol.json", metadata)
    shutil.copy2(__file__, folder / "check_direction_source.py")
    write(folder / "baseline.extxyz", atoms)
    for sign, name in ((-1, "minus_1deg"), (1, "plus_1deg")):
        case = folder / name
        case.mkdir()
        displaced = atoms.copy()
        rotated = unwrapped.copy()
        for chain in (0, 1):
            angle = np.deg2rad(sign * delta_deg * (1 if chain == 0 else -1))
            cosine, sine = np.cos(angle), np.sin(angle)
            rotation = np.array([[cosine, -sine, 0], [sine, cosine, 0], [0, 0, 1]])
            mask = chain_ids == chain
            rotated[mask] = (unwrapped[mask] - centers[chain]) @ rotation.T + centers[chain]
        displaced.positions = rotated - wrap @ atoms.cell.array
        displaced.info.update(source_candidate="cand_0001", rigid_setting_perturbation_deg=sign*delta_deg)
        audit = generator.audit_geometry(displaced, topology)
        bond_change = float(np.max(np.abs(np.array(
            [b["distance_A"] for b in audit["bonds"]]) - reference_bonds)))
        valid = (audit["passed"] and bond_change < 1e-10
                 and np.array_equal(displaced.cell.array, atoms.cell.array))
        dump(case / "geometry_validation.json",
             dict(passed=bool(valid), maximum_reference_bond_length_change_A=bond_change,
                  cell_unchanged=np.array_equal(displaced.cell.array, atoms.cell.array), audit=audit))
        if not valid:
            raise ValueError(f"Invalid rigid Pnam perturbation: {name}")
        write(case / "input.extxyz", displaced)
        start = time.monotonic()
        with (case / "native_output.log").open("w") as output:
            completed = subprocess.run(
                [sys.executable, str(WORKFLOW), "singlepoint", "--input", str(case/"input.extxyz"),
                 "--output", str(case/"result.json"), "--repeat", "2", "3", "8", "--accuracy", "0.01"],
                stdout=output, stderr=subprocess.STDOUT, check=False, timeout=600)
        if completed.returncode:
            raise RuntimeError(f"{name} singlepoint failed with {completed.returncode}")
        result = load(case / "result.json")
        result["total_process_seconds"] = time.monotonic() - start
        result["energy_change_eV_cell"] = result["energy_per_primitive_cell_eV"] - zero_energy
        result["energy_change_meV_CH2"] = result["energy_change_eV_cell"] * 1000 / 4
        result["geometry_validation_passed"] = bool(valid)
        result["maximum_reference_bond_length_change_A"] = bond_change
        if result["replica_force_max_deviation_eV_per_A"] > 1e-5:
            raise RuntimeError("Equivalent-replica forces failed")
        results[name] = result
        dump(folder / "progress.json", results)
        print(name, result["energy_change_meV_CH2"], "meV/CH2", flush=True)
    minus = results["minus_1deg"]["energy_per_primitive_cell_eV"]
    plus = results["plus_1deg"]["energy_per_primitive_cell_eV"]
    second_difference = plus + minus - 2*zero_energy
    delta_rad = np.deg2rad(delta_deg)
    numerical_margin = max(1e-7, 10*abs(zero_cold_energy-zero_energy))
    negative = second_difference < -2*numerical_margin and max(plus, minus) < zero_energy-numerical_margin
    summary = dict(
        **metadata, perturbations=results, second_energy_difference_eV_cell=second_difference,
        rigid_rotation_curvature_eV_cell_rad2=second_difference/delta_rad**2,
        central_angular_derivative_eV_cell_rad=(plus-minus)/(2*delta_rad),
        numerical_energy_margin_eV_cell=numerical_margin,
        negative_curvature_and_both_sides_lower=bool(negative),
        classification="unstable_stationary_along_tested_direction" if negative else "inconclusive_or_nonnegative_along_tested_direction",
        recommendation=("Keep original record; exclude from stable-minimum list; retain perturbed structures for a separately authorized reoptimization."
                        if negative else "Do not infer complete stability from one direction."),
        limitation="Finite +/-1 degree rigid rotation of a force-tolerance-converged structure; not a complete Hessian eigenvalue calculation.")
    dump(folder / "summary.json", summary)
    dump(AUDIT / "latest_result.json", dict(run_directory=str(folder), summary=str(folder/"summary.json"),
                                           classification=summary["classification"]))
    print(json.dumps({k: summary[k] for k in ("classification", "second_energy_difference_eV_cell",
          "rigid_rotation_curvature_eV_cell_rad2", "central_angular_derivative_eV_cell_rad")}), flush=True)


if __name__ == "__main__":
    main()
