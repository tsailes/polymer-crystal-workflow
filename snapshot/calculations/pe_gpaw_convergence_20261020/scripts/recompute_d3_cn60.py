"""Recombine saved GPAW DFT with D3 CN=60 Bohr; never rerun any SCF.

Original run folders are read-only. By default the complete 37-point PBE-D3
derivative set is required (one center plus 36 displacement/strain endpoints).
"""
from __future__ import annotations

import argparse
import copy
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
import time

import numpy as np
from ase.io import read
from ase.stress import voigt_6_to_full_3x3_stress
from ase.units import Bohr, GPa, Ry
from dftd3.ase import DFTD3


ROOT = Path(__file__).resolve().parents[1]
NEW_CUTOFF = {"disp2": 60 * Bohr, "disp3": 40 * Bohr, "cn": 60 * Bohr,
              "width2": 0.05 * Bohr, "width3": 0.05 * Bohr}
OLD_CUTOFF = {**NEW_CUTOFF, "cn": 40 * Bohr}


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


def exact_geometry_digest(atoms):
    fields = {
        "atomic_numbers": atoms.numbers.tolist(),
        "positions_A": atoms.positions.tolist(),
        "cell_A": atoms.cell.array.tolist(),
        "pbc": atoms.pbc.tolist(),
    }
    return hashlib.sha256(json.dumps(
        fields, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def inspect_run(run_dir):
    result_path = run_dir / "result.json"
    structure_path = run_dir / "input_structure.extxyz"
    source = json.loads(result_path.read_text())
    if source.get("status") != "passed_single_point" or not source.get("SCF_converged"):
        raise ValueError(f"{run_dir.name}: original SCF not complete and converged")
    if source.get("input", {}).get("method") != "PBE-D3(BJ)":
        raise ValueError(f"{run_dir.name}: unexpected method")
    prior = source["dispersion"]["calculator_parameters"]
    if prior.get("damping") != "d3bj":
        raise ValueError(f"{run_dir.name}: prior damping is not d3bj")
    tweaks = prior.get("params_tweaks", {})
    if tweaks.get("method", "").lower() != "pbe" or tweaks.get("atm") is not True:
        raise ValueError(f"{run_dir.name}: prior PBE/ATM is not explicit")
    cutoff = prior.get("realspace_cutoff", {})
    for key, value in OLD_CUTOFF.items():
        if key not in cutoff or not np.isclose(cutoff[key], value, rtol=0, atol=1e-12):
            raise ValueError(f"{run_dir.name}: prior {key} does not match CN40 smooth baseline")
    atoms = read(structure_path)
    if len(atoms) != 12 or atoms.get_chemical_formula() != "C4H8" or not atoms.pbc.all():
        raise ValueError(f"{run_dir.name}: unexpected periodic PE cell")
    exact_positions = np.asarray(source["positions_A"], float)
    exact_cell = np.asarray(source["cell_A"], float)
    position_rounding = float(np.max(np.abs(atoms.positions - exact_positions)))
    cell_rounding = float(np.max(np.abs(atoms.cell.array - exact_cell)))
    if max(position_rounding, cell_rounding) > 1e-7:
        raise ValueError(f"{run_dir.name}: saved geometry does not match original result")
    # ASE's extxyz writer rounds position columns. Use the exact double-precision
    # positions and cell archived before the original DFT, preserving atom order.
    atoms.set_cell(exact_cell, scale_atoms=False)
    atoms.set_positions(exact_positions)
    for key in ("energy_free_eV", "energy_extrapolated_eV", "forces_eV_A", "stress_eV_A3"):
        if not np.isfinite(np.asarray(source["dft"][key])).all():
            raise ValueError(f"{run_dir.name}: nonfinite original DFT {key}")
    return source, atoms, {
        "source_result": str(result_path.resolve()),
        "source_result_sha256": sha(result_path),
        "same_geometry_file": str(structure_path.resolve()),
        "same_geometry_sha256": sha(structure_path),
        "exact_geometry_sha256": exact_geometry_digest(atoms),
        "exact_geometry_source": "original result.json positions_A/cell_A; species and order from input_structure.extxyz",
        "extxyz_position_rounding_max_A": position_rounding,
        "extxyz_cell_rounding_max_A": cell_rounding,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", type=Path, default=ROOT / "runs")
    parser.add_argument("--output-dir", type=Path,
                        default=ROOT / "method_audit" / "derivative_d3_cn60")
    parser.add_argument("--pattern", default="deriv_smooth_pbe_d3bj_*")
    parser.add_argument("--expected-count", type=int, default=37)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if importlib.metadata.version("dftd3") != "1.6.0":
        raise RuntimeError("This revision was audited with dftd3-python 1.6.0.")
    runs = sorted(p for p in args.runs.glob(args.pattern) if p.is_dir())
    prepared = []
    invalid = []
    for run_dir in runs:
        try:
            source, atoms, provenance = inspect_run(run_dir)
            prepared.append((run_dir, source, atoms, provenance))
        except (OSError, ValueError, KeyError) as exc:
            invalid.append({"run_id": run_dir.name, "error": str(exc)})
    preflight = {
        "matched": len(runs), "ready": len(prepared), "expected": args.expected_count,
        "invalid": invalid, "output_dir": str(args.output_dir.resolve()),
        "no_SCF": True, "cutoff_A": NEW_CUTOFF,
    }
    print(json.dumps(preflight, indent=2), flush=True)
    if args.dry_run:
        return
    if invalid or len(prepared) != args.expected_count:
        raise RuntimeError("Complete valid derivative set required; no outputs written.")
    if args.output_dir.resolve() == args.runs.resolve() or args.runs.resolve() in args.output_dir.resolve().parents:
        raise ValueError("Revision output must be outside the original runs directory.")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    script_sha = sha(__file__)
    records = []
    for run_dir, original, atoms, provenance in prepared:
        output = args.output_dir / f"{run_dir.name}.json"
        if output.exists():
            prior = json.loads(output.read_text())
            if (prior.get("status") == "passed_dispersion_recombination"
                    and prior.get("source_result_sha256") == provenance["source_result_sha256"]
                    and prior.get("exact_geometry_sha256") == provenance["exact_geometry_sha256"]
                    and prior.get("revision_id") == "D3_CN60_smooth005bohr"
                    and prior.get("dispersion", {}).get("calculator_parameters", {}).get(
                        "realspace_cutoff") == NEW_CUTOFF):
                records.append({"run_id": run_dir.name, "file": str(output),
                                "sha256": sha(output), "reused": True})
                continue
            raise FileExistsError(f"Refusing to overwrite a different revision: {output}")
        start = time.time()
        calc = DFTD3(damping="d3bj", params_tweaks={"method": "PBE", "atm": True},
                     realspace_cutoff=NEW_CUTOFF)
        atoms.calc = calc
        de = float(atoms.get_potential_energy())
        df = atoms.get_forces()
        ds = atoms.get_stress()
        if not all(np.isfinite(x).all() for x in (de, df, ds)):
            raise ValueError(f"{run_dir.name}: nonfinite new dispersion")
        dft = copy.deepcopy(original["dft"])
        e_free = float(dft["energy_free_eV"]) + de
        e_extrap = float(dft["energy_extrapolated_eV"]) + de
        forces = np.asarray(dft["forces_eV_A"]) + df
        stress = np.asarray(dft["stress_eV_A3"]) + ds
        result = {
            "status": "passed_dispersion_recombination",
            "revision_id": "D3_CN60_smooth005bohr",
            "original_run_id": run_dir.name,
            "operation": "Recomputed only geometry-dependent D3; reused unchanged original GPAW DFT E/F/stress. No SCF rerun.",
            "DFT_recomputed": False,
            **provenance,
            "revision_script_sha256": script_sha,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "input": copy.deepcopy(original["input"]),
            "natoms": len(atoms), "formula": atoms.get_chemical_formula(),
            "cell_A": atoms.cell.array, "positions_A": atoms.positions,
            "volume_A3": atoms.get_volume(),
            "SCF_converged": original["SCF_converged"],
            "SCF_iterations": original["SCF_iterations"],
            "basis": copy.deepcopy(original.get("basis")),
            "datasets": copy.deepcopy(original.get("datasets")),
            "original_dispersion": copy.deepcopy(original["dispersion"]),
            "original_total": copy.deepcopy(original["total"]),
            "dft": dft,
            "dispersion": {
                "energy_eV": de, "forces_eV_A": df, "stress_eV_A3": ds,
                "stress_kbar": ds / GPa * 10,
                "calculator_class": str(type(calc)),
                "calculator_parameters": dict(calc.parameters),
                "software": {n: importlib.metadata.version(n) for n in ("ase", "dftd3")},
            },
            "total": {
                "energy_free_eV": e_free, "energy_extrapolated_eV": e_extrap,
                "energy_Ry": e_free / Ry, "energy_eV_atom": e_free / len(atoms),
                "energy_eV_CH2": e_free / 4,
                "forces_eV_A": forces, "stress_eV_A3": stress,
                "stress_tensor_eV_A3": voigt_6_to_full_3x3_stress(stress),
                "stress_kbar": stress / GPa * 10,
                "max_force_norm_eV_A": np.linalg.norm(forces, axis=1).max(),
                "net_force_eV_A": forces.sum(axis=0),
            },
            "difference_from_original_total": {
                "energy_free_eV": e_free - original["total"]["energy_free_eV"],
                "max_force_vector_change_eV_A": np.linalg.norm(
                    forces - original["total"]["forces_eV_A"], axis=1).max(),
                "max_stress_component_change_kbar": np.abs(
                    stress - original["total"]["stress_eV_A3"]).max() / GPa * 10,
            },
            "units": copy.deepcopy(original["units"]),
            "wall_seconds": time.time() - start,
        }
        # Exclusive creation: preserve every existing revision and original file.
        with output.open("x", encoding="utf8") as handle:
            json.dump(result, handle, indent=2, default=convert)
            handle.write("\n")
        records.append({"run_id": run_dir.name, "file": str(output),
                        "sha256": sha(output), "reused": False})
        print(json.dumps({"run_id": run_dir.name, "status": result["status"],
                          "wall_seconds": result["wall_seconds"]}), flush=True)
    manifest = {
        "status": "complete", "revision_id": "D3_CN60_smooth005bohr",
        "count": len(records), "DFT_recomputed": False, "records": records,
        "script_sha256": script_sha,
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    # Unique manifest names preserve successive resumptions.
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    with (args.output_dir / f"manifest_{stamp}.json").open("x", encoding="utf8") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
