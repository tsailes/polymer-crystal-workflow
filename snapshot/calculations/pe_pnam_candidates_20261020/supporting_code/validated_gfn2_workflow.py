"""Periodic Gamma-supercell convergence and fixed-cell primitive-coordinate relaxation.

All reported energies are per C4H8 cell (or explicitly per CH2). Optimization
uses E_super / number_of_replicas and the corresponding replica-averaged forces.
No literature coordinates or symmetry operations enter the relaxation.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase.io import read, write
from ase.neighborlist import neighbor_list
from ase.optimize import BFGS
from tblite.ase import TBLite

ROOT = Path(__file__).resolve().parent
PILOT = ROOT.parent
SETTINGS = dict(method="GFN2-xTB", charge=0, multiplicity=1, accuracy=0.01,
                max_iterations=500, electronic_temperature=300.0, verbosity=1)
THRESHOLDS = dict(energy_eV_per_CH2=0.001, force_vector_eV_per_A=0.005,
                  relaxation_fmax_eV_per_A=0.01, replica_force_eV_per_A=1e-5)


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                          encoding="utf-8")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fold_result(primitive, repeat, supercell):
    replicas = int(np.prod(repeat))
    forces = supercell.get_forces().reshape(replicas, len(primitive), 3)
    mean_forces = forces.mean(axis=0)
    charges = supercell.get_charges().reshape(replicas, len(primitive))
    energy = float(supercell.get_potential_energy()) / replicas
    return dict(
        repeats=list(repeat), natoms=len(supercell), primitive_natoms=len(primitive),
        energy_per_primitive_cell_eV=energy, energy_per_CH2_eV=energy / (len(primitive) / 3),
        forces_eV_per_A=mean_forces.tolist(),
        max_force_norm_eV_per_A=float(np.linalg.norm(mean_forces, axis=1).max()),
        replica_force_max_deviation_eV_per_A=float(np.linalg.norm(forces - mean_forces, axis=2).max()),
        net_force_eV_per_A=mean_forces.sum(axis=0).tolist(),
        stress_eV_per_A3_voigt=supercell.get_stress().tolist(),
        mean_charges_e=charges.mean(axis=0).tolist(),
        replica_charge_max_deviation_e=float(np.max(np.abs(charges - charges.mean(axis=0)))),
    )


def compare(first, second):
    de = second["energy_per_CH2_eV"] - first["energy_per_CH2_eV"]
    df = np.asarray(second["forces_eV_per_A"]) - np.asarray(first["forces_eV_per_A"])
    force_max = float(np.linalg.norm(df, axis=1).max())
    return dict(energy_change_meV_per_CH2=de * 1000,
                max_force_vector_change_eV_per_A=force_max,
                passed=bool(abs(de) < THRESHOLDS["energy_eV_per_CH2"] and
                            force_max < THRESHOLDS["force_vector_eV_per_A"]))


def singlepoint(inp, out, repeat, accuracy):
    atoms = read(inp)
    if not all(atoms.pbc):
        raise ValueError("Three periodic dimensions required")
    supercell = atoms.repeat(repeat)
    settings = dict(SETTINGS, accuracy=accuracy)
    supercell.calc = TBLite(**settings)
    start = time.perf_counter()
    energy = supercell.get_potential_energy()
    result = fold_result(atoms, repeat, supercell)
    if not all(np.isfinite(result[k]).all() for k in
               ("energy_per_CH2_eV", "forces_eV_per_A", "stress_eV_per_A3_voigt")):
        raise ValueError("Non-finite result")
    result.update(status="completed", settings=settings, input_sha256=digest(inp), script_sha256=digest(__file__),
                  supercell_energy_eV=float(energy), wall_seconds=time.perf_counter() - start)
    dump(out, result)
    print("RESULT", out, flush=True)


def launch_case(folder, atoms, repeat, accuracy=0.01):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=False)
    write(folder / "input.extxyz", atoms)
    start = time.perf_counter()
    with (folder / "native_output.log").open("w", encoding="utf-8") as log:
        completed = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "singlepoint",
             "--input", str(folder / "input.extxyz"), "--output", str(folder / "result.json"),
             "--repeat", *map(str, repeat), "--accuracy", str(accuracy)],
            stdout=log, stderr=subprocess.STDOUT, timeout=1800, check=False,
            env={**os.environ, "OMP_NUM_THREADS": "4", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"})
    if completed.returncode:
        raise RuntimeError(f"Calculation failed ({completed.returncode}): {folder}")
    result = json.loads((folder / "result.json").read_text())
    result["total_process_seconds"] = time.perf_counter() - start
    result["folder"] = str(folder)
    print(json.dumps({k: result[k] for k in
                      ["repeats", "energy_per_CH2_eV", "max_force_norm_eV_per_A",
                       "replica_force_max_deviation_eV_per_A", "wall_seconds"]}), flush=True)
    return result


def init_run():
    folder = ROOT / "runs" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    folder.mkdir(parents=True, exist_ok=False)
    model = PILOT / "model" / "pe_initial.extxyz"
    dump(folder / "protocol.json", dict(
        status="in_progress", initial_model=str(model), initial_model_sha256=digest(model),
        initial_cell_A=read(model).cell.array.tolist(),
        settings=SETTINGS, thresholds=THRESHOLDS,
        packages={p: importlib.metadata.version(p) for p in ["ase", "tblite", "numpy", "scipy"]},
        python=sys.executable, script_sha256_at_start=digest(__file__),
        geometry_scope="Fixed experimental-prior orthorhombic cell; initial constructed two-chain C4H8 period.",
        sampling_scope="Gamma-point calculations of exact replicas; no explicit k-point integration API.",
        optimization_space="12 primitive atoms, all coordinates free; translational cell period retained. Fixed lattice.",
        selection_rule="Compare c refinements, transverse refinements and combined refinements using energy AND force criteria; repeat checks after relaxation.",
        precision_rule="SCF accuracy=0.01; cross-check accuracy=0.001 before accepting results.",
        limitations=["Single starting packing", "No lattice relaxation", "No global minimum proof",
                     "No DFT calibration or experimental structure validation yet"],
    ))
    dump(ROOT / "active_run.json", {"run_directory": str(folder)})
    print(folder, flush=True)


def sample(args):
    folder = Path(args.run)
    atoms = read(args.input)
    suite = folder / args.name
    suite.mkdir(exist_ok=False)
    results = []
    for repeat_string in args.grids:
        repeat = tuple(map(int, repeat_string.split("x")))
        if len(repeat) != 3 or min(repeat) < 1:
            raise ValueError("Expected positive nx x ny x nz")
        result = launch_case(suite / repeat_string, atoms, repeat, args.accuracy)
        if results:
            result["change_from_previous"] = compare(results[-1], result)
        results.append(result)
        dump(suite / "summary.json", dict(input=str(args.input), input_sha256=digest(args.input),
                                        results=results, thresholds=THRESHOLDS))
    with (suite / "summary.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["repeat", "atoms", "energy_eV_per_CH2", "max_force_eV_A",
                         "replica_deviation_eV_A", "seconds"])
        for r in results:
            writer.writerow(["x".join(map(str, r["repeats"])), r["natoms"],
                             r["energy_per_CH2_eV"], r["max_force_norm_eV_per_A"],
                             r["replica_force_max_deviation_eV_per_A"], r["wall_seconds"]])


class ReplicatedGFN2(Calculator):
    """Derivative-consistent calculator for fixed-cell, primitive-period search."""
    implemented_properties = ["energy", "free_energy", "forces", "stress", "charges"]

    def __init__(self, repeat, accuracy=0.01, **kwargs):
        super().__init__(**kwargs)
        self.repeat = tuple(repeat)
        self.backend = TBLite(**dict(SETTINGS, accuracy=accuracy))
        self.last_folded = None

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        expanded = self.atoms.repeat(self.repeat)
        expanded.calc = self.backend
        result = fold_result(self.atoms, self.repeat, expanded)
        if result["replica_force_max_deviation_eV_per_A"] > THRESHOLDS["replica_force_eV_per_A"]:
            raise RuntimeError("Equivalent replicas have inequivalent forces; averaging is not accepted")
        self.last_folded = result
        self.results = dict(
            energy=result["energy_per_primitive_cell_eV"],
            free_energy=result["energy_per_primitive_cell_eV"],
            forces=np.asarray(result["forces_eV_per_A"]),
            stress=np.asarray(result["stress_eV_per_A3_voigt"]),
            charges=np.asarray(result["mean_charges_e"]),
        )


def geometry_audit(atoms):
    """Audit chemistry after relaxation without imposing the initial bond lengths."""
    ii, jj, shifts, lengths = neighbor_list("ijSd", atoms, 4.0)
    symbols = atoms.get_chemical_symbols()
    chain = atoms.arrays["chain_id"]
    graph = {i: [] for i in range(len(atoms))}
    edges = {}
    inter = []
    errors = []
    for i, j, shift, distance in zip(ii, jj, shifts, lengths):
        i, j = int(i), int(j)
        pair = "".join(sorted([symbols[i], symbols[j]]))
        if chain[i] != chain[j]:
            inter.append(float(distance))
        cutoff = {"CC": 1.80, "CH": 1.30}.get(pair, 0)
        if distance < cutoff:
            graph[i].append(j)
            if chain[i] != chain[j]:
                errors.append(f"Interchain covalent-distance contact: {i}, {j}")
            key = min((i, j, *map(int, shift)), (j, i, *map(int, -shift)))
            edges[key] = dict(i=key[0], j=key[1], shift=list(key[2:]), kind=pair,
                              distance_A=float(distance))
    for i, symbol in enumerate(symbols):
        counts = {s: sum(symbols[j] == s for j in graph[i]) for s in ["C", "H"]}
        expected = {"C": 2, "H": 2} if symbol == "C" else {"C": 1, "H": 0}
        if counts != expected:
            errors.append(f"Incorrect coordination for atom {i}: {counts}")
    cc = [v for v in edges.values() if v["kind"] == "CC"]
    ch = [v for v in edges.values() if v["kind"] == "CH"]
    if len(cc) != 4 or len(ch) != 8:
        errors.append("Expected 4 periodic C-C and 8 C-H edges")
    for n in [0, 1]:
        bonds = [v for v in cc if chain[v["i"]] == n and chain[v["j"]] == n]
        if len(bonds) != 2 or not np.array_equal(
                np.abs(np.asarray(bonds[0]["shift"]) - bonds[1]["shift"]), [0, 0, 1]):
            errors.append(f"Chain {n} lost its winding along c")
    for bond in cc + ch:
        lower, upper = (1.35, 1.75) if bond["kind"] == "CC" else (0.95, 1.22)
        if not lower < bond["distance_A"] < upper:
            errors.append(f"Unreasonable bond length: {bond}")
    if len(lengths) == 0 or not inter:
        errors.append("Expected periodic neighbours and interchain contacts")
    elif min(lengths) < 0.75 or min(inter) < 1.2:
        errors.append("Atomic close contact")
    def bounds(bonds):
        values = [v["distance_A"] for v in bonds]
        return [min(values), max(values)] if values else None
    return dict(passed=not errors, errors=errors, bonds=list(edges.values()),
                C_C_range_A=bounds(cc), C_H_range_A=bounds(ch),
                minimum_interchain_distance_A=min(inter) if inter else None,
                checks="Coordination, periodic c winding, bond-length ranges and close contacts")


def optimize(args):
    folder = Path(args.run) / args.name
    folder.mkdir(exist_ok=False)
    atoms = read(args.input)
    if (len(atoms) != 12 or atoms.get_chemical_formula() != "C4H8"
            or not all(atoms.pbc) or "chain_id" not in atoms.arrays
            or set(atoms.arrays["chain_id"]) != {0, 1}):
        raise ValueError("This optimizer is configured for the 12-atom, two-chain periodic PE model")
    cell_initial = atoms.cell.array.copy()
    atoms.calc = ReplicatedGFN2(args.repeat, args.accuracy)
    dump(folder / "configuration.json", dict(settings=dict(SETTINGS, accuracy=args.accuracy),
          repeat=args.repeat, input=str(args.input), input_sha256=digest(args.input),
          script_sha256=digest(__file__), optimizer="ASE BFGS", maxstep_A=0.06,
          max_steps=args.steps, fmax_eV_A=THRESHOLDS["relaxation_fmax_eV_per_A"]))
    history = []
    start = time.perf_counter()
    opt = BFGS(atoms, logfile=str(folder / "optimizer.log"),
               trajectory=str(folder / "optimization.traj"), maxstep=0.06)

    def record():
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        row = dict(step=opt.nsteps, energy_eV_per_cell=float(energy),
                   energy_eV_per_CH2=float(energy / 4),
                   max_force_eV_A=float(np.linalg.norm(forces, axis=1).max()),
                   replica_force_deviation_eV_A=atoms.calc.last_folded["replica_force_max_deviation_eV_per_A"],
                   elapsed_seconds=time.perf_counter() - start)
        history.append(row)
        dump(folder / "history.json", history)
        write(folder / "latest.extxyz", atoms)
        print("OPTIMIZATION_STEP", json.dumps(row), flush=True)

    opt.attach(record, interval=1)
    converged = bool(opt.run(fmax=THRESHOLDS["relaxation_fmax_eV_per_A"], steps=args.steps))
    if not np.array_equal(cell_initial, atoms.cell.array):
        raise RuntimeError("Cell changed in fixed-cell optimization")
    final = atoms.copy()
    final.wrap()
    final.info.update(model_role="GFN2_fixed_cell_relaxed_candidate",
                      sampling="x".join(map(str, args.repeat)), local_relaxation_converged=converged)
    write(folder / "pe_relaxed.extxyz", final)
    write(folder / "pe_relaxed.cif", final)
    audit = geometry_audit(final)
    dump(folder / "geometry_validation.json", audit)
    dump(folder / "final_result.json", dict(
        converged=converged, steps=opt.nsteps, cell_unchanged=True, cell_A=atoms.cell.array.tolist(),
        initial=history[0], final=history[-1], geometry_validation_passed=audit["passed"],
        folded_result=atoms.calc.last_folded, wall_seconds=time.perf_counter() - start,
        scope="Local stationary candidate in fixed cell and retained primitive period; final sampling recheck required"))
    with (folder / "atomic_forces.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["atom_1based", "element", "x_A", "y_A", "z_A", "Fx_eV_A", "Fy_eV_A", "Fz_eV_A", "Fnorm_eV_A"])
        for i, (atom, force) in enumerate(zip(final, atoms.get_forces()), 1):
            writer.writerow([i, atom.symbol, *atom.position, *force, np.linalg.norm(force)])
    if not converged or not audit["passed"]:
        raise RuntimeError("Optimization or final chemical audit did not pass")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("init")
    sp = commands.add_parser("singlepoint")
    sp.add_argument("--input", type=Path, required=True)
    sp.add_argument("--output", type=Path, required=True)
    sp.add_argument("--repeat", nargs=3, type=int, required=True)
    sp.add_argument("--accuracy", type=float, default=0.01)
    sa = commands.add_parser("sample")
    sa.add_argument("--run", required=True)
    sa.add_argument("--input", type=Path, required=True)
    sa.add_argument("--name", required=True)
    sa.add_argument("--grids", nargs="+", required=True)
    sa.add_argument("--accuracy", type=float, default=0.01)
    op = commands.add_parser("optimize")
    op.add_argument("--run", required=True)
    op.add_argument("--input", type=Path, required=True)
    op.add_argument("--name", required=True)
    op.add_argument("--repeat", nargs=3, type=int, required=True)
    op.add_argument("--accuracy", type=float, default=0.01)
    op.add_argument("--steps", type=int, default=300)
    args = parser.parse_args()
    if args.command == "init":
        init_run()
    elif args.command == "singlepoint":
        singlepoint(args.input, args.output, args.repeat, args.accuracy)
    elif args.command == "sample":
        sample(args)
    elif args.command == "optimize":
        optimize(args)


if __name__ == "__main__":
    main()
