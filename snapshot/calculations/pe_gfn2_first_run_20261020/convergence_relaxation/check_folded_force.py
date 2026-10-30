"""Independent finite-difference check of the normalized supercell objective."""
import argparse
from pathlib import Path
import numpy as np
from ase.io import read, write
from workflow import ReplicatedGFN2, dump, digest

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--repeat", type=int, nargs=3, required=True)
args = parser.parse_args()
args.output.mkdir(exist_ok=False)
atoms = read(args.input)
atoms.calc = ReplicatedGFN2(args.repeat)
energy = atoms.get_potential_energy()
forces = atoms.get_forces()
atom, component = np.unravel_index(np.abs(forces).argmax(), forces.shape)
analytic = float(forces[atom, component])
step = 1e-4
energies = {}
for name, sign in [("plus", 1), ("minus", -1)]:
    displaced = atoms.copy()
    displaced.positions[atom, component] += sign * step
    write(args.output / f"{name}.extxyz", displaced)
    displaced.calc = ReplicatedGFN2(args.repeat)
    energies[name] = float(displaced.get_potential_energy())
finite_difference = -(energies["plus"] - energies["minus"]) / (2 * step)
error = abs(finite_difference - analytic)
result = dict(
    input=str(args.input), input_sha256=digest(args.input), repeat=args.repeat,
    atom_1based=int(atom + 1), cartesian_component="xyz"[component], displacement_A=step,
    energy_eV_per_primitive_cell=float(energy), displaced_energies_eV_per_cell=energies,
    analytic_force_eV_A=analytic, finite_difference_force_eV_A=finite_difference,
    absolute_error_eV_A=error, threshold_eV_A=0.001, passed=bool(error < 0.001),
    explanation="Each perturbed primitive coordinate moves every equivalent replica. Energy is divided by replica count; forces are averaged, with no per-CH2 division.")
dump(args.output / "result.json", result)
print(result, flush=True)
if not result["passed"]:
    raise RuntimeError("Folded force finite difference failed")
