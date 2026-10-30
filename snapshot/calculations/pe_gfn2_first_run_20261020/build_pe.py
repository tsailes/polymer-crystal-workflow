"""Build an independent, fully hydrogenated periodic PE starting hypothesis.

No crystallographic coordinates are read.  The declared orthorhombic cell is
an experimental-cell prior from this project's Tashiro 1997 record; atomic coordinates follow
from C-C / C-H lengths, all-trans chain geometry and arbitrary packing choices.
This is an input for a first energy/force calculation, not a recovered crystal.

Usage:
    python build_pe.py --output-dir model

Dependencies: numpy and ASE only (plus Python's standard library).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import ase
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.neighborlist import neighbor_list


def angle_deg(first, second):
    first, second = np.asarray(first, float), np.asarray(second, float)
    cosine = np.dot(first, second) / (np.linalg.norm(first) * np.linalg.norm(second))
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def build_model(cell_lengths, cc, ch, rotations, phases, centers):
    a, b, c = (float(v) for v in cell_lengths)
    if min(a, b, c, cc, ch) <= 0:
        raise ValueError("Cell lengths and bond lengths must be positive")
    axial_step = c / 2.0
    if axial_step >= cc:
        raise ValueError("Two-carbon all-trans chain requires c/2 < C-C length")
    transverse_step = math.sqrt(cc * cc - axial_step * axial_step)
    # For carbon neighbours u1,u2, use their unit bisector b and a normal n.
    # h+/- = -q*b +/- sqrt(1-q*q)*n, with q=1/(3*cos(CCC/2)),
    # gives u1.h = u2.h = -1/3 (tetrahedral C-C-H angles).
    # The fixed chain period sets CCC, so H-C-H follows from this constraint;
    # it is not silently forced to an incompatible ideal tetrahedron.
    q = cc / (3.0 * transverse_step)
    if q >= 1:
        raise ValueError("Chain period leaves no real tetrahedral C-C-H construction")
    h_normal = math.sqrt(1.0 - q * q)
    cell = np.diag([a, b, c])
    positions, symbols, chain_ids, roles = [], [], [], []
    for chain, (rotation, phase, center) in enumerate(zip(rotations, phases, centers)):
        theta = math.radians(float(rotation))
        rot = np.array([[math.cos(theta), -math.sin(theta), 0.0],
                        [math.sin(theta), math.cos(theta), 0.0],
                        [0.0, 0.0, 1.0]])
        offset = np.array([center[0] * a, center[1] * b, phase * c])
        carbon_positions = [np.array([-transverse_step / 2, 0.0, 0.0]),
                            np.array([transverse_step / 2, 0.0, axial_step])]
        for carbon, position in enumerate(carbon_positions):
            positions.append(rot @ position + offset)
            symbols.append("C")
            chain_ids.append(chain)
            roles.append(f"chain{chain + 1}_C{carbon + 1}")
            bisector_sign = 1.0 if carbon == 0 else -1.0
            for hydrogen, sign in enumerate((1.0, -1.0)):
                h = np.array([-bisector_sign * q, sign * h_normal, 0.0])
                positions.append(rot @ (position + ch * h) + offset)
                symbols.append("H")
                chain_ids.append(chain)
                roles.append(f"chain{chain + 1}_C{carbon + 1}_H{hydrogen + 1}")
    atoms = Atoms(symbols, positions=positions, cell=cell, pbc=True)
    atoms.wrap()
    atoms.set_array("chain_id", np.array(chain_ids, dtype=int))
    atoms.info.update(polymer="PE", model_role="independent_initial_hypothesis",
                      chemical_repeat="CH2", chemical_repeat_count=4,
                      chain_count=2, electronic_charge=0, unpaired_electrons=0,
                      solved_coordinates_used=False)
    construction = {
        "carbon_axial_step_A": axial_step,
        "carbon_transverse_step_A": transverse_step,
        "CCC_angle_deg": angle_deg([transverse_step, 0, axial_step],
                                    [transverse_step, 0, -axial_step]),
        "CCH_angle_deg": math.degrees(math.acos(-1.0 / 3.0)),
        "HCH_angle_deg": angle_deg([-q, h_normal, 0], [-q, -h_normal, 0]),
        "hydrogen_rule": "All four C-C-H angles at each carbon equal acos(-1/3); "
                         "H-C-H adjusts to the C-C-C angle imposed by c and C-C length.",
        "atom_roles": roles,
    }
    return atoms, construction


def verify_geometry(atoms, cc, ch):
    """Inspect full periodic neighbour images, retaining both C-C links per C."""
    indices_i, indices_j, shifts, displacements, distances = neighbor_list(
        "ijSDd", atoms, 4.0, self_interaction=False)
    symbols = atoms.get_chemical_symbols()
    chain_ids = atoms.arrays["chain_id"]
    neighbours = {i: [] for i in range(len(atoms))}
    undirected = {}
    nonbonded, interchain = [], []
    for i, j, shift, vector, distance in zip(
            indices_i, indices_j, shifts, displacements, distances):
        i, j, distance = int(i), int(j), float(distance)
        pair = "".join(sorted([symbols[i], symbols[j]]))
        threshold = {"CC": 1.80, "CH": 1.30}.get(pair)
        bonded = threshold is not None and distance < threshold
        if chain_ids[i] != chain_ids[j]:
            interchain.append(distance)
        if bonded:
            record = {"index": j, "element": symbols[j], "shift": shift.tolist(),
                      "distance_A": distance, "vector_A": vector.tolist()}
            neighbours[i].append(record)
            forward = (i, j, *map(int, shift))
            backward = (j, i, *map(int, -shift))
            key = min(forward, backward)
            undirected[key] = {"i": key[0], "j": key[1], "shift": list(key[2:]),
                               "distance_A": distance, "kind": pair}
        else:
            nonbonded.append(distance)
    errors, atom_records, angles = [], [], []
    for i, symbol in enumerate(symbols):
        n = neighbours[i]
        counts = {element: sum(v["element"] == element for v in n) for element in ("C", "H")}
        expected = {"C": 2, "H": 2} if symbol == "C" else {"C": 1, "H": 0}
        if counts != expected:
            errors.append(f"Atom {i} {symbol}: neighbour counts {counts}; expected {expected}")
        if any(chain_ids[i] != chain_ids[v["index"]] for v in n):
            errors.append(f"Atom {i} has an unintended interchain covalent-distance neighbour")
        atom_records.append({"index": i, "element": symbol, "chain_id": int(chain_ids[i]),
                             "neighbour_counts": counts, "neighbours": n})
        if symbol == "C" and counts == expected:
            c_vectors = [v["vector_A"] for v in n if v["element"] == "C"]
            h_vectors = [v["vector_A"] for v in n if v["element"] == "H"]
            angles.append({"carbon_index": i, "CCC_deg": angle_deg(*c_vectors),
                           "HCH_deg": angle_deg(*h_vectors),
                           "CCH_deg": [angle_deg(v, w) for v in c_vectors for w in h_vectors]})
    bonds = list(undirected.values())
    if len(bonds) != 12:
        errors.append(f"Expected 12 periodic bond edges, found {len(bonds)}")
    cc_bonds = [v for v in bonds if v["kind"] == "CC"]
    ch_bonds = [v for v in bonds if v["kind"] == "CH"]
    if len(cc_bonds) != 4 or len(ch_bonds) != 8:
        errors.append("Expected four C-C and eight C-H periodic edges")
    if any(abs(v["distance_A"] - (cc if v["kind"] == "CC" else ch)) > 1e-8 for v in bonds):
        errors.append("Constructed bond length differs from requested length")
    crossing = [v for v in cc_bonds if v["shift"] != [0, 0, 0]]
    for chain in (0, 1):
        c_ids = [i for i in range(len(atoms)) if symbols[i] == "C" and chain_ids[i] == chain]
        two = [v for v in cc_bonds if v["i"] in c_ids and v["j"] in c_ids]
        if len(two) != 2:
            errors.append(f"Chain {chain} lacks its two periodic C-C edges")
        elif not np.array_equal(np.abs(np.asarray(two[0]["shift"]) - two[1]["shift"]), [0, 0, 1]):
            errors.append(f"Chain {chain} does not wind exactly once through the c boundary")
    minimum_distance = float(np.min(distances))
    minimum_nonbonded = float(min(nonbonded))
    minimum_interchain = float(min(interchain))
    if minimum_distance < 0.75 or minimum_nonbonded < 1.0 or minimum_interchain < 1.2:
        errors.append("Initial geometry fails the declared close-contact screen")
    volume = float(atoms.get_volume())
    return {
        "passed": not errors, "errors": errors,
        "formula": atoms.get_chemical_formula(), "atom_count": len(atoms),
        "composition": {element: symbols.count(element) for element in ("C", "H")},
        "pbc": atoms.pbc.tolist(), "cell_A": atoms.cell.array.tolist(), "volume_A3": volume,
        "density_g_cm3": float(atoms.get_masses().sum() * 1.66053906660 / volume),
        "chemical_repeat_unit": "CH2", "chemical_repeat_count": 4,
        "two_carbon_chain_periods_per_cell": 2,
        "periodic_CC_edge_count": len(cc_bonds), "CH_edge_count": len(ch_bonds),
        "CC_edges_crossing_cell_boundary": crossing,
        "minimum_pair_distance_A": minimum_distance,
        "minimum_nonbonded_distance_A": minimum_nonbonded,
        "minimum_interchain_distance_A": minimum_interchain,
        "screen_thresholds_A": {"CC_bond_cutoff": 1.80, "CH_bond_cutoff": 1.30,
                                "minimum_pair": 0.75, "minimum_nonbonded": 1.0,
                                "minimum_interchain": 1.2},
        "periodic_bonds": bonds, "atoms": atom_records, "angles": angles,
        "interpretation": "Connectivity and close-contact validation only; "
                          "no physical equilibrium, stable crystal or final phase is asserted.",
    }


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "model")
    parser.add_argument("--cell", type=float, nargs=3, default=[7.432, 4.929, 2.545], metavar=("A", "B", "C"))
    parser.add_argument("--cc", type=float, default=1.54)
    parser.add_argument("--ch", type=float, default=1.09)
    parser.add_argument("--rotations", type=float, nargs=2, default=[17.0, 43.0])
    parser.add_argument("--phases", type=float, nargs=2, default=[0.0, 0.17])
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    files = {name: output / name for name in ("pe_initial.cif", "pe_initial.extxyz", "geometry_validation.json", "provenance.json")}
    existing = [str(p) for p in files.values() if p.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite existing model artifacts: " + ", ".join(existing))
    centers = [[0.25, 0.25], [0.75, 0.75]]
    atoms, construction = build_model(args.cell, args.cc, args.ch, args.rotations, args.phases, centers)
    validation = verify_geometry(atoms, args.cc, args.ch)
    if not validation["passed"]:
        raise RuntimeError("Initial geometry rejected: " + "; ".join(validation["errors"]))
    write(files["pe_initial.cif"], atoms, format="cif")
    write(files["pe_initial.extxyz"], atoms, format="extxyz")
    roundtrip = {}
    for filename in ("pe_initial.cif", "pe_initial.extxyz"):
        reread = read(files[filename])
        if reread.get_chemical_symbols() != atoms.get_chemical_symbols():
            raise RuntimeError(f"Element order changed in {filename}")
        if not np.allclose(reread.cell.array, atoms.cell.array, atol=1e-8) or not np.all(reread.pbc):
            raise RuntimeError(f"Cell or PBC changed in {filename}")
        delta = reread.get_scaled_positions() - atoms.get_scaled_positions()
        delta -= np.rint(delta)
        maximum = float(np.linalg.norm(delta @ atoms.cell.array, axis=1).max())
        if maximum > 1e-7:
            raise RuntimeError(f"Coordinates changed on {filename} round trip")
        roundtrip[filename] = {"passed": True, "maximum_periodic_displacement_A": maximum}
    validation["file_roundtrip"] = roundtrip
    save_json(files["geometry_validation.json"], validation)
    provenance = {
        "schema_version": "pe_analytic_start/1.0",
        "purpose": "First periodic GFN2-xTB energy/force input; independent initial hypothesis",
        "builder": str(Path(__file__).resolve()),
        "builder_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "versions": {"numpy": np.__version__, "ase": ase.__version__},
        "solved_atomic_coordinates_read_or_copied": False,
        "input_structure_files": [],
        "cell_prior": {"lengths_A": args.cell, "angles_deg": [90.0, 90.0, 90.0],
                       "source": "Experimental PE cell prior from this project's Tashiro 1997 literature record "
                                 "and operation note: 7.432, 4.929, 2.545 A. These values were not supplied "
                                 "directly by the user. Only cell metadata was checked; no solved atomic "
                                 "coordinates are read or imported by this builder.",
                       "source_record": {
                           "path": "E:/PE_DFT/examples/PE_Tashiro_1997/sources/Tashiro_1997_literature_record.json",
                           "sha256": "ca1efa9c62a5cfc01f8b927ef4d532c9f6371e0edd854a1f623e19f3d73a1281",
                           "evidence_span_id": "evidence_pe_cell_symmetry"},
                       "source_note": {
                           "path": "E:/PE_DFT/examples/PE_Tashiro_1997/操作说明.txt",
                           "sha256": "2e3b17cc4f3f4daab0622ef99a4b4d3d464fe26c987de0390e3d2bd4e2d61bd6",
                           "lines": [10, 11]},
                       "role": "Declared experimental-cell prior; cell is not predicted in this first run."},
        "chemical_assumptions": {"polymer": "PE", "conformation": "all-trans infinite zigzag",
                                 "CC_length_A": args.cc, "CH_length_A": args.ch,
                                 "chain_axis": "c", "chain_count": 2,
                                 "formula_per_cell": "C4H8", "chemical_repeat": "CH2",
                                 "chemical_repeat_count": 4, "pbc": [True, True, True],
                                 "charge": 0, "unpaired_electrons": 0},
        "arbitrary_packing_choices": {"fractional_transverse_centers": centers,
                                      "chain_rotations_deg": args.rotations,
                                      "chain_axial_phases_fractional": args.phases,
                                      "selection_reason": "Distinct arbitrary rotations; 17/43 degrees replaces the "
                                                          "initial 17/103-degree geometry trial to avoid its short "
                                                          "1.77-A interchain contact. Selection uses only generated "
                                                          "geometries, with no energies, diffraction or solved coordinates.",
                                      "geometric_trials": {"first_chain_rotation_deg": 17.0,
                                                           "second_chain_rotations_deg": [43, 57, 73, 87, 103, 117, 137, 157],
                                                           "selection": "First tested geometry with minimum interchain distance > 2.1 A."}},
        "construction": construction,
        "scope": ["All H atoms are constructed analytically.",
                  "Cell metric is orthorhombic; no final space group or recovered packing is claimed.",
                  "Two chains and the two-carbon translation are explicit model assumptions.",
                  "A single initial model is an engine check, not a candidate search or method validation."],
        "artifacts_sha256": {name: hashlib.sha256(path.read_bytes()).hexdigest()
                             for name, path in files.items() if name != "provenance.json"},
    }
    save_json(files["provenance.json"], provenance)
    print(json.dumps({"output_dir": str(output), "passed": validation["passed"],
                      "formula": validation["formula"], "atom_count": len(atoms),
                      "chemical_repeat_count": 4,
                      "minimum_interchain_distance_A": validation["minimum_interchain_distance_A"],
                      "CCC_angle_deg": construction["CCC_angle_deg"],
                      "HCH_angle_deg": construction["HCH_angle_deg"]}, indent=2))


if __name__ == "__main__":
    main()
