"""Export PE as three fourfold Pnam sites with verified explicit symmetry.

This performs geometry/file-format work only, never energy calculations.
For trajectory.traj, a successful sibling relaxation result is required and
the last frame must equal its unrounded final_evaluation. Other structure
inputs (including the experimental reference POSCAR) are geometry-only inputs.
Existing source files, full-cell CIFs, and output files are never overwritten.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import importlib.metadata
import json
from pathlib import Path

import numpy as np
import spglib
from ase import Atoms
from ase.geometry import find_mic
from ase.io import read
from scipy.optimize import linear_sum_assignment


HALL_NUMBER = 297
SPACE_GROUP_NUMBER = 62


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def serialize(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(type(value).__name__)


def write_json_exclusive(path, value):
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False,
                  allow_nan=False, default=serialize)
        handle.write("\n")


def match_periodic(reference, other, tolerance):
    """Element-preserving one-to-one minimum-distance periodic assignment."""
    if (len(reference) != len(other)
            or sorted(reference.numbers) != sorted(other.numbers)):
        raise ValueError("Atom counts/elements differ")
    if not reference.pbc.all() or not other.pbc.all():
        raise ValueError("Three-dimensional periodicity is required")
    cell_error = float(np.max(np.abs(reference.cell.array - other.cell.array)))
    if cell_error > tolerance:
        raise ValueError(f"Cell changed by {cell_error:.9g} Angstrom")
    errors = []
    assignments = []
    for number in sorted(set(reference.numbers)):
        left = np.flatnonzero(reference.numbers == number)
        right = np.flatnonzero(other.numbers == number)
        delta = (reference.positions[left, None, :]
                 - other.positions[None, right, :])
        _, distances = find_mic(delta.reshape(-1, 3), reference.cell, pbc=True)
        distances = distances.reshape(len(left), len(right))
        rows, columns = linear_sum_assignment(distances)
        selected = distances[rows, columns]
        errors.extend(selected.tolist())
        assignments.extend(
            {"reference_atom": int(left[i]), "other_atom": int(right[j]),
             "atomic_number": int(number), "distance_A": float(distances[i, j])}
            for i, j in zip(rows, columns))
    maximum = max(errors, default=0.0)
    if maximum > tolerance:
        raise ValueError(f"Periodic atom mismatch {maximum:.9g} Angstrom")
    return {"maximum_periodic_distance_A": maximum,
            "rms_periodic_distance_A": float(np.sqrt(np.mean(np.square(errors)))),
            "maximum_cell_component_error_A": cell_error,
            "matching_tolerance_A": tolerance,
            "element_preserving_assignments": assignments}


def symmetry_summary(atoms, symprec, hall_number):
    dataset = spglib.get_symmetry_dataset(
        (atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers),
        symprec=symprec, hall_number=hall_number)
    if dataset is None:
        raise ValueError(f"Cannot identify requested Hall setting {hall_number}")
    summary = {
        "number": int(dataset.number), "hall_number": int(dataset.hall_number),
        "international": dataset.international, "hall": dataset.hall,
        "choice": dataset.choice, "operation_count": len(dataset.rotations),
        "transformation_matrix": dataset.transformation_matrix,
        "origin_shift": dataset.origin_shift,
        "equivalent_atoms": dataset.equivalent_atoms,
        "wyckoff_letters": list(dataset.wyckoffs),
    }
    return dataset, summary


def unique_orbit(fractional, rotations, translations, cell, tolerance):
    positions = []
    for rotation, translation in zip(rotations, translations):
        value = (rotation @ fractional + translation) % 1.0
        if positions:
            delta = (np.asarray(positions) - value) @ cell
            _, distances = find_mic(delta, cell, pbc=True)
            if np.min(distances) <= tolerance:
                continue
        positions.append(value)
    return np.array(positions)


def operation_xyz(rotation, translation):
    expressions = []
    for row, shift in zip(rotation, translation):
        terms = []
        for coefficient, name in zip(row, ("x", "y", "z")):
            if coefficient == 0:
                continue
            coefficient = int(coefficient)
            sign = "+" if coefficient > 0 else "-"
            magnitude = "" if abs(coefficient) == 1 else str(abs(coefficient)) + "*"
            terms.append(sign + magnitude + name)
        shift = float(shift % 1.0)
        if abs(shift - 1.0) < 1e-12:
            shift = 0.0
        fraction = Fraction(shift).limit_denominator(24)
        if abs(float(fraction) - shift) > 1e-12:
            raise ValueError("Symmetry translation cannot be represented exactly")
        if fraction:
            terms.append("+" + str(fraction))
        expression = "".join(terms).lstrip("+") or "0"
        expressions.append(expression)
    return ",".join(expressions)


def build_cif(atoms, symprec, tolerance):
    if (len(atoms) != 12 or atoms.get_chemical_formula() != "C4H8"
            or not atoms.pbc.all()):
        raise ValueError("Expected a fully periodic 12-atom C4H8 PE cell")
    if (not np.isfinite(atoms.positions).all()
            or not np.isfinite(atoms.cell.array).all()
            or np.linalg.det(atoms.cell.array) <= 0):
        raise ValueError("Nonfinite geometry or invalid cell orientation/volume")
    # A CIF specifies metric and fractional axes, not an arbitrary Cartesian
    # orientation. Reject such an input rather than silently rotate it.
    if (not np.allclose(atoms.cell.array, np.diag(atoms.cell.lengths()),
                        atol=tolerance, rtol=0)
            or not np.allclose(atoms.cell.angles(), 90.0, atol=1e-6, rtol=0)):
        raise ValueError("Expected orthorhombic a,b,c along the original positive Cartesian axes")
    dataset, requested = symmetry_summary(atoms, symprec, HALL_NUMBER)
    _, standard = symmetry_summary(atoms, symprec, 0)
    if (requested["number"] != SPACE_GROUP_NUMBER
            or requested["operation_count"] != 8
            or standard["number"] != SPACE_GROUP_NUMBER
            or standard["operation_count"] != 8
            or not np.allclose(dataset.transformation_matrix, np.eye(3),
                               atol=1e-8, rtol=0)):
        raise ValueError("Input is not exactly #62 in the original Pnam/Hall297 axes")

    operations = spglib.get_symmetry_from_database(HALL_NUMBER)
    rotations, translations = operations["rotations"], operations["translations"]
    if len(rotations) != 8:
        raise ValueError("Unexpected Hall297 operation count")
    fractional = atoms.get_scaled_positions(wrap=True)
    closure = []
    for index, (rotation, translation) in enumerate(zip(rotations, translations), 1):
        image = Atoms(numbers=atoms.numbers, cell=atoms.cell, pbc=True,
                      scaled_positions=fractional @ rotation.T + translation)
        closure.append({"operation_id": index,
                        **match_periodic(atoms, image, tolerance)})

    groups = {}
    for index, equivalent in enumerate(dataset.equivalent_atoms):
        groups.setdefault(int(equivalent), []).append(index)
    if len(groups) != 3 or sorted(map(len, groups.values())) != [4, 4, 4]:
        raise ValueError("Expected exactly three fourfold crystallographic orbits")
    groups = sorted(groups.values(), key=lambda indices: (0 if atoms[indices[0]].symbol == "C" else 1,
                                                           indices[0]))
    if [atoms[indices[0]].symbol for indices in groups] != ["C", "H", "H"]:
        raise ValueError("Expected one carbon and two hydrogen orbits")

    sites, expanded_positions, expanded_numbers = [], [], []
    for label, indices in zip(("C1", "H1", "H2"), groups):
        representative = min(indices)
        if len(set(atoms.numbers[indices])) != 1:
            raise ValueError("An equivalent-atom orbit mixes elements")
        orbit = unique_orbit(fractional[representative], rotations, translations,
                             atoms.cell.array, tolerance)
        if len(orbit) != 4:
            raise ValueError(f"{label} is not a fourfold site under Hall297 operations")
        orbit_atoms = Atoms(numbers=[int(atoms.numbers[representative])] * len(orbit),
                            cell=atoms.cell, pbc=True, scaled_positions=orbit)
        validation = match_periodic(atoms[indices], orbit_atoms, tolerance)
        sites.append({"label": label, "element": atoms[representative].symbol,
                      "representative_input_index": representative,
                      "equivalent_input_indices": indices, "multiplicity": 4,
                      "wyckoff_letter": str(dataset.wyckoffs[representative]),
                      "fractional_coordinates": fractional[representative],
                      "orbit_validation": validation})
        expanded_positions.extend(orbit)
        expanded_numbers.extend(orbit_atoms.numbers)
    expanded = Atoms(numbers=expanded_numbers, cell=atoms.cell, pbc=True,
                     scaled_positions=expanded_positions)
    reconstructed = match_periodic(atoms, expanded, tolerance)
    symbols = [operation_xyz(r, t) for r, t in zip(rotations, translations)]
    group_type = spglib.get_spacegroup_type(HALL_NUMBER)
    lengths, angles = atoms.cell.lengths(), atoms.cell.angles()
    lines = [
        "data_PE_Pnam", "# Geometry export only; no energy calculation or geometry refinement.",
        "# Hall number 297, original Pnam axes; c is the chain axis.",
        "_audit_creation_method 'PE Pnam symmetry-verified export'",
        "_chemical_formula_sum 'C4 H8'",
        "_space_group_IT_number 62",
        "_space_group_name_H-M_alt 'P n a m'",
        f"_space_group_name_Hall '{group_type.hall_symbol}'",
        f"_cell_length_a {lengths[0]:.16g}",
        f"_cell_length_b {lengths[1]:.16g}",
        f"_cell_length_c {lengths[2]:.16g}",
        f"_cell_angle_alpha {angles[0]:.16g}",
        f"_cell_angle_beta {angles[1]:.16g}",
        f"_cell_angle_gamma {angles[2]:.16g}",
        f"_cell_volume {atoms.get_volume():.16g}", "",
        "loop_", "_space_group_symop_id", "_space_group_symop_operation_xyz",
    ]
    lines.extend(f"{index} '{operation}'" for index, operation in enumerate(symbols, 1))
    lines += ["", "loop_", "_atom_site_label", "_atom_site_type_symbol",
              "_atom_site_symmetry_multiplicity", "_atom_site_Wyckoff_symbol",
              "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z",
              "_atom_site_occupancy"]
    for site in sites:
        coordinates = " ".join(f"{value:.16g}" for value in site["fractional_coordinates"])
        lines.append(f"{site['label']} {site['element']} 4 {site['wyckoff_letter']} "
                     f"{coordinates} 1.0")
    lines.append("")
    return "\n".join(lines), {
        "input_symmetry_Hall297": requested, "input_standard_symmetry": standard,
        "operations_xyz": symbols, "rotations": rotations, "translations": translations,
        "operation_count": 8, "independent_site_count": 3, "sites": sites,
        "full_cell_atom_count": 12, "formula": "C4H8",
        "occupancy": "All sites full occupancy; inherited PE model assumption.",
        "full_cell_operation_closure": closure,
        "reconstruction_from_three_sites": reconstructed,
    }


def load_structure(path):
    before_hash = sha(path)
    atoms = read(path, index=-1)
    if sha(path) != before_hash:
        raise ValueError("Source changed during reading; do not export a live trajectory")
    context = {"source_path": str(path), "source_sha256": before_hash,
               "frame_index": -1, "relaxation_status_checked": False}
    if path.name == "trajectory.traj":
        result_path = path.parent / "result.json"
        result = json.loads(result_path.read_text(encoding="utf-8"))
        if result.get("status") != "passed_relaxation" or result.get("exit_code") != 0:
            raise ValueError("Only a successfully completed relaxation trajectory may be exported")
        final = result["final_evaluation"]
        if (not np.array_equal(atoms.cell.array, final["cell_A"])
                or not np.array_equal(atoms.positions, final["positions_A"])):
            raise ValueError("Trajectory last frame differs from unrounded final_evaluation")
        context.update(relaxation_status_checked=True,
                       source_result_path=str(result_path),
                       source_result_sha256=sha(result_path))
    atoms.set_constraint([])
    atoms.calc = None
    return atoms, context


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("structure", type=Path,
                        help="Completed trajectory.traj (last frame), POSCAR, CIF or extxyz")
    parser.add_argument("--output", type=Path,
                        help="New CIF path; default: source folder/final_Pnam.cif")
    parser.add_argument("--symprec", type=float, default=1e-5)
    parser.add_argument("--tolerance", type=float, default=1e-7,
                        help="Maximum geometry reconstruction error in Angstrom")
    args = parser.parse_args()
    if not 0 < args.symprec <= 1e-4 or not 0 < args.tolerance <= 1e-6:
        raise ValueError("Use 0 < symprec <= 1e-4 and 0 < tolerance <= 1e-6 Angstrom")
    source = args.structure.resolve()
    output = (args.output or source.parent / "final_Pnam.cif").resolve()
    audit_path = output.with_suffix(".audit.json")
    if source == output or output.exists() or audit_path.exists():
        raise FileExistsError("Source/output or existing output would be overwritten")
    if not output.parent.is_dir():
        raise FileNotFoundError("Output parent must already exist")
    atoms, provenance = load_structure(source)
    text, geometry = build_cif(atoms, args.symprec, args.tolerance)
    audit = {
        "schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "pending_roundtrip", "quantum_calculations_started": False,
        "source": provenance, "exporter_path": str(Path(__file__).resolve()),
        "exporter_sha256": sha(__file__),
        "software": {name: importlib.metadata.version(name)
                     for name in ("ase", "spglib", "numpy", "scipy")},
        "settings": {"symprec_A": args.symprec, "matching_tolerance_A": args.tolerance,
                     "space_group_number": 62, "Hall_number": 297},
        "geometry": geometry, "output_path": str(output),
        "note": "Three inequivalent sites, not twelve full-cell atoms with symmetry added. "
                "No coordinate fitting, standardization, axis permutation or optimization was performed.",
    }
    with output.open("x", encoding="ascii", newline="\n") as handle:
        handle.write(text)
    try:
        reloaded = read(output, format="cif")
        audit["ASE_roundtrip"] = match_periodic(atoms, reloaded, args.tolerance)
        _, after = symmetry_summary(reloaded, args.symprec, HALL_NUMBER)
        if (after["number"] != 62 or after["operation_count"] != 8
                or not np.allclose(after["transformation_matrix"], np.eye(3),
                                   atol=1e-8, rtol=0)):
            raise ValueError("Roundtrip changed the original Pnam setting")
        audit["reloaded_symmetry_Hall297"] = after
        audit["reloaded_formula"] = reloaded.get_chemical_formula()
        audit["reloaded_atom_count"] = len(reloaded)
        audit["output_sha256"] = sha(output)
        if sha(source) != provenance["source_sha256"]:
            raise ValueError("Source changed during export")
        audit["status"] = "passed_symmetry_and_roundtrip"
    except Exception as exc:
        audit["status"] = "failed_roundtrip"
        audit["error"] = repr(exc)
        audit["output_sha256"] = sha(output)
        write_json_exclusive(audit_path, audit)
        raise
    write_json_exclusive(audit_path, audit)
    print(json.dumps({"status": audit["status"], "CIF": str(output), "audit": str(audit_path),
                      "independent_sites": 3, "symmetry_operations": 8,
                      "restored_atoms": len(reloaded),
                      "max_roundtrip_error_A": audit["ASE_roundtrip"]["maximum_periodic_distance_A"],
                      "CIF_sha256": audit["output_sha256"]}, indent=2))


if __name__ == "__main__":
    main()
