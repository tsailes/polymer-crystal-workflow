#!/usr/bin/env python3
"""Independent geometric checks for generated, full-occupancy PA6 O candidates.

Reads only generated structure JSON, never a paper or a reference coordinate file.
Cell vectors are rows; bond shift means r[j] + shift @ cell - r[i].
Physical admission uses chemistry, periodic connectivity, collision checks and
the declared branch symmetry. Hydrogen bonds are not evaluated by default.
Optional hydrogen-bond diagnostics never affect physical admission; their
distance/angle bounds are not experimental hydrogen-bond evidence.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict, deque
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

VERSION = "pa6-ordered-validator-0.3"
OPS = [
    (np.eye(3, dtype=int), np.zeros(3)),
    (np.diag([-1, -1, 1]), np.array([.5, .5, .5])),
    (-np.eye(3, dtype=int), np.zeros(3)),
    (np.diag([1, 1, -1]), np.array([.5, .5, .5])),
]
NUMBERS = {"H": 1, "C": 6, "N": 7, "O": 8}
# Independent broad operational bounds. No PA6 reference geometry was used.
BOND_LIMITS = {
    "CC": (1.25, 1.85), "CN": (1.15, 1.75), "CO": (1.05, 1.55),
    "CH": (.85, 1.30), "HN": (.80, 1.25),
}
COLLISION_FLOORS = {
    "HH": 1.10, "CH": 1.35, "HN": 1.25, "HO": 1.25,
    "CC": 2.00, "CN": 1.85, "CO": 1.85,
    "NN": 1.75, "NO": 1.75, "OO": 1.75,
}


def minimum_images(delta_fractional, cell):
    """Exact closest image by bounded enumeration, including skew cells.

    First find an upper bound R from componentwise rounding (NOT the answer).
    Any better Cartesian vector v has |f_i| <= R ||inv(cell)[:, i]||.
    Enumerating these bounded integer ranges around the initial image therefore
    includes every possible better image. This also handles unreduced bases.
    """
    delta = np.asarray(delta_fractional, float)
    shape = delta.shape
    flat = delta.reshape(-1, 3)
    cell = np.asarray(cell, float)
    if abs(np.linalg.det(cell)) < 1e-10:
        raise ValueError("Singular or nearly singular cell")
    initial = -np.rint(flat).astype(int)
    wrapped = flat + initial
    best = wrapped @ cell
    best2 = np.einsum("ij,ij->i", best, best)
    shifts = initial.copy()
    radius = np.sqrt(np.max(best2, initial=0.0))
    reach = np.ceil(radius * np.linalg.norm(np.linalg.inv(cell), axis=0) + .5).astype(int)
    count = int(np.prod(2 * reach + 1))
    if count > 200000:
        raise ValueError("Cell basis too ill-conditioned for bounded image enumeration")
    for offset in itertools.product(*(range(-int(k), int(k) + 1) for k in reach)):
        offset = np.array(offset, int)
        trial = (wrapped + offset) @ cell
        d2 = np.einsum("ij,ij->i", trial, trial)
        use = d2 < best2 - 1e-14
        best[use] = trial[use]
        best2[use] = d2[use]
        shifts[use] = initial[use] + offset
    return best.reshape(shape), shifts.reshape(shape)


def canonical_edge(i, j, shift):
    shift = tuple(int(x) for x in shift)
    return (int(i), int(j), *shift) if i < j else (int(j), int(i), *(-x for x in shift))


def structure_arrays(record):
    cell = np.asarray(record["cell"], float)
    symbols = np.asarray(record["symbols"], str)
    if "positions" in record:
        positions = np.asarray(record["positions"], float)
    else:
        positions = np.asarray(record["fractional_positions"], float) @ cell
    if positions.shape != (len(symbols), 3) or cell.shape != (3, 3):
        raise ValueError("Invalid cell or coordinate array shape")
    if not np.isfinite(positions).all() or not np.isfinite(cell).all():
        raise ValueError("Non-finite geometry")
    return cell, symbols, positions, positions @ np.linalg.inv(cell)


def assignment_residual(frac_a, frac_b, symbols_a, symbols_b, cell):
    """One global, species-preserving bijection; no chain-wise alignment."""
    mapping = np.full(len(symbols_a), -1, int)
    residuals = np.zeros(len(symbols_a))
    if Counter(symbols_a) != Counter(symbols_b):
        return None
    for element in sorted(set(symbols_a)):
        a = np.flatnonzero(symbols_a == element)
        b = np.flatnonzero(symbols_b == element)
        vectors, _ = minimum_images(frac_b[b][None, :, :] - frac_a[a][:, None, :], cell)
        distance = np.linalg.norm(vectors, axis=2)
        row, col = linear_sum_assignment(distance)
        mapping[a[row]] = b[col]
        residuals[a[row]] = distance[row, col]
    return {"mapping": mapping.tolist(), "rms_A": float(np.sqrt(np.mean(residuals**2))),
            "max_A": float(residuals.max(initial=0))}


def same_cell_equivalent(record_a, record_b, tolerance=1e-5):
    """Limited exact representation matcher: same cell, permutation + origin.

    This is NOT a general crystal matcher; basis changes, supercells, relaxed
    geometries, local defect matching and statistical ensembles are out of scope.
    """
    cell, sym_a, _, frac_a = structure_arrays(record_a)
    cell_b, sym_b, _, frac_b = structure_arrays(record_b)
    if not np.allclose(cell, cell_b, atol=1e-9, rtol=0) or Counter(sym_a) != Counter(sym_b):
        return {"equivalent": False, "reason": "cell_or_composition_mismatch"}
    rare = min(Counter(sym_a), key=lambda x: Counter(sym_a)[x])
    anchor = np.flatnonzero(sym_a == rare)[0]
    best = None
    for other in np.flatnonzero(sym_b == rare):
        origin = frac_b[other] - frac_a[anchor]
        result = assignment_residual(frac_a + origin, frac_b, sym_a, sym_b, cell)
        result["origin_fractional"] = origin.tolist()
        if best is None or result["max_A"] < best["max_A"]:
            best = result
        if result["max_A"] <= tolerance:
            return {"equivalent": True, **result}
    return {"equivalent": False, **best}


def infer_graph(symbols, fractional, cell):
    n = len(symbols)
    ia, ja = np.triu_indices(n, 1)
    vectors, shifts = minimum_images(fractional[ja] - fractional[ia], cell)
    distances = np.linalg.norm(vectors, axis=1)
    edges, errors = [], []
    for i, j, distance, shift in zip(ia, ja, distances, shifts):
        pair = "".join(sorted((symbols[i], symbols[j])))
        if pair in BOND_LIMITS and distance <= BOND_LIMITS[pair][1]:
            edges.append({"i": int(i), "j": int(j), "shift": shift.tolist(), "distance_A": float(distance)})
            if distance < BOND_LIMITS[pair][0]:
                errors.append({"i": int(i), "j": int(j), "pair": pair, "distance_A": float(distance), "reason": "bond_too_short"})
    adjacency = [[] for _ in range(n)]
    for edge in edges:
        i, j, shift = edge["i"], edge["j"], np.array(edge["shift"], int)
        adjacency[i].append((j, shift))
        adjacency[j].append((i, -shift))
    return edges, adjacency, errors, (ia, ja, distances)


def check_chemistry(symbols, adjacency):
    errors, carbonyls, amides = [], [], []
    for i, symbol in enumerate(symbols):
        count = Counter(symbols[j] for j, _ in adjacency[i])
        valid = False
        if symbol == "H":
            valid = len(adjacency[i]) == 1 and next(iter(count), None) in ("C", "N")
        elif symbol == "O":
            valid = count == {"C": 1}
        elif symbol == "N":
            valid = count == {"C": 2, "H": 1}
            if valid:
                n_carbonyl = sum(any(symbols[k] == "O" for k, _ in adjacency[j]) for j, _ in adjacency[i] if symbols[j] == "C")
                valid = n_carbonyl == 1
                if valid:
                    amides.append(i)
        elif symbol == "C":
            if count == {"O": 1, "N": 1, "C": 1}:
                valid = True
                carbonyls.append(i)
            else:
                valid = count in ({"C": 2, "H": 2}, {"C": 1, "N": 1, "H": 2})
        if not valid:
            errors.append({"index": i, "symbol": str(symbol), "neighbor_elements": dict(count)})
    # Each carbonyl must connect through exactly five CH2 atoms to the next N.
    for start in carbonyls:
        previous = start
        current = next(j for j, _ in adjacency[start] if symbols[j] == "C")
        visited = set()
        count = 0
        while symbols[current] == "C" and current not in visited and count <= 6:
            visited.add(current)
            count += 1
            following = [j for j, _ in adjacency[current] if j != previous and symbols[j] != "H"]
            if len(following) != 1:
                break
            previous, current = current, following[0]
        if count != 5 or symbols[current] != "N":
            errors.append({"carbonyl_index": start, "methylene_count": count, "reason": "not_PA6_backbone"})
    return {"passed": not errors and len(carbonyls) == 8 and len(amides) == 8,
            "carbonyl_count": len(carbonyls), "amide_count": len(amides), "errors": errors}


def periodic_components(adjacency):
    seen, components, membership = set(), [], {}
    for start in range(len(adjacency)):
        if start in seen:
            continue
        potential = {start: np.zeros(3, int)}
        queue = deque([start])
        members, windings = [], []
        while queue:
            i = queue.popleft()
            if i in seen:
                continue
            seen.add(i)
            members.append(i)
            membership[i] = len(components)
            for j, shift in adjacency[i]:
                implied = potential[i] + shift
                if j not in potential:
                    potential[j] = implied
                    queue.append(j)
                else:
                    winding = implied - potential[j]
                    if np.any(winding):
                        windings.append(tuple(int(x) for x in winding))
        unique = sorted(set(windings))
        rank = int(np.linalg.matrix_rank(np.array(unique, float))) if unique else 0
        only_c = bool(unique) and all(v[0] == v[1] == 0 for v in unique)
        gcd_c = math.gcd(*(abs(v[2]) for v in unique)) if unique else 0
        components.append({"indices": sorted(members), "atom_count": len(members), "winding_vectors": unique,
                           "winding_rank": rank, "primitive_c_chain": only_c and gcd_c == 1})
    passed = len(components) == 4 and all(c["atom_count"] == 38 and c["primitive_c_chain"] for c in components)
    return {"passed": passed, "components": components}, membership


def valence_geometry(symbols, fractional, cell, adjacency):
    """Broad independent chemical geometry screen, not a force-field fit."""
    errors, angles, planar = [], [], []
    for i, symbol in enumerate(symbols):
        neighbors = adjacency[i]
        if symbol not in ("C", "N") or len(neighbors) not in (3, 4):
            continue
        vectors = np.array([(fractional[j] + shift - fractional[i]) @ cell for j, shift in neighbors])
        unit = vectors / np.linalg.norm(vectors, axis=1)[:, None]
        local = []
        for a, b in itertools.combinations(range(len(neighbors)), 2):
            angle = float(np.degrees(np.arccos(np.clip(unit[a] @ unit[b], -1, 1))))
            local.append(angle)
        lower, upper = (90., 145.) if len(neighbors) == 3 else (80., 140.)
        angles.append({"center": i, "angles_deg": local, "allowed_range_deg": [lower, upper]})
        if min(local) < lower or max(local) > upper:
            errors.append({"index": i, "reason": "valence_angle_outside_broad_starting_range", "angles_deg": local})
        if len(neighbors) == 3:
            deviation = abs(float(sum(local)) - 360.)
            planar.append({"center": i, "angle_sum_deviation_deg": deviation})
            if deviation > 15.:
                errors.append({"index": i, "reason": "strongly_nonplanar_sp2_center", "angle_sum_deviation_deg": deviation})
    return {"passed": not errors, "errors": errors, "angles": angles, "sp2_planarity": planar,
            "role": "declared broad geometric screens, not literature-derived geometry or optimized minima"}


def validate_roles(record, symbols, adjacency, topology):
    roles = record.get("atom_roles", [])
    errors = []
    if len(roles) != len(symbols) or sorted(r.get("index", -1) for r in roles) != list(range(len(symbols))):
        return {"passed": False, "errors": ["Atom role list must cover every index exactly once"]}
    by_index = {r["index"]: r for r in roles}
    expected_prefix = {"N": "amide_N", "O": "carbonyl_O"}
    for i, element in enumerate(symbols):
        if element == "C":
            expected = "carbonyl_C" if any(symbols[j] == "O" for j, _ in adjacency[i]) else "methylene_C"
        elif element == "H":
            expected = "amide_H" if any(symbols[j] == "N" for j, _ in adjacency[i]) else "methylene_H"
        else:
            expected = expected_prefix.get(element, "invalid")
        if not by_index[i].get("role", "").startswith(expected):
            errors.append({"index": i, "reason": "role_does_not_match_inferred_chemistry", "expected_prefix": expected})
    distinct_chains = set()
    for component in topology["components"]:
        labels = {by_index[i].get("chain_id") for i in component["indices"]}
        distinct_chains.update(labels)
        repeats = Counter(by_index[i].get("repeat_index") for i in component["indices"])
        if len(labels) != 1 or None in labels or len(repeats) != 2 or sorted(repeats.values()) != [19, 19]:
            errors.append({"component": component["indices"][0], "reason": "invalid_chain_or_repeat_labels"})
    if len(distinct_chains) != 4:
        errors.append({"reason": "four_distinct_chain_labels_required"})
    return {"passed": not errors, "errors": errors}


def hydrogen_bonds(symbols, fractional, cell, adjacency, membership):
    bonds, donors = [], []
    oxygens = np.flatnonzero(symbols == "O")
    for n in np.flatnonzero(symbols == "N"):
        for h, nh_shift in adjacency[n]:
            if symbols[h] != "H":
                continue
            donors.append((int(n), int(h)))
            nh = (fractional[h] + nh_shift - fractional[n]) @ cell
            for o in oxygens:
                if membership[n] == membership[o]:
                    continue
                ho, image = minimum_images(fractional[o] - fractional[h], cell)
                ho_length = float(np.linalg.norm(ho))
                no_length = float(np.linalg.norm(ho + nh))
                angle = float(np.degrees(np.arccos(np.clip(np.dot(-nh, ho) / (np.linalg.norm(nh) * max(ho_length, 1e-15)), -1, 1))))
                if no_length <= 3.5 and ho_length <= 2.6 and angle >= 120:
                    normal_spacing = abs(float((fractional[o] + image - fractional[h])[1])) / np.linalg.norm(np.linalg.inv(cell)[:, 1])
                    bonds.append({"N": int(n), "H": int(h), "O": int(o), "oxygen_image_from_H": image.tolist(),
                                  "N_O_A": no_length, "H_O_A": ho_length, "N_H_O_deg": angle,
                                  "donor_chain": membership[n], "acceptor_chain": membership[o],
                                  "normal_to_ac_displacement_A": normal_spacing})
    covered_donors = {b["N"] for b in bonds}
    covered_acceptors = {b["O"] for b in bonds}
    chain_edges = sorted({tuple(sorted((b["donor_chain"], b["acceptor_chain"]))) for b in bonds})
    return {"criterion": {"max_N_O_A": 3.5, "max_H_O_A": 2.6, "min_N_H_O_deg": 120},
            "criterion_role": "optional geometric diagnostic only; never an admission rule or measured bond assignment",
            "donors": donors, "covered_donor_count": len(covered_donors),
            "covered_acceptor_count": len(covered_acceptors), "bonds": bonds, "chain_connections": chain_edges,
            "complete_donor_acceptor_coverage": len(covered_donors) == 8 and len(covered_acceptors) == 8}


def hydrogen_bond_sheet_topology(adjacency, hb):
    """Check infinite sheet connectivity, not just donor/acceptor coverage.

    H--O image shifts and the independently inferred covalent shifts are used
    together. Translation-cycle rank two in a,c demonstrates a periodic sheet;
    pairing chains without independent a winding would remain rank one.
    """
    combined = [[(j, shift.copy()) for j, shift in neighbors] for neighbors in adjacency]
    for bond in hb["bonds"]:
        h, o = bond["H"], bond["O"]
        shift = np.array(bond["oxygen_image_from_H"], int)
        combined[h].append((o, shift))
        combined[o].append((h, -shift))
    topology, _ = periodic_components(combined)
    components = []
    for component in topology["components"]:
        winding = component["winding_vectors"]
        ac_only = bool(winding) and all(v[1] == 0 for v in winding)
        sheet = component["winding_rank"] == 2 and ac_only
        components.append({"indices": component["indices"], "atom_count": component["atom_count"],
                           "winding_vectors": winding, "winding_rank": component["winding_rank"],
                           "periodic_ac_sheet": sheet})
    passed = len(components) == 2 and all(c["atom_count"] == 76 and c["periodic_ac_sheet"] for c in components)
    return {"passed": passed, "components": components,
            "interpretation": "Two infinite ac-connected sheets under the declared broad H-bond graph criterion; not measured or optimized connectivity"}


def validate(record, *, analyze_hbonds=False):
    cell, symbols, positions, fractional = structure_arrays(record)
    formula = dict(Counter(symbols))
    edges, adjacency, errors, pair_data = infer_graph(symbols, fractional, cell)
    chemistry = check_chemistry(symbols, adjacency)
    valence = valence_geometry(symbols, fractional, cell, adjacency)
    topology, membership = periodic_components(adjacency)
    roles = validate_roles(record, symbols, adjacency, topology)
    supplied = record.get("bonds", [])
    given_set = {canonical_edge(e["i"], e["j"], e["shift"]) for e in supplied}
    inferred_set = {canonical_edge(e["i"], e["j"], e["shift"]) for e in edges}
    graph_agreement = bool(supplied) and given_set == inferred_set and len(given_set) == len(supplied)
    collision = []
    for i, j, d in zip(*pair_data):
        near = {a for a, _ in adjacency[i]}
        near.update(b for a, _ in adjacency[i] for b, _ in adjacency[a])
        pair = "".join(sorted((symbols[i], symbols[j])))
        if j not in near and d < COLLISION_FLOORS[pair]:
            collision.append({"i": int(i), "j": int(j), "pair": pair, "distance_A": float(d), "floor_A": COLLISION_FLOORS[pair]})
    symmetry = []
    for rotation, translation in OPS:
        result = assignment_residual(fractional @ rotation.T + translation, fractional, symbols, symbols, cell)
        symmetry.append({"rotation": rotation.tolist(), "translation": translation.tolist(), **result})
    symmetry_pass = all(x["max_A"] < 1e-5 for x in symmetry)
    scan = []
    try:
        import spglib
        for tolerance in (1e-5, 1e-4, .001, .003, .01, .03, .05):
            dataset = spglib.get_symmetry_dataset((cell, fractional, [NUMBERS[x] for x in symbols]), symprec=tolerance)
            scan.append({"symprec_A": tolerance, "status": "identified" if dataset is not None else "failed",
                         "number": int(dataset.number) if dataset is not None else None,
                         "international": str(dataset.international) if dataset is not None else None,
                         "hall_number": int(dataset.hall_number) if dataset is not None else None})
        spglib_status = "completed"
    except ImportError:
        spglib_status = "unavailable"
    geometry_pass = formula == {"C": 48, "H": 88, "N": 8, "O": 8} and chemistry["passed"] and valence["passed"] and topology["passed"] and roles["passed"] and graph_agreement and not errors and not collision and symmetry_pass
    report = {
        "validator_version": VERSION, "scope": "O152 starting geometry, fixed declared basis and origin",
        "status": "passed" if geometry_pass else "failed",
        "passed": bool(geometry_pass),
        "geometry_topology_symmetry_passed": bool(geometry_pass),
        "hydrogen_bond_analysis": {"status": "not_requested", "affects_admission": False},
        "calculation_eligibility_scope": "physical starting-geometry check only; no optimization/DFT execution authorization, numerical validation, or search completeness implied",
        "formula": formula, "volume_A3": float(abs(np.linalg.det(cell))), "chemistry": chemistry,
        "bond_geometry_errors": errors, "valence_geometry": valence, "periodic_topology": topology, "atom_roles": roles,
        "declared_vs_inferred_graph": {"passed": graph_agreement, "declared_edge_count": len(supplied), "inferred_edge_count": len(edges),
                                       "missing_from_declared": sorted(inferred_set - given_set), "missing_from_inferred": sorted(given_set - inferred_set)},
        "collision_floors_A": COLLISION_FLOORS, "nonbonded_collisions": collision,
        "declared_symmetry": {"passed": symmetry_pass, "operations": symmetry,
                              "interpretation": "construction consistency, not independent experimental validation"},
        "spglib_status": spglib_status, "spglib_scan": scan,
        "identity_limitations": ["No arbitrary basis/supercell equivalence", "No relaxed endpoint clustering", "No defect dilution validation", "No experimental validation"],
        "provenance": "No reference coordinates read by validator; source independence requires a separate generation provenance audit",
    }
    if analyze_hbonds:
        hb = hydrogen_bonds(symbols, fractional, cell, adjacency, membership)
        sheets = hydrogen_bond_sheet_topology(adjacency, hb)
        report["hydrogen_bond_analysis"]["status"] = "completed"
        report["hydrogen_bonds"] = hb
        report["hydrogen_bond_sheet_topology"] = sheets
        report["hydrogen_bond_diagnostics_satisfied"] = bool(hb["complete_donor_acceptor_coverage"] and sheets["passed"])
    return report


def self_tests(record=None):
    skew = np.array([[4., 0., 0.], [3.6, 1., 0.], [0., 0., 5.]])
    delta = np.array([.49, .49, 0.])
    exact, _ = minimum_images(delta, skew)
    brute = min(np.linalg.norm((delta + s) @ skew) for s in itertools.product(range(-5, 6), repeat=3))
    tests = {"skew_MIC_vs_independent_brute_force": abs(np.linalg.norm(exact) - brute) < 1e-10,
             "skew_MIC_counterexample_to_fractional_rounding": np.linalg.norm(exact) + .1 < np.linalg.norm((delta - np.rint(delta)) @ skew)}
    if record is not None:
        cell, symbols, positions, fractional = structure_arrays(record)
        rng = np.random.default_rng(20261029)
        perm = rng.permutation(len(symbols))
        transformed = {"cell": cell.tolist(), "symbols": symbols[perm].tolist(),
                       "positions": (((fractional + [.217, -.131, .373])[perm] + rng.integers(-2, 3, size=(len(symbols), 3))) @ cell).tolist()}
        tests["global_origin_permutation_wrap_positive"] = same_cell_equivalent(record, transformed)["equivalent"]
        edges, adjacency, _, _ = infer_graph(symbols, fractional, cell)
        topology, _ = periodic_components(adjacency)
        displaced = fractional.copy()
        displaced[topology["components"][0]["indices"], 2] += .071
        negative = {"cell": cell.tolist(), "symbols": symbols.tolist(), "positions": (displaced @ cell).tolist()}
        result = same_cell_equivalent(record, negative)
        tests["single_chain_relative_c_shift_negative"] = not result["equivalent"]
        tests["negative_shift_max_residual_A"] = result.get("max_A")
    boolean_checks = {k: bool(v) for k, v in tests.items() if isinstance(v, (bool, np.bool_))}
    return {"passed": all(boolean_checks.values()), "checks": {**tests, **boolean_checks},
            "scope": "same-cell exact representation identity only; not complete M2 acceptance"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", type=Path, nargs="?")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--analyze-hbonds", action="store_true",
                        help="Explicit optional hydrogen-bond diagnostics; never changes physical admission")
    args = parser.parse_args()
    if not args.candidate and not args.self_test:
        parser.error("Supply a candidate JSON or --self-test")
    report = {}
    record = None
    if args.candidate:
        raw = args.candidate.read_bytes()
        record = json.loads(raw.decode("utf-8-sig"))
        report = validate(record, analyze_hbonds=args.analyze_hbonds)
        report["input_file"] = str(args.candidate)
        report["input_sha256"] = hashlib.sha256(raw).hexdigest()
        structure_file = args.candidate.with_suffix(".extxyz")
        report["structure_file"] = str(structure_file)
        report["structure_sha256"] = hashlib.sha256(structure_file.read_bytes()).hexdigest() if structure_file.exists() else None
        if structure_file.exists():
            from ase.io import read
            atoms = read(structure_file)
            cell, symbols, positions, _ = structure_arrays(record)
            match = atoms.get_chemical_symbols() == symbols.tolist() and np.allclose(atoms.cell.array, cell, atol=1e-7, rtol=0) and np.allclose(atoms.positions, positions, atol=1e-7, rtol=0) and bool(np.all(atoms.pbc))
            report["extxyz_matches_validated_JSON"] = bool(match)
            if not match:
                report["status"] = "failed"
                report["passed"] = False
    if args.self_test:
        report["self_tests"] = self_tests(record)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"validator": VERSION, "status": report.get("status"), "self_tests_passed": report.get("self_tests", {}).get("passed"), "output": str(args.output) if args.output else None}, ensure_ascii=False))
    return 0 if report.get("status", "passed") == "passed" and report.get("self_tests", {}).get("passed", True) else 2


if __name__ == "__main__":
    raise SystemExit(main())
