"""Cheap independent code/data audit. No calculator evaluation or optimization."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import spglib
from ase.io import read
from scipy.optimize import linear_sum_assignment

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def match(first, second):
    assert len(first) == len(second)
    assert sorted(first.numbers) == sorted(second.numbers)
    assert np.allclose(first.cell, second.cell, atol=1e-10, rtol=0)
    a = first.get_scaled_positions(wrap=False)
    b = second.get_scaled_positions(wrap=False)
    delta = a[:, None] - b[None, :]
    delta -= np.rint(delta)
    distance = np.linalg.norm(delta @ first.cell.array, axis=2)
    distance[first.numbers[:, None] != second.numbers[None, :]] = 1e30
    rows, cols = linear_sum_assignment(distance)
    mapping = np.empty(len(first), int)
    mapping[rows] = cols
    return mapping, float(distance[rows, cols].max())


def canonical_edge(i, j, shift):
    return min((int(i), int(j), *map(int, shift)),
               (int(j), int(i), *map(int, -np.asarray(shift))))


def signature(atoms):
    f = atoms.get_scaled_positions()
    d = f[:, None] - f[None, :]
    d -= np.rint(d)
    distance = np.linalg.norm(d @ atoms.cell.array, axis=2)
    signature_parts = {}
    for i in range(len(atoms)):
        for j in range(i):
            key = tuple(sorted((int(atoms.numbers[i]), int(atoms.numbers[j]))))
            signature_parts.setdefault(key, []).append(distance[i, j])
    return np.concatenate([np.sort(signature_parts[k]) for k in sorted(signature_parts)])


def main():
    independent = module("independent_construction", HERE / "independent_pnam_audit.py")
    gen = module("audited_generator", ROOT / "generate_candidates.py")
    relax = module("audited_relaxer", ROOT / "relax_candidates.py")
    manifest = json.loads((ROOT / "generation_manifest.json").read_text())
    rotations, translations, source = relax.target_operations(manifest)
    database = spglib.get_symmetry_from_database(297)

    def key(r, t):
        return tuple(np.rint(r).astype(int).ravel()), tuple(np.round(np.asarray(t) % 1, 12))

    assert {key(r, t) for r, t in zip(rotations, translations)} == {
        key(r, t) for r, t in zip(database["rotations"], database["translations"])}
    records, representatives = [], {}
    for entry in manifest["unique_survivors"]:
        cid = entry["candidate_id"]
        path = ROOT / entry["input_path"]
        atoms = read(path)
        p = entry["parameters"]
        independent_atoms = independent.make_pe(
            p["setting_angle_deg"], p["center_fractional"], p["carbon_phase_fractional"])
        mapping, construction_error = match(independent_atoms, atoms)
        assert construction_error < 1e-7
        # Transfer independently declared unwrapped chemical graph into the
        # candidate's wrapped atom order and image convention.
        offsets = np.rint(independent_atoms.get_scaled_positions(wrap=False) -
                          atoms.get_scaled_positions(wrap=False)[mapping]).astype(int)
        expected_edges = set()
        for i, neighbors in independent.bonded_graph().items():
            for j, shift in neighbors:
                expected_edges.add(canonical_edge(mapping[i], mapping[j],
                                                   shift + offsets[j] - offsets[i]))
        topology = json.loads((path.parent / "reference_topology.json").read_text())
        actual_edges = {canonical_edge(e["i"], e["j"], e["shift"])
                        for e in topology["periodic_bonds"]}
        assert expected_edges == actual_edges
        independent_contacts = independent.contact_check(independent_atoms)
        generated_contacts = gen.audit_geometry(atoms, topology)
        cc_delta = abs(independent_contacts["min_interchain_CC_A"] -
                       generated_contacts["minimum_interchain_CC"]["distance_A"])
        ratio_delta = abs(independent_contacts["min_nonbonded_r_over_vdw_sum"] / 0.85 -
                          generated_contacts["minimum_nonbonded_ratio"]["ratio_to_vdw_threshold"])
        assert cc_delta < 1e-7 and ratio_delta < 1e-7
        assert independent_contacts["passes_both"] == generated_contacts["passed"]
        assert independent_contacts["violating_nonbonded_pairs_per_cell"] == len(
            generated_contacts["vdw_clashes"])

        constraint = relax.TargetPnam(atoms, rotations, translations)
        basis = np.eye(3 * len(atoms))
        projector = np.column_stack([
            constraint.project(atoms, v.reshape(-1, 3)).ravel() for v in basis])
        rank = int(np.linalg.matrix_rank(projector, tol=1e-10))
        assert rank == 6
        assert np.max(np.abs(projector - projector.T)) < 1e-12
        assert np.max(np.abs(projector @ projector - projector)) < 1e-12
        vectors = np.random.default_rng(1973).normal(size=(len(atoms), 3))
        independent_projection = np.zeros_like(vectors)
        for rotation, mapping_op in zip(rotations, constraint.maps):
            cart_rotation = atoms.cell.array.T @ rotation @ np.linalg.inv(atoms.cell.array.T)
            independent_projection[mapping_op] += vectors @ cart_rotation.T / len(rotations)
        projector_error = float(np.max(np.abs(
            independent_projection - constraint.project(atoms, vectors))))
        assert projector_error < 1e-12
        relax.prepare_geometry_context(path)
        projection_test = relax.geometric_projection_test(atoms, constraint)
        assert projection_test["passed"]

        file_results = {}
        for filename in ("input.extxyz", "input.cif", "POSCAR"):
            file_path = path.parent / filename
            restored = read(file_path)
            _, maximum = match(atoms, restored)
            assert maximum < 1e-7
            assert all(restored.pbc)
            file_results[filename] = {
                "species_preserved": True, "periodic_positions_max_error_A": maximum,
                "sha256": hashlib.sha256(file_path.read_bytes()).hexdigest()}
        records.append({
            "candidate_id": cid, "setting_angle_deg": p["setting_angle_deg"],
            "independent_construction_max_error_A": construction_error,
            "periodic_graph_matches_independent_chemistry": True,
            "nonbonded_exclusions_and_contacts_match": True,
            "interchain_CC_comparison_error_A": cc_delta,
            "vdw_threshold_ratio_comparison_error": ratio_delta,
            "target_group_projector_rank": rank,
            "independent_cartesian_projection_max_error": projector_error,
            "projection_validation": projection_test, "file_roundtrips": file_results,
        })
        representatives[cid] = atoms

    invariant_differences = []
    keys = list(representatives)
    for index, first in enumerate(keys):
        for second in keys[:index]:
            delta = float(np.max(np.abs(signature(representatives[first]) -
                                        signature(representatives[second]))))
            invariant_differences.append(delta)
            assert delta > 1e-6  # Necessary invariant differs: no exact isometry can equate them.
    duplicate_error = 0.0
    duplicates = 0
    for entry in manifest["records"]:
        if entry["status"] != "exact_symmetry_duplicate":
            continue
        p, equivalence = entry["parameters"], entry["equivalence"]
        original, _ = gen.build(p["setting_angle_deg"], p["center_fractional"],
                                p["carbon_phase_fractional"])
        target = representatives[entry["representative_id"]]
        fractional = original.get_scaled_positions()
        transformed = (fractional * equivalence["signs"] +
                       equivalence["origin_shift_fractional"])
        mapping = np.array(equivalence["atom_mapping"])
        assert np.array_equal(original.numbers, target.numbers[mapping])
        difference = transformed - target.get_scaled_positions()[mapping]
        difference -= np.rint(difference)
        error = float(np.linalg.norm(difference @ original.cell.array, axis=1).max())
        assert error < 1e-7
        duplicate_error = max(duplicate_error, error)
        duplicates += 1

    result = {
        "passed": True, "quantum_calculations_performed": False,
        "experimental_refined_coordinates_read": False,
        "scope": "Generated candidate data and target-group algebra, before quantum optimization",
        "reviewed_scripts": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
                            (ROOT / "generate_candidates.py", ROOT / "relax_candidates.py")},
        "manifest_operations_match_spglib_hall_297": True,
        "symmetry_source": source, "candidate_count": len(records),
        "file_roundtrip_count": 3 * len(records),
        "validated_duplicate_mappings": duplicates,
        "duplicate_mapping_max_error_A": duplicate_error,
        "minimum_pair_signature_max_difference_A": min(invariant_differences),
        "all_distinct_representatives_have_different_isometry_invariants": True,
        "records": records,
        "nonblocking_scope_notes": [
            "Final preoptimized status tests the stationary point and chemical safety "
            "bounds, not original strict nonbonded filters. Consult the separate "
            "final_original_geometric_filters_passed field before calling the final "
            "structure a survivor of those filters.",
            "2x3x8 sampling was validated on the pilot. A new candidate's converged "
            "preoptimization is not by itself a completed dense-grid convergence proof.",
            "Stored wrap_image and immutable edge shifts remain valid because "
            "the optimizer does not independently wrap atom coordinates.",
        ],
    }
    destination = HERE / "generation_constraint_review.json"
    destination.write_text(json.dumps(result, ensure_ascii=False, indent=2,
                                     allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in
                      ("passed", "candidate_count", "file_roundtrip_count",
                       "validated_duplicate_mappings", "duplicate_mapping_max_error_A",
                       "minimum_pair_signature_max_difference_A")}))


if __name__ == "__main__":
    main()
