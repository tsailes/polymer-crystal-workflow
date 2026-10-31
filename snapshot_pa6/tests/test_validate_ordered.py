"""Independent geometric regression tests, without any reference coordinates."""
import importlib.util
import itertools
import json
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

SPEC = importlib.util.spec_from_file_location("validate_ordered", Path(__file__).parents[1] / "scripts/validate_ordered.py")
validator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(validator)


def frozen_grid_record(x, phase):
    """Recreate audited grid points by rigid placement of a frozen independent seed.

    Only whole-chain translations are applied. No template construction,
    minimization, energy evaluation, or literature coordinates are involved.
    Coordinates remain unwrapped so the original covalent image shifts apply.
    """
    source = Path(__file__).parents[1] / "branches/O/generated/O-c0001.json"
    saved = json.loads(source.read_text(encoding="utf-8-sig"))
    parameters = saved["generation_parameters"]
    assert saved["chemistry_template_id"] == "PA6-independent-chain-00"
    assert parameters["chain_rotation_about_c_deg"] == 0.0
    record = {key: saved[key] for key in ("cell", "symbols", "atom_roles", "bonds")}
    fractional = np.array(saved["fractional_positions"])
    shift = np.array([x - parameters["chain_x_fractional"], 0,
                      phase - parameters["chain_c_phase_fractional"]])
    for role in record["atom_roles"]:
        rotation = validator.OPS[role["symmetry_operation_index"]][0]
        fractional[role["index"]] += shift @ rotation.T
    record["fractional_positions"] = fractional.tolist()
    record["positions"] = (fractional @ np.array(record["cell"])).tolist()
    return record


class GeometricChecks(unittest.TestCase):
    def test_skew_cell_against_independent_brute_force(self):
        rng = np.random.default_rng(68134)
        cell = np.array([[4., 0., 0.], [3.6, 1., 0.], [.6, .4, 5.]])
        delta = rng.uniform(-2., 2., (50, 3))
        found, shifts = validator.minimum_images(delta, cell)
        offsets = np.array(list(itertools.product(range(-5, 6), repeat=3)))
        for d, result in zip(delta, found):
            independent = np.min(np.linalg.norm((d + offsets) @ cell, axis=1))
            self.assertAlmostEqual(np.linalg.norm(result), independent, places=10)
        np.testing.assert_allclose(found, (delta + shifts) @ cell, atol=1e-12)

    def test_rounding_is_not_a_skew_cell_distance_algorithm(self):
        self.assertTrue(validator.self_tests()["passed"])

    def test_four_independent_infinite_chains(self):
        adjacency = [[] for _ in range(152)]
        for chain in range(4):
            for k in range(38):
                i, j = chain * 38 + k, chain * 38 + (k + 1) % 38
                shift = np.array([0, 0, int(k == 37)])
                adjacency[i].append((j, shift))
                adjacency[j].append((i, -shift))
        result, _ = validator.periodic_components(adjacency)
        self.assertTrue(result["passed"])
        # A finite ring is not a translational polymer even with identical counts.
        adjacency[37][-1] = (0, np.zeros(3, int))
        adjacency[0][-1] = (37, np.zeros(3, int))
        result, _ = validator.periodic_components(adjacency)
        self.assertFalse(result["passed"])

    def test_permutation_origin_and_images_do_not_change_identity(self):
        cell = np.array([[9., 0., 0.], [2.7, 7., 0.], [0., 0., 17.]])
        a = {"cell": cell.tolist(), "symbols": ["N", "O", "C", "H"],
             "fractional_positions": [[.12, .3, .7], [.25, .6, .2], [.78, .21, .53], [.4, .8, .11]]}
        f = np.array(a["fractional_positions"])
        perm = [2, 0, 3, 1]
        b = {"cell": cell.tolist(), "symbols": [a["symbols"][i] for i in perm],
             "fractional_positions": (f[perm] + [.15, -.42, .31] + [[1, 0, 0], [0, 0, 1], [-1, 2, 0], [0, -1, 0]]).tolist()}
        self.assertTrue(validator.same_cell_equivalent(a, b)["equivalent"])
        b["fractional_positions"][0][2] += .071
        self.assertFalse(validator.same_cell_equivalent(a, b)["equivalent"])

    def test_hbond_pairing_is_not_an_infinite_sheet(self):
        adjacency = [[] for _ in range(152)]
        for chain in range(4):
            for k in range(38):
                i, j = chain * 38 + k, chain * 38 + (k + 1) % 38
                shift = np.array([0, 0, int(k == 37)])
                adjacency[i].append((j, shift))
                adjacency[j].append((i, -shift))
        hb = {"bonds": [{"H": 0, "O": 38, "oxygen_image_from_H": [0, 0, 0]},
                        {"H": 76, "O": 114, "oxygen_image_from_H": [0, 0, 0]}]}
        paired = validator.hydrogen_bond_sheet_topology(adjacency, hb)
        self.assertFalse(paired["passed"])
        self.assertEqual([c["winding_rank"] for c in paired["components"]], [1, 1])
        hb["bonds"].extend([{"H": 2, "O": 40, "oxygen_image_from_H": [1, 0, 0]},
                            {"H": 78, "O": 116, "oxygen_image_from_H": [1, 0, 0]}])
        sheets = validator.hydrogen_bond_sheet_topology(adjacency, hb)
        self.assertTrue(sheets["passed"])
        self.assertEqual([c["winding_rank"] for c in sheets["components"]], [2, 2])


class AdmissionWithoutHydrogenBonds(unittest.TestCase):
    def test_default_does_not_call_hbond_analysis(self):
        # grid-0329 was rejected by the historical full-coverage H-bond gate.
        record = frozen_grid_record(.225, 20 / 96)
        with patch.object(validator, "hydrogen_bonds", side_effect=AssertionError("HB analysis was called")) as bonds, \
             patch.object(validator, "hydrogen_bond_sheet_topology", side_effect=AssertionError("HB topology was called")) as sheets:
            report = validator.validate(record)
        bonds.assert_not_called()
        sheets.assert_not_called()
        self.assertTrue(report["passed"])
        self.assertEqual(report["hydrogen_bond_analysis"]["status"], "not_requested")
        self.assertFalse(report["hydrogen_bond_analysis"]["affects_admission"])
        self.assertNotIn("hydrogen_bonds", report)
        self.assertNotIn("hydrogen_bond_sheet_topology", report)

    def test_optional_incomplete_hbond_diagnostic_cannot_reject(self):
        report = validator.validate(frozen_grid_record(.225, 20 / 96), analyze_hbonds=True)
        self.assertTrue(report["passed"])
        self.assertEqual(report["hydrogen_bonds"]["covered_donor_count"], 4)
        self.assertFalse(report["hydrogen_bonds"]["complete_donor_acceptor_coverage"])
        self.assertFalse(report["hydrogen_bond_sheet_topology"]["passed"])
        self.assertFalse(report["hydrogen_bond_diagnostics_satisfied"])
        self.assertEqual(report["hydrogen_bond_analysis"]["status"], "completed")
        self.assertFalse(report["hydrogen_bond_analysis"]["affects_admission"])

    def test_severe_nonbonded_collision_still_rejects(self):
        # grid-0071 has four O...O overlaps (~1.645 A), independent of H bonds.
        report = validator.validate(frozen_grid_record(.25, 35 / 96))
        self.assertFalse(report["passed"])
        collisions = report["nonbonded_collisions"]
        self.assertEqual(len(collisions), 4)
        self.assertTrue(all(item["pair"] == "OO" for item in collisions))
        self.assertAlmostEqual(min(item["distance_A"] for item in collisions), 1.6448161227282667, places=8)


if __name__ == "__main__":
    unittest.main()
