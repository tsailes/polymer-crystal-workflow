"""Audit approximate periodic matching with an optimized global origin.

Only reads completed structures and performs coordinate arithmetic. No quantum
calculator is imported or evaluated. Outputs are confined to this audit folder.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import sys

import numpy as np
from ase.io import read
from scipy.optimize import linear_sum_assignment

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def assignment(first_fractional, second_fractional, cell, first_numbers, second_numbers):
    delta = first_fractional[:, None, :] - second_fractional[None, :, :]
    delta -= np.rint(delta)
    squared_distance = np.sum((delta @ cell)**2, axis=2)
    squared_distance[first_numbers[:, None] != second_numbers[None, :]] = 1e30
    rows, cols = linear_sum_assignment(squared_distance)
    mapping = np.empty(len(first_fractional), int)
    mapping[rows] = cols
    assert np.array_equal(first_numbers, second_numbers[mapping])
    residual_fractional = first_fractional - second_fractional[mapping]
    residual_fractional -= np.rint(residual_fractional)
    return mapping, residual_fractional


def record(mapping, residual, cell, signs, shift, steps, start):
    distance = np.linalg.norm(residual @ cell, axis=1)
    return {
        "rms_A": float(np.sqrt(np.mean(distance**2))),
        "max_A": float(distance.max()),
        "signs": list(signs),
        "origin_shift_fractional": (shift % 1).tolist(),
        "atom_mapping": mapping.tolist(),
        "mean_residual_cartesian_A": np.mean(residual @ cell, axis=0).tolist(),
        "origin_refinement_steps": steps,
        "anchor_target_index": start,
    }


def compare(first, second, rms_tol=0.02, max_tol=0.05):
    if len(first) != len(second) or sorted(first.numbers) != sorted(second.numbers):
        raise ValueError("Atom count or composition differ")
    if not all(first.pbc) or not all(second.pbc):
        raise ValueError("Three-dimensional PBC required")
    cell = first.cell.array
    if not np.allclose(cell, second.cell.array, atol=1e-10, rtol=0):
        raise ValueError("Cells differ")
    if not np.allclose(cell, np.diag(np.diag(cell)), atol=1e-12, rtol=0):
        raise ValueError("This image-distance implementation requires orthorhombic cells")
    a = first.get_scaled_positions(wrap=False)
    b = second.get_scaled_positions(wrap=False)
    anchor = int(np.flatnonzero(first.numbers == 6)[0])
    original_records, optimized_records = [], []
    for signs in itertools.product((-1, 1), repeat=3):
        transformed = a * signs
        for target in np.flatnonzero(second.numbers == 6):
            shift = b[target] - transformed[anchor]
            mapping, residual = assignment(
                transformed + shift, b, cell, first.numbers, second.numbers)
            original_records.append(record(mapping, residual, cell, signs, shift, 0, int(target)))
            previous_objective = float(np.sum((residual @ cell)**2))
            for step in range(1, 101):
                correction = residual.mean(axis=0)
                shift -= correction
                next_mapping, next_residual = assignment(
                    transformed + shift, b, cell, first.numbers, second.numbers)
                objective = float(np.sum((next_residual @ cell)**2))
                assert objective <= previous_objective + 1e-10
                converged = (np.array_equal(mapping, next_mapping) and
                             np.linalg.norm(correction @ cell) < 1e-12)
                mapping, residual = next_mapping, next_residual
                previous_objective = objective
                if converged:
                    break
            else:
                raise RuntimeError("Origin/assignment refinement did not converge")
            optimized_records.append(record(
                mapping, residual, cell, signs, shift, step, int(target)))

    def best(records):
        return min(records, key=lambda r: (r["rms_A"], r["max_A"]))

    best_anchor = best(original_records)
    best_optimized = best(optimized_records)
    # Both bounds must hold for the same mapping. An RMS-only winner can in
    # principle have a larger maximum error than another valid mapping.
    feasible = [r for r in optimized_records
                if r["rms_A"] <= rms_tol and r["max_A"] <= max_tol]
    return {
        "thresholds_A": {"rms": rms_tol, "maximum": max_tol},
        "anchor_only_squared_assignment_best": best_anchor,
        "origin_optimized_best": best_optimized,
        "origin_optimized_match": best(feasible) if feasible else None,
        "origin_optimized_passes": bool(feasible),
        "anchor_only_passes": best_anchor["rms_A"] <= rms_tol
                              and best_anchor["max_A"] <= max_tol,
        "rms_improvement_A": best_anchor["rms_A"] - best_optimized["rms_A"],
        "optimized_rms_margin_A": rms_tol - best_optimized["rms_A"],
        "optimized_max_margin_A": max_tol - best_optimized["max_A"],
        "number_of_anchor_isometry_starts": len(optimized_records),
        "max_origin_refinement_steps": max(r["origin_refinement_steps"]
                                          for r in optimized_records),
    }


def load_candidate(candidate):
    status = json.loads((ROOT / "candidates" / candidate / "relaxation/status.json").read_text())
    if status["status"] != "preoptimized":
        raise ValueError(f"{candidate} is not completed")
    result_path = Path(status["result"])
    result = json.loads(result_path.read_text())
    if result.get("passed") is not True:
        raise ValueError(f"{candidate} failed final validation")
    model_path = result_path.parent / "final.extxyz"
    return read(model_path), result, {
        "candidate_id": candidate, "model_path": str(model_path),
        "model_sha256": hashlib.sha256(model_path.read_bytes()).hexdigest(),
        "energy_eV_per_CH2": result["energy_eV_per_CH2"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--first", default="cand_0025")
    parser.add_argument("--second", default="cand_0041")
    args = parser.parse_args()
    a, result_a, meta_a = load_candidate(args.first)
    b, result_b, meta_b = load_candidate(args.second)
    spec = importlib.util.spec_from_file_location("main_generator_for_comparison",
                                                 ROOT / "generate_candidates.py")
    gen = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gen)
    actual_old = gen.equivalent(a, b, rms_tol=100, max_tol=100)
    actual_old_pass = gen.equivalent(a, b, rms_tol=0.02, max_tol=0.05)
    report = compare(a, b)
    difference_meV = abs(result_a["energy_eV_per_CH2"] - result_b["energy_eV_per_CH2"]) * 1000
    report.update({
        "inputs": [meta_a, meta_b],
        "energy_difference_meV_per_CH2": difference_meV,
        "energy_threshold_passed": difference_meV <= 1.0,
        "actual_generator_anchor_only_best": actual_old,
        "actual_generator_anchor_only_passes": actual_old_pass is not None,
        "quantum_calculations_performed": False,
        "algorithm": "For each orthorhombic sign isometry and carbon-anchor seed, "
                     "alternate element-preserving minimum squared-distance assignment "
                     "with subtracting the mean minimum-image residual from the global "
                     "translation. Keep the mapping/image-local optimum and test RMS "
                     "and maximum distance for the same map.",
        "limitation": "Multistart local translation/assignment refinement, not a general "
                      "global optimality proof for arbitrary molecular crystals.",
    })
    destination = HERE / f"approx_match_{args.first}_{args.second}.json"
    destination.write_text(json.dumps(report, ensure_ascii=False, indent=2,
                                     allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "audit_file": str(destination), "energy_difference_meV_CH2": difference_meV,
        "original_anchor_rms_A": actual_old["rms_A"], "original_anchor_max_A": actual_old["max_A"],
        "original_anchor_pass": actual_old_pass is not None,
        "optimized_rms_A": report["origin_optimized_best"]["rms_A"],
        "optimized_max_A": report["origin_optimized_best"]["max_A"],
        "optimized_match_pass": report["origin_optimized_passes"],
        "rms_margin_A": report["optimized_rms_margin_A"],
        "max_margin_A": report["optimized_max_margin_A"],
    }))


if __name__ == "__main__":
    main()
