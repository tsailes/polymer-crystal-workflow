"""Read-only delivery validation; writes only audits/final_delivery_audit.json.

Run in the project WSL Python environment. Exit codes:
0 = passed; 1 = inconsistent/failed evidence; 2 = incomplete delivery.
No quantum backend, optimizer, subprocess, or source-coordinate template is used.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import sys

sys.dont_write_bytecode = True

import gemmi
import numpy as np
from ase.io import read
from ase.units import Hartree
from scipy.optimize import linear_sum_assignment

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "final_delivery_audit.json"
CELL = np.diag([7.432, 4.929, 2.545])


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def close(actual, expected, tolerance, message):
    a, b = np.asarray(actual, float), np.asarray(expected, float)
    require(a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all()
            and np.max(np.abs(a - b), initial=0.0) <= tolerance, message)


def path_of(value):
    value = str(value).replace("\\", "/")
    if value.startswith("E:/") and ROOT.as_posix().startswith("/mnt/e/"):
        value = "/mnt/e/" + value[3:]
    p = Path(value)
    return p if p.is_absolute() else ROOT / p


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


GENERATOR = None


def periodic_match(first, second):
    """Independent element-wise bijection; no origin changes for file roundtrips."""
    require(len(first) == len(second), "Roundtrip atom count changed")
    require(sorted(first.numbers) == sorted(second.numbers), "Roundtrip composition changed")
    close(first.cell.array, second.cell.array, 1e-10, "Roundtrip cell changed")
    delta = first.get_scaled_positions()[:, None] - second.get_scaled_positions()[None, :]
    delta -= np.rint(delta)
    cost = np.sum((delta @ first.cell.array)**2, axis=2)
    cost[first.numbers[:, None] != second.numbers[None, :]] = 1e30
    i, j = linear_sum_assignment(cost)
    return float(np.sqrt(cost[i, j]).max())


def structure_check(atoms):
    require(len(atoms) == 12 and np.sum(atoms.numbers == 6) == 4
            and np.sum(atoms.numbers == 1) == 8, "Expected C4H8 / 12 atoms")
    require(all(atoms.pbc), "Three-dimensional PBC missing")
    require(np.isfinite(atoms.positions).all(), "Nonfinite atom coordinates")
    close(atoms.cell.array, CELL, 1e-10, "Fixed cell changed")
    f = atoms.get_scaled_positions()
    maximum = 0.0
    for operation in gemmi.SpaceGroup("P n a m").operations():
        r, t = np.array(operation.rot) / operation.DEN, np.array(operation.tran) / operation.DEN
        transformed = atoms.copy()
        transformed.set_scaled_positions(f @ r.T + t)
        maximum = max(maximum, periodic_match(transformed, atoms))
    require(maximum < 1e-5, "Target Pnam operations not preserved")
    return {"formula": "C4H8", "atom_count": 12, "fixed_cell": True,
            "full_pbc": True, "target_Pnam_max_error_A": maximum}


def file_family(folder, stem):
    reference = read(folder / f"{stem}.extxyz")
    structure = structure_check(reference)
    validation = load(folder / "file_validation.json")
    files = []
    for name in (f"{stem}.extxyz", f"{stem}.cif", "POSCAR"):
        path = folder / name
        restored = read(path)
        structure_check(restored)
        error = periodic_match(reference, restored)
        require(error <= 2e-7, f"File roundtrip coordinate error: {path}")
        require(validation[name]["passed"] is True, f"Stored file validation failed: {path}")
        require(validation[name]["sha256"] == sha(path), f"File validation hash mismatch: {path}")
        files.append({"path": str(path), "sha256": sha(path), "maximum_roundtrip_error_A": error})
    return reference, {"structure": structure, "files": files}


def native_final_energy(path, repeats, energy_per_cell):
    content = Path(path).read_text(encoding="utf-8", errors="replace")
    require("Traceback (most recent call last)" not in content, f"Traceback in {path}")
    require(not re.search(r"SCF\s+(?:failed|not\s+converged)|did\s+not\s+converge",
                          content, flags=re.I), f"SCF failure in {path}")
    pattern = re.compile(r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)"
                         r"\s+[-+\d.Ee]+\s+[-+\d.Ee]+\s*$", re.M)
    rows = list(pattern.finditer(content))
    require(bool(rows), f"Native SCF cycle table missing: {path}")
    last = rows[-1]
    require("total:" in content[last.end():], f"Last native SCF timing not completed: {path}")
    cycle, hartree_energy = int(last.group(1)), float(last.group(2))
    require(cycle < 500, f"Native SCF reached iteration ceiling: {path}")
    folded = hartree_energy * Hartree / np.prod(repeats)
    close(folded, energy_per_cell, 1e-7, f"Native final energy inconsistent: {path}")
    return {"path": str(path), "sha256": sha(path), "last_scf_cycle": cycle,
            "native_folded_energy_eV_per_cell": float(folded)}


def completed_candidate(entry, status, row):
    cid = entry["candidate_id"]
    initial_path = ROOT / entry["input_path"]
    initial, initial_family = file_family(initial_path.parent, "input")
    require(sha(initial_path) == status["input_sha256"], f"{cid}: input hash changed")
    topology_path = initial_path.parent / "reference_topology.json"
    topology = load(topology_path)
    require(sha(topology_path) == status["reference_topology_sha256"], f"{cid}: topology hash changed")
    initial_audit = GENERATOR.audit_geometry(initial, topology)
    require(initial_audit["passed"], f"{cid}: generated survivor now fails geometry")
    result_path = path_of(status["result"])
    result = load(result_path)
    require(result["passed"] is True and result["status"] == "preoptimized",
            f"{cid}: result contradicts successful status")
    require(path_of(row["final_result_path"]).resolve() == result_path.resolve(),
            f"{cid}: summary points to another result")
    require(row["final_result_sha256"] == sha(result_path), f"{cid}: result hash mismatch")
    final_path = result_path.parent / "final.extxyz"
    final, final_family = file_family(result_path.parent, "final")
    require(path_of(row["final_model_path"]).resolve() == final_path.resolve(),
            f"{cid}: summary points to another model")
    require(row["final_model_sha256"] == sha(final_path), f"{cid}: model hash mismatch")
    require(np.array_equal(final.numbers, initial.numbers), f"{cid}: atom order changed")
    for name in ("chain_id", "wrap_image"):
        require(np.array_equal(final.arrays[name], initial.arrays[name]),
                f"{cid}: immutable {name} changed")
    raw = np.asarray(result["raw_forces_eV_A"], float)
    projected = np.asarray(result["projected_forces_eV_A"], float)
    require(raw.shape == projected.shape == (12, 3), f"{cid}: invalid force shape")
    require(np.isfinite(raw).all() and np.isfinite(projected).all(), f"{cid}: nonfinite forces")
    fmax, pmax = float(np.linalg.norm(raw, axis=1).max()), float(np.linalg.norm(projected, axis=1).max())
    require(fmax < 0.01 and pmax < 0.01, f"{cid}: force tolerance not met")
    close(fmax, result["raw_fmax_eV_A"], 1e-10, f"{cid}: raw force norm mismatch")
    close(pmax, result["projected_fmax_eV_A"], 1e-10, f"{cid}: projected force norm mismatch")
    close(raw, final.arrays["raw_forces"], 1e-7, f"{cid}: exported raw forces differ")
    close(projected, final.arrays["projected_forces"], 1e-7, f"{cid}: exported projected forces differ")
    energy = result["energy_eV_per_CH2"]
    close(energy, final.info["energy_eV_per_CH2"], 1e-10, f"{cid}: geometry energy metadata differs")
    close(energy, row["energy_eV_per_CH2"], 1e-10, f"{cid}: summary energy differs")
    close(4 * energy, result["energy_eV_per_cell"], 1e-10, f"{cid}: energy units inconsistent")
    folded = result["folded_result"]
    require(folded["repeats"] == [2, 3, 8] and folded["natoms"] == 576,
            f"{cid}: unexpected fine sampling")
    close(folded["energy_per_CH2_eV"], energy, 1e-10, f"{cid}: folded energy mismatch")
    close(folded["forces_eV_per_A"], raw, 1e-10, f"{cid}: folded force mismatch")
    require(folded["replica_force_max_deviation_eV_per_A"] < 1e-6,
            f"{cid}: non-equivalent replica forces")
    fine = load(result_path.parent / "fine/result.json")
    require(fine["converged"] is True, f"{cid}: fine optimizer did not converge")
    de = energy - fine["final"]["energy_eV_per_CH2"]
    df = float(np.linalg.norm(raw - np.asarray(fine["final"]["raw_forces_eV_A"]), axis=1).max())
    close(de, result["accuracy_energy_change_eV_per_CH2"], 1e-10, f"{cid}: accuracy delta E mismatch")
    close(df, result["accuracy_force_change_eV_A"], 1e-10, f"{cid}: accuracy delta F mismatch")
    require(abs(de) < 0.001 and df < 0.005, f"{cid}: tighter accuracy test failed")
    fine_geometry = read(result_path.parent / "fine/final.extxyz")
    close(fine_geometry.positions, final.positions, 1e-10, f"{cid}: cold check changed geometry")
    geometry = GENERATOR.audit_geometry(final, topology)
    require(bool(row["geometry_passed"]) == bool(geometry["passed"]),
            f"{cid}: final geometry classification mismatch")
    require(bool(result["final_original_geometric_filters_passed"]) == bool(geometry["passed"]),
            f"{cid}: stored final geometric filters inconsistent")
    require(bool(row["final_CC_strict_1p53_1p54"]) == bool(geometry["CC_strict_1p53_1p54_passed"]),
            f"{cid}: strict C-C interval misreported")
    close([row["final_CC_min_A"], row["final_CC_max_A"]], geometry["CC_range_A"], 1e-8,
          f"{cid}: C-C range misreported")
    close(row["final_min_interchain_CC_A"], geometry["minimum_interchain_CC"]["distance_A"],
          1e-8, f"{cid}: interchain C-C misreported")
    require(row["final_detected_SG"] == geometry["symmetry"]["detected_number"],
            f"{cid}: detected space group misreported")
    angles = []
    for edge in topology["periodic_bonds"]:
        if edge["kind"] != "CC":
            continue
        vector = (final.positions[edge["j"]] + np.asarray(edge["shift"]) @ final.cell.array
                  - final.positions[edge["i"]])
        angle = float(np.degrees(np.arctan2(vector[1], vector[0])) % 180)
        angles.append(min(angle, 180 - angle))
    close(angles, [row["final_chain_plane_angle_deg"]] * len(angles), 1e-7,
          f"{cid}: catalog chain angle disagrees with periodic bonds")
    return final, result, {
        "candidate_id": cid, "passed": True, "initial_files": initial_family,
        "final_files": final_family, "energy_eV_per_CH2": energy, "raw_fmax_eV_A": fmax,
        "geometry_filters_passed": bool(geometry["passed"]),
        "strict_CC_interval_passed": bool(geometry["CC_strict_1p53_1p54_passed"]),
        "detected_spacegroup": geometry["symmetry"]["detected_number"],
        "target_Pnam_contained": geometry["symmetry"]["contains_target_group"],
        "native_log": native_final_energy(result_path.parent / "native_output.log",
                                          (2, 3, 8), 4 * energy)}


def directional_exclusion(cid, row, model):
    path = path_of(row["directional_audit"])
    evidence = load(path)
    require(sha(path) == row["directional_audit_sha256"], f"{cid}: instability audit hash mismatch")
    require(evidence["candidate"] == cid and evidence["negative_curvature_and_both_sides_lower"],
            f"{cid}: invalid instability classification")
    require(evidence["source_sha256"] == row["final_model_sha256"], f"{cid}: instability model mismatch")
    require(sha(path_of(evidence["source_path"])) == evidence["source_sha256"],
            f"{cid}: instability source changed")
    for prefix in ("topology", "baseline_energy", "baseline_cold"):
        require(sha(path_of(evidence[f"{prefix}_path"])) == evidence[f"{prefix}_sha256"],
                f"{cid}: instability {prefix} hash changed")
    baseline = load(path_of(evidence["baseline_energy_path"]))["final"]["energy_eV_per_cell"]
    close(baseline, evidence["baseline_energy_accuracy_0p01_eV_cell"], 1e-10,
          f"{cid}: instability baseline energy differs")
    topology = load(path_of(evidence["topology_path"]))
    base_bonds = np.array([e["distance_A"] for e in GENERATOR.audit_geometry(model, topology)["bonds"]])
    energies = []
    side_records = []
    for name in ("minus_1deg", "plus_1deg"):
        side = evidence["perturbations"][name]
        folder = path.parent / name
        raw_result = load(folder / "result.json")
        displaced = read(folder / "input.extxyz")
        structure_check(displaced)
        require(raw_result["status"] == "completed" and side["status"] == "completed",
                f"{cid}: incomplete angular perturbation")
        require(raw_result["input_sha256"] == sha(folder / "input.extxyz")
                == side["input_sha256"], f"{cid}: perturbation input hash mismatch")
        geometry = GENERATOR.audit_geometry(displaced, topology)
        require(geometry["passed"], f"{cid}: invalid angular perturbation geometry")
        close([e["distance_A"] for e in geometry["bonds"]], base_bonds, 2e-8,
              f"{cid}: angular perturbation changed chemical bond lengths")
        energy = raw_result["energy_per_primitive_cell_eV"]
        close(energy, side["energy_per_primitive_cell_eV"], 1e-10,
              f"{cid}: perturbation energy record differs")
        require(energy < baseline - evidence["numerical_energy_margin_eV_cell"],
                f"{cid}: perturbed side did not lower energy")
        energies.append(energy)
        side_records.append(native_final_energy(folder / "native_output.log", (2, 3, 8), energy))
    second = sum(energies) - 2 * baseline
    curvature = second / np.deg2rad(evidence["delta_degrees"])**2
    require(curvature < 0, f"{cid}: tested direction does not have negative curvature")
    close(second, evidence["second_energy_difference_eV_cell"], 1e-9,
          f"{cid}: central second difference mismatch")
    close(curvature, evidence["rigid_rotation_curvature_eV_cell_rad2"], 1e-6,
          f"{cid}: directional curvature mismatch")
    return {"candidate_id": cid, "audit_sha256": sha(path), "source_matches": True,
            "second_energy_difference_eV_cell": second,
            "directional_curvature_eV_cell_rad2": curvature, "native_logs": side_records,
            "scope": "One tested finite angular direction; no full-Hessian claim"}


def group_check(group, rows, models, results):
    gid, representative = group["cluster_id"], group["representative_id"]
    folder = ROOT / "representatives" / gid
    delivered, family = file_family(folder, "representative")
    source = models[representative]
    require(np.array_equal(delivered.numbers, source.numbers), f"{gid}: representative atom order changed")
    close(delivered.positions, source.positions, 1e-10, f"{gid}: representative geometry changed")
    stored_members = load(folder / "members.json")
    require(stored_members == group, f"{gid}: stored membership differs from summary")
    require(representative in group["members"], f"{gid}: representative missing from its own group")
    member_records = []
    for cid in group["members"]:
        row, atom = rows[cid], models[cid]
        require(row["status"] == "preoptimized" and row["geometry_passed"],
                f"{gid}: excluded member was retained: {cid}")
        require(row["cluster_id"] == gid and row["representative_id"] == representative,
                f"{cid}: row/group membership mismatch")
        energy_delta = abs(results[cid]["energy_eV_per_CH2"] -
                           results[representative]["energy_eV_per_CH2"])
        require(energy_delta <= 0.001 + 1e-12, f"{cid}: group energy threshold failed")
        close(energy_delta * 1000, row["dedup_delta_energy_meV_CH2"], 1e-7,
              f"{cid}: recorded dedup energy delta differs")
        if cid == representative:
            mapping, signs, shift = np.arange(len(atom)), np.ones(3), np.zeros(3)
        else:
            mapping = np.asarray(row["equivalent_atom_map"])
            signs = np.asarray(row["equivalence_isometry_signs"], float)
            shift = np.asarray(row["equivalence_origin_shift"], float)
        require(mapping.shape == (12,) and np.array_equal(np.sort(mapping), np.arange(12)),
                f"{cid}: atom map is not a bijection")
        require(np.issubdtype(mapping.dtype, np.integer), f"{cid}: atom map is not integer")
        require(signs.shape == shift.shape == (3,) and set(signs).issubset({-1.0, 1.0})
                and np.isfinite(shift).all(), f"{cid}: invalid isometry/origin shift")
        require(np.array_equal(atom.numbers, source.numbers[mapping]), f"{cid}: map changes elements")
        difference = (atom.get_scaled_positions() * signs + shift
                      - source.get_scaled_positions()[mapping])
        difference -= np.rint(difference)
        distances = np.linalg.norm(difference @ atom.cell.array, axis=1)
        rms, maximum = float(np.sqrt(np.mean(distances**2))), float(distances.max())
        require(rms <= 0.02 + 1e-8 and maximum <= 0.05 + 1e-8, f"{cid}: actual atom map exceeds tolerance")
        close(rms, row["dedup_rms_A"], 1e-8, f"{cid}: recorded RMS differs from actual mapping")
        close(maximum, row["dedup_max_A"], 1e-8, f"{cid}: recorded maximum differs from actual mapping")
        require(results[representative]["energy_eV_per_CH2"]
                <= results[cid]["energy_eV_per_CH2"] + 1e-10, f"{gid}: representative is not lowest energy member")
        member_records.append({"candidate_id": cid, "rms_A": rms, "max_A": maximum,
                               "energy_delta_meV_per_CH2": 1000 * energy_delta})
    return {"group_id": gid, "representative_id": representative, "files": family,
            "member_count": len(member_records), "members": member_records}


def convergence_check(group, comparison, rows, models, results):
    gid, cid = group["cluster_id"], group["representative_id"]
    folder = ROOT / "representative_convergence" / gid
    require(load(folder / "comparison.json") == comparison, f"{gid}: convergence summary differs from comparison file")
    require(comparison["candidate_id"] == cid and comparison["group_id"] == gid,
            f"{gid}: convergence representative identity mismatch")
    representative_path = ROOT / "representatives" / gid / "representative.extxyz"
    require(sha(representative_path) == comparison["representative_sha256"],
            f"{gid}: convergence bound to different representative hash")
    require(path_of(comparison["representative_file"]).resolve() == representative_path.resolve(),
            f"{gid}: convergence representative path mismatch")
    require(sha(path_of(rows[cid]["final_result_path"])) == comparison["source_final_result_sha256"],
            f"{gid}: source fine result hash mismatch")
    case = folder / "3x4x10"
    dense = load(case / "result.json")
    input_path = case / "input.extxyz"
    require(sha(input_path) == comparison["dense_input_sha256"] == dense["input_sha256"],
            f"{gid}: dense input hash mismatch")
    atoms = read(input_path)
    structure_check(atoms)
    require(np.array_equal(atoms.numbers, models[cid].numbers), f"{gid}: dense atom order differs")
    close(atoms.positions, models[cid].positions, 1e-10, f"{gid}: dense input is another geometry")
    require(dense["status"] == "completed" and dense["repeats"] == [3, 4, 10]
            and dense["natoms"] == 1440, f"{gid}: dense run incomplete or wrong sampling")
    require(dense["settings"]["accuracy"] == comparison["comparison_accuracy"] == 0.001,
            f"{gid}: dense/fine cold-check accuracy differs")
    fine = results[cid]
    raw = np.asarray(dense["forces_eV_per_A"], float)
    require(raw.shape == (12, 3) and np.isfinite(raw).all(), f"{gid}: invalid dense forces")
    de = (dense["energy_per_CH2_eV"] - fine["energy_eV_per_CH2"]) * 1000
    df = float(np.linalg.norm(raw - np.asarray(fine["raw_forces_eV_A"]), axis=1).max())
    fmax = float(np.linalg.norm(raw, axis=1).max())
    close(de, comparison["energy_change_meV_per_CH2"], 1e-8, f"{gid}: convergence delta E mismatch")
    close(df, comparison["max_force_vector_change_eV_A"], 1e-10, f"{gid}: convergence delta F mismatch")
    close(fmax, comparison["dense_fmax_eV_A"], 1e-10, f"{gid}: dense Fmax mismatch")
    close(dense["energy_per_primitive_cell_eV"], 4 * dense["energy_per_CH2_eV"], 1e-10,
          f"{gid}: dense energy units mismatch")
    require(abs(de) < 1 and df < 0.005 and fmax < 0.01, f"{gid}: convergence/force thresholds failed")
    require(dense["replica_force_max_deviation_eV_per_A"]
            < comparison["thresholds"]["replica_force_eV_A"],
            f"{gid}: equivalent replicas disagree")
    require(comparison["passed"] and comparison["sampling_change_passed"]
            and comparison["dense_force_target_passed"] and comparison["replica_equivalence_passed"],
            f"{gid}: convergence report did not pass")
    return {"group_id": gid, "candidate_id": cid, "representative_sha256": sha(representative_path),
            "dense_input_sha256": sha(input_path), "same_order_same_coordinates": True,
            "delta_energy_meV_per_CH2": de, "max_delta_force_eV_A": df, "dense_fmax_eV_A": fmax,
            "native_log": native_final_energy(case / "native_output.log", (3, 4, 10),
                                              dense["energy_per_primitive_cell_eV"])}


def xrd_check(summary, rows):
    folder = ROOT / "xrd_benchmark_final"
    if not folder.exists():
        return {"status": "not_present_optional", "checked": False}
    xsummary = load(folder / "summary.json")
    require(xsummary["status"] == "position_only_diagnostic_completed", "Final XRD incomplete")
    require(xsummary["packing_ranking_produced"] is False and xsummary["energy_calculations_run"] == 0,
            "Final XRD incorrectly claims packing ranking/energy calculation")
    require(xsummary["candidate_count"] == len(xsummary["candidates"]), "Final XRD count mismatch")
    known_final = {path_of(row["final_model_path"]).resolve() for row in rows.values()}
    known_representatives = {
        (ROOT / "representatives" / g["cluster_id"] / "representative.extxyz").resolve()
        for g in summary["groups"]}
    covered, reports = set(), []
    for xrow in xsummary["candidates"]:
        path = path_of(xrow["result_path"])
        data = load(path)
        model_path = path_of(data["candidate"]["path"]).resolve()
        require(model_path in known_final | known_representatives, "XRD is bound to a non-delivered model")
        require(sha(model_path) == data["candidate"]["sha256"], "XRD geometry hash mismatch")
        require(data["target_symmetry_check"]["contains_supplied_target_operations"] is True,
                "Final XRD target group closure failed")
        require(data["target_space_group"]["number"] == 62
                and data["target_space_group"]["hall"] == "-P 2c 2n", "Final XRD target group changed")
        close(data["cell_A"], CELL, 1e-10, "Final XRD cell mismatch")
        require(sha(path_of(data["observed_input"]["path"])) == data["observed_input"]["sha256"],
                "Final XRD observation hash mismatch")
        with (path.parent / "observed_nearest_allowed.csv").open(encoding="utf-8-sig", newline="") as stream:
            observed = list(csv.DictReader(stream))
        with (path.parent / "allowed_positions.csv").open(encoding="utf-8-sig", newline="") as stream:
            allowed = list(csv.DictReader(stream))
        require(bool(allowed), "Final XRD allowed-position table is empty")
        allowed_by_id = {peak["position_id"]: peak for peak in allowed}
        require(len(allowed_by_id) == len(allowed), "Final XRD allowed-position IDs repeat")
        allowed_q = np.array([float(peak["q_A_inv"]) for peak in allowed])
        require(len(observed) == data["observed_input"]["count"] == xsummary["observed_count"],
                "Final XRD observed peak count mismatch")
        differences = []
        for peak in observed:
            qobs, qpred = float(peak["q_A_inv"]), float(peak["nearest_allowed_q_A_inv"])
            close(qobs, 2 * math.pi / float(peak["d_obs_A"]), 1e-10, "XRD Qobs/d units mismatch")
            close(qpred, 2 * math.pi / float(peak["nearest_allowed_d_A"]), 1e-10,
                  "XRD Qpred/d units mismatch")
            close(qpred, float(allowed_by_id[peak["allowed_position_id"]]["q_A_inv"]), 1e-10,
                  "XRD selected peak disagrees with allowed-position table")
            difference = qpred - qobs
            close(abs(difference), np.min(np.abs(allowed_q - qobs)), 1e-10,
                  "XRD selected position is not nearest in Q")
            close(difference, float(peak["delta_q_pred_minus_obs_A_inv"]), 1e-10,
                  "XRD signed Q residual mismatch")
            differences.append(difference)
        dq = np.asarray(differences)
        rms = float(np.sqrt(np.mean(dq**2)))
        within = int(np.sum(np.abs(dq) <= xsummary["q_tolerance_A_inv"]))
        close(rms, data["descriptive_statistics"]["rms_delta_q_A_inv"], 1e-10, "XRD RMS mismatch")
        close(rms, xrow["rms_delta_q_A_inv"], 1e-10, "XRD summary RMS mismatch")
        require(within == data["descriptive_statistics"]["observations_within_tolerance"],
                "XRD tolerance count mismatch")
        require(within == xrow["observations_within_tolerance"], "XRD summary tolerance count mismatch")
        require(data["descriptive_statistics"]["eligible_for_packing_ranking"] is False,
                "XRD result incorrectly ranks packings")
        covered.add(model_path)
        reports.append({"model": str(model_path), "sha256": sha(model_path),
                        "observed_count": len(observed), "within_tolerance": within,
                        "rms_delta_Q_A_inv": rms})
    require(known_representatives <= covered or known_final <= covered,
            "Final XRD does not cover all delivered representatives or all final candidates")
    return {"status": "passed", "checked": True, "summary_sha256": sha(folder / "summary.json"),
            "packing_ranking_produced": False, "records": reports}


def main():
    global GENERATOR
    report = {"status": "incomplete", "passed": False,
              "generated_at_utc": datetime.now(timezone.utc).isoformat(),
              "script_sha256": sha(__file__), "quantum_calculations_performed": False,
              "experimental_refined_coordinates_read": False,
              "errors": [], "pending_requirements": [], "candidates": [],
              "directional_exclusions": [], "groups": [], "representative_convergence": []}

    def checked(label, action):
        try:
            return action()
        except Exception as exc:
            report["errors"].append({"check": label, "error": f"{type(exc).__name__}: {exc}"})
            return None

    def inspect():
        manifest = load(ROOT / "generation_manifest.json")
        entries = manifest["unique_survivors"]
        ids = [e["candidate_id"] for e in entries]
        require(len(ids) == len(set(ids)) == manifest["unique_initial_count"], "Manifest unique count mismatch")
        summary_path = ROOT / "optimization_summary.json"
        summary = load(summary_path) if summary_path.exists() else {}
        statuses = {}
        for cid in ids:
            path = ROOT / "candidates" / cid / "relaxation/status.json"
            statuses[cid] = load(path) if path.exists() else {"status": "not_started"}
            if statuses[cid].get("status") != "preoptimized":
                report["pending_requirements"].append(f"{cid}: {statuses[cid].get('status')}")
        if summary.get("status") != "completed" or not summary.get("all_preoptimization_attempts_completed"):
            report["pending_requirements"].append("optimization_summary is not completed")
        report["expected_initial_count"] = len(ids)
        report["preoptimized_status_count"] = sum(s.get("status") == "preoptimized" for s in statuses.values())
        report["optimization_summary_status"] = summary.get("status", "not_present")
        if report["pending_requirements"]:
            return
        report["optimization_summary_sha256"] = sha(summary_path)
        require(summary["generation_manifest_sha256"] == sha(ROOT / "generation_manifest.json"),
                "Summary generation manifest hash mismatch")
        require(summary["matching_helper_sha256"] == sha(HERE / "approx_match_check.py"),
                "Summary matching helper hash mismatch")
        require(summary["program_sha256"] == sha(ROOT / "assess_candidates.py"),
                "Summary assessment program hash mismatch")
        rows = {row["candidate_id"]: row for row in summary["candidates"]}
        require(len(rows) == len(summary["candidates"]) and set(rows) == set(ids), "Summary candidate coverage differs")
        require(summary["preoptimized_count"] == len(ids) and not summary["pending_or_failed"],
                "Completed summary contains unfinished candidates")
        models, results, geometry_flags = {}, {}, {}
        global GENERATOR
        GENERATOR = module("delivery_geometry_reader", ROOT / "generate_candidates.py")
        for entry in entries:
            cid = entry["candidate_id"]
            item = checked(f"candidate:{cid}", lambda e=entry: completed_candidate(e, statuses[e["candidate_id"]],
                                                                                 rows[e["candidate_id"]]))
            if item is not None:
                models[cid], results[cid], item_report = item
                report["candidates"].append(item_report)
                geometry_flags[cid] = item_report["geometry_filters_passed"]
        if set(models) != set(ids):
            return
        unstable = set()
        for entry in entries:
            cid, row = entry["candidate_id"], rows[entry["candidate_id"]]
            if row["status"] == "unstable_stationary":
                evidence = checked(f"direction:{cid}", lambda c=cid, r=row: directional_exclusion(c, r, models[c]))
                if evidence is not None:
                    report["directional_exclusions"].append(evidence)
                    unstable.add(cid)
                require(row.get("excluded_from_retained_pool") is True, f"{cid}: unstable row not excluded")
            elif not geometry_flags[cid]:
                require(row["status"] == "preoptimized_but_geometric_filter_failed",
                        f"{cid}: geometric failure not classified")
            else:
                require(row["status"] == "preoptimized", f"{cid}: unexplained final classification")
            if entry["parameters"]["setting_angle_deg"] == 0:
                require(row["status"] == "unstable_stationary" and cid in unstable,
                        "Known zero-angle directional exclusion is missing or invalid")
        retained = {cid for cid in ids if geometry_flags[cid] and cid not in unstable}
        require(summary["accepted_geometry_count"] == sum(geometry_flags.values()), "Geometric survivor count mismatch")
        require(summary["retained_count"] == len(retained), "Retained count mismatch")
        require(summary["directionally_unstable_count"] == len(unstable), "Instability count mismatch")
        groups = summary["groups"]
        require(summary["geometrically_distinct_final_count"] == len(groups), "Group count mismatch")
        require(len({g["cluster_id"] for g in groups}) == len(groups), "Duplicate group IDs")
        members = [cid for group in groups for cid in group["members"]]
        require(len(members) == len(set(members)) and set(members) == retained,
                "Group membership does not cover every retained candidate exactly once")
        report["retained_candidate_count"] = len(retained)
        report["excluded_geometry_count"] = sum(not passed for passed in geometry_flags.values())
        report["excluded_direction_count"] = len(unstable)
        report["retained_candidate_ids"] = sorted(retained)
        if retained:
            emin = min(results[c]["energy_eV_per_CH2"] for c in retained)
            for cid in retained:
                close(rows[cid]["relative_energy_meV_CH2"],
                      (results[cid]["energy_eV_per_CH2"] - emin) * 1000, 1e-7,
                      f"{cid}: relative energy misreported")
        for group in groups:
            record = checked(f"group:{group['cluster_id']}", lambda g=group: group_check(g, rows, models, results))
            if record is not None:
                report["groups"].append(record)
        convergence_path = ROOT / "representative_convergence/summary.json"
        if not convergence_path.exists():
            report["pending_requirements"].append("representative_convergence/summary.json is missing")
        else:
            convergence = load(convergence_path)
            comparisons = {r["group_id"]: r for r in convergence["results"]}
            require(len(comparisons) == len(convergence["results"])
                    and set(comparisons) == {g["cluster_id"] for g in groups},
                    "Convergence summary group coverage differs from final representatives")
            for group in groups:
                record = checked(f"convergence:{group['cluster_id']}",
                                 lambda g=group: convergence_check(g, comparisons[g["cluster_id"]], rows, models, results))
                if record is not None:
                    report["representative_convergence"].append(record)
            require(convergence["passed"] is True, "Representative convergence summary did not pass")
            report["representative_convergence_summary_sha256"] = sha(convergence_path)
        record = checked("final_xrd", lambda: xrd_check(summary, rows))
        if record is not None:
            report["final_xrd"] = record
        report["structure_families_checked"] = 2 * len(report["candidates"]) + len(report["groups"])
        report["structure_files_checked"] = 3 * report["structure_families_checked"]

    checked("delivery", inspect)
    report["passed"] = not report["errors"] and not report["pending_requirements"]
    report["status"] = ("passed" if report["passed"] else
                        "failed" if report["errors"] else "incomplete")
    report["scope"] = ("Files, arithmetic, geometry and recorded numerical evidence only; "
                       "no phonon/full-Hessian stability claim; dense sampling validated "
                       "only for the exact delivered representatives.")
    OUTPUT.write_text(json.dumps(report, ensure_ascii=False, indent=2,
                                 allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "passed": report["passed"],
                      "errors": report["errors"], "pending_requirements": report["pending_requirements"],
                      "audit_file": str(OUTPUT)}, ensure_ascii=False))
    return 0 if report["passed"] else 1 if report["errors"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
