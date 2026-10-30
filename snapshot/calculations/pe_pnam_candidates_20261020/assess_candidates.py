"""Collect every preoptimized candidate; confirm geometry before energy dedup."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read

from generate_candidates import ROOT, audit_geometry, dump, export_structure, sha, symmetry_fingerprint
from audits.approx_match_check import compare as compare_periodic


def local_path(value):
    return Path(value)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()
    manifest = json.loads((ROOT/"generation_manifest.json").read_text())
    direction_audits = {}
    for audit_file in sorted((ROOT/"audits/zero_angle_curvature").glob("*/summary.json")):
        evidence = json.loads(audit_file.read_text())
        if evidence.get("negative_curvature_and_both_sides_lower"):
            direction_audits[evidence["candidate"]] = (audit_file, evidence)
    rows, finished, pending = [], [], []
    for entry in manifest["unique_survivors"]:
        cid = entry["candidate_id"]
        status_file = ROOT/"candidates"/cid/"relaxation/status.json"
        status = json.loads(status_file.read_text()) if status_file.exists() else {"status": "not_started"}
        row = dict(candidate_id=cid, initial_angle_deg=entry["parameters"]["setting_angle_deg"],
                   status=status["status"], stage=status.get("stage"), step=status.get("step"),
                   preoptimization_completed=status["status"] == "preoptimized")
        if status["status"] != "preoptimized":
            pending.append(cid)
            rows.append(row)
            continue
        result_path = local_path(status["result"])
        result = json.loads(result_path.read_text())
        if result.get("passed") is not True or result.get("status") != "preoptimized":
            raise ValueError(f"{cid}: status file contradicts final validation")
        model_path = result_path.parent/"final.extxyz"
        atoms = read(model_path)
        reference_atoms = read(ROOT/entry["input_path"])
        if not all(atoms.pbc) or not np.allclose(atoms.cell, reference_atoms.cell, rtol=0, atol=1e-12):
            raise ValueError(f"{cid}: cell or PBC changed")
        force = np.asarray(result["raw_forces_eV_A"])
        if force.shape != (12, 3) or not np.isfinite(force).all() or not np.isfinite(result["energy_eV_per_CH2"]):
            raise ValueError(f"{cid}: invalid energy or forces")
        if abs(np.linalg.norm(force, axis=1).max()-result["raw_fmax_eV_A"]) > 1e-10:
            raise ValueError(f"{cid}: inconsistent force norm")
        if np.linalg.norm(force, axis=1).max() >= .01:
            raise ValueError(f"{cid}: force threshold not met")
        if abs(atoms.info["energy_eV_per_CH2"]-result["energy_eV_per_CH2"]) > 1e-10:
            raise ValueError(f"{cid}: final geometry metadata contradicts energy record")
        topology = json.loads((ROOT/"candidates"/cid/"reference_topology.json").read_text())
        audit = audit_geometry(atoms, topology)
        row.update(energy_eV_per_CH2=result["energy_eV_per_CH2"],
                   max_raw_force_eV_A=result["raw_fmax_eV_A"],
                   geometry_passed=audit["passed"], final_model_path=str(model_path),
                   final_model_sha256=sha(model_path),
                   final_result_path=str(result_path), final_result_sha256=sha(result_path),
                   final_min_interchain_CC_A=audit["minimum_interchain_CC"]["distance_A"],
                   final_CC_min_A=audit["CC_range_A"][0], final_CC_max_A=audit["CC_range_A"][1],
                   final_CC_strict_1p53_1p54=audit["CC_strict_1p53_1p54_passed"],
                   final_detected_SG=audit["symmetry"]["detected_number"])
        # Same-carbon sequence around C at z=1/4, use immutable topology vectors.
        cc_edge = next(e for e in topology["periodic_bonds"] if e["kind"] == "CC" and e["i"] == 0)
        v = atoms.positions[cc_edge["j"]] + np.asarray(cc_edge["shift"])@atoms.cell.array-atoms.positions[0]
        angle = float(np.degrees(np.arctan2(v[1], v[0])) % 180)
        row["final_chain_plane_angle_deg"] = min(angle, 180-angle)
        rows.append(row)
        if not audit["passed"]:
            row["status"] = "preoptimized_but_geometric_filter_failed"
            continue
        if cid in direction_audits:
            audit_file, evidence = direction_audits[cid]
            if evidence["source_sha256"] != sha(model_path):
                raise ValueError(f"{cid}: directional-instability audit belongs to another model")
            row.update(status="unstable_stationary", excluded_from_retained_pool=True,
                       directional_audit=str(audit_file),
                       directional_audit_sha256=sha(audit_file),
                       stability_note="Both +/-1 degree Pnam-preserving rotations lower energy; full Hessian not evaluated")
            continue
        finished.append(dict(row=row, atoms=atoms, result=result))
    finished.sort(key=lambda x: x["row"]["energy_eV_per_CH2"])
    groups = []
    for entry in finished:
        row, atoms = entry["row"], entry["atoms"]
        fingerprint = symmetry_fingerprint(atoms)
        duplicate = None
        for group in groups:
            delta_e = abs(row["energy_eV_per_CH2"]-group["energy_eV_per_CH2"])
            if delta_e > .001:
                continue
            other_fp = group["fingerprint"]
            if (fingerprint["composition"] != other_fp["composition"] or
                    not fingerprint["target_group_contained"] or not other_fp["target_group_contained"]):
                continue
            # Necessary descriptor condition for a maximum atom displacement of
            # 0.05 A: each matched pair distance can change by at most 0.10 A.
            if any(np.max(np.abs(np.array(values)-other_fp["pair_distributions_A"][kind])) > .10000001
                   for kind, values in fingerprint["pair_distributions_A"].items()):
                continue
            match = compare_periodic(atoms, group["atoms"], rms_tol=.02, max_tol=.05)["origin_optimized_match"]
            if match:
                duplicate = (group, match, delta_e)
                break
        if duplicate:
            group, match, delta_e = duplicate
            row.update(cluster_id=group["cluster_id"], representative_id=group["representative_id"],
                       dedup_rms_A=match["rms_A"], dedup_max_A=match["max_A"],
                       dedup_delta_energy_meV_CH2=1000*delta_e,
                       equivalent_atom_map=match["atom_mapping"],
                       equivalence_origin_shift=match["origin_shift_fractional"],
                       equivalence_isometry_signs=match["signs"])
            group["members"].append(row["candidate_id"])
        else:
            gid = f"group_{len(groups)+1:03d}"
            row.update(cluster_id=gid, representative_id=row["candidate_id"],
                       dedup_rms_A=0, dedup_max_A=0, dedup_delta_energy_meV_CH2=0)
            groups.append(dict(cluster_id=gid, representative_id=row["candidate_id"],
                               energy_eV_per_CH2=row["energy_eV_per_CH2"],
                               members=[row["candidate_id"]], atoms=atoms,
                               fingerprint=fingerprint))
    if finished:
        emin = finished[0]["row"]["energy_eV_per_CH2"]
        for item in finished:
            item["row"]["relative_energy_meV_CH2"] = (item["row"]["energy_eV_per_CH2"]-emin)*1000
    for group in groups:
        folder = ROOT/("representatives_preview" if pending else "representatives")/group["cluster_id"]
        export_structure(group["atoms"], folder, stem="representative")
        dump(folder/"symmetry_fingerprint.json", group["fingerprint"])
        dump(folder/"members.json", {k: v for k, v in group.items() if k not in ("atoms", "fingerprint")})
    summary = dict(
        status="completed" if not pending else "in_progress",
        all_preoptimization_attempts_completed=not pending,
        generation_counts={k: manifest[k] for k in ("enumerated_count", "physical_survivor_count",
                           "unique_initial_count", "rejected_count", "exact_duplicates_count")},
        preoptimized_count=sum(r["preoptimization_completed"] for r in rows),
        accepted_geometry_count=sum(bool(r.get("geometry_passed")) for r in rows),
        retained_count=len(finished),
        directionally_unstable_count=sum(r["status"] == "unstable_stationary" for r in rows),
        pending_or_failed=pending,
        geometrically_distinct_final_count=len(groups),
        deduplication_criteria=dict(energy_tolerance_meV_per_CH2=1, periodic_rms_tolerance_A=.02,
                                   maximum_displacement_tolerance_A=.05,
                                   requires="Target-symmetry fingerprint AND explicit element-preserving atom mapping",
                                   lattice_isometries="All eight sign inversions of unequal orthorhombic axes",
                                   origins="Every carbon anchor, followed by alternating least-squares global-translation and element-preserving atom-assignment refinement",
                                   matching_limitation="Multistart local refinement, not a general global matching proof",
                                   note="Same energy or same spacegroup alone never establishes equivalence"),
        groups=[{k:v for k,v in g.items() if k not in ("atoms", "fingerprint")} for g in groups],
        candidates=rows,
        scientific_scope="Stationary preoptimized structures within fixed-cell Pnam manifold; not phonon-verified minima",
        sampling_status="Fine 2x3x8; denser-grid representative checks recorded separately",
        program_sha256=sha(__file__),
        matching_helper_sha256=sha(ROOT/"audits/approx_match_check.py"),
        generation_manifest_sha256=sha(ROOT/"generation_manifest.json"))
    dump(ROOT/"optimization_summary.json", summary)
    fields = ["candidate_id", "initial_angle_deg", "status", "stage", "step", "energy_eV_per_CH2",
              "relative_energy_meV_CH2", "max_raw_force_eV_A", "geometry_passed",
              "final_chain_plane_angle_deg", "final_min_interchain_CC_A",
              "final_CC_min_A", "final_CC_max_A", "final_CC_strict_1p53_1p54", "final_detected_SG",
              "cluster_id", "representative_id", "dedup_rms_A", "dedup_max_A",
              "dedup_delta_energy_meV_CH2", "final_model_path"]
    with (ROOT/"optimization_catalog.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({k:summary[k] for k in ("status", "preoptimized_count",
                     "accepted_geometry_count", "retained_count", "directionally_unstable_count",
                     "geometrically_distinct_final_count", "pending_or_failed")}))
    if args.require_complete and pending:
        raise SystemExit("Incomplete candidates: "+", ".join(pending))


if __name__ == "__main__":
    main()
