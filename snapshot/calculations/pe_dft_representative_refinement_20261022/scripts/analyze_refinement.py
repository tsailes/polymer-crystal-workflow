"""Read-only scientific analysis of four fixed-cell PE DFT refinements.

Writes only derived tables and results_summary.md under --root.  No quantum
calculator is imported, attached or evaluated.  Running/failed points are never
ranked.  Complete-link groups require geometry AND energy for every member pair.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import os
from pathlib import Path
import sys
from datetime import datetime, timezone

sys.dont_write_bytecode = True
import numpy as np
from ase.io import read
from ase.units import GPa
from structure_match import compare

DEFAULT_ROOT = Path(__file__).resolve().parents[1]
EXPECTED_MERGE = {"rms_A": 0.02, "max_A": 0.05, "energy_meV_CH2": 1.0}
EXPECTED_SETTINGS = {"method": "PBE-D3(BJ)", "ecut_eV": 1800,
                     "kmesh": [4, 8, 12], "gpts": [112, 72, 44],
                     "xc_stencil": 4, "nbands": 16,
                     "scf_energy_eV_electron": 1e-9,
                     "scf_density_electron_electron": 1e-9,
                     "scf_eigenstates_eV2_electron": 1e-11}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def local_path(value):
    value = str(value)
    if os.name == "nt" and value.startswith("/mnt/") and len(value) > 7:
        value = value[5].upper() + ":/" + value[7:]
    return Path(value)


def load_json(path, warnings=None):
    path = Path(path)
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, ValueError) as exc:
        if warnings is None:
            raise
        warnings.append(f"Could not read {path}: {exc}")
        return None


def jsonl(path, warnings):
    if not path.is_file():
        return []
    rows = []
    for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        try:
            rows.append(json.loads(line))
        except ValueError as exc:
            warnings.append(f"Unreadable JSONL record {path}:{lineno}: {exc}")
    return rows


def write_json(path, value):
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)
                   + "\n", encoding="utf-8")
    tmp.replace(path)


def write_csv(path, rows, required=()):
    fields = list(dict.fromkeys([*required, *(k for r in rows for k in r)]))
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (dict, list))
                             else v for k, v in row.items()})
    tmp.replace(path)


def validate_protocol(protocol):
    for key, expected in EXPECTED_SETTINGS.items():
        if protocol.get(key) != expected:
            raise ValueError(f"Protocol {key} differs from audited settings: {protocol.get(key)!r}")
    for key, expected in EXPECTED_MERGE.items():
        if protocol.get("merge", {}).get(key) != expected:
            raise ValueError(f"Protocol merge.{key} differs from audited criterion")
    if protocol.get("cell_fixed") is not True or protocol.get("N_CH2") != 4:
        raise ValueError("Expected fixed cell with four CH2 units")
    if protocol.get("cell_lengths_A") != [7.432, 4.929, 2.545]:
        raise ValueError("Expected frozen Tashiro cell")
    if protocol.get("optimizer", {}).get("fmax_eV_A") != 0.002:
        raise ValueError("Expected raw force convergence of 0.002 eV/Angstrom")
    if protocol.get("cold_final_verification") is not True:
        raise ValueError("Independent cold-start final check is required")


def minimum_image_delta(a, b, cell):
    delta = (np.asarray(a) - np.asarray(b)) @ np.linalg.inv(cell)
    delta -= np.rint(delta)
    return delta @ cell


def setting_angle(atoms, sample, root):
    topology = local_path(sample.get("topology_path", root.parent /
        "pe_pnam_candidates_20261020" / "candidates" /
        sample["source_candidate_id"] / "reference_topology.json"))
    data = load_json(topology)
    edge = next(e for e in data["periodic_bonds"] if e["kind"] == "CC" and e["i"] == 0)
    # Keep the original bond identities, choosing a periodic image after export.
    # Either of the +/-c/2 C-C images gives the same chain-plane orientation.
    vector = minimum_image_delta(atoms.positions[edge["j"]], atoms.positions[edge["i"]], atoms.cell.array)
    angle = float(np.degrees(np.arctan2(vector[1], vector[0])) % 180)
    return {"setting_angle_deg": min(angle, 180 - angle),
            "oriented_angle_mod180_deg": angle,
            "CC_bond_indices": [edge["i"], edge["j"]],
            "CC_bond_minimum_image_length_A": float(np.linalg.norm(vector)),
            "topology_path": str(topology), "topology_sha256": sha(topology)}


def evaluation_row(evaluation, sample, attempt, source):
    total, dft, dispersion = (evaluation.get(k, {}) for k in ("total", "dft", "dispersion"))
    row = {"sample_id": sample["id"], "source_candidate_id": sample["source_candidate_id"],
           "attempt": attempt.name, "source_file": str(source),
           **{k: evaluation.get(k) for k in ("evaluation_id", "phase", "optimizer_step", "wall_seconds",
               "SCF_converged", "SCF_iterations", "raw_max_force_norm_eV_A",
               "raw_max_stress_component_kbar", "Pnam_operation_max_residual_A", "fixed_cell_max_delta_A")}}
    for label, data, keys in (("total", total, ("energy_extrapolated_eV", "energy_free_eV")),
                             ("dft", dft, ("energy_extrapolated_eV", "energy_free_eV")),
                             ("dispersion", dispersion, ("energy_eV", "free_energy_eV"))):
        for key in keys:
            row[label + "_" + key + "_cell"] = data.get(key)
    if total.get("energy_extrapolated_eV") is not None:
        row["total_energy_extrapolated_eV_CH2"] = total["energy_extrapolated_eV"] / sample["N_CH2"]
    for label, data in (("total", total), ("projected", evaluation.get("projected", {})),
                        ("dft", dft), ("dispersion", dispersion)):
        row[label + "_forces_eV_A"] = data.get("forces_eV_A")
        row[label + "_stress_eV_A3"] = data.get("stress_eV_A3")
    row["positions_A"] = evaluation.get("positions_A")
    row["cell_A"] = evaluation.get("cell_A")
    return row


def endpoint(sample, root, protocol, warnings):
    result = {"sample_id": sample["id"], "source_candidate_id": sample["source_candidate_id"],
              "initial_angle_deg": sample.get("initial_angle_deg"), "accepted": False,
              "analysis_status": "not_started", "reasons": [], "rank": None,
              "relative_energy_meV_CH2": None, "N_CH2": sample["N_CH2"]}
    source = local_path(sample["structure"])
    initial = read(source)
    if sha(source) != sample["structure_sha256"]:
        raise ValueError(f"Frozen input hash changed: {source}")
    result["source_sha256"] = sha(source)
    result["source_structure"] = str(source)
    result["starting_setting_angle"] = setting_angle(initial, sample, root)
    attempts = sorted((root / "runs" / sample["id"]).glob("attempt_*"))
    attempts = [p for p in attempts if p.is_dir()]
    if not attempts:
        return result, None, [], []
    trajectories, scf_rows = [], []
    for attempt in attempts:
        evalfile = attempt / "evaluations.jsonl"
        for item in jsonl(evalfile, warnings):
            trajectories.append(evaluation_row(item, sample, attempt, evalfile))
        for log in ("scf_iterations.jsonl", "final_cold_scf.jsonl"):
            for item in jsonl(attempt / log, warnings):
                scf_rows.append({"sample_id": sample["id"], "attempt": attempt.name,
                                 "source_file": str(attempt / log), **item})
    # A later attempt supersedes earlier endpoints; older data are still exported.
    work = attempts[-1]
    result["attempt"] = work.name
    result["attempt_path"] = str(work)
    state = load_json(work / "result.json", warnings)
    process = load_json(work / "process.json", warnings)
    status = load_json(work / "status.json", warnings) or {}
    result["runner_status"] = (state or status).get("status", "pending")
    result["process_exit_code"] = (process or {}).get("exit_code")
    result["analysis_status"] = "pending_or_running"
    current = [r for r in trajectories if r["attempt"] == work.name]
    if current:
        result["last_completed_evaluation"] = {k: current[-1].get(k) for k in
            ("evaluation_id", "phase", "optimizer_step", "raw_max_force_norm_eV_A", "wall_seconds")}
    if state is None or process is None or process.get("exit_code") is None:
        return result, None, trajectories, scf_rows
    result["analysis_status"] = "not_accepted"
    if state.get("status") != "passed_relaxation" or process.get("exit_code") != 0 or state.get("exit_code") != 0:
        result["reasons"].append("Runner and process must both finish successfully")
        return result, None, trajectories, scf_rows
    try:
        def require(condition, message):
            if not condition:
                raise ValueError(message)
        require(state.get("protocol_sha256") == sha(root / "protocol.json"), "Protocol hash mismatch")
        require(state.get("source_geometry_sha256") == sample["structure_sha256"], "Source hash mismatch")
        require(state.get("final_cold_verification", {}).get("passed") is True, "Cold final check did not pass")
        require(state.get("final_symmetry", {}).get("passed") is True, "Pnam final check did not pass")
        for key, expected in EXPECTED_SETTINGS.items():
            require(state.get("input", {}).get(key) == expected, f"Final input mismatch: {key}")
        ev = state["final_evaluation"]
        require(ev.get("phase") == "final_cold" and ev.get("SCF_converged") is True, "Final evaluation is not converged cold SCF")
        atoms = read(work / "final_structure.extxyz")
        require(sha(work / "final_structure.extxyz") == state.get("final_geometry_sha256"), "Final geometry hash mismatch")
        require(np.array_equal(atoms.numbers, initial.numbers) and all(atoms.pbc), "Composition/order/PBC changed")
        require(len(atoms) == 12 and sum(atoms.numbers == 6) == 4 and sum(atoms.numbers == 1) == 8, "Not C4H8")
        cell = np.diag(protocol["cell_lengths_A"])
        require(np.allclose(atoms.cell.array, cell, atol=protocol["cell_tolerance_A"], rtol=0), "Cell was not fixed")
        require(np.allclose(ev["cell_A"], cell, atol=protocol["cell_tolerance_A"], rtol=0), "Evaluation cell changed")
        require(np.max(np.abs(minimum_image_delta(atoms.positions, ev["positions_A"], cell))) < 1e-7,
                "Final file and final evaluation positions differ")
        total = ev["total"]
        energy = float(total["energy_extrapolated_eV"])
        forces = np.asarray(total["forces_eV_A"], float)
        stress = np.asarray(total["stress_eV_A3"], float)
        require(forces.shape == (12, 3) and stress.shape == (6,), "Property array shape mismatch")
        require(np.isfinite([energy, total["energy_free_eV"]]).all() and np.isfinite(forces).all()
                and np.isfinite(stress).all() and np.isfinite(atoms.positions).all(), "Nonfinite final properties")
        fmax = float(np.linalg.norm(forces, axis=1).max())
        require(fmax <= protocol["optimizer"]["fmax_eV_A"], "Raw final forces exceed tolerance")
        require(abs(energy - ev["dft"]["energy_extrapolated_eV"] - ev["dispersion"]["energy_eV"]) < 1e-9,
                "Energy component sum mismatch")
        require(abs(total["energy_free_eV"] - ev["dft"]["energy_free_eV"] - ev["dispersion"]["free_energy_eV"]) < 1e-9,
                "Free energy component sum mismatch")
        result.update({"accepted": True, "analysis_status": "accepted", "energy_eV_cell": energy,
                       "free_energy_eV_cell": total["energy_free_eV"], "energy_eV_CH2": energy / 4,
                       "dft_energy_eV_cell": ev["dft"]["energy_extrapolated_eV"],
                       "dispersion_energy_eV_cell": ev["dispersion"]["energy_eV"],
                       "raw_fmax_eV_A": fmax, "stress_kbar": (stress / GPa * 10).tolist(),
                       "optimizer_steps": state.get("optimizer_steps"), "wall_seconds": state.get("wall_seconds"),
                       "final_structure": str(work / "final_structure.extxyz"),
                       "final_geometry_sha256": state["final_geometry_sha256"],
                       "final_setting_angle": setting_angle(atoms, sample, root),
                       "cold_verification": state["final_cold_verification"],
                       "result_sha256": sha(work / "result.json"), "process_sha256": sha(work / "process.json")})
        return result, atoms, trajectories, scf_rows
    except (KeyError, ValueError, OSError, TypeError) as exc:
        result["reasons"].append(str(exc))
        return result, None, trajectories, scf_rows


def compare_pair(first_id, second_id, first, second, first_energy, second_energy, criteria):
    match = compare(first, second, rms_tol=criteria["rms_A"], max_tol=criteria["max_A"])
    mapping = match["origin_optimized_match"] or match["origin_optimized_best"]
    delta = abs(first_energy - second_energy) * 1000 / 4
    geometry_pass = match["origin_optimized_passes"]
    return {"first": first_id, "second": second_id, "status": "compared",
            "rms_A": mapping["rms_A"], "max_A": mapping["max_A"],
            "energy_difference_meV_CH2": delta, "geometry_pass": geometry_pass,
            "energy_pass": delta <= criteria["energy_meV_CH2"],
            "merge_pass": bool(geometry_pass and delta <= criteria["energy_meV_CH2"]),
            "selected_mapping_role": "feasible_mapping" if geometry_pass else "best_RMS_mapping",
            "atom_mapping": mapping["atom_mapping"], "signs": mapping["signs"],
            "origin_shift_fractional": mapping["origin_shift_fractional"], "matching_details": match}


def complete_link_partitions(ids, pair_rows):
    """Exhaustive clique partitions for four nodes; never connected components."""
    passed = {frozenset((r["first"], r["second"])) for r in pair_rows if r.get("merge_pass") is True}
    def fits(node, group):
        return all(frozenset((node, other)) in passed for other in group)
    def partitions(left, groups):
        if not left:
            yield [list(g) for g in groups]
            return
        node, *rest = left
        for index, group in enumerate(groups):
            if fits(node, group):
                updated = [list(g) for g in groups]
                updated[index].append(node)
                yield from partitions(rest, updated)
        yield from partitions(rest, groups + [[node]])
    feasible = list(partitions(list(ids), []))
    minimum = min(map(len, feasible))
    best = [p for p in feasible if len(p) == minimum]
    triples = []
    for a, b, c in itertools.combinations(ids, 3):
        edges = [frozenset(p) in passed for p in ((a, b), (a, c), (b, c))]
        if sum(edges) == 2:
            triples.append([a, b, c])
    return {"minimum_complete_link_group_count": minimum,
            "minimum_group_partitions": best, "selected_partition": best[0],
            "partition_ambiguous": len(best) > 1, "nontransitive_triples": triples,
            "all_four_merge": len(ids) == 4 and len(passed) == 6,
            "algorithm": "Enumerate all set partitions; every pair inside every group must pass both bounds. Select fewest groups; retain all tied partitions."}


def analyze(root):
    protocol = load_json(root / "protocol.json")
    validate_protocol(protocol)
    samples = load_json(root / "sample_manifest.json")["samples"]
    if [s["id"] for s in samples] != [f"group_{i:03d}" for i in range(1, 5)]:
        raise ValueError("Expected the four frozen representative sample IDs in order")
    warnings, endpoints, geometries, trajectories, scfs = [], [], {}, [], []
    for sample in samples:
        result, atoms, trajectory, scf = endpoint(sample, root, protocol, warnings)
        endpoints.append(result)
        if atoms is not None:
            geometries[sample["id"]] = atoms
        trajectories.extend(trajectory)
        scfs.extend(scf)
    lookup = {r["sample_id"]: r for r in endpoints}
    pair_rows, mapping_rows = [], []
    for first, second in itertools.combinations(lookup, 2):
        if first not in geometries or second not in geometries:
            pair_rows.append({"first": first, "second": second, "status": "awaiting_accepted_endpoints",
                              "merge_pass": None})
            continue
        pair = compare_pair(first, second, geometries[first], geometries[second],
                            lookup[first]["energy_eV_cell"], lookup[second]["energy_eV_cell"], protocol["merge"])
        pair_rows.append(pair)
        a, b = geometries[first], geometries[second]
        x = a.get_scaled_positions(wrap=False) * pair["signs"] + pair["origin_shift_fractional"]
        y = b.get_scaled_positions(wrap=False)[pair["atom_mapping"]]
        residual = x - y
        residual -= np.rint(residual)
        residual = residual @ a.cell.array
        for i, j in enumerate(pair["atom_mapping"]):
            mapping_rows.append({"first": first, "second": second, "first_atom_index": i,
                                 "second_atom_index": j, "element": a[i].symbol,
                                 "dx_A": float(residual[i, 0]), "dy_A": float(residual[i, 1]),
                                 "dz_A": float(residual[i, 2]), "distance_A": float(np.linalg.norm(residual[i]))})
    complete = len(geometries) == 4
    grouping = None
    if complete:
        ordered = sorted(endpoints, key=lambda r: (r["energy_eV_cell"], r["sample_id"]))
        for rank, row in enumerate(ordered, 1):
            row["rank"] = rank
            row["relative_energy_meV_CH2"] = (row["energy_eV_cell"] - ordered[0]["energy_eV_cell"]) * 1000 / 4
        grouping = complete_link_partitions([r["sample_id"] for r in ordered], pair_rows)
        grouping["selected_groups"] = [{"id": f"dft_group_{i:03d}", "members": members,
            "representative": min(members, key=lambda s: (lookup[s]["energy_eV_cell"], s))}
            for i, members in enumerate(grouping["selected_partition"], 1)]
    summary = {"schema_version": 1, "generated_utc": datetime.now(timezone.utc).isoformat(),
               "status": "complete" if complete else "partial", "complete": complete,
               "accepted_count": len(geometries), "expected_count": 4,
               "method": protocol["method"], "ecut_eV": protocol["ecut_eV"],
               "kmesh": protocol["kmesh"], "gpts": protocol["gpts"],
               "fmax_eV_A": protocol["optimizer"]["fmax_eV_A"],
               "cell_lengths_A": protocol["cell_lengths_A"], "N_CH2": 4,
               "merge_criteria": protocol["merge"], "endpoints": endpoints,
               "pair_comparisons": pair_rows, "grouping": grouping,
               "warnings": warnings, "energy_zero": "Lowest accepted cold-final energy, defined only after all four pass",
               "scope": "PBE-D3(BJ), fixed Tashiro cell and Pnam-constrained atomic relaxation; numerical tolerance-based grouping, not Hessian-verified minima or polymorph identification",
               "quantum_calculations_performed_by_analysis": False,
               "provenance": {"protocol_sha256": sha(root / "protocol.json"),
                   "manifest_sha256": sha(root / "sample_manifest.json"),
                   "analysis_script_sha256": sha(__file__),
                   "structure_match_sha256": sha(Path(__file__).with_name("structure_match.py"))}}
    table = root / "tables"
    table.mkdir(exist_ok=True)
    write_json(table / "refinement_summary.json", summary)
    write_json(table / "pairwise_geometry_energy.json", pair_rows)
    write_csv(table / "final_energy_results.csv", endpoints, ("sample_id", "analysis_status", "accepted", "rank", "relative_energy_meV_CH2"))
    write_csv(table / "all_evaluation_energy_trajectory.csv", trajectories, ("sample_id", "attempt", "evaluation_id", "phase"))
    write_csv(table / "all_scf_iteration_energies.csv", scfs, ("sample_id", "attempt", "evaluation_id", "phase", "SCF_iteration"))
    write_csv(table / "pairwise_geometry_energy.csv", [{k: v for k, v in r.items() if k != "matching_details"} for r in pair_rows], ("first", "second", "status", "merge_pass"))
    write_csv(table / "pairwise_atom_mappings.csv", mapping_rows, ("first", "second", "first_atom_index", "second_atom_index", "distance_A"))
    lines = ["# 四个 PE 代表结构的固定晶胞 DFT 精修", "",
             f"状态：{'全部完成并通过验收' if complete else '计算未完整结束；不作最终排序或合并结论'}。已接受 {len(geometries)}/4。", "",
             "方法：PBE-D3(BJ)，截断能 **1800 eV**，Γ 中心 **k = 4 × 8 × 12**，实空间网格 **112 × 72 × 44**。",
             "固定 Tashiro 晶胞 a/b/c = **7.432/4.929/2.545 Å**；Pnam 对称约束下优化原子位置；原始最大力 ≤ **0.002 eV/Å**，并通过独立冷启动 SCF 复核。应力保留但不作为固定晶胞收敛目标。", "",
             "|代表|来源候选|验收状态|最终 E (eV/晶胞)|ΔE (meV/CH₂)|最终 setting angle (°)|", "|---|---|---|---:|---:|---:|"]
    for row in endpoints:
        e = f"{row['energy_eV_cell']:.10f}" if row.get("accepted") else "—"
        d = f"{row['relative_energy_meV_CH2']:.6f}" if complete else "—"
        angle = f"{row['final_setting_angle']['setting_angle_deg']:.6f}" if row.get("accepted") else "—"
        lines.append(f"|{row['sample_id']}|{row['source_candidate_id']}|{row['analysis_status']}|{e}|{d}|{angle}|")
    lines += ["", "合并要求同一原子映射同时满足 RMS ≤ **0.02 Å**、最大位移 ≤ **0.05 Å**，并且能差 ≤ **1.0 meV/CH₂**。使用周期边界、全局原点平移、正交晶胞反射与同元素原子置换；完整比较六对。每组内部所有配对都须通过，不能仅凭近能或传递连通合并。", ""]
    if complete:
        lines.append(f"最少完整相容分组数：**{grouping['minimum_complete_link_group_count']}**。")
        for group in grouping["selected_groups"]:
            lines.append(f"- {group['id']}：{', '.join(group['members'])}；代表 {group['representative']}。")
        if grouping["partition_ambiguous"] or grouping["nontransitive_triples"]:
            lines.append("\n存在非传递/分组不唯一的容差匹配；JSON 保留所有并列最少分组和相关三元组，不能解释为唯一结构类别。")
        if grouping["all_four_merge"]:
            lines.append("\n四个终态的全部六对均满足几何与能量条件，可在本轮容差内合并为一组。")
    for row in endpoints:
        if row["reasons"]:
            lines.append(f"\n{row['sample_id']} 未接受原因：{'；'.join(row['reasons'])}")
    if warnings:
        lines += ["", "读取提示（完整保存在 JSON）：", *[f"- {w}" for w in warnings]]
    lines += ["", "输出：`tables/final_energy_results.csv`（终态绝对能量/组分/力/应力/角度）；`all_evaluation_energy_trajectory.csv`（所有尝试的每次几何评估、总能量/DFT/色散分量及原子力应力）；`all_scf_iteration_energies.csv`（SCF 每步能量）；`pairwise_geometry_energy.csv/.json` 与 `pairwise_atom_mappings.csv`（六对比较及周期映射）；`refinement_summary.json`（完整可追溯结果）。原始运行目录不改写。", "",
              "范围：这是固定晶胞、Pnam 约束下的局部精修与容差去重。没有进行声子/Hessian 验证，不据此声称全局最低能、已证实动力学稳定或独立多晶型；旧 16 初态 σ 不用于本轮终态淘汰。", ""]
    tmp = root / "results_summary.md.tmp"
    tmp.write_text("\n".join(lines), encoding="utf-8")
    tmp.replace(root / "results_summary.md")
    return summary


def self_test(root):
    samples = load_json(root / "sample_manifest.json")["samples"]
    atoms = [read(local_path(s["structure"])) for s in samples]
    output = {"quantum_calculations_performed": False, "self_matches": [], "original_four_pairs": []}
    for sample, a in zip(samples, atoms):
        m = compare(a, a)
        assert m["origin_optimized_passes"] and m["origin_optimized_best"]["max_A"] < 1e-10
        b = a.copy()
        b.set_scaled_positions(a.get_scaled_positions(wrap=False) * [-1, 1, -1] + [0.137, 0.271, 0.419])
        b = b[list(reversed(range(len(b))))]
        m = compare(a, b)
        assert m["origin_optimized_passes"] and m["origin_optimized_best"]["max_A"] < 1e-10
        output["self_matches"].append({"id": sample["id"], "reflection_translation_permutation_max_A": m["origin_optimized_best"]["max_A"]})
    for (i, a), (j, b) in itertools.combinations(enumerate(atoms), 2):
        m = compare_pair(samples[i]["id"], samples[j]["id"], a, b, 0.0, 0.0, EXPECTED_MERGE)
        output["original_four_pairs"].append({k: m[k] for k in ("first", "second", "rms_A", "max_A", "geometry_pass", "merge_pass")})
    distorted = atoms[0].copy()
    distorted.positions[4] += [0.4, 0.0, 0.0]
    equal_energy = compare_pair("a", "b", atoms[0], distorted, -100, -100, EXPECTED_MERGE)
    assert equal_energy["energy_pass"] and not equal_energy["merge_pass"]
    high_energy = compare_pair("a", "b", atoms[0], atoms[0], -100, -99.992, EXPECTED_MERGE)
    assert high_energy["geometry_pass"] and not high_energy["merge_pass"]
    ids = ["a", "b", "c", "d"]
    pairs = [{"first": a, "second": b, "merge_pass": (a, b) in [("a", "b"), ("b", "c")]} for a, b in itertools.combinations(ids, 2)]
    grouping = complete_link_partitions(ids, pairs)
    assert grouping["minimum_complete_link_group_count"] == 3 and grouping["partition_ambiguous"]
    assert not grouping["all_four_merge"] and grouping["nontransitive_triples"] == [["a", "b", "c"]]
    for p in pairs:
        p["merge_pass"] = True
    assert complete_link_partitions(ids, pairs)["all_four_merge"]
    output.update({"equal_energy_distorted_geometry_rejected": True,
                   "identical_geometry_large_energy_difference_rejected": True,
                   "nontransitive_complete_link_test": grouping, "status": "passed"})
    # Exercise endpoint acceptance entirely in memory; never fabricate an on-disk
    # result in the scientific run directory.
    from unittest.mock import patch
    sample, a = samples[0], atoms[0]
    protocol = load_json(root / "protocol.json")
    work = root / "runs" / sample["id"] / "attempt_01"
    fixture = {"status": "passed_relaxation", "exit_code": 0,
        "protocol_sha256": sha(root / "protocol.json"),
        "source_geometry_sha256": sample["structure_sha256"],
        "final_geometry_sha256": "synthetic-in-memory-only",
        "final_cold_verification": {"passed": True}, "final_symmetry": {"passed": True},
        "input": dict(EXPECTED_SETTINGS),
        "final_evaluation": {"phase": "final_cold", "SCF_converged": True,
            "cell_A": a.cell.array.tolist(), "positions_A": a.positions.tolist(),
            "total": {"energy_extrapolated_eV": -100.0, "energy_free_eV": -100.0,
                      "forces_eV_A": np.zeros((12, 3)).tolist(), "stress_eV_A3": [0.] * 6},
            "dft": {"energy_extrapolated_eV": -99.0, "energy_free_eV": -99.0},
            "dispersion": {"energy_eV": -1.0, "free_energy_eV": -1.0}}}
    process = {"exit_code": 0}
    real_load, real_sha, real_read, real_glob, real_is_dir = load_json, sha, read, Path.glob, Path.is_dir
    def fake_load(path, warnings=None):
        path = Path(path)
        if path == work / "result.json":
            return fixture
        if path == work / "process.json":
            return process
        if path == work / "status.json":
            return {"status": fixture["status"]}
        return real_load(path, warnings)
    def fake_sha(path):
        if Path(path) in (work / "final_structure.extxyz", work / "result.json", work / "process.json"):
            return "synthetic-in-memory-only"
        return real_sha(path)
    def fake_read(path, *args, **kwargs):
        if Path(path) == work / "final_structure.extxyz":
            return a.copy()
        return real_read(path, *args, **kwargs)
    module = sys.modules[__name__]
    with patch.object(module, "load_json", fake_load), patch.object(module, "sha", fake_sha), \
         patch.object(module, "read", fake_read), patch.object(module, "jsonl", lambda *args: []), \
         patch.object(Path, "glob", lambda p, pat: iter([work]) if p == work.parent else real_glob(p, pat)), \
         patch.object(Path, "is_dir", lambda p: True if p == work else real_is_dir(p)):
        result, _, _, _ = endpoint(sample, root, protocol, [])
        assert result["accepted"], result["reasons"]
        process["exit_code"] = 1
        assert not endpoint(sample, root, protocol, [])[0]["accepted"]
        process["exit_code"] = 0
        fixture["final_evaluation"]["phase"] = "initial"
        assert not endpoint(sample, root, protocol, [])[0]["accepted"]
        fixture["final_evaluation"]["phase"] = "final_cold"
        fixture["final_evaluation"]["total"]["forces_eV_A"][0][0] = 0.003
        assert not endpoint(sample, root, protocol, [])[0]["accepted"]
    output["in_memory_endpoint_acceptance_tests"] = {
        "valid_cold_final_accepted": True, "nonzero_process_exit_rejected": True,
        "noncold_final_rejected": True, "excess_final_force_rejected": True,
        "synthetic_scientific_result_files_written": False}
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps(self_test(args.root.resolve()), indent=2, allow_nan=False))
    else:
        summary = analyze(args.root.resolve())
        print(json.dumps({"status": summary["status"], "accepted_count": summary["accepted_count"],
                          "expected_count": summary["expected_count"], "grouping": summary["grouping"],
                          "summary_path": str(args.root.resolve() / "tables/refinement_summary.json")}))


if __name__ == "__main__":
    main()
