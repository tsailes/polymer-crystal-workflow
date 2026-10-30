"""Post-hoc threshold sensitivity only; no quantum calculations or production edits."""
from pathlib import Path
import csv
import hashlib
import json
import sys
import time

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from ase.io import read
import numpy as np
from generate_candidates import symmetry_fingerprint
from audits.approx_match_check import compare


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    started = time.perf_counter()
    source = ROOT / "optimization_summary.json"
    source_hash = sha(source)
    data = json.loads(source.read_text())
    original = [r for r in data["candidates"] if r.get("cluster_id")]
    rows = sorted(original, key=lambda r: r["energy_eV_per_CH2"])
    atoms, fingerprints, provenance = {}, {}, []
    for row in rows:
        ident = row["candidate_id"]
        model, result_path = Path(row["final_model_path"]), Path(row["final_result_path"])
        assert sha(model) == row["final_model_sha256"]
        assert sha(result_path) == row["final_result_sha256"]
        result = json.loads(result_path.read_text())
        assert result["energy_eV_per_CH2"] == row["energy_eV_per_CH2"]
        atoms[ident] = read(model)
        fingerprints[ident] = symmetry_fingerprint(atoms[ident])
        provenance.append({"candidate_id": ident, "structure_path": str(model),
            "structure_sha256": sha(model), "result_path": str(result_path),
            "result_sha256": sha(result_path), "energy_eV_CH2": row["energy_eV_per_CH2"],
            "initial_angle_deg": row["initial_angle_deg"],
            "final_angle_deg": row["final_chain_plane_angle_deg"]})
    assert len(rows) == 15
    cache = {}

    def grouping(rms, maximum, energy):
        groups, comparisons = [], []
        for row in rows:
            ident, fp = row["candidate_id"], fingerprints[row["candidate_id"]]
            found = False
            for group in groups:
                rep = group["representative"]
                delta = abs(row["energy_eV_per_CH2"] - group["energy_eV_CH2"]) * 1000
                if delta > energy:
                    continue
                other = fingerprints[rep]
                if fp["composition"] != other["composition"] or not fp["target_group_contained"] or not other["target_group_contained"]:
                    continue
                # Same necessary descriptor filter as the production algorithm,
                # generalized to 2*max displacement for each sensitivity case.
                if any(np.max(np.abs(np.array(v) - other["pair_distributions_A"][kind])) > 2 * maximum + 1e-8
                       for kind, v in fp["pair_distributions_A"].items()):
                    continue
                key = (ident, rep, rms, maximum)
                if key not in cache:
                    cache[key] = compare(atoms[ident], atoms[rep], rms_tol=rms, max_tol=maximum)
                match = cache[key]["origin_optimized_match"]
                best = match or cache[key]["origin_optimized_best"]
                comparisons.append({"candidate_id": ident, "representative": rep,
                    "delta_energy_meV_CH2": delta, "geometry_pass": match is not None,
                    "rms_A": best["rms_A"], "max_A": best["max_A"],
                    "atom_mapping": best["atom_mapping"], "signs": best["signs"],
                    "origin_shift_fractional": best["origin_shift_fractional"]})
                if match:
                    group["members"].append(ident)
                    found = True
                    break
            if not found:
                groups.append({"representative": ident, "energy_eV_CH2": row["energy_eV_per_CH2"],
                               "members": [ident]})
        return {"rms_A": rms, "max_A": maximum, "energy_meV_CH2": energy,
                "group_count": len(groups), "groups": groups, "comparisons": comparisons}

    geometry = [grouping(rms, maximum, 1.0) for rms, maximum in
                ((.01, .025), (.02, .05), (.03, .075), (.04, .10))]
    energy = [grouping(.02, .05, threshold) for threshold in (.1, .5, 1.0, 2.0)]
    expected = [g["members"] for g in data["groups"]]
    assert [g["members"] for g in geometry[1]["groups"]] == expected
    assert geometry[1]["group_count"] == 4
    assert sha(source) == source_hash
    for row in provenance:
        assert sha(row["structure_path"]) == row["structure_sha256"]
        assert sha(row["result_path"]) == row["result_sha256"]
    report = {"status": "posthoc_sensitivity_not_production", "quantum_calculations": False,
        "production_inputs_thresholds_and_results_modified": False,
        "purpose": "Describe sensitivity of prior representative-based grouping; not validate a physical tolerance or choose a new threshold",
        "algorithm": "Same energy-sorted greedy representative matching as assess_candidates.py; every member compared with representative, not complete linkage",
        "descriptor_filter": "Necessary pair-distance prefilter uses 2*maximum_displacement (the production .10 A filter at max=.05 A)",
        "geometries_and_energies": "Existing GFN2 optimized endpoints, 300 K electronic occupations; no new optimization and not the subsequent 0 K singlepoint calibration",
        "excluded_candidate": "cand_0001, original directional-instability exclusion retained",
        "default_4_groups_and_exact_members_reproduced": True,
        "source_summary_path": str(source), "source_summary_sha256": source_hash,
        "analysis_script_sha256": sha(__file__),
        "production_grouping_script_sha256": sha(ROOT / "assess_candidates.py"),
        "matching_helper_sha256": sha(ROOT / "audits/approx_match_check.py"),
        "source_records": provenance, "geometry_sensitivity": geometry,
        "energy_sensitivity": energy, "wall_seconds": time.perf_counter() - started}
    (HERE / "results.json").write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    with (HERE / "group_counts.csv").open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.writer(stream)
        writer.writerow(["axis", "rms_A", "max_A", "energy_meV_CH2", "group_count", "members"])
        for axis, cases in (("geometry", geometry), ("energy", energy)):
            for case in cases:
                writer.writerow([axis, case["rms_A"], case["max_A"], case["energy_meV_CH2"],
                    case["group_count"], json.dumps([g["members"] for g in case["groups"]])])
    lines = ["# 事后去重阈值敏感性检查（非正式生产判据）", "",
        "只读取已有 15 个 GFN2 终态坐标与能量；保持原先能量排序、成员与代表匹配的算法。不启动量化、不改正式协议或结果。", "",
        "|RMS Å|最大位移 Å|能差 meV/CH₂|组数|", "|---:|---:|---:|---:|"]
    for case in geometry + energy:
        lines.append(f"|{case['rms_A']}|{case['max_A']}|{case['energy_meV_CH2']}|{case['group_count']}|")
    lines += ["", "默认 0.02 Å / 0.05 Å / 1 meV/CH₂ 精确复现原来的 4 组及全部成员。详细成员、映射、原始坐标与能量文件哈希见 results.json。", "",
        "这项事后检查只说明分组对所选容差是否敏感，不提供这三个数值的文献依据，也不证明应采用哪组新阈值。这里使用的是原 GFN2 预优化末态（300 K 电子占据），不是随后独立计算的 0 K 标定单点。", ""]
    (HERE / "README.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({"default_reproduced": True, "wall_seconds": report["wall_seconds"],
        "geometry": [{k: c[k] for k in ("rms_A", "max_A", "energy_meV_CH2", "group_count", "groups")} for c in geometry],
        "energy": [{k: c[k] for k in ("rms_A", "max_A", "energy_meV_CH2", "group_count")} for c in energy]}))


if __name__ == "__main__":
    main()
