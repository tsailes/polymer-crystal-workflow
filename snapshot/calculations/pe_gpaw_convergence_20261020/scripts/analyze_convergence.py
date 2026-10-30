"""Audit fixed-geometry convergence; reads archived results, never runs GPAW.

Only explicitly allowed stages are considered.  Unknown input fields remain
part of the group key.  A group needs three distinct, monotonically refined
settings and two stable terminal intervals before its endpoint may serve as a
provisional reference.  The endpoint itself is never recommended.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import math
from pathlib import Path
import re
import tempfile

STAGES = {
    "numerical_cutoff": ["ecut_eV"],
    "numerical_kmesh": ["kmesh"],
    "numerical_grid": ["gpts", "h_A"],
    "numerical_joint": ["ecut_eV", "kmesh", "gpts", "h_A"],
    "functional_transfer": ["ecut_eV", "kmesh", "gpts", "h_A"],
    "geometry_transfer": ["ecut_eV", "kmesh", "gpts", "h_A"],
}
SCAN_FIELDS = {"ecut_eV", "kmesh", "gpts", "h_A"}
METADATA = {
    "id", "stage", "structure", "checkpoint", "protocol_sha256", "queue_sha256",
    "runner_sha256_at_submission", "scan_parameters", "scan_label", "scan_index",
    "description", "notes",
}
LIMITS = {"energy_meV_atom": 0.5, "force_vector_eV_A": 0.005,
          "stress_component_kbar": 0.1}
VOIGT = ["xx", "yy", "zz", "yz", "xz", "xy"]
KBAR_PER_EV_A3 = 1602.1766208  # ASE CODATA-2014 conversion used by run_point.py.


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False, allow_nan=False)


def fingerprint(obj):
    return hashlib.sha256(canonical(obj).encode()).hexdigest()


def file_sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def numeric(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name}: expected finite numeric value")
    if not math.isfinite(value):
        raise ValueError(f"{name}: nonfinite value")
    return float(value)


def vector(value, count, name):
    if not isinstance(value, list) or len(value) != count:
        raise ValueError(f"{name}: expected length {count}")
    return [numeric(v, name) for v in value]


def matrix(value, count, name):
    if not isinstance(value, list) or len(value) != count:
        raise ValueError(f"{name}: expected {count} rows")
    return [vector(row, 3, name) for row in value]


def positive_triplet(value, name):
    values = vector(value, 3, name)
    if any(v <= 0 or int(v) != v for v in values):
        raise ValueError(f"{name}: expected positive integers")
    return [int(v) for v in values]


def source_identity(result):
    datasets = result.get("datasets")
    if not isinstance(datasets, list) or len(datasets) != result["natoms"]:
        raise ValueError("missing per-atom PAW dataset identities")
    identity = []
    for item in datasets:
        if not item.get("file_sha256") or not item.get("symbol"):
            raise ValueError("PAW dataset SHA256 or atom identity missing")
        identity.append({k: item.get(k) for k in
                         ("symbol", "type", "setupname", "fingerprint", "file_sha256")})
    return identity


def extract_record(directory, spec, result):
    """Reject unsuitable results before deriving a strict physics/config key."""
    if result.get("status") != "passed_single_point" or result.get("SCF_converged") is not True:
        raise ValueError("result is not a successful, SCF-converged single point")
    if result.get("input") != spec:
        raise ValueError("result.input differs from archived input.json")
    if result.get("input_sha256") != file_sha(directory / "input.json"):
        raise ValueError("archived input SHA256 mismatch")
    runner = result.get("runner_sha256")
    if not runner or runner != spec.get("runner_sha256_at_submission"):
        raise ValueError("runner provenance is missing or inconsistent")
    snapshot = directory / "run_point_snapshot.py"
    if snapshot.exists() and file_sha(snapshot) != runner:
        raise ValueError("runner snapshot SHA256 mismatch")
    units = result.get("units", {})
    if (units.get("energy") != "eV/cell" or units.get("forces") != "eV/Angstrom"
            or units.get("stress_voigt_order") != VOIGT
            or not units.get("stress", "").startswith("eV/Angstrom^3")):
        raise ValueError("unknown or incomplete physical units")
    n = result.get("natoms")
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise ValueError("invalid atom count")
    total = result["total"]
    energy = numeric(total["energy_free_eV"], "free energy")
    forces = matrix(total["forces_eV_A"], n, "forces")
    stress = vector(total["stress_kbar"], 6, "stress kbar")
    stress_native = vector(total["stress_eV_A3"], 6, "stress eV/A3")
    for x, y in zip(stress, stress_native):
        if not math.isclose(x, y * KBAR_PER_EV_A3, rel_tol=2e-8, abs_tol=1e-7):
            raise ValueError("stress unit conversion is inconsistent")
    datasets = source_identity(result)
    geometry = {
        "source_geometry_sha256": result.get("source_geometry_sha256"),
        "natoms": n, "symbols": [x["symbol"] for x in datasets],
        "cell_A": matrix(result["cell_A"], 3, "cell"),
        "positions_A": matrix(result["positions_A"], n, "positions"),
    }
    if not geometry["source_geometry_sha256"]:
        raise ValueError("source geometry SHA256 missing")
    method = spec["method"]
    dispersion = result.get("dispersion", {})
    dispersion_parameters = dispersion.get("calculator_parameters")
    if method in {"PBE-D3(BJ)", "r2SCAN-D4"}:
        if not isinstance(dispersion_parameters, dict):
            raise ValueError("dispersion calculator settings missing")
        cutoffs = dispersion_parameters.get("realspace_cutoff", {})
        if any(numeric(cutoffs.get(k, 0), k) <= 0 for k in ("width2", "width3")):
            raise ValueError("hard-cutoff dispersion excluded from production convergence")
    scan = spec.get("scan_parameters", STAGES[spec["stage"]])
    if not isinstance(scan, list) or not scan or len(scan) != len(set(scan)):
        raise ValueError("scan_parameters must be a nonempty unique list")
    if not set(scan) <= SCAN_FIELDS:
        raise ValueError("unsupported scan parameter")
    if spec["stage"] in {"numerical_cutoff", "numerical_kmesh", "numerical_grid"}:
        if not set(scan) <= set(STAGES[spec["stage"]]):
            raise ValueError("scan_parameters conflicts with stage")
    scan = sorted(scan)
    settings = {k: v for k, v in spec.items() if k not in METADATA and k not in scan}
    # GPAW input defaults are protected by the exact runner hash.  Unknown
    # input fields are deliberately retained, so unrelated conditions split.
    identity = {
        "stage": spec["stage"], "method": method,
        "geometry": geometry, "fixed_input_parameters": settings,
        "scan_parameters": scan, "runner_sha256": runner,
        "software": result.get("software", {}),
        "units": units, "datasets": datasets,
        "dispersion_calculator_class": dispersion.get("calculator_class"),
        "dispersion_parameters": dispersion_parameters,
    }
    coarse = positive_triplet(result["basis"]["coarse_grid"]["shape"], "coarse grid")
    fine = positive_triplet(result["basis"]["fine_grid"]["shape"], "fine grid")
    ecut = numeric(spec["ecut_eV"], "ecut")
    if ecut <= 0:
        raise ValueError("nonpositive cutoff")
    kmesh = positive_triplet(spec["kmesh"], "k mesh")
    coordinates = []
    if "ecut_eV" in scan:
        coordinates.append(ecut)
    if "kmesh" in scan:
        coordinates.extend(kmesh)
    if {"gpts", "h_A"} & set(scan):
        # Distinct requested h values mapping to the same actual grid do not
        # constitute independent refinement points.
        coordinates.extend(coarse)
        coordinates.extend(fine)
    if not coordinates:
        raise ValueError("no effective refinement coordinate")
    actual_setting = {
        "ecut_eV": ecut, "kmesh": kmesh,
        "coarse_grid": coarse, "fine_grid": fine,
        "fixed_input": settings,
    }
    attempt_match = re.search(r"__attempt(\d+)$", directory.name)
    attempt = int(attempt_match.group(1)) if attempt_match else 1
    return {
        "run_id": directory.name, "stage": spec["stage"], "method": method,
        "result_path": str(directory / "result.json"),
        "result_sha256": file_sha(directory / "result.json"),
        "input_sha256": result["input_sha256"],
        "runner_sha256": runner,
        "geometry_sha256": fingerprint(geometry),
        "source_geometry_sha256": geometry["source_geometry_sha256"],
        "group_id": fingerprint(identity)[:20], "group_identity": identity,
        "setting_id": fingerprint(actual_setting),
        "attempt": attempt, "started_unix": result.get("start_unix", 0),
        "result_mtime_ns": (directory / "result.json").stat().st_mtime_ns,
        "coordinates": coordinates, "scan_parameters": scan,
        "scan_input": {k: spec.get(k) for k in scan},
        "ecut_eV": ecut, "kmesh": kmesh,
        "coarse_grid": coarse, "fine_grid": fine,
        "natoms": n, "energy_free_eV_cell": energy,
        "energy_extrapolated_eV_cell": numeric(
            total["energy_extrapolated_eV"], "extrapolated energy"),
        "forces_eV_A": forces, "stress_kbar": stress,
    }


def compare(left, right):
    if left["natoms"] != right["natoms"] or left["geometry_sha256"] != right["geometry_sha256"]:
        raise ValueError("comparisons require the same ordered atoms and exact geometry")
    delta_e = (left["energy_free_eV_cell"] - right["energy_free_eV_cell"]) * 1000 / left["natoms"]
    delta_f = max(math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))
                  for a, b in zip(left["forces_eV_A"], right["forces_eV_A"]))
    delta_s = max(abs(x - y) for x, y in zip(left["stress_kbar"], right["stress_kbar"]))
    passed = {
        "energy": abs(delta_e) <= LIMITS["energy_meV_atom"],
        "forces": delta_f <= LIMITS["force_vector_eV_A"],
        "stress": delta_s <= LIMITS["stress_component_kbar"],
    }
    return {
        "left_run": left["run_id"], "right_run": right["run_id"],
        "delta_energy_signed_meV_atom": delta_e,
        "delta_energy_abs_meV_atom": abs(delta_e),
        "max_atom_force_vector_difference_eV_A": delta_f,
        "max_stress_component_difference_kbar": delta_s,
        "threshold_pass": passed, "all_pass": all(passed.values()),
    }


def analyze_group(records):
    identity = records[0]["group_identity"]
    points = sorted(records, key=lambda r: (tuple(r["coordinates"]), r["run_id"]))
    group = {
        "group_id": points[0]["group_id"], "stage": points[0]["stage"],
        "method": points[0]["method"], "identity": identity,
        "n_distinct_settings": len(points), "terminal_run": points[-1]["run_id"],
        "terminal_is_provisional_reference": False,
        "recommended_run": None, "points": [],
        "terminal_stability": None,
    }
    comparable = all(
        len(a["coordinates"]) == len(b["coordinates"])
        and all(x <= y for x, y in zip(a["coordinates"], b["coordinates"]))
        and any(x < y for x, y in zip(a["coordinates"], b["coordinates"]))
        for a, b in zip(points, points[1:]))
    if len(points) < 3:
        group["status"] = "not_enough"
    elif not comparable:
        group["status"] = "nonmonotonic_refinement"
    else:
        intervals = [compare(points[-3], points[-2]), compare(points[-2], points[-1])]
        stable = all(item["all_pass"] for item in intervals)
        group["terminal_stability"] = {"last_two_intervals": intervals, "all_pass": stable}
        group["terminal_is_provisional_reference"] = stable
        group["status"] = "terminal_unstable" if not stable else "no_candidate_passed"
    for i, point in enumerate(points):
        row = {key: value for key, value in point.items() if key not in
               {"group_identity", "forces_eV_A", "stress_kbar", "result_mtime_ns"}}
        row["against_terminal"] = compare(point, points[-1])
        row["against_previous"] = compare(point, points[i - 1]) if i else None
        # Require this and every intervening higher setting to agree with the
        # terminal reference; a fortuitous isolated agreement is insufficient.
        row["all_higher_points_pass_terminal"] = all(
            compare(other, points[-1])["all_pass"] for other in points[i:])
        row["eligible_as_recommended"] = (
            group["terminal_is_provisional_reference"]
            and i < len(points) - 1 and row["all_higher_points_pass_terminal"])
        row["is_terminal_reference_only"] = i == len(points) - 1
        if row["eligible_as_recommended"] and group["recommended_run"] is None:
            group["recommended_run"] = point["run_id"]
            group["status"] = "converged_fixed_geometry"
        group["points"].append(row)
    return group


def collect(root):
    excluded, candidates = [], []
    for directory in sorted((root / "runs").glob("*")):
        if not directory.is_dir() or not (directory / "input.json").exists():
            continue
        try:
            spec = read_json(directory / "input.json")
            if spec.get("stage") not in STAGES:
                excluded.append({"run_id": directory.name, "reason": "stage_not_allowed",
                                 "stage": spec.get("stage")})
                continue
            if not (directory / "result.json").exists():
                excluded.append({"run_id": directory.name, "reason": "pending_no_result"})
                continue
            result = read_json(directory / "result.json")
            candidates.append(extract_record(directory, spec, result))
        except (KeyError, ValueError, TypeError, OSError) as exc:
            excluded.append({"run_id": directory.name, "reason": "invalid_or_unsuccessful",
                             "detail": str(exc)})
    # The latest successful physical setting wins, irrespective of directory
    # suffix.  A newer failure never erases the previous successful result.
    selected = {}
    for record in sorted(candidates, key=lambda r: (
            r["started_unix"], r["attempt"], r["result_mtime_ns"], r["run_id"])):
        key = (record["group_id"], record["setting_id"])
        if key in selected:
            excluded.append({"run_id": selected[key]["run_id"],
                             "reason": "superseded_successful_attempt",
                             "selected_run_id": record["run_id"]})
        selected[key] = record
    grouped = {}
    for record in selected.values():
        grouped.setdefault(record["group_id"], []).append(record)
    groups = [analyze_group(records) for _, records in sorted(grouped.items())]
    return groups, excluded


def render_report(summary):
    groups = summary["groups"]
    lines = [
        "# 固定几何收敛进度",
        "",
        f"生成时间：{summary['generated_utc']}。",
        "",
        "仅分析数值收敛与方法/几何转移的指定阶段；pilot、导数验证、"
        "旧硬截断色散、失败或未完成记录均不参加收敛判定，原始文件全部保留。",
        "",
        "判据为自由能差 ≤ 0.5 meV/atom、逐原子力矢量差的最大范数 "
        "≤ 0.005 eV/Å、六个应力分量绝对差的最大值 ≤ 0.1 kbar。"
        "力判据不是两个最大力数值相减。",
        "",
        "按方法、实际晶胞/有序原子坐标、PAW 文件、软件与运行程序哈希、"
        "色散参数及全部非扫描输入分别分组。每组至少需要三个不同的实际设置；"
        "最后三点中的两段相邻差均须通过，末点才作为暂定参照。"
        "推荐点以及它之后所有点还须对末点通过全部门槛，最高点本身不被自动推荐。",
        "",
        "**这里只验证固定几何的能量、力和应力数值稳定性。"
        "晶胞长度差 ≤ 0.5% 尚未由此验证，须在接受参数与更高精度下分别优化后比较。"
        "泛函物理误差、PAW 数据集误差及所有其他几何的可靠性也不能由本表消除。**",
        "",
    ]
    if not groups:
        lines += ["当前没有符合范围且成功完成的数值扫描数据；状态为 `no_data`。",
                  "没有锁定任何生产参数。", ""]
    else:
        lines += [
            "| 阶段 | 方法 | 组 | 独立点数 | 状态 | 暂定末点 | 推荐点 |",
            "|---|---|---|---:|---|---|---|",
        ]
        for group in groups:
            lines.append("| " + " | ".join(str(x) for x in [
                group["stage"], group["method"], group["group_id"],
                group["n_distinct_settings"], group["status"],
                group["terminal_run"] if group["terminal_is_provisional_reference"] else "未确认",
                group["recommended_run"] or "无",
            ]) + " |")
        lines += ["", "逐点自由能、与末点及相邻点的差值、分组条件和来源哈希见"
                  " `tables/convergence_summary.json` 与 `.csv`。", ""]
        for group in groups:
            fixed = group["identity"]["fixed_input_parameters"]
            lines += [
                f"## {group['stage']} / {group['method']} / {group['group_id']}",
                "",
                "固定输入：" + json.dumps(fixed, ensure_ascii=False, sort_keys=True) + "。",
                "",
                f"当前最高设置：`{group['terminal_run']}`。"
                + ("末端两段已通过，仅作为本组暂定参照。"
                   if group["terminal_is_provisional_reference"] else
                   "末端稳定性尚未确认；下表差值只用于诊断，不代表已知真实误差。"),
                "",
                "| 设置 | 自由能 (eV/cell) | 对末点 ΔE (meV/atom) | 对末点 ΔF (eV/Å) | 对末点 Δσ (kbar) | 三项均通过 |",
                "|---|---:|---:|---:|---:|---|",
            ]
            for point in group["points"]:
                diff = point["against_terminal"]
                setting = (f"{point['ecut_eV']:g} eV; "
                           f"k={'×'.join(map(str, point['kmesh']))}; "
                           f"grid={'×'.join(map(str, point['coarse_grid']))}")
                label = ("末点自身" if point["is_terminal_reference_only"]
                         else "是" if diff["all_pass"] else "否")
                lines.append(
                    f"| {setting} | {point['energy_free_eV_cell']:.12f} | "
                    f"{diff['delta_energy_abs_meV_atom']:.7g} | "
                    f"{diff['max_atom_force_vector_difference_eV_A']:.7g} | "
                    f"{diff['max_stress_component_difference_kbar']:.7g} | {label} |")
            stability = group["terminal_stability"]
            if stability:
                lines += ["", "末端相邻差（依次为 ΔE、ΔF、Δσ，单位同上）：", ""]
                for interval in stability["last_two_intervals"]:
                    lines.append(
                        f"- `{interval['left_run']}` → `{interval['right_run']}`："
                        f"{interval['delta_energy_abs_meV_atom']:.7g}、"
                        f"{interval['max_atom_force_vector_difference_eV_A']:.7g}、"
                        f"{interval['max_stress_component_difference_kbar']:.7g}；"
                        + ("通过。" if interval["all_pass"] else "未通过。"))
            lines.append("")
    lines.append(f"已排除或尚待完成的记录数：{len(summary['excluded_runs'])}。")
    if summary["warnings"]:
        lines += ["", "额外提示："]
        lines += [f"- {warning}" for warning in summary["warnings"]]
    return "\n".join(lines) + "\n"


def write_outputs(root, groups, excluded):
    hashes_by_stage = {}
    for group in groups:
        key = (group["stage"], group["method"])
        hashes_by_stage.setdefault(key, set()).add(group["identity"]["runner_sha256"])
    warning_lines = [
        f"{stage} / {method} 包含多个 runner 哈希，严格分开；合并前需独立审查。"
        for (stage, method), hashes in sorted(hashes_by_stage.items()) if len(hashes) > 1]
    summary = {
        "schema_version": 1,
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "analyzer_sha256": file_sha(__file__),
        "root": str(root), "status": "no_data" if not groups else "partial_fixed_geometry_analysis",
        "allowed_stages": list(STAGES), "thresholds": LIMITS,
        "energy_definition": "total.energy_free_eV (eV per cell), divided by natoms",
        "force_metric": "max_i norm(F_i(candidate)-F_i(reference))",
        "stress_metric": "max over six ASE Voigt components in kbar",
        "minimum_distinct_settings": 3,
        "endpoint_policy": "Both terminal adjacent intervals must pass; endpoint is never a candidate.",
        "cell_length_relative_tolerance": 0.005,
        "cell_length_validation_status": "not_verified_by_fixed_geometry_data",
        "scientific_model_accuracy_status": "not_established_by_numerical_convergence",
        "groups": groups, "excluded_runs": excluded, "warnings": warning_lines,
    }
    columns = [
        "group_id", "stage", "method", "group_status", "run_id", "source_geometry_sha256",
        "ecut_eV", "kmesh", "coarse_grid", "fine_grid", "natoms",
        "energy_free_eV_cell", "energy_extrapolated_eV_cell",
        "terminal_run", "terminal_confirmed",
        "delta_energy_meV_atom", "delta_force_vector_eV_A", "delta_stress_kbar",
        "energy_pass", "force_pass", "stress_pass", "all_pass_terminal",
        "eligible_as_recommended", "recommended_run", "result_path", "result_sha256",
    ]
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=columns)
    writer.writeheader()
    for group in groups:
        for point in group["points"]:
            diff = point["against_terminal"]
            row = {key: point[key] for key in columns if key in point}
            for key in ("kmesh", "coarse_grid", "fine_grid"):
                row[key] = "x".join(map(str, point[key]))
            row.update(
                group_status=group["status"], terminal_run=group["terminal_run"],
                terminal_confirmed=group["terminal_is_provisional_reference"],
                delta_energy_meV_atom=diff["delta_energy_signed_meV_atom"],
                delta_force_vector_eV_A=diff["max_atom_force_vector_difference_eV_A"],
                delta_stress_kbar=diff["max_stress_component_difference_kbar"],
                energy_pass=diff["threshold_pass"]["energy"],
                force_pass=diff["threshold_pass"]["forces"],
                stress_pass=diff["threshold_pass"]["stress"],
                all_pass_terminal=diff["all_pass"],
                recommended_run=group["recommended_run"])
            writer.writerow(row)
    outputs = {
        root / "tables/convergence_summary.json": json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        root / "tables/convergence_summary.csv": "\ufeff" + buffer.getvalue(),
        root / "reports/固定几何收敛进度.md": render_report(summary),
    }
    for path, content in outputs.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_text(content, encoding="utf-8")
        temporary.replace(path)
    return summary


def self_test():
    """Small synthetic decision checks; no quantum calculations."""
    def record(i, energy, force=0, stress=0):
        return {
            "group_id": "test", "group_identity": {}, "stage": "numerical_cutoff",
            "method": "PBE", "run_id": f"p{i}", "coordinates": [i],
            "natoms": 1, "geometry_sha256": "same", "energy_free_eV_cell": energy,
            "forces_eV_A": [[force, 0, 0]], "stress_kbar": [stress, 0, 0, 0, 0, 0],
        }
    group = analyze_group([record(1, 0), record(2, .0001), record(3, .0002)])
    assert group["status"] == "converged_fixed_geometry"
    assert group["recommended_run"] == "p1"
    assert not group["points"][-1]["eligible_as_recommended"]
    assert analyze_group([record(1, 0), record(2, 0)])["status"] == "not_enough"
    assert analyze_group([record(1, 0), record(2, .005), record(3, .0051)])[
        "status"] == "terminal_unstable"
    assert analyze_group([record(1, 0), record(2, .0001), record(3, .0002, stress=.2)])[
        "status"] == "terminal_unstable"
    # Equal force magnitudes can have a large vector difference.
    a, b = record(1, 0, force=1), record(2, 0, force=-1)
    assert compare(a, b)["max_atom_force_vector_difference_eV_A"] == 2
    a["forces_eV_A"] = [[.004, .004, 0]]
    b["forces_eV_A"] = [[0, 0, 0]]
    assert not compare(a, b)["threshold_pass"]["forces"]
    a, b, c = record(1, 0), record(2, 0), record(3, 0)
    a["coordinates"], b["coordinates"], c["coordinates"] = [2, 6], [3, 4], [4, 8]
    assert analyze_group([a, b, c])["status"] == "nonmonotonic_refinement"
    with tempfile.TemporaryDirectory(prefix="gpaw-convergence-empty-") as tmp:
        groups, excluded = collect(Path(tmp))
        summary = write_outputs(Path(tmp), groups, excluded)
        assert summary["status"] == "no_data"
        assert (Path(tmp) / "tables/convergence_summary.csv").exists()
        root = Path(tmp)
        def fixture(name, ecut, when, *, nbands=16, shift=0, status="passed_single_point",
                    method="PBE", stage="numerical_cutoff", width=.02):
            directory = root / "runs" / name
            directory.mkdir(parents=True)
            spec = {
                "id": name, "stage": stage, "method": method, "ecut_eV": ecut,
                "kmesh": [2, 4, 6], "gpts": [64, 40, 24], "nbands": nbands,
                "scan_parameters": ["ecut_eV"], "runner_sha256_at_submission": "a" * 64,
            }
            (directory / "input.json").write_text(json.dumps(spec))
            result = {
                "input": spec, "input_sha256": file_sha(directory / "input.json"),
                "runner_sha256": "a" * 64, "status": status, "SCF_converged": True,
                "natoms": 1, "start_unix": when, "source_geometry_sha256": "b" * 64,
                "cell_A": [[4, 0, 0], [0, 4, 0], [0, 0, 4]], "positions_A": [[shift, 0, 0]],
                "units": {"energy": "eV/cell", "forces": "eV/Angstrom",
                          "stress": "eV/Angstrom^3; ASE tensile-positive",
                          "stress_voigt_order": VOIGT},
                "datasets": [{"symbol": "H", "file_sha256": "c" * 64}],
                "total": {"energy_free_eV": ecut * 1e-6, "energy_extrapolated_eV": ecut * 1e-6,
                          "forces_eV_A": [[0, 0, 0]], "stress_kbar": [0] * 6,
                          "stress_eV_A3": [0] * 6},
                "basis": {"coarse_grid": {"shape": [64, 40, 24]},
                          "fine_grid": {"shape": [128, 80, 48]}},
            }
            if method == "PBE-D3(BJ)":
                result["dispersion"] = {"calculator_parameters": {
                    "realspace_cutoff": {"width2": width, "width3": width}}}
            (directory / "result.json").write_text(json.dumps(result))
        fixture("ecut400", 400, 1)
        fixture("ecut500", 500, 2)
        fixture("ecut600", 600, 3)
        fixture("ecut500__attempt2", 500, 4)
        fixture("ecut500__attempt3", 500, 5, status="failed")
        fixture("different_nbands", 400, 6, nbands=32)
        fixture("different_geometry", 400, 7, shift=.01)
        fixture("pilot", 400, 8, stage="pilot")
        fixture("hard_cutoff", 400, 9, method="PBE-D3(BJ)", width=0)
        groups, excluded = collect(root)
        assert len(groups) == 3, groups
        converged = next(g for g in groups if g["status"] == "converged_fixed_geometry")
        assert converged["n_distinct_settings"] == 3
        assert "ecut500__attempt2" in [p["run_id"] for p in converged["points"]]
        assert "ecut500" not in [p["run_id"] for p in converged["points"]]
        assert sum(g["status"] == "not_enough" for g in groups) == 2
        assert any(x["reason"] == "superseded_successful_attempt" for x in excluded)
        assert any("hard-cutoff" in x.get("detail", "") for x in excluded)
        result = write_outputs(root, groups, excluded)
        assert result["cell_length_validation_status"] == "not_verified_by_fixed_geometry_data"
    print("Synthetic checks passed: empty input, minimum sample size, both tail intervals, "
          "endpoint exclusion, vector-force metric, stress threshold, monotonic refinement, "
          "latest successful attempt, strict fixed-parameter/geometry isolation, "
          "pilot/hard-cutoff exclusion, and JSON/CSV/report outputs.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    root = args.root.resolve()
    groups, excluded = collect(root)
    result = write_outputs(root, groups, excluded)
    print(json.dumps({"status": result["status"], "groups": len(groups),
                      "converged_groups": sum(g["status"] == "converged_fixed_geometry"
                                              for g in groups),
                      "excluded_or_pending_runs": len(excluded),
                      "cell_length_validation": result["cell_length_validation_status"]}))


if __name__ == "__main__":
    main()
