"""Read-only comparison of saved base/tighter-SCF PE stress endpoints.

Writes validation_audit/scf_noise_validation.{json,csv,md}.  Missing successful
endpoints remain pending.  No calculator or quantum engine is imported.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import numpy as np
import analyze_derivatives as derivative

ROOT = Path(__file__).resolve().parents[1]
SCF_KEYS = ("scf_energy_eV_electron", "scf_density_electron_electron",
            "scf_eigenstates_eV2_electron")
ENERGY_TOL_EV_CELL = 1e-8
FD_STRESS_TOL_KBAR = .02


def scientific_settings(spec):
    settings = derivative.effective_settings(spec)
    settings = {key: val for key, val in settings.items()
                if key not in {*SCF_KEYS, "original_id"}}
    # Unlike center versus strained endpoints, this comparison is between
    # identical geometries and must preserve even the explicit grid inputs.
    settings.update({key: spec.get(key) for key in ("gpts", "h_A", "mpi_ranks")})
    return settings


def validate_pair(base, tight, checks):
    """Reject any change except explicitly tightened SCF criteria and metadata."""
    a, b = base["result"], tight["result"]
    sa, sb = base["input"], tight["input"]
    ca, cb = checks[base["folder"].name], checks[tight["folder"].name]
    if derivative.logical_id(sb.get("original_id", "")) != base["id"]:
        raise ValueError("Tighter endpoint original_id does not identify the selected base point")
    if scientific_settings(sa) != scientific_settings(sb):
        different = sorted(k for k in scientific_settings(sa).keys() |
                           scientific_settings(sb).keys()
                           if scientific_settings(sa).get(k) != scientific_settings(sb).get(k))
        raise ValueError("Non-SCF numerical/method settings changed: " + ", ".join(different))
    for key in SCF_KEYS:
        original = float(sa.get(key, derivative.DEFAULT_SETTINGS[key]))
        tightened = float(sb.get(key, derivative.DEFAULT_SETTINGS[key]))
        if original <= 0 or not np.isclose(tightened / original, .1, atol=0, rtol=1e-12):
            raise ValueError("SCF criterion is not tightened by exactly ten: " + key)
    # These metadata are omitted from effective_settings but are essential for
    # identifying the physical endpoint and avoiding swapped +/- comparisons.
    for key in ("center_id", "derivative_kind", "component", "step", "sign",
                "displacement", "strain"):
        if sa.get(key) != sb.get(key):
            raise ValueError("Endpoint transformation/label changed: " + key)
    if a["source_geometry_sha256"] != b["source_geometry_sha256"]:
        raise ValueError("Base and tighter point have different original probe sources")
    if ca["symbols"] != cb["symbols"] or ca["datasets"] != cb["datasets"]:
        raise ValueError("Element order or PAW dataset fingerprints/files changed")
    if a["electrons"] != b["electrons"]:
        raise ValueError("Electron count changed")
    if a["mpi_world_size"] != b["mpi_world_size"]:
        raise ValueError("MPI size changed, confounding the SCF-only comparison")
    if (not np.allclose(a["cell_A"], b["cell_A"], atol=1e-12, rtol=0)
            or not np.allclose(a["positions_A"], b["positions_A"], atol=1e-12, rtol=0)
            or not np.isclose(a["volume_A3"], b["volume_A3"], atol=1e-12, rtol=0)):
        raise ValueError("Unrounded result geometry changed")
    if a["software"] != b["software"] or a["runner_sha256"] != b["runner_sha256"]:
        raise ValueError("Software or saved runner changed")
    for key in ("xc_name", "xc_class"):
        if a["basis"][key] != b["basis"][key]:
            raise ValueError("Actual XC implementation changed")
    if a["dispersion"].get("calculator_parameters") != b["dispersion"].get("calculator_parameters"):
        raise ValueError("Dispersion damping, CN, cutoffs or smoothing changed")
    if not (ca["basis"]["density_hamiltonian_basis_recorded"] and
            cb["basis"]["density_hamiltonian_basis_recorded"]):
        raise ValueError("Full wavefunction/density/Hamiltonian basis provenance is required")
    if ca["basis"]["merged_sha256"] != cb["basis"]["merged_sha256"]:
        raise ValueError("Actual wavefunction/density/Hamiltonian G sets, grids or k points changed")
    return {
        "non_SCF_settings_identical": True,
        "geometry_identical_tolerance_A": 1e-12,
        "basis_merged_sha256": ca["basis"]["merged_sha256"],
        "SCF_ratio_tight_over_base": {
            key: float(sb.get(key, derivative.DEFAULT_SETTINGS[key])) /
            float(sa.get(key, derivative.DEFAULT_SETTINGS[key])) for key in SCF_KEYS},
    }


def endpoint_changes(base, tight):
    values = {}
    for level in derivative.LEVELS:
        old, new = base["result"][level], tight["result"][level]
        energy_difference = derivative.energy(tight["result"], level) - derivative.energy(base["result"], level)
        forces = np.asarray(new["forces_eV_A"]) - np.asarray(old["forces_eV_A"])
        stress = (np.asarray(new["stress_eV_A3"]) - np.asarray(old["stress_eV_A3"])) * derivative.KBAR
        values[level] = {
            "base_free_energy_eV_cell": derivative.energy(base["result"], level),
            "tight_free_energy_eV_cell": derivative.energy(tight["result"], level),
            "energy_change_eV_cell": energy_difference,
            "absolute_energy_change_eV_cell": abs(energy_difference),
            "energy_change_passed": abs(energy_difference) <= ENERGY_TOL_EV_CELL,
            "force_change_eV_A": forces.tolist(),
            "max_abs_force_component_change_eV_A": float(np.max(np.abs(forces))),
            "max_force_vector_change_eV_A": float(np.max(np.linalg.norm(forces, axis=1))),
            "stress_change_kbar": stress.tolist(),
            "max_abs_stress_component_change_kbar": float(np.max(np.abs(stress))),
        }
    return values


def run_provenance(run):
    return {
        "logical_id": run["id"], "selected_run": run["folder"].name,
        "input_sha256": derivative.digest(run["folder"] / "input.json"),
        "result_sha256": derivative.digest(run["folder"] / "result.json"),
        "runner_sha256": run["result"].get("runner_sha256"),
        "SCF_converged": run["result"].get("SCF_converged"),
        "SCF_iterations": run["result"].get("SCF_iterations"),
        "source_geometry_sha256": run["result"].get("source_geometry_sha256"),
    }


def analyze(root, queue_path, original_prefix, tight_prefix):
    queue = derivative.load(queue_path)
    requested = [{**queue.get("defaults", {}), **point} for point in queue["points"]]
    if len(requested) != 6 or len({spec["id"] for spec in requested}) != 6:
        raise ValueError("Protocol requires six unique tighter-SCF endpoints")
    expected_methods = {method for _, method in derivative.METHODS}
    grouped = {}
    for spec in requested:
        if (spec.get("method") not in expected_methods or spec.get("derivative_kind") != "stress"
                or spec.get("component") != 5 or spec.get("step") != 1e-5
                or spec.get("sign") not in (-1, 1) or not spec.get("original_id")):
            raise ValueError("Queue differs from the agreed three-method xy +/- 1e-5 check")
        grouped.setdefault(spec["method"], []).append(spec)
    if set(grouped) != expected_methods or any(
            len(points) != 2 or {point["sign"] for point in points} != {-1, 1}
            for points in grouped.values()):
        raise ValueError("Each method must have exactly one negative and one positive endpoint")

    original_attempts, unreadable_original = derivative.discover(root, original_prefix)
    tight_attempts, unreadable_tight = derivative.discover(root, tight_prefix)
    checks, invalid, selected, points, by_method = {}, {}, {}, [], {}

    def get_run(attempts, ident):
        run = derivative.select(attempts, ident)
        if run is not None and run["folder"].name not in selected:
            name = run["folder"].name
            selected[name] = run
            try:
                checks[name] = derivative.validate_run(run)
            except (OSError, ValueError, KeyError, TypeError, IndexError) as exc:
                invalid[name] = str(exc)
        return run

    for spec in requested:
        base = get_run(original_attempts, spec["original_id"])
        tight = get_run(tight_attempts, spec["id"])
        center = get_run(original_attempts, spec["center_id"])
        point = {
            "method": spec["method"], "component": "xy", "step": spec["step"],
            "sign": spec["sign"], "original_id": spec["original_id"],
            "tight_id": spec["id"], "center_id": spec["center_id"],
            "status": "pending", "coverage_complete": bool(base and tight and center),
            "reasons": [], "levels": {},
            "base": run_provenance(base) if base else None,
            "tight": run_provenance(tight) if tight else None,
            "original_attempts": derivative.history(original_attempts, spec["original_id"]),
            "tight_attempts": derivative.history(tight_attempts, spec["id"]),
        }
        for run in (base, tight, center):
            if run and run["folder"].name in invalid:
                point["reasons"].append(run["folder"].name + ": " + invalid[run["folder"].name])
        if tight is not None:
            # Ensure the executed job is the declared queue job, allowing only
            # run_queue's provenance fields to be added.
            for key, value in spec.items():
                if tight["input"].get(key) != value:
                    point["reasons"].append("Executed tighter input differs from queue: " + key)
        if point["coverage_complete"] and not point["reasons"]:
            try:
                point["original_endpoint_validation"] = derivative.compare_endpoint(
                    center, base, checks, "stress", 5, 1e-5, spec["sign"])
                if not point["original_endpoint_validation"]["basis_unchanged"]:
                    raise ValueError("Base endpoint basis differs from its center")
                point["pair_validation"] = validate_pair(base, tight, checks)
                point["levels"] = endpoint_changes(base, tight)
                point["center_volume_A3"] = center["result"]["volume_A3"]
                point["status"] = ("passed" if all(
                    point["levels"][level]["energy_change_passed"]
                    for level in ("dft", "total")) else "failed")
            except (OSError, ValueError, KeyError, TypeError, IndexError) as exc:
                point["reasons"].append(str(exc))
        if point["reasons"]:
            point["status"] = "invalid_comparison"
        elif not point["coverage_complete"]:
            point["reasons"].append("Base/tighter endpoint or base center is not yet successfully completed")
        points.append(point)
        by_method.setdefault(spec["method"], {})[spec["sign"]] = point

    groups = []
    for method in sorted(expected_methods):
        pair = by_method[method]
        group = {"method": method, "component": "xy", "step": 1e-5,
                 "status": "pending", "coverage_complete": all(
                     point["coverage_complete"] for point in pair.values()),
                 "reasons": [], "levels": {}}
        if any(point["status"] == "invalid_comparison" for point in pair.values()):
            group["status"] = "invalid_comparison"
            group["reasons"] = [reason for point in pair.values() for reason in point["reasons"]]
        elif group["coverage_complete"]:
            minus, plus = pair[-1], pair[1]
            if minus["center_id"] != plus["center_id"]:
                group["status"] = "invalid_comparison"
                group["reasons"].append("Positive and negative endpoints reference different centers")
            else:
                volume = minus["center_volume_A3"]
                group["center_volume_A3"] = volume
                for level in derivative.LEVELS:
                    m, p = minus["levels"][level], plus["levels"][level]
                    old_fd = ((p["base_free_energy_eV_cell"] - m["base_free_energy_eV_cell"])
                              / (2e-5 * volume) * derivative.KBAR)
                    new_fd = ((p["tight_free_energy_eV_cell"] - m["tight_free_energy_eV_cell"])
                              / (2e-5 * volume) * derivative.KBAR)
                    delta = new_fd - old_fd
                    group["levels"][level] = {
                        "base_FD_stress_kbar": old_fd, "tight_FD_stress_kbar": new_fd,
                        "FD_stress_change_kbar": delta, "absolute_FD_stress_change_kbar": abs(delta),
                        "FD_change_passed": abs(delta) <= FD_STRESS_TOL_KBAR}
                group["status"] = ("passed" if
                    all(point["status"] == "passed" for point in pair.values()) and
                    all(group["levels"][level]["FD_change_passed"] for level in ("dft", "total"))
                    else "failed")
        else:
            group["reasons"].append("Both successful base/tighter endpoints and a center are required")
        groups.append(group)

    complete = all(group["coverage_complete"] for group in groups)
    passed = complete and all(group["status"] == "passed" for group in groups)
    counts = {status: sum(group["status"] == status for group in groups)
              for status in ("passed", "pending", "failed", "invalid_comparison")}
    return {
        "schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": ("SCF_noise_checks_passed" if passed else "needs_attention"
                   if counts["invalid_comparison"] or counts["failed"] else "pending"),
        "coverage_complete": complete, "passed": passed,
        "quantum_calculations_started": False, "counts": counts,
        "queue_path": str(queue_path), "queue_sha256": derivative.digest(queue_path),
        "analyzer_sha256": derivative.digest(__file__),
        "derivative_helper_sha256": derivative.digest(Path(derivative.__file__)),
        "thresholds": {"absolute_energy_change_eV_cell_per_endpoint": ENERGY_TOL_EV_CELL,
                       "absolute_FD_stress_change_kbar": FD_STRESS_TOL_KBAR},
        "scope": "SCF-threshold sensitivity of the same three-method xy +/- 1e-5 endpoints. "
                 "This is a stability check, not a bound on the exact unknown SCF error. "
                 "All reported derivatives use free energies and the base probe-center volume. "
                 "Forces and endpoint stress changes are reported without extra acceptance thresholds.",
        "expected_tighter_points": 6, "expected_method_pairs": 3,
        "unreadable_in_progress_inputs": unreadable_original + unreadable_tight,
        "selected_runs": {name: {
            **run_provenance(run), "validation_error": invalid.get(name),
            "merged_basis": checks.get(name, {}).get("basis")}
            for name, run in sorted(selected.items())},
        "points": points, "groups": groups,
    }


def output(root, report):
    folder = root / "validation_audit"
    folder.mkdir(parents=True, exist_ok=True)
    derivative.write_json(folder / "scf_noise_validation.json", report)
    rows = []
    for point in report["points"]:
        for level in derivative.LEVELS:
            val = point["levels"].get(level, {})
            rows.append({
                "method": point["method"], "record_type": "endpoint", "contribution": level,
                "component": "xy", "step": point["step"], "sign": point["sign"],
                "base_run": point["base"]["selected_run"] if point["base"] else "",
                "tight_run": point["tight"]["selected_run"] if point["tight"] else "",
                "energy_change_eV_cell": val.get("energy_change_eV_cell"),
                "max_force_vector_change_eV_A": val.get("max_force_vector_change_eV_A"),
                "max_abs_stress_component_change_kbar": val.get("max_abs_stress_component_change_kbar"),
                "FD_stress_change_kbar": None, "status": point["status"],
                "reasons": " | ".join(point["reasons"])})
    for group in report["groups"]:
        for level in derivative.LEVELS:
            rows.append({
                "method": group["method"], "record_type": "FD_pair", "contribution": level,
                "component": "xy", "step": group["step"], "sign": "",
                "base_run": "", "tight_run": "", "energy_change_eV_cell": None,
                "max_force_vector_change_eV_A": None,
                "max_abs_stress_component_change_kbar": None,
                "FD_stress_change_kbar": group["levels"].get(level, {}).get("FD_stress_change_kbar"),
                "status": group["status"], "reasons": " | ".join(group["reasons"])})
    with (folder / "scf_noise_validation.csv").open("w", encoding="utf-8-sig", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    lines = [
        "# 电子自洽阈值敏感性复核", "",
        f"状态：**{report['status']}**；六端点覆盖完成：{report['coverage_complete']}；通过：{report['passed']}。", "",
        "三个方法分别在同一探针的 xy 工程剪切 ±1e−5 端点，将三项 SCF 阈值收紧10倍。"
        "比较使用未舍入几何、相同PAW及电子数、泛函与色散参数、软件/运行器，"
        "并要求全部rank的波函数/密度/Hamiltonian实际G集合、FFT网格与k点一致。", "",
        "门槛：每端点自由能变化≤1e−8 eV/晶胞；正负端点的差分应力变化≤0.02 kbar。"
        "差分分母为2hV，V使用原探针中心体积；DFT与total均需满足门槛。"
        "全部原子力和六个应力分量的变化另行报告，不引入未经约定的新阈值。", "",
        "这是对SCF阈值收紧的稳定性证据，不是精确SCF误差的真值；"
        "也不能替代最终批的完整力/应力导数检查或cutoff/k网格收敛。", "",
        "| 方法 | 分项 | 两端最大能量变化 (eV/胞) | 差分应力变化 (kbar) | 状态 |",
        "|---|---|---:|---:|---|",
    ]
    for group in report["groups"]:
        for level in ("dft", "total"):
            changes = [point["levels"][level]["absolute_energy_change_eV_cell"]
                       for point in report["points"] if point["method"] == group["method"]
                       and level in point["levels"]]
            energy = f"{max(changes):.9g}" if len(changes) == 2 else "待完成"
            delta = group["levels"].get(level, {}).get("FD_stress_change_kbar")
            stress = f"{delta:.9g}" if delta is not None else "待完成"
            lines.append(f"| {group['method']} | {level} | {energy} | {stress} | {group['status']} |")
    reasons = [(point["method"], point["sign"], reason)
               for point in report["points"] for reason in point["reasons"]]
    if reasons:
        lines += ["", "## 待完成或拒绝原因", ""]
        lines += [f"- {method}，符号{sign:+d}：{reason}" for method, sign, reason in reasons]
    lines += ["", f"生成时间：{report['created_utc']}", ""]
    (folder / "scf_noise_validation.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--queue", type=Path)
    parser.add_argument("--original-prefix", default="deriv_final")
    parser.add_argument("--tight-prefix", default="scfnoise_final")
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    queue = args.queue.resolve() if args.queue else root / "scf_noise_queue.json"
    report = analyze(root, queue, args.original_prefix, args.tight_prefix)
    output(root, report)
    print(json.dumps({key: report[key] for key in
                      ("status", "coverage_complete", "passed", "counts")}))
    if args.require_complete and not report["passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
