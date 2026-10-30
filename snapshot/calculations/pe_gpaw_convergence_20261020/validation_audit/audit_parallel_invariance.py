"""Read-only numerical audit of the explicitly scheduled four/eight MPI pair."""
from pathlib import Path
import datetime as dt
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parents[1]
IDS = ["parallel4_pbe_E1400_k3x6x9_g96", "grid_pbe_E1400_k3x6x9_g96x60x36"]
LIMITS = {"energy_eV_cell": 1e-6, "force_eV_A": 1e-6, "stress_kbar": 1e-5}
METADATA = {"id", "stage", "scan_parameters", "scan_label", "scan_index",
            "queue_sha256", "description", "notes", "checkpoint"}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def norm(v):
    return math.sqrt(sum(x * x for x in v))


def read(run_id):
    directory = ROOT / "runs" / run_id
    result = json.loads((directory / "result.json").read_text())
    spec = json.loads((directory / "input.json").read_text())
    errors = []
    if result.get("status") != "passed_single_point" or result.get("SCF_converged") is not True:
        errors.append("Single point did not pass SCF")
    if result.get("input") != spec or result.get("input_sha256") != sha(directory / "input.json"):
        errors.append("Input JSON/content SHA mismatch")
    if (result.get("runner_sha256") != sha(directory / "run_point_snapshot.py")
            or result.get("runner_sha256") != spec.get("runner_sha256_at_submission")):
        errors.append("Runner snapshot/submission SHA mismatch")
    if result.get("mpi_world_size") != spec.get("mpi_ranks"):
        errors.append("Actual MPI world size differs from requested size")
    source = Path(spec["structure"])
    if not source.is_absolute():
        source = directory / source
    if not source.is_file() or sha(source) != result.get("source_geometry_sha256"):
        errors.append("Source geometry file SHA mismatch")
    units = result.get("units", {})
    if (units.get("energy") != "eV/cell" or units.get("forces") != "eV/Angstrom"
            or units.get("stress_voigt_order") != ["xx", "yy", "zz", "yz", "xz", "xy"]):
        errors.append("Physical units mismatch")
    data = result["total"]
    natoms = result["natoms"]
    if len(data["forces_eV_A"]) != natoms or any(len(row) != 3 for row in data["forces_eV_A"]):
        errors.append("Force array shape mismatch")
    if len(data["stress_kbar"]) != 6 or len(data["stress_eV_A3"]) != 6:
        errors.append("Stress array shape mismatch")
    values = [data["energy_free_eV"], data["energy_extrapolated_eV"],
              *data["stress_kbar"], *data["stress_eV_A3"],
              *[value for row in data["forces_eV_A"] for value in row]]
    if not all(math.isfinite(value) for value in values):
        errors.append("Nonfinite energy/force/stress")
    if any(abs(a - b * 1602.1766208) > 1e-7
           for a, b in zip(data["stress_kbar"], data["stress_eV_A3"])):
        errors.append("Stress conversion mismatch")
    paw = [{key: atom.get(key) for key in ("symbol", "type", "setupname", "fingerprint", "file_sha256")}
           for atom in result["datasets"]]
    if len(paw) != natoms or any(not atom["file_sha256"] for atom in paw):
        errors.append("Incomplete PAW evidence")
    for atom in result["datasets"]:
        path = Path(atom["filename"])
        if not path.is_file() or sha(path) != atom["file_sha256"]:
            errors.append("Actual installed PAW file hash mismatch")
    identity = {
        "input_except_mpi_metadata": {k: v for k, v in spec.items()
                                      if k not in METADATA | {"mpi_ranks"}},
        "natoms": natoms, "formula": result["formula"], "cell_A": result["cell_A"],
        "positions_A": result["positions_A"], "source_geometry_sha256": result["source_geometry_sha256"],
        "datasets": paw, "runner_sha256": result["runner_sha256"], "software": result["software"],
        "units": units,
        "physical_basis": {key: result["basis"].get(key)
                           for key in ("coarse_grid", "fine_grid", "ngmin", "ngmax",
                                       "bz_kpoints", "ibz_kpoints", "xc_class", "xc_name")},
        "dispersion_class": result.get("dispersion", {}).get("calculator_class"),
        "dispersion_parameters": result.get("dispersion", {}).get("calculator_parameters"),
    }
    return {"run_id": run_id, "result_path": str(directory / "result.json"),
            "result_sha256": sha(directory / "result.json"), "input_sha256": sha(directory / "input.json"),
            "mpi_ranks": result["mpi_world_size"], "SCF_iterations": result["SCF_iterations"],
            "wall_seconds": result["wall_seconds"], "identity": identity,
            "total": data, "errors": errors}


def main():
    if not all((ROOT / "runs" / name / "result.json").exists() for name in IDS):
        print(json.dumps({"status": "pending", "reason": "Both scheduled MPI point results are required"}))
        return
    try:
        a, b = map(read, IDS)
        differing = [key for key in a["identity"] if a["identity"][key] != b["identity"][key]]
        errors = a["errors"] + b["errors"]
        if sorted([a["mpi_ranks"], b["mpi_ranks"]]) != [4, 8]:
            errors.append("Expected four/eight MPI pair")
        left, right = a["total"], b["total"]
        differences = {
            "energy_eV_cell": abs(left["energy_free_eV"] - right["energy_free_eV"]),
            "energy_extrapolated_eV_cell": abs(left["energy_extrapolated_eV"] - right["energy_extrapolated_eV"]),
            "force_eV_A": max(norm([x - y for x, y in zip(u, v)])
                             for u, v in zip(left["forces_eV_A"], right["forces_eV_A"])),
            "stress_kbar": max(abs(x - y) for x, y in zip(left["stress_kbar"], right["stress_kbar"])),
        }
        flags = {key: differences[key] <= limit for key, limit in LIMITS.items()}
        passed = not errors and not differing and all(flags.values())
        output = {
            "status": "passed" if passed else "failed",
            "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "audit_script_sha256": sha(Path(__file__)),
            "method": a["identity"]["input_except_mpi_metadata"]["method"],
            "checks": {"same_geometry_PAW_runner_settings_except_MPI": not differing,
                       "differing_identity_fields": differing, "evidence_errors": errors},
            "limits": LIMITS, "differences": differences, "threshold_pass": flags,
            "runs": [{key: value for key, value in run.items() if key not in ("total", "identity", "errors")}
                     for run in (a, b)],
            "common_identity": a["identity"],
            "note": "Technical MPI invariance at this one geometry/configuration; not an independent physical accuracy or complete transfer test.",
        }
    except (ValueError, KeyError, TypeError, OSError) as exc:
        output = {"status": "failed", "error": repr(exc),
                  "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat()}
    destination = ROOT / "validation_audit"
    (destination / "parallel_invariance.json").write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n")
    lines = ["# 4 与 8 MPI 进程不变性验收", "", f"状态：**{output['status']}**。", ""]
    if "differences" in output:
        lines += [
            "检查同一 PE 几何、PAW、软件和程序快照，除 MPI 进程数和行政元数据外输入完全一致。",
            "", "| 指标 | 实际绝对差 | 阈值 | 通过 |", "|---|---:|---:|---|"]
        for key, label in (("energy_eV_cell", "自由能 eV/cell"),
                           ("force_eV_A", "最大原子力矢量差 eV/Å"),
                           ("stress_kbar", "最大应力分量差 kbar")):
            lines.append(f"| {label} | {output['differences'][key]:.12g} | "
                         f"{LIMITS[key]:.12g} | {output['threshold_pass'][key]} |")
        lines += ["", "完整几何、PAW 与文件哈希、实际进程数及输入一致性证据见同名 JSON。",
                  "这里只验证本几何与本参数的并行数值不变性，不代替方法、其他几何或晶胞优化的验收。"]
        if output["checks"]["evidence_errors"] or output["checks"]["differing_identity_fields"]:
            lines += ["", "证据问题：" + json.dumps(output["checks"], ensure_ascii=False)]
    else:
        lines += [str(output.get("error"))]
    (destination / "parallel_invariance.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({key: output[key] for key in ("status", "differences", "checks") if key in output},
                     ensure_ascii=False))


if __name__ == "__main__":
    main()
