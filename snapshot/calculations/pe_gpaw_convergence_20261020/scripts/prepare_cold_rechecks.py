"""Prepare or verify final-geometry cold starts; never starts quantum jobs."""
import argparse
import json
from pathlib import Path

import analyze_relaxations as audit

ROOT = Path(__file__).resolve().parents[1]


def write_new(path, content):
    with path.open("x", encoding="utf-8") as handle:
        json.dump(content, handle, indent=2)
        handle.write("\n")


def inspect_pair(candidate_name, reference_name):
    if candidate_name == reference_name:
        raise ValueError("Candidate and reference must be distinct relaxation runs")
    pair = [audit.inspect_run(ROOT / "relaxations" / name)
            for name in (candidate_name, reference_name)]
    if (pair[0]["run_id"] == pair[1]["run_id"]
            or Path(pair[0]["directory"]).resolve() == Path(pair[1]["directory"]).resolve()):
        raise ValueError("Candidate and reference resolve to the same run or run identifier")
    for run in pair:
        if run["acceptance"] != "verified_relaxation":
            raise ValueError(f"{run['run_id']}: relaxation not independently accepted")
        if not run.get("full_precision_trajectory_sha256"):
            raise ValueError("Full-precision trajectory is required")
    if audit.canonical(pair[0]["comparison_identity"]) != audit.canonical(
            pair[1]["comparison_identity"]):
        raise ValueError("Relaxation methods, provenance or fixed settings differ")
    candidate_precision, reference_precision = map(audit.precision_coordinates, pair)
    if (not all(a <= b for a, b in zip(candidate_precision, reference_precision))
            or candidate_precision == reference_precision
            or pair[0]["precision"]["xc_stencil"] != pair[1]["precision"]["xc_stencil"]):
        raise ValueError("Reference must refine ecut/kmesh/grid monotonically at the same XC stencil")
    return pair


def prepare(pair):
    candidate, reference = pair
    defaults = {
        **reference["precision"],
        **reference["comparison_identity"]["fixed_numerical_parameters"],
        "method": reference["input"]["method"],
        "mpi_ranks": reference["input"]["mpi_ranks"],
        "stage": "optimized_geometry_recheck",
        "checkpoint": True,
    }
    points = [{
        "id": f"cold_{run['run_id']}",
        "structure": str(Path(run["directory"]) / "trajectory.traj"),
        "notes": "Fresh GPAW object at the full-precision final trajectory frame; no density or wavefunction restart.",
    } for run in pair]
    path = ROOT / "cold_final_geometries_queue.json"
    write_new(path, {"stop_on_failure": True, "defaults": defaults, "points": points})
    print(json.dumps({"status": "prepared_no_quantum_calculation",
                      "queue": str(path)}))


def verify(pair, candidate_result=None, reference_result=None):
    candidate, reference = pair
    report = {
        "status": "verified",
        "status_meaning": "Provenance/evidence verified; scientific pass is separately recomputed from the differences.",
        "method": reference["input"]["method"],
        "reference_precision": reference["precision"],
        "thresholds": audit.CROSS_LIMITS,
        "checks": [],
    }
    for role, run, selected_path in zip(
            ("candidate", "reference"), pair, (candidate_result, reference_result)):
        report[role] = {key: run[key] for key in
                        ("run_id", "result_sha256", "final_geometry_sha256")}
        point_path = (Path(selected_path) if selected_path is not None else
                      ROOT / "runs" / f"cold_{run['run_id']}" / "result.json")
        if not point_path.is_absolute():
            point_path = ROOT / point_path
        point_path = point_path.resolve()
        check = {
            "relaxation_run_id": run["run_id"],
            "status": "verified",
            "geometry_sha256": run["final_geometry_sha256"],
            "source_geometry_sha256": run["full_precision_trajectory_sha256"],
            "result_path": str(point_path),
            "point_result_sha256": audit.file_sha(point_path),
        }
        # This reads real point/input/snapshot/geometry files, checks every
        # provenance constraint, and recomputes E/F/stress differences.
        check["audit"] = audit.verify_cross_point(check, run, reference, ROOT)
        report["checks"].append(check)
    report["all_numerical_checks_pass"] = all(
        check["audit"]["all_pass"] for check in report["checks"])
    path = ROOT / "relaxations" / "cross_precision_verification.json"
    write_new(path, report)
    print(json.dumps({"evidence_status": report["status"],
                      "all_numerical_checks_pass": report["all_numerical_checks_pass"],
                      "report": str(path)}))
    if not report["all_numerical_checks_pass"]:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("prepare", "verify"))
    parser.add_argument("--candidate", default="cell_candidate")
    parser.add_argument("--reference", default="cell_reference")
    parser.add_argument("--candidate-result", type=Path,
                        help="verify only: explicit successful cold-point result.json, including an __attemptN retry")
    parser.add_argument("--reference-result", type=Path,
                        help="verify only: explicit successful reference cold-point result.json; relative paths use project root")
    args = parser.parse_args()
    if args.action == "prepare" and (args.candidate_result is not None or args.reference_result is not None):
        parser.error("--candidate-result/--reference-result apply only to verify")
    pair = inspect_pair(args.candidate, args.reference)
    if args.action == "prepare":
        prepare(pair)
    else:
        verify(pair, args.candidate_result, args.reference_result)


if __name__ == "__main__":
    main()
