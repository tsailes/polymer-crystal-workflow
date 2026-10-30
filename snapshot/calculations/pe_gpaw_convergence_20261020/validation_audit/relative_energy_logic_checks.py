"""Small synthetic checks for relative-energy decisions; no quantum engine."""
from copy import deepcopy
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from unittest.mock import patch

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import analyze_relative_energies as analysis


def point(index, relative=10.0, kmesh=None):
    return {
        "pair_id": str(index), "relative_energy_meV_atom": relative,
        "precision": {
            "ecut_eV": 600 + 100 * index, "kmesh": kmesh or [4, 8, 12],
            "gpts": [96, 60, 36], "fine_grid": [192, 120, 72],
            "actual_density_cutoffs_eV": {
                "density_pd2": 2 * (600 + 100 * index),
                "density_pd3": 8 * (600 + 100 * index),
                "hamiltonian_pd2": 2 * (600 + 100 * index),
                "hamiltonian_pd3": 8 * (600 + 100 * index)}}}


def records_for(points):
    records = []
    for p in points:
        for role in analysis.ROLES:
            energy = -60 + (12 * p["relative_energy_meV_atom"] / 1000
                            if role == "GFN2" else 0)
            records.append({
                "pair_id": p["pair_id"], "family_id": "synthetic-family",
                "family_identity": {"method": "PBE-D3(BJ)"},
                "role": role, "method": "PBE-D3(BJ)", "energy_free_eV_cell": energy,
                "precision": p["precision"], "started_unix": 1, "mtime_ns": 1,
                "run_id": role + "_" + p["pair_id"], "result_sha256": "synthetic"})
    return records


def analyze_records(records, pending=None, rejected=None, unreadable=None):
    refs = {role: {"source_sha256": "synthetic_" + role} for role in analysis.ROLES}
    with patch.object(analysis, "references", return_value=refs), patch.object(
            analysis, "scan", return_value=(records, pending or [], rejected or [],
                                            {}, unreadable or [])):
        return analysis.analyze(ROOT, required_methods=["PBE-D3(BJ)"])


def main():
    checks = []

    def check(name, condition):
        if not condition:
            raise AssertionError(name)
        checks.append({"name": name, "passed": True})

    good = [point(0, 10), point(1, 10.1), point(2, 10.15)]
    result = analysis.assess_family(good)
    check("three stable refinements accept the lowest qualified point",
          result["passed"] and result["minimal_accepted_candidate_pair_ids"] == ["0"])
    result = analysis.assess_family(good[:2])
    check("two precisions remain pending",
          not result["passed"] and result["status"] == "pending")
    result = analysis.assess_family([point(0, 0), point(1, .4), point(2, .8)])
    check("small adjacent changes cannot hide cumulative drift",
          result["passed"] and result["minimal_accepted_candidate_pair_ids"] == ["1"])
    result = analysis.assess_family([point(0, 0), point(1, .2), point(2, .9)])
    check("an unstable final interval fails",
          not result["passed"] and result["status"] == "needs_attention")
    incomparable = [point(0), point(1), point(2, kmesh=[2, 4, 6])]
    result = analysis.assess_family(incomparable)
    check("cutoff growth cannot compensate for lower kmesh",
          not analysis.refines(incomparable[-1], incomparable[-2])
          and result["status"] == "pending")
    equal_precision = deepcopy(good[2])
    equal_precision["pair_id"] = "different-BZ-same-mesh"
    result = analysis.assess_family(good + [equal_precision])
    check("distinct points at equal precision do not fabricate a unique highest",
          result["status"] == "pending")

    records = records_for(good)
    report = analyze_records(records)
    check("complete paired series passes", report["passed"])
    check("cell-energy difference is divided by twelve and converted to meV",
          abs(report["pairs"][0]["relative_energy_meV_atom"] - 10) < 1e-10)
    check("no available calculations remain pending",
          analyze_records([])["status"] == "pending")
    report = analyze_records(records[:-1])
    check("one missing structure cannot complete the series",
          not report["passed"] and report["unmatched_pairs"]
          and not report["coverage_complete"])
    incompatible_records = deepcopy(records)
    incompatible_records[-1]["pair_id"] += "_different_provenance"
    report = analyze_records(incompatible_records)
    check("different settings or provenance keys are never paired",
          not report["passed"] and len(report["unmatched_pairs"]) == 2)
    report = analyze_records(records, pending=[{"logical_id": "planned_extra"}])
    check("unstarted planned point prevents global completion",
          not report["passed"] and not report["coverage_complete"])
    report = analyze_records(records, unreadable=[{"path": "corrupt-input.json"}])
    check("unreadable input prevents global pass",
          not report["passed"] and report["status"] == "needs_attention")
    other_family = deepcopy(records[:2])
    for record in other_family:
        record["pair_id"] = "new-settings-one-point"
        record["family_id"] = "different-runner-family"
    report = analyze_records(records + other_family)
    check("passing family cannot hide an incomplete incompatible family",
          not report["passed"])
    duplicates = deepcopy(records[:2])
    duplicates[0]["energy_free_eV_cell"] += .4 * 12 / 1000
    duplicates[1]["energy_free_eV_cell"] -= .4 * 12 / 1000
    for row in duplicates:
        row["run_id"] += "_duplicate"
        row["started_unix"] = 2
    report = analyze_records(records + duplicates)
    check("opposite duplicate drifts are bounded on the relative energy itself",
          not report["passed"] and report["rejected_records"]
          and abs(report["rejected_records"][0]["relative_duplicate_spread_meV_atom"]
                  - .8) < 1e-10)
    with TemporaryDirectory(prefix="relative_logic_", dir=ROOT / "validation_audit") as temp:
        directory = Path(temp)
        analysis.save(directory / "planned_queue.json", {
            "defaults": {"stage": "numerical_joint", "method": "PBE"},
            "points": [{"id": "planned_without_explicit_structure", "ecut_eV": 600}]})
        _, pending, rejected, _, unreadable = analysis.scan(
            directory, {}, analysis.STAGES)
        check("unstarted queue uses the scheduler default reference and stays pending",
              len(pending) == 1 and not rejected and not unreadable
              and Path(pending[0]["structure"]) == directory / "reference" / "POSCAR")

    # Parse one existing native result as a schema/units regression.  This is
    # excluded from transfer output and is not a newly calculated energy.
    refs = analysis.references(ROOT)
    native = ROOT / "runs" / "cutoff_grid96_pbe_E1000_k4x8x12"
    record = analysis.extract(native, analysis.load(native / "input.json"),
                              analysis.load(native / "result.json"), refs, ROOT)
    check("existing native output passes independent geometry and units checks",
          record["role"] == "Tashiro" and record["method"] == "PBE"
          and record["precision"]["actual_density_cutoffs_eV"]["density_pd2"] == 2000)
    evidence = {
        "passed": True, "checks": checks,
        "quantum_calculations_started": False,
        "synthetic_energy_values_are_not_physical_results": True,
        "analyzer_sha256": analysis.sha(analysis.__file__),
        "test_sha256": analysis.sha(__file__),
        "native_regression_result_path": str(native / "result.json"),
        "native_regression_result_sha256": analysis.sha(native / "result.json")}
    analysis.save(ROOT / "validation_audit" / "relative_energy_logic_checks.json", evidence)
    print(f"{len(checks)} relative-energy logic/schema checks passed; no quantum calculations")


if __name__ == "__main__":
    main()
