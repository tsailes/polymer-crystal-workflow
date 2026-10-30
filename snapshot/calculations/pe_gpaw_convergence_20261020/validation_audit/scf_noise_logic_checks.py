"""Cheap SCF comparison logic checks, with no quantum calculations."""
import copy
import json
from pathlib import Path
import sys
import tempfile

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import analyze_scf_noise as noise

a = noise.derivative
results = []


def check(name, function):
    try:
        function()
        results.append({"test": name, "passed": True})
    except Exception as exc:
        results.append({"test": name, "passed": False, "error": repr(exc)})


def must_raise(function):
    try:
        function()
    except ValueError:
        return
    raise AssertionError("Expected rejection")


attempts, _ = a.discover(ROOT, "scfnoise_final")
tight = a.select(attempts, "scfnoise_final_pbe_d3bj_S5_h1e-05_m")
tight_check = a.validate_run(tight)
base = copy.deepcopy(tight)
base["id"] = tight["input"]["original_id"]
base["input"]["id"] = base["id"]
base["input"].pop("original_id")
base["folder"] = Path("synthetic_base")
for key in noise.SCF_KEYS:
    base["input"][key] *= 10
checks = {base["folder"].name: copy.deepcopy(tight_check),
          tight["folder"].name: tight_check}


def compare(changed=None, changed_checks=None):
    return noise.validate_pair(base, changed or tight, changed_checks or checks)


def changed_input(key, value):
    altered = copy.deepcopy(tight)
    altered["input"][key] = value
    must_raise(lambda: compare(altered))


def geometry():
    altered = copy.deepcopy(tight)
    altered["result"]["positions_A"][0][0] += 1e-7
    must_raise(lambda: compare(altered))


def basis():
    altered = copy.deepcopy(checks)
    altered[tight["folder"].name]["basis"]["merged_sha256"] = "f" * 64
    must_raise(lambda: compare(changed_checks=altered))


def endpoint_values():
    altered = copy.deepcopy(tight)
    for level in ("dft", "total"):
        altered["result"][level]["energy_free_eV"] += 2e-8
    values = noise.endpoint_changes(base, altered)
    assert not values["dft"]["energy_change_passed"]
    assert not values["total"]["energy_change_passed"]
    assert abs(values["total"]["energy_change_eV_cell"] - 2e-8) < 1e-13


def pending():
    with tempfile.TemporaryDirectory(dir=ROOT / "validation_audit", prefix="noise_fixture_") as temp:
        report = noise.analyze(Path(temp), ROOT / "scf_noise_queue.json", "deriv_final", "scfnoise_final")
        assert report["status"] == "pending"
        assert report["counts"]["pending"] == 3
        assert not report["passed"] and not report["coverage_complete"]


def incomplete_provenance():
    altered = copy.deepcopy(tight)
    altered["result"] = {}
    value = noise.run_provenance(altered)
    assert value["SCF_converged"] is None and value["runner_sha256"] is None


check("identical geometry and tenfold tighter SCF accepted", compare)
check("changed XC stencil rejected", lambda: changed_input("xc_stencil", 2))
check("changed explicit grid rejected", lambda: changed_input("gpts", [32, 20, 12]))
check("reversed strain sign rejected", lambda: changed_input("sign", 1))
check("insufficient SCF tightening rejected", lambda: changed_input(
    "scf_energy_eV_electron", base["input"]["scf_energy_eV_electron"]))
check("unrounded geometry mismatch rejected", geometry)
check("actual G-set mismatch rejected", basis)
check("energy change threshold applies to eV per cell", endpoint_values)
check("missing jobs stay pending", pending)
check("incomplete result metadata can still be reported as invalid", incomplete_provenance)
passed = all(result["passed"] for result in results)
a.write_json(ROOT / "validation_audit/scf_noise_analyzer_validation.json", {
    "quantum_calculations_started": False, "passed": passed,
    "analyzer_sha256": a.digest(ROOT / "scripts/analyze_scf_noise.py"),
    "derivative_helper_sha256": a.digest(ROOT / "scripts/analyze_derivatives.py"),
    "test_script_sha256": a.digest(__file__), "tests": results})
print(json.dumps({"passed": passed, "tests": results}, indent=2))
if not passed:
    raise SystemExit(1)
