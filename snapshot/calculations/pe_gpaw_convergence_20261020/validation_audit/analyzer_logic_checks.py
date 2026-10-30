"""Cheap checks of the read-only analyzer. Never imports a quantum engine."""
import copy
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import tempfile

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("analyzer", ROOT / "scripts/analyze_derivatives.py")
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)
RESULTS = []


def check(name, function):
    try:
        function()
        RESULTS.append({"test": name, "passed": True})
    except Exception as exc:
        RESULTS.append({"test": name, "passed": False, "error": repr(exc)})


def must_raise(function):
    try:
        function()
    except ValueError:
        return
    raise AssertionError("Expected rejection")


def selection():
    values = [{"state": status, "attempt": number} for number, status in
              enumerate(["success", "failed", "success", "pending"], 1)]
    assert a.select({"point": values}, "point")["attempt"] == 3


def strain():
    stress = a.np.array([[1., 4., 5.], [4., 2., 6.], [5., 6., 3.]])
    for index, (i, j) in enumerate(a.PAIRS):
        eps = a.strain_matrix(index, 1e-5)
        assert abs(a.np.sum(stress * eps) / 1e-5 - stress[i, j]) < 1e-12


attempts, _ = a.discover(ROOT, "deriv_dense")
center = a.select(attempts, "deriv_dense_pbe_d3bj_center_E600")
endpoint = a.select(attempts, "deriv_dense_pbe_d3bj_S1_h1e-05_m")


def wrong_geometry():
    changed = copy.deepcopy(center)
    changed["result"]["positions_A"][0][0] += .01
    must_raise(lambda: a.validate_run(changed))


def actual_basis_change():
    checks = {run["folder"].name: a.validate_run(run) for run in (center, endpoint)}
    result = a.compare_endpoint(center, endpoint, checks, "stress", 1, 1e-5, -1)
    assert result["wavefunction_basis_unchanged"] is True
    assert result["density_hamiltonian_basis_unchanged"] is False
    assert result["basis_unchanged"] is False


def absent_basis_metadata():
    checks = {run["folder"].name: a.validate_run(run) for run in (center, endpoint)}
    checks[center["folder"].name]["basis"]["density_hamiltonian_basis_recorded"] = False
    result = a.compare_endpoint(center, endpoint, checks, "stress", 1, 1e-5, -1)
    assert result["density_hamiltonian_basis_unchanged"] is None
    assert result["density_hamiltonian_basis_recorded"] is False


def missing_rank():
    with tempfile.TemporaryDirectory(dir=ROOT / "validation_audit", prefix="analysis_fixture_") as temp:
        folder = Path(temp)
        shutil.copy(center["folder"] / "basis_rank0000.json", folder)
        changed = copy.deepcopy(center)
        changed["folder"] = folder
        must_raise(lambda: a.merged_basis(changed))


def pending():
    with tempfile.TemporaryDirectory(dir=ROOT / "validation_audit", prefix="analysis_fixture_") as temp:
        report, _ = a.analyze(Path(temp), "deriv_future", None)
        assert report["status"] == "pending"
        assert not report["passed"]
        assert report["counts"]["pending"] == 27


def auto_center():
    changed = copy.deepcopy(center)
    changed["input"]["ecut_eV"] = 600.123
    changed["id"] = "deriv_future_pbe_d3bj_center_E600.123"
    changed["input"]["id"] = changed["id"]
    centers, _ = a.expected_points(ROOT, "deriv_future", None, {changed["id"]: [changed]})
    assert changed["id"] in centers


def cn60():
    with tempfile.TemporaryDirectory(dir=ROOT / "validation_audit", prefix="analysis_fixture_") as temp:
        folder = Path(temp) / "method_audit/derivative_d3_cn60"
        folder.mkdir(parents=True)
        revision = copy.deepcopy(center["result"])
        atoms = a.read(center["folder"] / "input_structure.extxyz")
        geometry = {"atomic_numbers": atoms.numbers.tolist(),
                    "positions_A": revision["positions_A"], "cell_A": revision["cell_A"],
                    "pbc": atoms.pbc.tolist()}
        revision.update(
            status="passed_dispersion_recombination", revision_id="D3_CN60_smooth005bohr",
            DFT_recomputed=False,
            source_result_sha256=a.digest(center["folder"] / "result.json"),
            same_geometry_sha256=a.digest(center["folder"] / "input_structure.extxyz"),
            exact_geometry_sha256=a.object_digest(geometry), revision_script_sha256="a" * 64)
        revision["dispersion"]["calculator_parameters"]["realspace_cutoff"]["cn"] = 60 * a.Bohr
        revision["dispersion"]["energy_eV"] += .01
        revision["total"]["energy_free_eV"] += .01
        path = folder / (center["folder"].name + ".json")
        a.write_json(path, revision)
        overlay = a.cn60_view(Path(temp), center)
        assert overlay["result"]["total"]["energy_free_eV"] == revision["total"]["energy_free_eV"]
        assert overlay["result"]["dft"] == center["result"]["dft"]
        revision["source_result_sha256"] = "0" * 64
        a.write_json(path, revision)
        must_raise(lambda: a.cn60_view(Path(temp), center))


check("latest successful attempt skips newer failed/pending", selection)
check("engineering shear has half offdiagonal and no factor-two error", strain)
check("actual completed dense center validates", lambda: a.validate_run(center))
check("mutated result geometry rejected", wrong_geometry)
check("density shell crossing detected despite unchanged wavefunction basis", actual_basis_change)
check("missing metadata remains unknown rather than falsely claiming a shell crossing", absent_basis_metadata)
check("missing MPI rank basis rejected", missing_rank)
check("absent jobs remain pending and never pass", pending)
check("decimal cutoff discovered from center input", auto_center)
check("CN60 overlay preserves DFT, no double addition, rejects wrong source hash", cn60)
passed = all(result["passed"] for result in RESULTS)
a.write_json(ROOT / "validation_audit/analyzer_validation.json", {
    "quantum_calculations_started": False,
    "analyzer_sha256": a.digest(ROOT / "scripts/analyze_derivatives.py"),
    "test_script_sha256": a.digest(__file__), "passed": passed, "tests": RESULTS})
print(json.dumps({"passed": passed, "tests": RESULTS}, indent=2))
if not passed:
    raise SystemExit(1)
