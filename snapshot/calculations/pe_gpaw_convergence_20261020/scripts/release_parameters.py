"""Fail-closed release gate for this PE/GPAW experiment; never runs a calculation.

Default execution writes a new release_checks/<UTC>/ directory only.
--publish writes an immutable release only after every gate passes.
Only read-only analyzer functions are imported: their output/main functions
are never called. Historical unsuccessful scan groups are retained, not gated.
Run with the audited WSL Python; absolute Linux paths in evidence are intentional.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib
import importlib.metadata
import json
import math
from pathlib import Path
import shutil
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
METHODS = ("PBE-D3(BJ)", "optB88-vdW", "r2SCAN-D4")
ROLES = ("Tashiro", "GFN2")
PRECISION = {
    "L": {"ecut_eV": 1400, "kmesh": [3, 6, 9], "gpts": [96, 60, 36]},
    "M": {"ecut_eV": 1600, "kmesh": [4, 8, 12], "gpts": [112, 72, 44]},
    "H": {"ecut_eV": 1800, "kmesh": [4, 8, 12], "gpts": [112, 72, 44]},
}
FIXED = {"mpi_ranks": 4, "point_group": True, "xc_stencil": 4, "nbands": 16,
         "scf_energy_eV_electron": 1e-9, "scf_density_electron_electron": 1e-9,
         "scf_eigenstates_eV2_electron": 1e-11}
LIMITS = {"energy_meV_atom": .5, "force_vector_eV_A": .005,
          "stress_component_kbar": .1}
PAW = {
    "C": "5a769dd14bb5d462c57fadc4dd8d3b7e70428b90701eca0d20b6f0666495cdfd",
    "H": "43da6533cea08dec1ee7861140be503c9945863d6d8c7791e01e55e2161a88b8",
}
SOURCE = {
    "Tashiro": ("reference/POSCAR",
                "e54173a343c8a1bac20d048fdeaa99241986ee699e9f4591d65fc0985b76ef73"),
    "GFN2": ("reference/POSCAR_GFN2",
             "f3db63cbf993da085a1c5bb4f772eb30d2e4debf4ae80fb474e738db79599bd0"),
}
# These audited code versions are a trust boundary. A changed implementation
# requires review and an explicit revision here, never --force or auto-accept.
CODE = {
    "run_point": "cc4cd9d563a1d5eb69bc72346b096473ea785832487a26d984c441c72e2f57ae",
    "run_queue": "544c7814c0ba0f7f2deabbab79791335ad8ba78515c7326430ab0e9ba729d544",
    "run_relax": "b5009e3d38d422643a39a04a10a7f22caad954b4509d204d75e5e1824898bea6",
    "analyze_convergence": "d2541f7032612d1d2f92bd4dd1db27d37c2d1792757cf9b9e7d8e70d43a686b6",
    "analyze_relative_energies": "bca2274b7f91b2d54719d733a2d6dc201a625c0924bc95549dd3ba75b05df321",
    "analyze_relaxations": "a1f4f44266e8ce9f118cd68c68d7f7f6fc57e61c0b769a35ea11e558b383b694",
    "analyze_derivatives": "a9a2c14bb7c32f81d9b2f6573ab19267d1cd24cd723c0b08490fee148bc0d689",
    "analyze_scf_noise": "40e9936e045c8dc5c21affc337db415d315288550c41e96816be1a3cc09f4a35",
}
PBE_GROUPS = {
    "cutoff": ("7955d8aca0855caa0359", "cutoff_grid96_pbe_E1400_k4x8x12"),
    "kmesh": ("1e586d4b70a001501ffc", "kmesh_grid96_pbe_E1200_k3x6x9"),
    "grid": ("564c252252407663531b", "grid_pbe_E1400_k3x6x9_g96x60x36"),
}
BOHR = 0.5291772105638411


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False,
                      separators=(",", ":"), allow_nan=False)


def save_new(path, data):
    with Path(path).open("x", encoding="utf-8") as handle:
        json.dump(data, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")


def precision(spec):
    return {key: spec[key] for key in ("ecut_eV", "kmesh", "gpts")}


def scientific_software(value):
    return {key: value[key] for key in ("gpaw", "ase", "numpy", "scipy", "dftd3", "dftd4")}


def point_settings(spec, excluded=()):
    metadata = {
        "id", "stage", "structure", "checkpoint", "protocol_sha256", "queue_sha256",
        "runner_sha256_at_submission", "scan_parameters", "scan_label", "scan_index",
        "description", "notes", "geometry_role", "structure_role", "precision_label",
    }
    result = {"nbands": 16, "point_group": True, "scf_energy_eV_electron": 1e-8,
              "scf_density_electron_electron": 1e-8, "scf_eigenstates_eV2_electron": 1e-10}
    result.update({k: v for k, v in spec.items() if k not in metadata})
    return {k: v for k, v in result.items() if k not in excluded}


def same_geometry_and_settings(left, right, excluded):
    require(all(left[key] == right[key] for key in (
        "source_geometry_sha256", "natoms", "formula", "cell_A", "positions_A")),
        "Comparison does not have identical source and actual geometry")
    require(point_settings(left["input"], excluded) == point_settings(right["input"], excluded),
            "Nonexcluded scientific input settings differ")
    require(left["software"] == right["software"], "Comparison software differs")
    fields = ("symbol", "type", "setupname", "file_sha256")
    require([{k: d[k] for k in fields} for d in left["datasets"]]
            == [{k: d[k] for k in fields} for d in right["datasets"]],
            "Comparison ordered PAW datasets differ")


class Evidence:
    def __init__(self, root):
        self.root, self.files = root, {}

    def path(self, path):
        p = Path(path)
        return (p if p.is_absolute() else self.root / p).resolve()

    def track(self, path, expected=None, runtime=False):
        p = self.path(path)
        actual = digest(p)
        require(expected is None or actual == expected, f"SHA mismatch: {p}")
        prior = self.files.get(str(p))
        require(prior is None or prior["sha256"] == actual,
                f"Evidence changed during checking: {p}")
        self.files[str(p)] = {
            "path": str(p), "sha256": actual, "bytes": p.stat().st_size,
            "runtime_dependency": runtime or bool(prior and prior["runtime_dependency"])}
        return p

    def load(self, path, expected=None):
        p = self.track(path, expected)
        value = json.loads(p.read_text(encoding="utf-8-sig"))
        self.track(p, self.files[str(p)]["sha256"])
        return value

    def unchanged(self):
        for item in list(self.files.values()):
            require(digest(item["path"]) == item["sha256"],
                    f"Evidence changed; refresh and rerun: {item['path']}")


class Gate:
    def __init__(self, root, candidate, reference):
        self.root, self.ev = root, Evidence(root)
        self.candidate, self.reference = candidate, reference
        self.checks, self.details, self.modules = [], {}, {}
        self.points = {}
        self.relative = None
        self.convergence_saved, self.convergence_fresh = {}, {}

    def check(self, name, function):
        try:
            detail = function()
            self.details[name] = detail
            self.checks.append({"gate": name, "status": "passed"})
        except (ValueError, KeyError, TypeError, OSError, ImportError) as exc:
            self.checks.append({"gate": name, "status": "blocked", "reason": str(exc)})

    def modules_and_environment(self):
        for name, expected in CODE.items():
            self.ev.track(f"scripts/{name}.py", expected,
                          runtime=name in ("run_point", "run_queue", "run_relax"))
        self.ev.track(__file__, runtime=True)
        sys.path.insert(0, str(self.root / "scripts"))
        for name in CODE:
            if name.startswith("analyze_"):
                self.modules[name] = importlib.import_module(name)
        env = self.ev.load("environment/environment_verified.json")
        require(Path(sys.executable).resolve() == Path(env["python"]).resolve(),
                "Use the audited WSL Python environment")
        require(env["libvdwxc_stress"] is True, "libvdwxc stress unavailable")
        for name, version in env["versions"].items():
            require(importlib.metadata.version(name) == version,
                    f"Installed version changed: {name}")
        for item in env["important_files"] + env["setup_files"]:
            self.ev.track(item["path"], item["sha256"], runtime=True)
        # The editable source is verified again against the official archive
        # by a read-only file comparison, not by rebuilding/importing GPAW.
        build = self.ev.load("environment/gpaw_build_audit.json")
        self.ev.track(build["official_source_archive"], build["source_sha256"])
        self.ev.track("environment/conda_explicit_lock.txt")
        self.ev.track("environment/conda_packages.json")
        self.ev.track("environment/downloaded_package_hashes.json")
        self.ev.track("environment/siteconfig.py", build["configuration_sha256"])
        self.ev.track("environment/runtime_environment.sh", runtime=True)
        import tarfile
        source_root = Path(build["source_directory"])
        compared = 0
        with tarfile.open(build["official_source_archive"]) as archive:
            for member in archive.getmembers():
                if not member.isfile():
                    continue
                relative = Path(*Path(member.name).parts[1:])
                if relative.suffix not in (".py", ".c", ".h", ".cpp", ".hpp"):
                    continue
                actual = source_root / relative
                expected = hashlib.sha256(archive.extractfile(member).read()).hexdigest()
                self.ev.track(actual, expected, runtime=True)
                compared += 1
        require(compared > 2000, "Official source verification coverage unexpectedly small")
        self.method = self.ev.load("method_audit/method_specification.json")
        require(self.method["paw_sha256"] == PAW, "PAW method specification changed")
        require(set(self.method["methods"]) == set(METHODS), "Method set changed")
        binding = self.ev.load("method_audit/gpaw_binding_check.json")
        self.ev.track(binding["extension_path"], binding["extension_sha256"], runtime=True)
        for role, (path, expected) in SOURCE.items():
            self.ev.track(path, expected, runtime=True)
        self.protocol = self.ev.load("protocol.json")
        require(self.protocol["schema_version"] == 5, "Protocol revision needs review")
        require(set(self.protocol["methods"]) == set(METHODS), "Protocol methods changed")
        numerical = self.protocol["numerical_convergence"]
        require(numerical["energy_tolerance_meV_atom"] == .5
                and numerical["force_max_vector_difference_eV_A"] == .005
                and numerical["stress_max_component_difference_kbar"] == .1
                and numerical["relative_structure_energy_tolerance_meV_atom"] == .5,
                "Protocol tolerances changed")
        self.software = scientific_software(env["versions"])
        return {"versions": env["versions"], "important_binary_files": env["important_files"],
                "PAW": env["setup_files"], "official_scientific_source_files_checked": compared,
                "method_definition": self.method, "engine": self.protocol["engine"],
                "note": "Historical environment interface status is not a quantum validation gate."}

    def point(self, path, expected=None):
        path = self.ev.path(path)
        result = self.ev.load(path, expected)
        folder = path.parent
        require(result["status"] == "passed_single_point" and result["SCF_converged"] is True,
                f"Unsuccessful point: {folder.name}")
        process = self.ev.load(folder / "process.json")
        require(process.get("exit_code") == 0, f"Point process did not finish successfully: {folder.name}")
        spec = self.ev.load(folder / "input.json", result["input_sha256"])
        require(spec == result["input"], f"Input contents differ: {folder.name}")
        require(result["runner_sha256"] == CODE["run_point"]
                and spec["runner_sha256_at_submission"] == CODE["run_point"],
                f"Unaudited point runner: {folder.name}")
        self.ev.track(folder / "run_point_snapshot.py", CODE["run_point"])
        self.ev.track(folder / "protocol_snapshot.json", spec["protocol_sha256"])
        self.ev.track(spec["structure"], result["source_geometry_sha256"])
        require(result["natoms"] == 12 and result["formula"] == "C4H8", "Unsupported formula")
        require(result["mpi_world_size"] == spec.get("mpi_ranks", 4), "Actual MPI mismatch")
        require(scientific_software(result["software"]) == self.software, "Software differs")
        require(result["units"]["energy"] == "eV/cell", "Unexpected energy units")
        counts = {"C": 0, "H": 0}
        for item in result["datasets"]:
            symbol = item["symbol"]
            require(symbol in PAW and item["file_sha256"] == PAW[symbol]
                    and item["type"] == "paw" and item["setupname"] == "PBE",
                    "PAW dataset differs")
            self.ev.track(item["filename"], PAW[symbol], runtime=True)
            counts[symbol] += 1
        require(counts == {"C": 4, "H": 8}, "PAW element multiplicities differ")
        total, dft, disp = result["total"], result["dft"], result["dispersion"]
        require(abs(total["energy_free_eV"] - dft["energy_free_eV"] - disp["energy_eV"]) < 1e-9,
                "Total free energy does not equal DFT plus dispersion")
        for key, shape in (("forces_eV_A", (12, 3)), ("stress_eV_A3", (6,))):
            import numpy as np
            a, b, c = (np.asarray(item[key], dtype=float) for item in (total, dft, disp))
            require(a.shape == b.shape == c.shape == shape
                    and np.isfinite(a).all() and np.isfinite(b).all() and np.isfinite(c).all()
                    and float(np.max(np.abs(a - b - c))) < 1e-9,
                    f"Invalid or inconsistent {key}")
        method = spec["method"]
        params = disp.get("calculator_parameters")
        if method in ("PBE-D3(BJ)", "r2SCAN-D4"):
            require(params["params_tweaks"]["atm"] is True, "ATM is not enabled")
            target = self.method["methods"][method]
            for key, radius in target["cutoff_bohr"].items():
                require(math.isclose(params["realspace_cutoff"][key], radius * BOHR,
                                     rel_tol=1e-12, abs_tol=1e-12), "Dispersion cutoff differs")
            if method == "PBE-D3(BJ)":
                require(params["params_tweaks"]["method"] == "PBE"
                        and params["damping"] == "d3bj", "D3 parameter set differs")
            else:
                require(params["method"] == "r2scan" and params.get("model", "d4") == "d4",
                        "D4 is not r2SCAN-D4")
        else:
            require(method in ("PBE", "optB88-vdW") and disp["energy_eV"] == 0,
                    "Unexpected additive dispersion")
        for name in ("input_structure.extxyz", "gpaw.txt", "scf_iterations.jsonl",
                     "stdout.log", "stderr.log", "process.json"):
            p = folder / name
            if p.exists():
                self.ev.track(p)
        basis = list(folder.glob("basis_rank*.json"))
        require(len(basis) == result["mpi_world_size"], "Missing rank basis metadata")
        for p in basis:
            self.ev.track(p)
        self.points[str(path)] = result
        return result

    def validated_reports(self):
        derivative = self.modules["analyze_derivatives"]
        saved = self.ev.load("validation_audit/derivative_results_deriv_final.json")
        fresh, _ = derivative.analyze(self.root, "deriv_final", None, "original")
        for report in (saved, fresh):
            require(report["passed"] is True and report["coverage_complete"] is True
                    and report["status"] == "representative_derivative_checks_passed"
                    and report["counts"]["passed"] == 27
                    and report["selected_successful_points"] == 111,
                    "Derivative coverage incomplete or failed")
            require(report["analyzer_sha256"] == CODE["analyze_derivatives"],
                    "Derivative report code mismatch")
        self.bind_selected(saved, fresh, derivative=True)
        queue = self.ev.track("scf_noise_queue.json")
        noise = self.modules["analyze_scf_noise"]
        saved_noise = self.ev.load("validation_audit/scf_noise_validation.json")
        fresh_noise = noise.analyze(self.root, queue, "deriv_final", "scfnoise_final")
        for report in (saved_noise, fresh_noise):
            require(report["passed"] is True and report["coverage_complete"] is True
                    and report["status"] == "SCF_noise_checks_passed"
                    and report["counts"]["passed"] == 3, "SCF noise check incomplete/failed")
            require(report["analyzer_sha256"] == CODE["analyze_scf_noise"],
                    "SCF analyzer mismatch")
        self.bind_selected(saved_noise, fresh_noise, derivative=False)
        return {"derivative_counts": fresh["counts"], "noise_counts": fresh_noise["counts"],
                "coverage": "3 representative forces and all 6 stress components per method"}

    def bind_selected(self, saved, fresh, derivative):
        require(set(saved["selected_runs"]) == set(fresh["selected_runs"]), "Stale report run coverage")
        for key, value in fresh["selected_runs"].items():
            old = saved["selected_runs"][key]
            require(old["result_sha256"] == value["result_sha256"]
                    and old["input_sha256"] == value["input_sha256"], "Stale report point SHA")
            folder = (Path(value["folder"]) if derivative else
                      self.root / "runs" / value["selected_run"])
            self.point(folder / "result.json", value["result_sha256"])

    def fixed_geometry(self):
        saved = self.ev.load("tables/convergence_summary.json")
        require(saved["analyzer_sha256"] == CODE["analyze_convergence"],
                "Fixed-geometry analyzer mismatch")
        module = self.modules["analyze_convergence"]
        groups, _ = module.collect(self.root)
        by_id = {g["group_id"]: g for g in groups}
        stored = {g["group_id"]: g for g in saved["groups"]}
        accepted = {}
        for name, (group_id, candidate) in PBE_GROUPS.items():
            group, archived = by_id[group_id], stored[group_id]
            require(group["status"] == archived["status"] == "converged_fixed_geometry"
                    and group["method"] == "PBE" and group["terminal_is_provisional_reference"],
                    f"PBE {name} group incomplete/failed")
            require(group["recommended_run"] == archived["recommended_run"] == candidate,
                    f"PBE {name} candidate changed; review required")
            old = {p["run_id"]: p["result_sha256"] for p in archived["points"]}
            for point in group["points"]:
                require(old.get(point["run_id"]) == point["result_sha256"], "Stale convergence table")
                self.point(point["result_path"], point["result_sha256"])
            accepted[name] = group
        self.convergence_saved, self.convergence_fresh = stored, by_id
        return {"selected_groups": accepted,
                "historical_failed_groups": "Retained and excluded from this explicit release selection"}

    def relative_energies(self):
        saved = self.ev.load("tables/relative_structure_energies.json")
        module = self.modules["analyze_relative_energies"]
        for path in ("transfer_tashiro_queue.json", "transfer_gfn2_queue.json"):
            self.ev.track(path)
        fresh = module.analyze(self.root)
        self.relative = fresh  # Allows useful independent transfer checks even while pending.
        for record in fresh["records"]:
            self.point(record["result_path"], record["result_sha256"])
        for report in (saved, fresh):
            require(report["passed"] is True and report["coverage_complete"] is True
                    and report["status"] == "relative_structure_energy_checks_passed",
                    f"Relative-structure check {report['status']}; "
                    f"pending={len(report['pending_points'])}, pairs={report['complete_pairs']}")
            require(set(report["required_methods"]) == set(METHODS)
                    and report["analyzer_sha256"] == CODE["analyze_relative_energies"],
                    "Relative-energy method coverage/code mismatch")
            require(not any(report[k] for k in (
                "missing_or_unpassed_required_methods", "pending_points",
                "unmatched_pairs", "rejected_records", "unreadable_inputs")),
                "Relative-energy report contains unfinished or rejected evidence")
            require(report["families"] and all(f["passed"] is True for f in report["families"]),
                    "An included relative-energy family failed")
        old = {r["run_id"]: r["result_sha256"] for r in saved["records"]}
        require(old == {r["run_id"]: r["result_sha256"] for r in fresh["records"]},
                "Relative-energy summary must be refreshed")
        by_pair = {p["pair_id"]: p for p in fresh["pairs"]}
        accepted = {}
        for method in METHODS:
            families = [f for f in fresh["families"] if f["method"] == method]
            require(len(families) == 1, "Multiple provenance families require explicit release review")
            family = families[0]
            candidates = [by_pair[i] for i in family["accepted_candidate_pair_ids"]
                          if precision(by_pair[i]["precision"]) == PRECISION["L"]]
            require(len(candidates) == 1
                    and precision(by_pair[family["highest_pair_id"]]["precision"]) == PRECISION["H"],
                    f"L/H relative-energy acceptance missing: {method}")
            accepted[method] = {"candidate": candidates[0]["pair_id"],
                                "reference": family["highest_pair_id"],
                                "identity": family["identity"]}
        return {"accepted": accepted, "complete_pairs": fresh["complete_pairs"]}

    def transfers(self):
        require(self.relative is not None, "Relative records unavailable")
        module = self.modules["analyze_convergence"]
        accepted = {}
        for method in METHODS:
            for role in ROLES:
                records = []
                for level, target in PRECISION.items():
                    found = [r for r in self.relative["records"]
                             if r["method"] == method and r["role"] == role
                             and precision(r["precision"]) == target]
                    require(len(found) == 1, f"Missing or ambiguous {method}/{role}/{level}")
                    r = found[0]
                    result = self.point(r["result_path"], r["result_sha256"])
                    spec = result["input"]
                    for key, expected in FIXED.items():
                        require(spec.get(key, 16 if key == "nbands" else None) == expected,
                                f"Release fixed parameter differs: {key}")
                    records.append(module.extract_record(Path(r["result_path"]).parent, spec, result))
                require(len({r["group_id"] for r in records}) == 1, "Transfer inputs incompatible")
                group = module.analyze_group(records)
                require(group["status"] == "converged_fixed_geometry"
                        and group["recommended_run"] == records[0]["run_id"],
                        f"L/M/H E/F/stress transfer failed: {method}/{role}")
                archived = self.convergence_saved.get(group["group_id"])
                require(archived and archived["status"] == "converged_fixed_geometry",
                        f"Refresh persisted convergence summary: {method}/{role}")
                old = {p["run_id"]: p["result_sha256"] for p in archived["points"]}
                require(all(old.get(p["run_id"]) == p["result_sha256"] for p in records),
                        "Persisted transfer table does not bind actual results")
                accepted[f"{method}/{role}"] = group
        return accepted

    def bridge_and_parallel(self):
        parallel = self.ev.load("validation_audit/parallel_invariance.json")
        require(parallel["status"] == "passed" and all(parallel["threshold_pass"].values()),
                "4/8 MPI check did not pass")
        self.ev.track("validation_audit/audit_parallel_invariance.py", parallel["audit_script_sha256"])
        results = [self.point(x["result_path"], x["result_sha256"]) for x in parallel["runs"]]
        require({x["mpi_world_size"] for x in results} == {4, 8}, "MPI coverage is not 4/8")
        same_geometry_and_settings(results[0], results[1], {"mpi_ranks"})
        require(all(r["input"]["method"] == "PBE" for r in results),
                "Parallel invariance is no longer the audited PBE comparison")
        p_diff = self.property_delta(results[0]["total"], results[1]["total"])
        require(p_diff["energy_meV_atom"] * 12 / 1000 <= 1e-6
                and p_diff["force_vector_eV_A"] <= 1e-6
                and p_diff["stress_component_kbar"] <= 1e-5, "MPI delta no longer passes")
        bridge = self.ev.load("validation_audit/pbe_reference_grid_bridge.json")
        require(bridge["status"] == "passed" and not bridge["evidence_errors"]
                and not bridge["identity_mismatches"]
                and all(bridge["threshold_pass"].values()), "Reference grid bridge failed")
        self.ev.track("validation_audit/audit_pbe_reference_grid_bridge.py",
                      bridge["audit_script_sha256"])
        b_results = [self.point(x["result_path"], x["result_sha256"]) for x in bridge["runs"]]
        require({r["input"]["method"] for r in b_results} == {"PBE", "PBE-D3(BJ)"},
                "Bridge is not the intended common PBE DFT comparison")
        require(all(r["input"]["ecut_eV"] == 1800 and r["input"]["kmesh"] == [4, 8, 12]
                    and r["source_geometry_sha256"] == SOURCE["Tashiro"][1] for r in b_results),
                "Bridge geometry/cutoff/k differs")
        require({tuple(r["input"]["gpts"]) for r in b_results}
                == {(96, 60, 36), (112, 72, 44)}, "Bridge grid coverage changed")
        same_geometry_and_settings(b_results[0], b_results[1], {"mpi_ranks", "gpts", "method"})
        require(all(r["basis"]["xc_name"] == "PBE" for r in b_results),
                "Bridge is not the same PBE DFT functional")
        b_diff = self.property_delta(b_results[0]["dft"], b_results[1]["dft"])
        require(all(b_diff[k] <= LIMITS[k] for k in LIMITS), "Bridge actual differences failed")
        decision = self.ev.load("relaxations/cell_reference/grid_selection_decision.json")
        self.ev.track("relaxations/cell_reference/input.json", decision["current_input_sha256"])
        self.ev.track("relaxations/cell_reference/input_before_grid_bridge.json")
        return {"parallel": parallel, "parallel_recomputed": p_diff,
                "reference_grid_bridge": bridge, "bridge_recomputed_DFT": b_diff,
                "selection_decision": decision}

    @staticmethod
    def property_delta(left, right):
        import numpy as np
        return {
            "energy_meV_atom": abs(left["energy_free_eV"] - right["energy_free_eV"]) * 1000 / 12,
            "force_vector_eV_A": float(np.linalg.norm(
                np.asarray(left["forces_eV_A"]) - right["forces_eV_A"], axis=1).max()),
            "stress_component_kbar": float(np.abs(
                np.asarray(left["stress_kbar"]) - right["stress_kbar"]).max())}

    def track_relaxation_lineage(self, run, module):
        """Freeze every verified ancestor, including for a still-pending child."""
        folder = Path(run["directory"])
        evidence = self.ev.load(folder / "provenance.json")
        spec = self.ev.load(folder / "input.json", evidence["input_sha256"])
        require(spec == evidence["input"], "Relaxation input differs from its provenance")
        lineage = module.source_lineage(folder, evidence)
        if run.get("source_lineage") is not None:
            require(canonical(run["source_lineage"]) == canonical(lineage),
                    "Relaxation lineage changed after result inspection")
        require(lineage["root_source_geometry_sha256"] == SOURCE["Tashiro"][1],
                "Cell validation does not descend from the frozen Tashiro structure")
        require(self.ev.path(lineage["root_source_geometry_path"])
                == self.ev.path(SOURCE["Tashiro"][0]),
                "Cell validation root source path differs from this experiment")
        self.ev.track(lineage["root_source_geometry_path"], SOURCE["Tashiro"][1])
        source = Path(evidence["source_geometry_path"])
        self.ev.track(source if source.is_absolute() else folder / source,
                      evidence["source_geometry_sha256"])
        for continuation in lineage["continuations"]:
            manifest_path = self.ev.path(continuation["manifest_path"])
            manifest = self.ev.load(manifest_path, continuation["manifest_sha256"])
            parent = self.ev.path(manifest["parent_directory"])
            for name, expected in manifest["parent_evidence_sha256"].items():
                path = self.ev.path(parent / name)
                require(path.parent == parent, "Continuation parent evidence escapes its attempt")
                self.ev.track(path, expected)
            self.ev.track(manifest_path.parent / manifest["source_file"],
                          manifest["source_sha256"])
        return lineage

    def relaxations(self):
        saved = self.ev.load("tables/relaxation_energies.json")
        require(saved["analyzer_sha256"] == CODE["analyze_relaxations"], "Relaxation analyzer changed")
        module = self.modules["analyze_relaxations"]
        runs = [module.inspect_run(self.root / "relaxations" / name)
                for name in (self.candidate, self.reference)]
        self.details["relaxation_source_lineages"] = {
            run["run_id"]: self.track_relaxation_lineage(run, module) for run in runs}
        require(all(r["acceptance"] == "verified_relaxation" for r in runs),
                "Both selected relaxation attempts must succeed independently")
        cross_path = self.root / "relaxations/cross_precision_verification.json"
        cross = self.ev.load(cross_path)
        cross["_report_directory"] = str(cross_path.parent)
        fresh = module.compare_runs(runs, self.candidate, self.reference, cross)
        require(saved["status"] == fresh["status"] == "passed"
                and fresh["cell_length_status"] == "passed"
                and fresh["cross_validation"]["status"] == "verified",
                "Cell comparison or common-precision cold-start evidence not passed")
        require(saved["comparison"]["candidate_run"] == self.candidate
                and saved["comparison"]["reference_run"] == self.reference,
                "Saved relaxation pair differs")
        candidate, reference = runs
        require(precision(candidate["precision"]) == PRECISION["L"], "Candidate cell precision differs from L")
        high = {**PRECISION["H"], "gpts": [96, 60, 36]}
        require(precision(reference["precision"]) == high, "Reference differs from approved grid bridge")
        for run in runs:
            require(run["comparison_identity"]["method"] == "PBE-D3(BJ)", "Unexpected cell method")
            require(run["point_runner_sha256"] == CODE["run_point"]
                    and run["runner_sha256"] == CODE["run_relax"], "Relaxation runner mismatch")
            for key, expected in FIXED.items():
                value = run["input"].get(key, 16 if key == "nbands" else None)
                require(value == expected, f"Cell fixed parameter differs from transfer: {key}")
            folder = Path(run["directory"])
            process = self.ev.load(folder / "process.json")
            require(process.get("exit_code") == 0,
                    f"Relaxation process has not exited successfully: {run['run_id']}")
            for path in folder.iterdir():
                if path.is_file() and path.suffix != ".gpw":
                    self.ev.track(path)
            self.ev.track(folder / "result.json", run["result_sha256"])
            stored = next(r for r in saved["runs"] if r["run_id"] == run["run_id"])
            require(stored["result_sha256"] == run["result_sha256"], "Stale relaxation summary")
        for item in fresh["cross_validation"]["checks"]:
            self.point(item["result_path"], item["point_result_sha256"])
        return {"comparison": fresh, "runs": [{k: v for k, v in r.items() if k != "rows"} for r in runs]}

    def run(self):
        self.check("environment_method_and_code", self.modules_and_environment)
        if self.checks[-1]["status"] == "passed":
            self.check("derivatives_and_SCF_noise", self.validated_reports)
            self.check("PBE_selected_cutoff_k_grid", self.fixed_geometry)
            self.check("relative_structure_energies", self.relative_energies)
            self.check("three_methods_two_geometries_LMH_E_F_stress", self.transfers)
            self.check("MPI_and_reference_grid_bridge", self.bridge_and_parallel)
            self.check("two_relaxations_and_cold_start", self.relaxations)
        self.check("evidence_unchanged_during_check", self.ev.unchanged)
        passed = len(self.checks) == 8 and all(c["status"] == "passed" for c in self.checks)
        return {"schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
                "status": "ready_for_explicit_publish" if passed else "release_refused",
                "passed": passed, "quantum_calculations_started": False,
                "checks": self.checks, "details": self.details,
                "candidate_precision": PRECISION["L"], "fixed_parameters": FIXED,
                "required_methods": METHODS, "required_structures": SOURCE,
                "release_script_sha256": digest(__file__),
                "source_manifest": list(self.ev.files.values()),
                "scope": "Validated two frozen C4H8 geometries, and a PBE-D3(BJ) Pnam cell comparison. "
                         "Reuse on PE candidates of similar composition/size requires sampled geometry "
                         "transfer checks; new elements, substantially changed cells, other phases, "
                         "or untested methods are not automatically validated."}


def publish(gate, report, destination, release_id):
    require(report["passed"] is True, "Refusing release: one or more gates incomplete")
    require(not destination.exists(), "Release directory already exists")
    gate.ev.unchanged()
    destination.mkdir(parents=True, exist_ok=False)
    archive = destination / "evidence"
    archive.mkdir()
    manifest = []
    for item in report["source_manifest"]:
        row = dict(item)
        # Freeze evidence documents and scientific scripts. Large runtime
        # binaries remain hash-bound in their installed/archive locations.
        path = Path(item["path"])
        if path.suffix in (".json", ".jsonl", ".py", ".md", ".txt", ".sh", ".csv") \
                or path.name.startswith("POSCAR") or path.suffix in (".extxyz", ".cif", ".traj"):
            target = archive / (item["sha256"] + "_" + path.name)
            if not target.exists():
                shutil.copyfile(path, target)
            require(digest(target) == item["sha256"], "Evidence copy changed")
            row["archived_path"] = str(target.relative_to(destination))
        manifest.append(row)
    save_new(destination / "release_check.json", report)
    save_new(destination / "evidence_manifest.json", manifest)
    templates = {}
    for method, slug in zip(METHODS, ("pbe_d3bj", "optb88_vdw", "r2scan_d4")):
        points = []
        for role, (source, _) in SOURCE.items():
            points.append({"id": f"released_{release_id}_{slug}_{role.lower()}",
                           "method": method, "structure": str(gate.root / source)})
        queue = {
            "stop_on_failure": True,
            "setup_path": str(Path(report["details"]["environment_method_and_code"]["PAW"][0]["path"]).parent),
            "defaults": {**FIXED, **PRECISION["L"], "stage": "released_reference_replay",
                         "checkpoint": True},
            "points": points,
            "release_note": "Runnable reference examples; no calculation started by release. "
                            "Changing geometry requires sampling at a higher tested precision; "
                            "do not infer every PE candidate is validated.",
        }
        name = f"queue_{slug}.json"
        save_new(destination / name, queue)
        templates[name] = digest(destination / name)
    gate.ev.unchanged()
    lock = {
        "schema_version": 1, "status": "locked_after_all_required_checks_passed",
        "created_utc": datetime.now(timezone.utc).isoformat(), "release_id": release_id,
        "scope": report["scope"], "methods": METHODS, "candidate_precision": PRECISION["L"],
        "fixed_parameters": FIXED, "physics": gate.method,
        "environment": report["details"]["environment_method_and_code"],
        "code_sha256": CODE, "source_geometry": SOURCE,
        "checks": report["checks"], "template_sha256": templates,
        "evidence_manifest_sha256": digest(destination / "evidence_manifest.json"),
        "release_check_sha256": digest(destination / "release_check.json"),
        "runtime_instructions": {
            "verify": f"{sys.executable} {Path(__file__).resolve()} --verify-release {destination}",
            "run_queue": f"{sys.executable} {gate.root / 'scripts/run_queue.py'} <queue_PATH>",
            "required_order": "Verify release first; only then explicitly invoke run_queue. "
                              "This release tool never launches the queue."}}
    # Written last: absence of this file means an incomplete staging directory.
    save_new(destination / "parameters_lock.json", lock)


def verify_release(path):
    path = path.resolve()
    lock = json.loads((path / "parameters_lock.json").read_text())
    require(lock["status"] == "locked_after_all_required_checks_passed", "Not a completed release")
    require(digest(path / "evidence_manifest.json") == lock["evidence_manifest_sha256"]
            and digest(path / "release_check.json") == lock["release_check_sha256"],
            "Release evidence manifest/check changed")
    for name, expected in lock["environment"]["versions"].items():
        require(importlib.metadata.version(name) == expected,
                f"Live installed package version changed: {name}")
    for name, expected in lock["template_sha256"].items():
        require(digest(path / name) == expected, "Queue template changed")
    manifest = json.loads((path / "evidence_manifest.json").read_text())
    for item in manifest:
        if "archived_path" in item:
            require(digest(path / item["archived_path"]) == item["sha256"], "Archived evidence changed")
        if item["runtime_dependency"]:
            require(digest(item["path"]) == item["sha256"], "Live runtime/geometry changed")
    return {"status": "release_files_and_runtime_verified", "quantum_calculations_started": False,
            "release": str(path)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--candidate", default="cell_candidate")
    parser.add_argument("--reference", default="cell_reference")
    parser.add_argument("--publish", action="store_true")
    parser.add_argument("--verify-release", type=Path)
    args = parser.parse_args()
    if args.verify_release:
        require(not args.publish, "Do not combine verify and publish")
        print(json.dumps(verify_release(args.verify_release), ensure_ascii=False))
        return
    root = args.root.resolve()
    require(root == ROOT, "This audited release policy is scoped to its own project directory")
    for value in (args.candidate, args.reference):
        require(Path(value).name == value and value not in (".", ".."), "Invalid attempt name")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    gate = Gate(root, args.candidate, args.reference)
    report = gate.run()
    output = root / "release_checks" / stamp
    output.mkdir(parents=True, exist_ok=False)
    save_new(output / "check.json", report)
    save_new(output / "release_plan.json", {
        "status": "plan_only_not_a_parameter_lock",
        "candidate": PRECISION["L"], "fixed_parameters": FIXED,
        "transfer_levels": PRECISION, "PBE_selected_groups": PBE_GROUPS,
        "required_methods": METHODS, "reference_files": SOURCE,
        "relaxation_pair": [args.candidate, args.reference],
        "proposed_output": str(root / "releases" / stamp), "scope": report["scope"],
        "publish_condition": "All eight gates pass in the same invocation, then explicit --publish.",
        "incomplete_evidence": [c for c in report["checks"] if c["status"] != "passed"]})
    if args.publish:
        publish(gate, report, root / "releases" / stamp, stamp)
    print(json.dumps({"status": report["status"], "passed": report["passed"],
                      "published": bool(args.publish and report["passed"]),
                      "check_report": str(output / "check.json"), "checks": report["checks"]},
                     ensure_ascii=False, indent=2))
    if not report["passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
