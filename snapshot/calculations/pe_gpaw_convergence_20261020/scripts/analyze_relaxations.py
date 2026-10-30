"""Read relaxation evidence and summarize energies; never starts a calculation.

Default: inspect every relaxations/*/result.json and evaluations.jsonl.
Use --candidate/--reference when more than two accepted attempts exist.
The optional cross-validation report must bind both result and final-geometry
hashes. Missing cross-validation remains pending even when cell lengths agree.
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
import tempfile


FORCE_LIMIT = 0.002
STRESS_LIMIT = 0.05
CELL_RELATIVE_LIMIT = 0.005
KBAR_PER_EV_A3 = 1602.1766208
PRECISION_KEYS = ("ecut_eV", "kmesh", "gpts", "xc_stencil")
IDENTITY = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
# Independently inspected run_point.py: each invocation constructs a new GPAW
# object with no restart/density input. New runner versions need a new audit.
COLD_START_RUNNER_SHA256 = "cc4cd9d563a1d5eb69bc72346b096473ea785832487a26d984c441c72e2f57ae"
CROSS_LIMITS = {"energy_meV_atom": 0.5, "force_eV_A": 0.005, "stress_kbar": 0.1}
FIXED_DEFAULTS = {"nbands": 16, "point_group": True,
                  "scf_energy_eV_electron": 1e-8,
                  "scf_density_electron_electron": 1e-8,
                  "scf_eigenstates_eV2_electron": 1e-10}


def file_sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False)


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def number(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("Expected finite numeric value")
    if not math.isfinite(value):
        raise ValueError("Nonfinite numeric value")
    return float(value)


def vector(value, size):
    if len(value) != size:
        raise ValueError(f"Expected vector of size {size}")
    return [number(x) for x in value]


def norm(value):
    return math.sqrt(sum(x * x for x in value))


def close(left, right, tolerance=1e-9):
    return abs(number(left) - number(right)) <= tolerance


def evaluation_metrics(item):
    """Recompute raw physical metrics; projected forces never enter acceptance."""
    total, dft, disp = item["total"], item["dft"], item["dispersion"]
    natoms = int(item["natoms"])
    if natoms <= 0 or len(total["forces_eV_A"]) != natoms:
        raise ValueError("Invalid atom count/raw force shape")
    forces = [vector(row, 3) for row in total["forces_eV_A"]]
    stress = vector(total["stress_kbar"], 6)
    stress_ev = vector(total["stress_eV_A3"], 6)
    if any(not close(a, b * KBAR_PER_EV_A3, 1e-7) for a, b in zip(stress, stress_ev)):
        raise ValueError("Stress units are inconsistent")
    component_forces = []
    component_stress = []
    for component in (dft, disp):
        if len(component["forces_eV_A"]) != natoms:
            raise ValueError("DFT/dispersion force shape mismatch")
        component_forces.append([vector(row, 3) for row in component["forces_eV_A"]])
        component_stress.append(vector(component["stress_eV_A3"], 6))
    if any(not close(total_value, dft_value + disp_value, 1e-9)
           for total_row, dft_row, disp_row in zip(forces, *component_forces)
           for total_value, dft_value, disp_value in zip(total_row, dft_row, disp_row)):
        raise ValueError("Raw total forces do not equal DFT plus dispersion")
    if any(not close(a, b + c, 1e-11) for a, b, c in zip(stress_ev, *component_stress)):
        raise ValueError("Raw total stress does not equal DFT plus dispersion")
    cell = [vector(row, 3) for row in item["cell_A"]]
    if len(cell) != 3:
        raise ValueError("Expected 3x3 cell")
    lengths = [norm(row) for row in cell]
    if min(lengths) <= 0:
        raise ValueError("Degenerate cell")
    reported_lengths = vector(item["cell_lengths_A"], 3)
    if any(not close(a, b) for a, b in zip(lengths, reported_lengths)):
        raise ValueError("Stored cell lengths do not match cell vectors")
    fmax, smax = max(map(norm, forces)), max(map(abs, stress))
    if not close(fmax, item["raw_max_force_norm_eV_A"]):
        raise ValueError("Stored raw fmax does not match raw forces")
    if not close(smax, item["raw_max_stress_component_kbar"]):
        raise ValueError("Stored raw stress maximum does not match raw stress")
    energy_fields = {
        "total_energy_extrapolated_eV": total["energy_extrapolated_eV"],
        "total_energy_free_eV": total["energy_free_eV"],
        "dft_energy_extrapolated_eV": dft["energy_extrapolated_eV"],
        "dft_energy_free_eV": dft["energy_free_eV"],
        "dispersion_energy_eV": disp["energy_eV"],
        "dispersion_free_energy_eV": disp["free_energy_eV"],
    }
    energies = {key: number(value) for key, value in energy_fields.items()}
    if (not close(energies["total_energy_extrapolated_eV"],
                  energies["dft_energy_extrapolated_eV"] + energies["dispersion_energy_eV"], 1e-8)
            or not close(energies["total_energy_free_eV"],
                         energies["dft_energy_free_eV"] + energies["dispersion_free_energy_eV"], 1e-8)
            or not close(energies["dispersion_energy_eV"],
                         energies["dispersion_free_energy_eV"], 1e-12)):
        raise ValueError("Total/DFT/dispersion energy definitions are inconsistent")
    if item.get("SCF_converged") is not True:
        raise ValueError("Evaluation SCF is not converged")
    return {**energies, "natoms": natoms, "raw_fmax_eV_A": fmax,
            "raw_stress_max_kbar": smax, "cell_A": cell,
            "cell_lengths_A": lengths, "stress_kbar": stress}


def source_lineage(directory, evidence, seen=None):
    """Verify a geometry-only continuation; never infer ancestry from labels."""
    directory = Path(directory).resolve()
    seen = set() if seen is None else set(seen)
    if directory in seen or len(seen) >= 16:
        raise ValueError("Cyclic or excessive continuation ancestry")
    seen.add(directory)
    spec = evidence["input"]
    attempt_maxsteps = spec.get("maxsteps", 200)
    if type(attempt_maxsteps) is not int or not 1 <= attempt_maxsteps <= 200:
        raise ValueError("Invalid integer optimization-step budget")
    source = Path(evidence["source_geometry_path"])
    if not source.is_absolute():
        source = directory / source
    source = source.resolve()
    if not source.is_file() or file_sha(source) != evidence["source_geometry_sha256"]:
        raise ValueError("Continuation/source geometry hash mismatch")
    if spec.get("structure"):
        requested = Path(spec["structure"])
        if not requested.is_absolute():
            requested = directory / requested
        if requested.resolve() != source:
            raise ValueError("Source geometry differs from the actual input structure")
    direct = {"root_source_geometry_path": str(source),
              "root_source_geometry_sha256": evidence["source_geometry_sha256"],
              "original_maxsteps": attempt_maxsteps,
              "completed_ancestor_steps": 0, "ancestor_steps": [],
              "continuations": []}
    manifest_name = spec.get("continuation_manifest")
    if not manifest_name:
        if spec.get("continuation_manifest_sha256"):
            raise ValueError("Continuation manifest hash has no manifest")
        return direct
    manifest_path = (directory / manifest_name).resolve()
    if (manifest_path.parent != directory or not manifest_path.is_file()
            or file_sha(manifest_path) != spec.get("continuation_manifest_sha256")):
        raise ValueError("Missing, external, or mismatched continuation manifest")
    manifest = read_json(manifest_path)
    if (manifest.get("schema_version") != 1 or manifest.get("restart_mode") !=
            "last_completed_geometry_electronic_cold_LBFGS_reinitialized"):
        raise ValueError("Unsupported continuation mode/schema")
    parent = Path(manifest["parent_directory"]).resolve()
    if parent.parent != directory.parent or parent == directory:
        raise ValueError("Continuation parent must be another attempt in this batch")
    if manifest.get("parent_run_id") != parent.name:
        raise ValueError("Continuation parent ID/path mismatch")
    last_filename = manifest["parent_last_evaluation_file"]
    hashes = manifest["parent_evidence_sha256"]
    required = {"input.json", "provenance.json", "trajectory.traj", "evaluations.jsonl",
                "run_relax_snapshot.py", "run_point_snapshot.py", "datasets.json",
                last_filename}
    if not required.issubset(hashes):
        raise ValueError("Continuation manifest omits required parent evidence")
    for name, digest in hashes.items():
        actual = (parent / name).resolve()
        if actual.parent != parent or not actual.is_file() or file_sha(actual) != digest:
            raise ValueError(f"Missing or changed continuation parent evidence: {name}")
    parent_evidence = read_json(parent / "provenance.json")
    parent_input = read_json(parent / "input.json")
    if (parent_evidence["input"] != parent_input
            or parent_evidence["input_sha256"] != hashes["input.json"]):
        raise ValueError("Continuation parent input/provenance mismatch")
    for filename, key in (("run_relax_snapshot.py", "runner_sha256"),
                          ("run_point_snapshot.py", "point_runner_sha256")):
        if parent_evidence[key] != hashes[filename] or evidence[key] != hashes[filename]:
            raise ValueError("Continuation changed a scientific runner")
    ignored = {"structure", "continuation_manifest", "continuation_manifest_sha256"}
    budget = manifest.get("step_budget")
    if budget is not None:
        # A changed maxsteps is permitted only after the ancestry-derived
        # remaining budget is independently verified below. Other scientific
        # settings retain the original strict equality rule.
        if not isinstance(budget, dict):
            raise ValueError("Continuation step_budget must be an explicit mapping")
        ignored.add("maxsteps")
    if ({k: v for k, v in spec.items() if k not in ignored} !=
            {k: v for k, v in parent_input.items() if k not in ignored}):
        raise ValueError("Continuation changed non-structural input parameters")
    for key in ("software", "mpi_world_size"):
        if key in evidence and evidence[key] != parent_evidence[key]:
            raise ValueError(f"Continuation changed {key}")
    if "datasets" in evidence and evidence["datasets"] != read_json(parent / "datasets.json"):
        raise ValueError("Continuation changed PAW datasets")
    history = [json.loads(line) for line in (parent / "evaluations.jsonl").read_text().splitlines()
               if line.strip()]
    last = read_json(parent / last_filename)
    if (not history or canonical(history[-1]) != canonical(last)
            or [item["evaluation_id"] for item in history] != list(range(len(history)))
            or last_filename != f"evaluation_{last['evaluation_id']:04d}.json"
            or manifest["parent_last_evaluation_id"] != last["evaluation_id"]
            or manifest["parent_optimizer_step"] != last["optimizer_step"]):
        raise ValueError("Continuation does not bind the last completed evaluation")
    evaluation_metrics(last)  # Requires finite, converged, consistent raw properties.
    parent_steps = last["optimizer_step"]
    parent_maxsteps = parent_input.get("maxsteps", 200)
    if (type(parent_steps) is not int or type(parent_maxsteps) is not int
            or not 0 <= parent_steps <= parent_maxsteps <= 200):
        raise ValueError("Completed parent steps exceed the verified parent budget")
    if (source != (directory / manifest["source_file"]).resolve()
            or file_sha(source) != manifest["source_sha256"]
            or manifest["source_sha256"] != hashes["trajectory.traj"]):
        raise ValueError("Continuation source is not the frozen parent trajectory")
    from ase.io.trajectory import Trajectory
    with Trajectory(str(source), "r") as trajectory:
        if (len(trajectory) != len(history) or len(trajectory) != manifest["source_frame_count"]
                or manifest["selected_frame_index"] != len(trajectory) - 1):
            raise ValueError("Continuation trajectory/frame count mismatch")
        frame = trajectory[-1]
    if (frame.cell.array.tolist() != last["cell_A"]
            or frame.positions.tolist() != last["positions_A"]
            or len(frame) != last["natoms"] or not frame.pbc.all()
            or frame.numbers.tolist().count(6) * 2 != frame.numbers.tolist().count(1)
            or set(frame.numbers) != {1, 6}):
        raise ValueError("Continuation binary last-frame geometry does not match its evaluation")
    for property_name, evaluation_name in (("energy", "energy_extrapolated_eV"),
                                           ("free_energy", "energy_free_eV"),
                                           ("forces", "forces_eV_A"),
                                           ("stress", "stress_eV_A3")):
        value = frame.calc.results[property_name]
        if hasattr(value, "tolist"):
            value = value.tolist()
        if value != last["total"][evaluation_name]:
            raise ValueError(f"Continuation cached raw {property_name} differs from its evaluation")
    ancestry = source_lineage(parent, parent_evidence, seen)
    if (manifest["original_source_geometry_sha256"] != ancestry["root_source_geometry_sha256"]
            or Path(manifest["original_source_geometry_path"]).resolve() !=
            Path(ancestry["root_source_geometry_path"]).resolve()):
        raise ValueError("Continuation claimed original source differs from verified ancestry")
    ancestor_steps = ancestry["ancestor_steps"] + [{"run_id": parent.name,
                                                   "steps": parent_steps}]
    completed_steps = sum(item["steps"] for item in ancestor_steps)
    if budget is not None:
        required_budget_keys = {"scope", "original_maxsteps", "completed_ancestor_steps",
                                "remaining_maxsteps", "ancestor_steps"}
        if set(budget) != required_budget_keys or budget["scope"] != "entire_reference_ancestry":
            raise ValueError("Unsupported or incomplete ancestry step-budget declaration")
        original_maxsteps = ancestry["original_maxsteps"]
        remaining = original_maxsteps - completed_steps
        declared_steps = budget["ancestor_steps"]
        if (not isinstance(declared_steps, list)
                or any(not isinstance(item, dict) or set(item) != {"run_id", "steps"}
                       or type(item["steps"]) is not int for item in declared_steps)
                or declared_steps != ancestor_steps):
            raise ValueError("Declared ancestor step counts differ from completed evaluations")
        expected = {"original_maxsteps": original_maxsteps,
                    "completed_ancestor_steps": completed_steps,
                    "remaining_maxsteps": remaining}
        if any(type(budget[key]) is not int or budget[key] != value
               for key, value in expected.items()):
            raise ValueError("Continuation ancestry step-budget arithmetic does not match evidence")
        if (remaining < 1 or attempt_maxsteps != remaining
                or attempt_maxsteps > parent_maxsteps):
            raise ValueError("Continuation must use exactly the remaining original step budget")
    ancestry["ancestor_steps"] = ancestor_steps
    ancestry["completed_ancestor_steps"] = completed_steps
    ancestry["continuations"].append({
        "run_id": directory.name, "parent_run_id": parent.name,
        "manifest_path": str(manifest_path), "manifest_sha256": file_sha(manifest_path),
        "direct_source_geometry_sha256": evidence["source_geometry_sha256"],
        "parent_last_evaluation_id": last["evaluation_id"],
        "parent_optimizer_step": last["optimizer_step"],
        "step_budget": budget,
        "restart_mode": manifest["restart_mode"]})
    return ancestry


def inspect_run(directory):
    result_path, history_path = directory / "result.json", directory / "evaluations.jsonl"
    run = {"run_id": directory.name, "directory": str(directory),
           "status": "pending", "acceptance": "pending", "errors": [],
           "history_complete": True, "evaluations_count": 0, "rows": []}
    result = None
    if result_path.exists():
        try:
            result = read_json(result_path)
            run.update(result_sha256=file_sha(result_path),
                       result_path=str(result_path), status=result.get("status", "unknown"))
        except (ValueError, OSError) as exc:
            run["errors"].append(f"result: {exc}")
    history = []
    if history_path.exists():
        raw_bytes = history_path.read_bytes()
        run["evaluations_sha256"] = hashlib.sha256(raw_bytes).hexdigest()
        for line_number, line in enumerate(raw_bytes.decode("utf-8-sig").splitlines(), 1):
            if not line.strip():
                continue
            try:
                item = json.loads(line)
                history.append(item)
                metrics = evaluation_metrics(item)
                row = {"run_id": directory.name, "line_number": line_number,
                       "evaluation_id": item["evaluation_id"],
                       "optimizer_step": item["optimizer_step"], "phase": item["phase"],
                       "wall_seconds": item.get("wall_seconds"),
                       "SCF_iterations": item.get("SCF_iterations"),
                       "row_status": "valid", "source": str(history_path), **metrics}
            except (ValueError, KeyError, TypeError) as exc:
                run["errors"].append(f"history line {line_number}: {exc}")
                run["history_complete"] = False
                # Retain the malformed/partial record as text; never silently drop it.
                row = {"run_id": directory.name, "line_number": line_number,
                       "row_status": "invalid_or_incomplete", "error": str(exc),
                       "raw_line": line, "source": str(history_path)}
            run["rows"].append(row)
    run["evaluations_count"] = len(history)
    if history:
        ids = [item.get("evaluation_id") for item in history]
        if ids != list(range(len(history))):
            run["errors"].append("Evaluation IDs must be consecutive from zero")
        steps = [item.get("optimizer_step") for item in history]
        if any(not isinstance(x, int) or x < 0 for x in steps) or steps != sorted(steps):
            run["errors"].append("Optimizer step history is invalid")
    if result is None:
        return run
    run["input"] = result.get("input", {})
    run["optimizer_steps"] = result.get("optimizer_steps")
    run["final_geometry_sha256"] = result.get("final_geometry_sha256")
    run["runner_sha256"] = result.get("runner_sha256")
    run["point_runner_sha256"] = result.get("point_runner_sha256")
    run["source_geometry_sha256"] = result.get("source_geometry_sha256")
    for row in run["rows"]:
        row.update(run_status=run["status"], method=run["input"].get("method"),
                   result_sha256=run["result_sha256"],
                   final_optimizer_steps=run["optimizer_steps"])
        row.update({key: run["input"].get(key) for key in PRECISION_KEYS})
    if result.get("status") != "passed_relaxation":
        run["acceptance"] = "not_passed"
        return run
    try:
        if result.get("exit_code") != 0:
            raise ValueError("Passed status has nonzero/missing exit code")
        provenance_files = {
            "input.json": "input_sha256",
            "run_relax_snapshot.py": "runner_sha256",
            "run_point_snapshot.py": "point_runner_sha256",
            "final_structure.extxyz": "final_geometry_sha256",
        }
        for filename, field in provenance_files.items():
            if not (directory / filename).is_file() or file_sha(directory / filename) != result.get(field):
                raise ValueError(f"Missing or mismatched actual provenance file: {filename}")
        if read_json(directory / "input.json") != result["input"]:
            raise ValueError("Archived input differs from result.input")
        source_path = Path(result["source_geometry_path"])
        if not source_path.is_absolute():
            source_path = directory / source_path
        if not source_path.is_file() or file_sha(source_path) != result["source_geometry_sha256"]:
            raise ValueError("Actual source geometry file/hash mismatch")
        units = result["units"]
        if units.get("energy") != "eV/cell" or units.get("forces") != "eV/Angstrom":
            raise ValueError("Unexpected energy/force units")
        if units.get("stress_voigt_order") != ["xx", "yy", "zz", "yz", "xz", "xy"]:
            raise ValueError("Unexpected stress component ordering")
        final = result["final_evaluation"]
        if not history or canonical(final) != canonical(history[-1]):
            raise ValueError("Final evaluation differs from the last history geometry/properties")
        metrics = evaluation_metrics(final)
        run["final_metrics"] = metrics
        run["final_geometry"] = {"cell_A": final["cell_A"], "positions_A": final["positions_A"]}
        run["final_raw_properties"] = final["total"]
        if (directory / "trajectory.traj").is_file():
            run["full_precision_trajectory_sha256"] = file_sha(directory / "trajectory.traj")
        acceptance = result["acceptance"]
        force_target = min(number(acceptance["raw_atomic_fmax_eV_A"]), FORCE_LIMIT)
        stress_target = min(number(acceptance["raw_stress_max_kbar"]), STRESS_LIMIT)
        if metrics["raw_fmax_eV_A"] > force_target:
            raise ValueError("Raw atomic forces exceed acceptance")
        if metrics["raw_stress_max_kbar"] > stress_target:
            raise ValueError("Raw stress exceeds acceptance")
        if force_target <= 0 or stress_target <= 0:
            raise ValueError("Invalid physical tolerance")
        sym = result["final_symmetry"]
        hall, standard = sym["Hall297_input_setting"], sym["default_standard_setting"]
        transform = [vector(row, 3) for row in hall["transformation_matrix"]]
        if (sym.get("passed") is not True or not hall.get("valid")
                or hall.get("number") != 62 or hall.get("operations") != 8
                or hall.get("hall_number") != 297 or len(transform) != 3
                or not standard.get("valid") or standard.get("number") != 62
                or standard.get("operations") != 8
                or any(not close(a, b, 1e-8)
                       for row, expected in zip(transform, IDENTITY)
                       for a, b in zip(row, expected))):
            raise ValueError("Final symmetry is not #62 Pnam in the original axes")
        symprec = number(run["input"].get("symprec_A", 1e-5))
        if number(final["Pnam_operation_max_residual_A"]) > symprec:
            raise ValueError("Final geometry fails target Pnam operations")
        if any(abs(number(a) - 90) > 1e-6 for a in vector(final["cell_angles_deg"], 3)):
            raise ValueError("Final cell is not orthorhombic")
        phases = result["optimization_phases"]
        if (not phases or phases[-1].get("optimizer_converged") is not True
                or phases[-1].get("physical_acceptance_pass") is not True):
            raise ValueError("No accepted final optimizer phase")
        if (not isinstance(run["optimizer_steps"], int) or run["optimizer_steps"] < 0
                or run["optimizer_steps"] > 200):
            raise ValueError("Invalid optimizer step count")
        datasets = result["datasets"]
        if len(datasets) != metrics["natoms"]:
            raise ValueError("PAW evidence does not cover every atom")
        paw = [{key: item[key] for key in ("symbol", "type", "setupname", "file_sha256")}
               for item in datasets]
        if any(not item["file_sha256"] for item in paw):
            raise ValueError("Missing PAW file hash")
        for key in ("runner_sha256", "point_runner_sha256", "source_geometry_sha256",
                    "final_geometry_sha256"):
            if not run.get(key):
                raise ValueError(f"Missing {key}")
        run["source_lineage"] = source_lineage(directory, result)
        run["comparison_identity"] = {
            "method": run["input"]["method"], "paw": paw,
            "runner_sha256": run["runner_sha256"],
            "point_runner_sha256": run["point_runner_sha256"],
            "software": result["software"],
            "source_geometry_sha256": run["source_lineage"]["root_source_geometry_sha256"],
            "dispersion_parameters": final["dispersion"]["calculator_parameters"],
            "fixed_numerical_parameters": {key: run["input"].get(key, default)
                                            for key, default in FIXED_DEFAULTS.items()},
            "mpi_ranks": result.get("mpi_world_size"),
            "natoms": metrics["natoms"], "symmetry": "#62 Pnam / Hall297 original axes"}
        run["precision"] = {key: run["input"][key] for key in PRECISION_KEYS}
    except (ValueError, KeyError, TypeError, OSError) as exc:
        run["errors"].append(str(exc))
    run["acceptance"] = "verified_relaxation" if not run["errors"] else "invalid_passed_claim"
    return run


def verify_cross_point(check, run, reference, base):
    """Verify actual archived point evidence and recompute all three differences."""
    path = Path(check["result_path"])
    if not path.is_absolute():
        path = base / path
    if not path.is_file() or file_sha(path) != check.get("point_result_sha256"):
        raise ValueError("Actual cross-check result file/hash mismatch")
    result = read_json(path)
    spec = result["input"]
    folder = path.parent
    if result.get("status") != "passed_single_point" or result.get("SCF_converged") is not True:
        raise ValueError("Cross-check single point is not SCF-converged/passed")
    if result.get("runner_sha256") != reference["point_runner_sha256"]:
        raise ValueError("Cross-check runner differs from the accepted relaxation runner")
    if result["runner_sha256"] != COLD_START_RUNNER_SHA256:
        raise ValueError("Runner version lacks the independent cold-start implementation audit")
    if (not (folder / "run_point_snapshot.py").is_file()
            or file_sha(folder / "run_point_snapshot.py") != result["runner_sha256"]):
        raise ValueError("Actual cross-check runner snapshot/hash mismatch")
    if (not (folder / "input.json").is_file()
            or file_sha(folder / "input.json") != result.get("input_sha256")
            or read_json(folder / "input.json") != spec):
        raise ValueError("Actual cross-check input/hash mismatch")
    if any(key in spec for key in ("restart", "restart_file", "density_restart", "gpw_restart")):
        raise ValueError("Cross-check unexpectedly requests a restart")
    if spec.get("method") != reference["comparison_identity"]["method"]:
        raise ValueError("Cross-check method mismatch")
    if any(spec.get(key) != reference["precision"][key] for key in PRECISION_KEYS):
        raise ValueError("Cross-check is not at the common reference precision")
    fixed = {key: spec.get(key, default) for key, default in FIXED_DEFAULTS.items()}
    if fixed != reference["comparison_identity"]["fixed_numerical_parameters"]:
        raise ValueError("Cross-check fixed numerical parameters differ")
    point_software = dict(result.get("software", {}))
    relaxation_software = dict(reference["comparison_identity"]["software"])
    relaxation_software.pop("spglib", None)  # Geometry tool is only recorded by run_relax.
    if point_software != relaxation_software:
        raise ValueError("Cross-check scientific software versions differ")
    source = Path(spec["structure"])
    if not source.is_absolute():
        source = folder / source
    if not source.is_file() or file_sha(source) != result.get("source_geometry_sha256"):
        raise ValueError("Cross-check actual geometry source file/hash mismatch")
    source_hash = result["source_geometry_sha256"]
    allowed_sources = {run["final_geometry_sha256"]: "final_structure.extxyz"}
    if run.get("full_precision_trajectory_sha256"):
        allowed_sources[run["full_precision_trajectory_sha256"]] = "trajectory.traj last frame (full precision)"
    if source_hash not in allowed_sources:
        raise ValueError("Cross-check source is not the bound final structure/trajectory")
    if check.get("source_geometry_sha256", source_hash) != source_hash:
        raise ValueError("Cross-check report geometry-source SHA is stale")
    if result.get("natoms") != run["final_metrics"]["natoms"]:
        raise ValueError("Cross-check atom count mismatch")
    if any(result.get(key) != run["final_geometry"][key] for key in ("cell_A", "positions_A")):
        raise ValueError("Cross-check is not at the identical final cell/coordinates")
    paw = [{key: item[key] for key in ("symbol", "type", "setupname", "file_sha256")}
           for item in result["datasets"]]
    if paw != reference["comparison_identity"]["paw"]:
        raise ValueError("Cross-check PAW hashes/atom order differ")
    if result["dispersion"].get("calculator_parameters") != reference["comparison_identity"]["dispersion_parameters"]:
        raise ValueError("Cross-check dispersion parameters differ")
    units = result.get("units", {})
    if (units.get("energy") != "eV/cell" or units.get("forces") != "eV/Angstrom"
            or units.get("stress_voigt_order") != ["xx", "yy", "zz", "yz", "xz", "xy"]):
        raise ValueError("Cross-check units differ")
    point, final = result["total"], run["final_raw_properties"]
    force = [vector(row, 3) for row in point["forces_eV_A"]]
    if len(force) != result["natoms"]:
        raise ValueError("Cross-check force shape mismatch")
    stress = vector(point["stress_kbar"], 6)
    stress_ev = vector(point["stress_eV_A3"], 6)
    if any(not close(a, b * KBAR_PER_EV_A3, 1e-7) for a, b in zip(stress, stress_ev)):
        raise ValueError("Cross-check stress conversion mismatch")
    delta = {
        "energy_meV_atom": abs(number(point["energy_free_eV"]) - number(final["energy_free_eV"]))
                            * 1000 / result["natoms"],
        "force_eV_A": max(norm([a - b for a, b in zip(row, old)])
                         for row, old in zip(force, final["forces_eV_A"])),
        "stress_kbar": max(abs(a - b) for a, b in zip(stress, final["stress_kbar"])),
    }
    flags = {key: value <= CROSS_LIMITS[key] for key, value in delta.items()}
    return {"relaxation_run_id": run["run_id"], "result_path": str(path),
            "point_result_sha256": file_sha(path), "differences": delta,
            "source_geometry_sha256": source_hash,
            "geometry_source": allowed_sources[source_hash],
            "threshold_pass": flags, "all_pass": all(flags.values()),
            "cold_start_evidence": "Actual frozen snapshot matches the independently audited fresh-GPAW runner SHA256"}


def cross_validation(report, candidate, reference):
    if report is None:
        return {"status": "pending", "reason": "High-precision final-geometry cross-check is missing"}
    try:
        if report.get("status") != "verified":
            return {"status": "pending", "reason": "Cross-check status is not verified"}
        for name, run in (("candidate", candidate), ("reference", reference)):
            binding = report[name]
            for key in ("run_id", "result_sha256", "final_geometry_sha256"):
                if binding.get(key) != run.get(key):
                    raise ValueError(f"Cross-check {name}.{key} is stale or mismatched")
        if report.get("method") != reference["comparison_identity"]["method"]:
            raise ValueError("Cross-check method mismatch")
        if report.get("reference_precision") != reference["precision"]:
            raise ValueError("Cross-check precision does not match the accepted reference run")
        checks = report["checks"]
        by_run = {item["relaxation_run_id"]: item for item in checks}
        if len(by_run) != 2 or len(checks) != 2:
            raise ValueError("Exactly two distinct final-geometry checks are required")
        for run in (candidate, reference):
            check = by_run[run["run_id"]]
            if (check.get("status") != "verified"
                    or check.get("geometry_sha256") != run["final_geometry_sha256"]
                    or not check.get("point_result_sha256") or not check.get("result_path")):
                raise ValueError("Missing verified, geometry-bound single-point evidence")
        base = Path(report.get("_report_directory", reference["directory"])).resolve()
        verified = [verify_cross_point(by_run[run["run_id"]], run, reference, base)
                    for run in (candidate, reference)]
        return {"status": "verified" if all(x["all_pass"] for x in verified) else "failed",
                "source_status": report["status"], "limits": CROSS_LIMITS,
                "note": "Actual point/input/snapshot files, exact geometry and all three differences independently verified; no SCF run by this analyzer.",
                "checks": verified}
    except (ValueError, KeyError, TypeError, OSError) as exc:
        return {"status": "pending", "reason": str(exc)}


def precision_coordinates(run):
    p = run["precision"]
    return [number(p["ecut_eV"]), *vector(p["kmesh"], 3), *vector(p["gpts"], 3)]


def compare_runs(runs, candidate_id=None, reference_id=None, cross_report=None):
    comparison = {"status": "pending", "cell_length_status": "pending",
                  "relative_tolerance": CELL_RELATIVE_LIMIT,
                  "relative_definition": "abs(candidate-reference)/reference for a,b,c",
                  "strict_inequality": True, "cross_validation": {"status": "pending"}}
    accepted = [run for run in runs if run["acceptance"] == "verified_relaxation"]
    if bool(candidate_id) != bool(reference_id):
        raise ValueError("--candidate and --reference must be specified together")
    if candidate_id:
        mapping = {run["run_id"]: run for run in runs}
        if candidate_id == reference_id or candidate_id not in mapping or reference_id not in mapping:
            comparison["reason"] = "Select two distinct, existing relaxation attempts"
            return comparison
        candidate, reference = mapping[candidate_id], mapping[reference_id]
        if any(run["acceptance"] != "verified_relaxation" for run in (candidate, reference)):
            comparison["reason"] = "Both selected relaxations must independently pass"
            return comparison
    elif len(accepted) != 2:
        comparison["reason"] = "Exactly two accepted runs or an explicit pair are required"
        return comparison
    else:
        candidate, reference = sorted(accepted, key=precision_coordinates)
    comparison.update(candidate_run=candidate["run_id"], reference_run=reference["run_id"])
    if canonical(candidate["comparison_identity"]) != canonical(reference["comparison_identity"]):
        comparison["reason"] = "Method, PAW, symmetry, runner, software, source geometry, or dispersion differ"
        return comparison
    left, right = precision_coordinates(candidate), precision_coordinates(reference)
    if (not all(a <= b for a, b in zip(left, right)) or left == right
            or candidate["precision"]["xc_stencil"] != reference["precision"]["xc_stencil"]):
        comparison["reason"] = "Reference must refine ecut/kmesh/grid without changing XC stencil"
        return comparison
    ca, cb = candidate["final_metrics"]["cell_lengths_A"], reference["final_metrics"]["cell_lengths_A"]
    relative = [abs(a - b) / b for a, b in zip(ca, cb)]
    # Do not let binary rounding turn an exactly-on-threshold value into a pass.
    cell_pass = all(value < CELL_RELATIVE_LIMIT - 1e-12 for value in relative)
    comparison.update(candidate_cell_lengths_A=ca, reference_cell_lengths_A=cb,
                      relative_differences=dict(zip(("a", "b", "c"), relative)),
                      max_relative_difference=max(relative),
                      cell_length_status="passed" if cell_pass else "failed")
    comparison["cross_validation"] = cross_validation(cross_report, candidate, reference)
    if not cell_pass:
        comparison.update(status="failed", reason="At least one cell length differs by 0.5% or more")
    elif comparison["cross_validation"]["status"] == "verified":
        comparison.update(status="passed", reason="Cell agreement and bound cross-check evidence both pass")
    elif comparison["cross_validation"]["status"] == "failed":
        comparison.update(status="failed", reason="High-precision same-geometry differences exceed acceptance")
    else:
        comparison["reason"] = "Cell lengths agree, but high-precision cross-check remains pending"
    return comparison


def analyze(root, candidate=None, reference=None, cross_path=None):
    relaxation_root = root / "relaxations"
    runs = [inspect_run(path) for path in sorted(relaxation_root.glob("*"))
            if path.is_dir() and ((path / "result.json").exists()
                                 or (path / "evaluations.jsonl").exists()
                                 or (path / "input.json").exists())]
    cross_path = cross_path or relaxation_root / "cross_precision_verification.json"
    cross_report, cross_read_error = None, None
    if cross_path.exists():
        try:
            cross_report = read_json(cross_path)
            cross_report["_report_directory"] = str(cross_path.resolve().parent)
        except (ValueError, OSError) as exc:
            cross_read_error = str(exc)
    comparison = compare_runs(runs, candidate, reference, cross_report)
    if cross_read_error:
        comparison["cross_validation"] = {"status": "pending", "reason": cross_read_error}
    if cross_path.exists():
        comparison["cross_validation"]["report_path"] = str(cross_path)
        comparison["cross_validation"]["report_sha256"] = file_sha(cross_path)
    rows = [row for run in runs for row in run["rows"]]
    summary = {
        "schema_version": 1, "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "analyzer_sha256": file_sha(__file__), "root": str(root),
        "status": comparison["status"], "data_status": "available" if runs else "not_started",
        "limits": {"raw_force_eV_A": FORCE_LIMIT, "raw_stress_kbar": STRESS_LIMIT,
                   "cell_relative_strict_upper": CELL_RELATIVE_LIMIT},
        "comparison": comparison,
        "runs": [{key: value for key, value in run.items() if key != "rows"} for run in runs],
        "evaluations": rows,
        "note": "All attempts and evaluations are retained. No total success without verified high-precision cross-check evidence."
    }
    columns = ["run_id", "run_status", "method", *PRECISION_KEYS, "result_sha256",
               "line_number", "evaluation_id", "optimizer_step", "final_optimizer_steps", "phase",
               "wall_seconds", "SCF_iterations", "row_status", "natoms",
               "total_energy_extrapolated_eV", "total_energy_free_eV",
               "dft_energy_extrapolated_eV", "dft_energy_free_eV",
               "dispersion_energy_eV", "dispersion_free_energy_eV",
               "raw_fmax_eV_A", "raw_stress_max_kbar",
               "cell_lengths_A", "cell_A", "stress_kbar", "source", "error", "raw_line"]
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=columns, extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow({key: canonical(value) if isinstance(value, (list, dict)) else value
                         for key, value in row.items()})
    tables = root / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    for name, content in (("relaxation_energies.json", json.dumps(summary, ensure_ascii=False, indent=2,
                                                                allow_nan=False) + "\n"),
                          ("relaxation_energies.csv", output.getvalue())):
        destination = tables / name
        temporary = destination.with_suffix(destination.suffix + ".tmp")
        temporary.write_text(content, encoding="utf-8")
        temporary.replace(destination)
    return summary


def self_test():
    global COLD_START_RUNNER_SHA256
    import copy
    item = {
        "evaluation_id": 0, "optimizer_step": 0, "phase": "initial", "SCF_converged": True,
        "SCF_iterations": 8, "natoms": 1, "cell_A": [[7., 0, 0], [0, 5., 0], [0, 0, 2.5]],
        "positions_A": [[1., 1., 1.]],
        "cell_lengths_A": [7., 5., 2.5], "cell_angles_deg": [90, 90, 90],
        "total": {"energy_extrapolated_eV": -11., "energy_free_eV": -11.001,
                  "forces_eV_A": [[.001, 0, 0]], "stress_kbar": [.01, 0, 0, 0, 0, 0],
                  "stress_eV_A3": [.01 / KBAR_PER_EV_A3, 0, 0, 0, 0, 0]},
        "dft": {"energy_extrapolated_eV": -10., "energy_free_eV": -10.001,
                "forces_eV_A": [[.0007, 0, 0]],
                "stress_eV_A3": [.007 / KBAR_PER_EV_A3, 0, 0, 0, 0, 0]},
        "dispersion": {"energy_eV": -1., "free_energy_eV": -1., "calculator_parameters": {"atm": True},
                       "forces_eV_A": [[.0003, 0, 0]],
                       "stress_eV_A3": [.003 / KBAR_PER_EV_A3, 0, 0, 0, 0, 0]},
        "raw_max_force_norm_eV_A": .001, "raw_max_stress_component_kbar": .01,
        "Pnam_operation_max_residual_A": 1e-14}
    spec = {"method": "PBE-D3(BJ)", "ecut_eV": 1400, "kmesh": [3, 6, 9],
            "gpts": [96, 60, 36], "xc_stencil": 4}
    hall = {"valid": True, "number": 62, "hall_number": 297,
            "operations": 8, "transformation_matrix": IDENTITY}
    result = {
        "status": "passed_relaxation", "exit_code": 0, "input": spec, "optimizer_steps": 0,
        "units": {"energy": "eV/cell", "forces": "eV/Angstrom",
                  "stress_voigt_order": ["xx", "yy", "zz", "yz", "xz", "xy"]},
        "final_evaluation": item, "acceptance": {"raw_atomic_fmax_eV_A": .002,
                                                "raw_stress_max_kbar": .05},
        "final_symmetry": {"passed": True, "Hall297_input_setting": hall,
                           "default_standard_setting": hall},
        "optimization_phases": [{"optimizer_converged": True, "physical_acceptance_pass": True}],
        "datasets": [{"symbol": "C", "type": "paw", "setupname": "PBE", "file_sha256": "paw"}],
        "runner_sha256": "relax", "point_runner_sha256": "point", "source_geometry_sha256": "source",
        "final_geometry_sha256": "geometry", "software": {"gpaw": "test"}}
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        assert analyze(root)["data_status"] == "not_started"
        source = root / "source.extxyz"
        source.write_text("controlled synthetic source, never used for calculation")
        for name, cutoff in (("candidate", 1400), ("reference", 1600)):
            folder = root / "relaxations" / name
            folder.mkdir(parents=True)
            value = copy.deepcopy(result)
            value["input"]["ecut_eV"] = cutoff
            (folder / "input.json").write_text(json.dumps(value["input"]))
            (folder / "run_relax_snapshot.py").write_text("# synthetic relaxation fixture\n")
            (folder / "run_point_snapshot.py").write_text("# synthetic cold-start fixture\n")
            (folder / "final_structure.extxyz").write_text("controlled synthetic final structure")
            for filename, field in (("input.json", "input_sha256"),
                                    ("run_relax_snapshot.py", "runner_sha256"),
                                    ("run_point_snapshot.py", "point_runner_sha256"),
                                    ("final_structure.extxyz", "final_geometry_sha256")):
                value[field] = file_sha(folder / filename)
            value["source_geometry_path"] = str(source)
            value["source_geometry_sha256"] = file_sha(source)
            (folder / "result.json").write_text(json.dumps(value))
            (folder / "evaluations.jsonl").write_text(json.dumps(item) + "\n")
        summary = analyze(root)
        assert summary["comparison"]["cell_length_status"] == "passed"
        assert summary["status"] == "pending"
        accepted = [inspect_run(root / "relaxations" / x) for x in ("candidate", "reference")]
        cross = {"status": "verified", "method": "PBE-D3(BJ)",
                 "reference_precision": accepted[1]["precision"], "checks": []}
        audited_hash = COLD_START_RUNNER_SHA256
        COLD_START_RUNNER_SHA256 = accepted[1]["point_runner_sha256"]
        for role, run in zip(("candidate", "reference"), accepted):
            cross[role] = {key: run[key] for key in ("run_id", "result_sha256", "final_geometry_sha256")}
            folder = root / "cross" / role
            folder.mkdir(parents=True)
            point_spec = dict(accepted[1]["input"])
            point_spec["structure"] = str(Path(run["directory"]) / "final_structure.extxyz")
            (folder / "input.json").write_text(json.dumps(point_spec))
            (folder / "run_point_snapshot.py").write_text("# synthetic cold-start fixture\n")
            point = {"status": "passed_single_point", "SCF_converged": True,
                     "input": point_spec, "input_sha256": file_sha(folder / "input.json"),
                     "runner_sha256": COLD_START_RUNNER_SHA256,
                     "source_geometry_sha256": run["final_geometry_sha256"],
                     "natoms": 1, "software": {"gpaw": "test"},
                     **run["final_geometry"], "datasets": result["datasets"],
                     "dispersion": item["dispersion"], "units": result["units"],
                     "total": item["total"]}
            (folder / "result.json").write_text(json.dumps(point))
            cross["checks"].append({"relaxation_run_id": run["run_id"], "status": "verified",
                                    "geometry_sha256": run["final_geometry_sha256"],
                                    "point_result_sha256": file_sha(folder / "result.json"),
                                    "result_path": str(folder / "result.json")})
        assert compare_runs(accepted, cross_report=cross)["status"] == "passed"
        stale = copy.deepcopy(cross)
        stale["candidate"]["result_sha256"] = "old"
        assert compare_runs(accepted, cross_report=stale)["status"] == "pending"
        changed = copy.deepcopy(accepted)
        changed[0]["final_metrics"]["cell_lengths_A"][0] *= 1.006
        assert compare_runs(changed, cross_report=cross)["status"] == "failed"
        changed = copy.deepcopy(accepted)
        changed[0]["final_metrics"]["cell_lengths_A"][0] *= 1.005
        assert compare_runs(changed, cross_report=cross)["status"] == "failed"
        changed = copy.deepcopy(accepted)
        changed[0]["comparison_identity"]["paw"][0]["file_sha256"] = "other"
        assert compare_runs(changed, cross_report=cross)["status"] == "pending"
        missing = copy.deepcopy(cross)
        missing["checks"][0]["result_path"] = str(root / "nonexistent.json")
        assert compare_runs(accepted, cross_report=missing)["status"] == "pending"
        COLD_START_RUNNER_SHA256 = audited_hash
        bad = read_json(root / "relaxations/candidate/result.json")
        bad["final_evaluation"]["total"]["forces_eV_A"] = [[.003, 0, 0]]
        bad["final_evaluation"]["dft"]["forces_eV_A"] = [[.0027, 0, 0]]
        bad["final_evaluation"]["raw_max_force_norm_eV_A"] = .003
        folder = root / "relaxations/candidate"
        (folder / "result.json").write_text(json.dumps(bad))
        (folder / "evaluations.jsonl").write_text(json.dumps(bad["final_evaluation"]) + "\n")
        assert inspect_run(folder)["acceptance"] == "invalid_passed_claim"
        (folder / "evaluations.jsonl").write_text('{"evaluation_id":')
        assert inspect_run(folder)["history_complete"] is False
        assert (root / "tables/relaxation_energies.csv").exists()
    print("Synthetic checks passed: empty data, raw-force rejection, exact final/history match, "
          "partial history, cell threshold, PAW isolation, stale cross-check rejection, and pending/verified states.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--candidate")
    parser.add_argument("--reference")
    parser.add_argument("--cross-validation", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    summary = analyze(args.root.resolve(), args.candidate, args.reference, args.cross_validation)
    print(json.dumps({"status": summary["status"], "data_status": summary["data_status"],
                      "runs": len(summary["runs"]), "evaluations": len(summary["evaluations"]),
                      "comparison": summary["comparison"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
