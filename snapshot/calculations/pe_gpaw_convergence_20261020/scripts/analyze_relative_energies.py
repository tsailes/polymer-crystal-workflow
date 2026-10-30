"""Read-only, provenance-checked relative energies of two frozen PE structures.

No quantum engine is imported or called.  Pairing spans the allowed transfer/
joint stages, but never spans methods, PAW datasets, runners or software.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys

sys.dont_write_bytecode = True
import numpy as np
from ase.io import read
from ase.units import Hartree

ROOT = Path(__file__).resolve().parents[1]
STAGES = {"geometry_transfer", "functional_transfer", "numerical_joint"}
METHODS = ("PBE-D3(BJ)", "optB88-vdW", "r2SCAN-D4")
ROLES = ("Tashiro", "GFN2")
SCAN = {"ecut_eV", "kmesh", "gpts"}
META = {
    "id", "stage", "structure", "checkpoint", "protocol_sha256", "queue_sha256",
    "runner_sha256_at_submission", "scan_parameters", "scan_label", "scan_index",
    "description", "notes", "geometry_role", "structure_role", "precision_label",
}
DEFAULTS = {"nbands": 16, "point_group": True, "scf_energy_eV_electron": 1e-8,
            "scf_density_electron_electron": 1e-8, "scf_eigenstates_eV2_electron": 1e-10}
TOL = .5  # meV per atom, for changes in the two-structure energy difference.


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def fingerprint(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def save(path, value):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                    encoding="utf-8")
    temp.replace(path)


def local_path(value, root):
    value = str(value)
    if os.name == "nt" and re.match(r"^/mnt/[a-zA-Z]/", value):
        value = value[5].upper() + ":/" + value[7:]
    elif os.name != "nt" and re.match(r"^[a-zA-Z]:[/\\]", value):
        value = "/mnt/" + value[0].lower() + "/" + value[3:].replace("\\", "/")
    path = Path(value)
    return path if path.is_absolute() else root / path


def logical_id(value):
    return re.sub(r"__attempt\d+$", "", value)


def number(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError("Expected a finite number")
    return float(value)


def triplet(value):
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ValueError("Expected a positive integer triplet")
    values = [number(x) for x in value]
    if any(x <= 0 or x != int(x) for x in values):
        raise ValueError("Expected positive integer mesh dimensions")
    return [int(x) for x in values]


def settings(spec):
    value = dict(DEFAULTS)
    value.update({key: val for key, val in spec.items() if key not in META})
    return value


def references(root):
    out = {}
    protocol = load(root / "protocol.json") if (root / "protocol.json").exists() else {}
    for role, name, key in (("Tashiro", "POSCAR", "primary_sha256"),
                            ("GFN2", "POSCAR_GFN2", "transfer_sha256")):
        path = root / "reference" / name
        atoms = read(path)
        digest = sha(path)
        if protocol.get("geometry", {}).get(key, digest) != digest:
            raise ValueError("Reference file differs from protocol SHA256: " + role)
        if len(atoms) != 12 or atoms.get_chemical_formula() != "C4H8" or not atoms.pbc.all():
            raise ValueError("Expected periodic C4H8 reference: " + role)
        out[role] = {"path": str(path), "source_sha256": digest, "atoms": atoms,
                     "exact_geometry_sha256": fingerprint({
                         "numbers": atoms.numbers.tolist(), "pbc": atoms.pbc.tolist(),
                         "cell_A": atoms.cell.array.tolist(), "positions_A": atoms.positions.tolist()})}
    if out["Tashiro"]["source_sha256"] == out["GFN2"]["source_sha256"]:
        raise ValueError("The two reference sources are unexpectedly identical")
    return out


def role_of(spec, refs, root):
    path = local_path(spec["structure"], root)
    digest = sha(path)
    matches = [role for role in ROLES if refs[role]["source_sha256"] == digest]
    if len(matches) != 1:
        raise ValueError("Structure source SHA does not identify either frozen reference")
    return matches[0], path, digest


def extract(folder, spec, result, refs, root):
    """Validate a completed point independently of the main convergence analyzer."""
    if result.get("status") != "passed_single_point" or result.get("SCF_converged") is not True:
        raise ValueError("A successful native SCF single point is required")
    if result.get("input") != spec or result.get("input_sha256") != sha(folder / "input.json"):
        raise ValueError("Saved input/result input hash mismatch")
    runner = sha(folder / "run_point_snapshot.py")
    if runner != result.get("runner_sha256") or runner != spec.get("runner_sha256_at_submission"):
        raise ValueError("Runner snapshot, submission and result hashes disagree")
    if spec.get("displacement") or spec.get("strain") or spec.get("reference_cell_A"):
        raise ValueError("Only the unmodified frozen structures and unscaled cutoff are accepted")
    if result.get("natoms") != 12 or result.get("formula") != "C4H8" or result.get("electrons") != 24:
        raise ValueError("Expected 12-atom C4H8 with 24 valence electrons")
    if result.get("units", {}).get("energy") != "eV/cell":
        raise ValueError("Free energy must be expressed in eV per cell")
    role, source, source_sha = role_of(spec, refs, root)
    if result.get("source_geometry_sha256") != source_sha:
        raise ValueError("Result refers to a different source-file hash")
    atoms = refs[role]["atoms"]
    positions = np.asarray(result["positions_A"], float)
    cell = np.asarray(result["cell_A"], float)
    if (positions.shape != (12, 3) or cell.shape != (3, 3)
            or not np.isfinite(positions).all() or not np.isfinite(cell).all()
            or not np.allclose(positions, atoms.positions, atol=1e-12, rtol=0)
            or not np.allclose(cell, atoms.cell.array, atol=1e-12, rtol=0)):
        raise ValueError("Unrounded geometry differs from the frozen reference")
    if not np.isclose(number(result["volume_A3"]), atoms.get_volume(), atol=1e-10, rtol=0):
        raise ValueError("Volume differs from the frozen cell")
    saved = read(folder / "input_structure.extxyz")
    if (saved.numbers.tolist() != atoms.numbers.tolist() or not saved.pbc.all()
            or not np.allclose(saved.positions, positions, atol=2e-8, rtol=0)
            or not np.allclose(saved.cell.array, cell, atol=1e-10, rtol=0)):
        raise ValueError("Archived input structure differs from the exact result geometry")
    datasets = result["datasets"]
    if len(datasets) != 12 or [x["symbol"] for x in datasets] != atoms.get_chemical_symbols():
        raise ValueError("PAW element order differs from reference atom order")
    fields = ("symbol", "type", "setupname", "fingerprint", "file_sha256")
    paw = [{key: data.get(key) for key in fields} for data in datasets]
    if any(not row["file_sha256"] or row["type"] != "paw" for row in paw):
        raise ValueError("Incomplete native PAW identity")
    # Atom ordering may differ between the two structures; energies require
    # identical element/dataset multiplicities, not a coordinate assignment.
    paw = sorted(paw, key=canonical)
    software = result.get("software")
    if not isinstance(software, dict) or not software.get("gpaw"):
        raise ValueError("Software identity is missing")
    basis = result["basis"]
    coarse = triplet(basis["coarse_grid"]["shape"])
    fine = triplet(basis["fine_grid"]["shape"])
    requested_grid = triplet(spec["gpts"])
    if coarse != requested_grid:
        raise ValueError("Actual coarse grid differs from explicit gpts")
    for name in ("coarse_grid", "fine_grid"):
        if not np.allclose(basis[name]["cell_A"], cell, atol=1e-8, rtol=0):
            raise ValueError("Native grid cell differs from the frozen geometry")
    ecut = number(spec["ecut_eV"])
    kmesh = triplet(spec["kmesh"])
    if ecut <= 0:
        raise ValueError("Cutoff must be positive")
    descriptors = basis.get("reciprocal_descriptors", {})
    expected_descriptors = {"density_pd2", "density_pd3", "hamiltonian_pd2", "hamiltonian_pd3"}
    if set(descriptors) != expected_descriptors:
        raise ValueError("Effective density/Hamiltonian cutoffs were not archived")
    density_cutoffs = {key: number(val["ecut_Ha"]) * Hartree
                       for key, val in sorted(descriptors.items())}
    energy = number(result["total"]["energy_free_eV"])
    dft = number(result["dft"]["energy_free_eV"])
    dispersion = number(result["dispersion"]["energy_eV"])
    if abs(energy - dft - dispersion) > 1e-10:
        raise ValueError("Total free energy differs from DFT plus dispersion")
    dispersion_identity = {key: result["dispersion"].get(key)
                           for key in ("calculator_class", "calculator_parameters")}
    if spec["method"] in ("PBE-D3(BJ)", "r2SCAN-D4"):
        if not isinstance(dispersion_identity["calculator_parameters"], dict):
            raise ValueError("Additive dispersion parameters are missing")
    physics = {
        "method": spec["method"], "fixed_settings": {
            key: val for key, val in settings(spec).items() if key not in SCAN},
        "PAW_multiset": paw, "software": software, "runner_sha256": runner,
        "native_XC": {key: basis.get(key) for key in
                      ("xc_name", "xc_class", "actual_XC_stencil")},
        "dispersion": dispersion_identity, "electrons": 24,
        "reference_sources": {key: refs[key]["source_sha256"] for key in ROLES},
    }
    precision = {"ecut_eV": ecut, "kmesh": kmesh, "gpts": coarse, "fine_grid": fine,
                 "actual_density_cutoffs_eV": density_cutoffs}
    pair_identity = {**physics, "all_input_settings": settings(spec), "precision": precision,
                     "BZ_kpoints": basis["bz_kpoints"]}
    return {
        "logical_id": logical_id(spec["id"]), "run_id": folder.name, "stage": spec["stage"],
        "role": role, "method": spec["method"], "energy_free_eV_cell": energy,
        "DFT_free_energy_eV_cell": dft, "dispersion_energy_eV_cell": dispersion,
        "source_path": str(source), "source_sha256": source_sha,
        "exact_geometry_sha256": refs[role]["exact_geometry_sha256"],
        "result_path": str(folder / "result.json"), "result_sha256": sha(folder / "result.json"),
        "input_path": str(folder / "input.json"), "input_sha256": sha(folder / "input.json"),
        "saved_geometry_path": str(folder / "input_structure.extxyz"),
        "saved_geometry_sha256": sha(folder / "input_structure.extxyz"),
        "runner_path": str(folder / "run_point_snapshot.py"), "runner_sha256": runner,
        "started_unix": result.get("start_unix", 0),
        "mtime_ns": (folder / "result.json").stat().st_mtime_ns,
        "family_id": fingerprint(physics), "family_identity": physics,
        "pair_id": fingerprint(pair_identity), "pair_identity": pair_identity,
        "precision": precision,
    }


def refinement_vector(point):
    p = point["precision"]
    return (p["ecut_eV"], *p["kmesh"], *p["gpts"], *p["fine_grid"],
            *p["actual_density_cutoffs_eV"].values())


def refines(high, low):
    a, b = refinement_vector(high), refinement_vector(low)
    return len(a) == len(b) and all(x >= y for x, y in zip(a, b)) and any(
        x > y for x, y in zip(a, b))


def delta(left, right):
    value = left["relative_energy_meV_atom"] - right["relative_energy_meV_atom"]
    return {"left_pair_id": left["pair_id"], "right_pair_id": right["pair_id"],
            "signed_change_meV_atom": value, "absolute_change_meV_atom": abs(value),
            "passed": abs(value) <= TOL}


def assess_family(points):
    """Use the refinement partial order, never a label or energy ordering."""
    points = sorted(points, key=lambda x: refinement_vector(x))
    out = {"status": "pending", "passed": False, "complete_pair_count": len(points),
           "reasons": [], "chains": [], "accepted_candidate_pair_ids": [],
           "minimal_accepted_candidate_pair_ids": []}
    if len(points) < 3:
        out["reasons"].append("At least three distinct compatible paired precisions are required")
        return out
    highest = [p for p in points if all(p["pair_id"] == q["pair_id"] or refines(p, q)
                                      for q in points)]
    if len(highest) != 1:
        out["reasons"].append("No unique common highest precision; incomparable branches remain pending")
        return out
    reference = highest[0]
    out["highest_pair_id"] = reference["pair_id"]
    # Hasse edges omit transitive shortcuts, so one totally ordered series
    # yields one full chain instead of exponentially many subsequences.
    edges = {p["pair_id"]: [] for p in points}
    by_id = {p["pair_id"]: p for p in points}
    incoming = {key: 0 for key in edges}
    for low in points:
        for high in points:
            if not refines(high, low):
                continue
            if any(refines(mid, low) and refines(high, mid) for mid in points):
                continue
            edges[low["pair_id"]].append(high["pair_id"])
            incoming[high["pair_id"]] += 1
    paths = []

    def visit(path):
        if len(paths) >= 256:
            raise ValueError("Too many incomparable refinement branches; use separate fixed-setting series")
        children = edges[path[-1]]
        if not children:
            paths.append(path)
        for child in children:
            visit(path + [child])

    try:
        for key, count in incoming.items():
            if count == 0:
                visit([key])
    except ValueError as exc:
        out["reasons"].append(str(exc))
        return out
    accepted = set()
    for path in paths:
        chain = [by_id[key] for key in path]
        intervals = [delta(b, a) for a, b in zip(chain, chain[1:])]
        comparisons = [delta(p, reference) for p in chain]
        plateau = len(chain) >= 3 and all(x["passed"] for x in intervals[-2:])
        candidates = []
        if plateau:
            for p in chain[:-1]:
                # All measured settings componentwise above the candidate,
                # including other branches, must agree with the same high
                # reference.  A favorable path cannot hide a contradictory one.
                above = [q for q in points if q["pair_id"] == p["pair_id"] or refines(q, p)]
                if all(delta(q, reference)["passed"] for q in above):
                    candidates.append(p["pair_id"])
                    accepted.add(p["pair_id"])
        out["chains"].append({
            "pair_ids": path, "status": "passed" if candidates else
                "pending" if len(chain) < 3 else "failed",
            "terminal_two_intervals_stable": plateau,
            "adjacent_comparisons": intervals, "comparisons_to_highest": comparisons,
            "candidate_pair_ids": candidates})
    out["accepted_candidate_pair_ids"] = sorted(accepted)
    out["minimal_accepted_candidate_pair_ids"] = sorted(
        key for key in accepted if not any(
            other != key and refines(by_id[key], by_id[other]) for other in accepted))
    out["passed"] = bool(accepted)
    out["status"] = "relative_energy_convergence_passed" if accepted else "needs_attention"
    if not accepted:
        out["reasons"].append("No candidate has a three-point stable refinement chain and highest-reference agreement")
    out["comparisons_to_highest"] = [delta(p, reference) for p in points]
    return out


def scan(root, refs, stages):
    attempts, expected, unreadable, conflicts = {}, {}, [], []

    def remember(spec, path):
        ident = logical_id(spec["id"])
        previous = expected.get(ident)
        if previous is not None and (
                settings(previous) != settings(spec)
                or local_path(previous["structure"], root).resolve()
                != local_path(spec["structure"], root).resolve()):
            conflicts.append({"logical_id": ident, "path": str(path),
                              "reason": "A logical run ID has conflicting planned or archived inputs"})
        else:
            expected.setdefault(ident, spec)

    # Queue coverage includes planned endpoints whose folders do not yet exist.
    for path in sorted(root.glob("*queue*.json")):
        try:
            queue = load(path)
            if not isinstance(queue, dict):
                raise ValueError("Queue is not a JSON object")
            for item in queue.get("points", []):
                spec = {**queue.get("defaults", {}), **item}
                if spec.get("stage") in stages and spec.get("id"):
                    # run_queue.py uses this same default before writing input.
                    spec.setdefault("structure", str(root / "reference" / "POSCAR"))
                    remember(spec, path)
        except (OSError, ValueError, TypeError) as exc:
            unreadable.append({"path": str(path), "reason": str(exc)})
    for folder in sorted((root / "runs").iterdir()) if (root / "runs").exists() else []:
        if not folder.is_dir() or not (folder / "input.json").exists():
            continue
        try:
            spec = load(folder / "input.json")
        except (OSError, ValueError) as exc:
            unreadable.append({"path": str(folder / "input.json"), "reason": str(exc)})
            continue
        if spec.get("stage") not in stages:
            continue
        ident = logical_id(spec.get("id", folder.name))
        if not spec.get("id") or not spec.get("structure"):
            unreadable.append({"path": str(folder / "input.json"),
                               "reason": "Stage input is missing id or structure"})
            continue
        remember(spec, folder / "input.json")
        result, process, reason = {}, {}, None
        try:
            if (folder / "result.json").exists():
                result = load(folder / "result.json")
            if (folder / "process.json").exists():
                process = load(folder / "process.json")
        except (OSError, ValueError) as exc:
            reason = str(exc)
        success = result.get("status") == "passed_single_point" and process.get("exit_code") == 0
        state = "success" if success else "failed" if (
            result.get("status") == "failed" or process.get("exit_code") not in (None, 0)) else "pending"
        match = re.search(r"__attempt(\d+)$", folder.name)
        attempts.setdefault(ident, []).append({
            "folder": folder, "spec": spec, "result": result, "state": state, "reason": reason,
            "attempt": int(match[1]) if match else 1,
            "started_unix": result.get("start_unix", process.get("start_unix", 0))})
    records, pending, rejected, history = [], [], list(conflicts), {}
    for ident, spec in expected.items():
        items = sorted(attempts.get(ident, []),
                       key=lambda x: (x["attempt"], x["started_unix"], x["folder"].name))
        history[ident] = [{"run_id": x["folder"].name, "state": x["state"], "reason": x["reason"]}
                          for x in items]
        successful = [item for item in items if item["state"] == "success"]
        if not successful:
            pending.append({"logical_id": ident, "method": spec.get("method"),
                            "stage": spec.get("stage"), "structure": spec.get("structure"),
                            "reason": "No completed successful attempt; failed attempts remain archived"})
            continue
        chosen = successful[-1]
        try:
            record = extract(chosen["folder"], chosen["spec"], chosen["result"], refs, root)
            record["attempts"] = history[ident]
            records.append(record)
        except (OSError, ValueError, KeyError, TypeError, IndexError) as exc:
            rejected.append({"logical_id": ident, "run_id": chosen["folder"].name,
                             "reason": str(exc), "attempts": history[ident]})
    return records, pending, rejected, history, unreadable


def analyze(root, stages=STAGES, required_methods=METHODS):
    refs = references(root)
    records, pending, rejected, history, unreadable = scan(root, refs, stages)
    buckets = {}
    for record in records:
        buckets.setdefault(record["pair_id"], {}).setdefault(record["role"], []).append(record)
    pairs, unmatched = [], []
    for pair_id, roles in sorted(buckets.items()):
        if set(roles) != set(ROLES):
            unmatched.append({
                "pair_id": pair_id, "method": next(iter(roles.values()))[0]["method"],
                "roles_present": sorted(roles), "missing_roles": sorted(set(ROLES) - set(roles)),
                "reason": "The two structures lack an exactly matching method/settings/provenance point"})
            continue
        chosen, duplicate_spread = {}, {}
        for role, values in roles.items():
            values.sort(key=lambda x: (x["started_unix"], x["mtime_ns"], x["run_id"]))
            chosen[role] = values[-1]
            energies = [x["energy_free_eV_cell"] for x in values]
            duplicate_spread[role] = (max(energies) - min(energies)) * 1000 / 12
        a, b = chosen["Tashiro"], chosen["GFN2"]
        # All compatible duplicate combinations are evidence.  The possible
        # range of E_GFN2 - E_Tashiro is the sum of the two energy ranges, not
        # their maximum; opposite 0.4 meV/atom drifts would otherwise hide 0.8.
        relative_duplicate_spread = sum(duplicate_spread.values())
        if relative_duplicate_spread > TOL:
            rejected.append({"pair_id": pair_id,
                             "reason": "Repeated identical settings permit a relative-energy spread beyond tolerance",
                             "duplicate_spread_meV_atom": duplicate_spread,
                             "relative_duplicate_spread_meV_atom": relative_duplicate_spread})
            continue
        value_eV = (b["energy_free_eV_cell"] - a["energy_free_eV_cell"]) / 12
        pairs.append({
            "pair_id": pair_id, "family_id": a["family_id"], "method": a["method"],
            "family_identity": a["family_identity"], "precision": a["precision"],
            "relative_energy_eV_atom": value_eV, "relative_energy_meV_atom": value_eV * 1000,
            "definition": "(E_GFN2 - E_Tashiro)/12; negative favors GFN2 at this fixed geometry",
            "records": chosen, "duplicate_spread_meV_atom": duplicate_spread,
            "relative_duplicate_spread_meV_atom": relative_duplicate_spread,
            "compatible_duplicate_runs": {role: [
                {"run_id": x["run_id"], "result_sha256": x["result_sha256"],
                 "energy_free_eV_cell": x["energy_free_eV_cell"]} for x in values]
                for role, values in roles.items()}})
    families = []
    for family_id in sorted({pair["family_id"] for pair in pairs}):
        members = [pair for pair in pairs if pair["family_id"] == family_id]
        assessment = assess_family(members)
        families.append({
            "family_id": family_id, "method": members[0]["method"],
            "identity": members[0]["family_identity"], "paired_precisions": [
                pair["pair_id"] for pair in sorted(members, key=refinement_vector)],
            **assessment})
    passed_methods = {family["method"] for family in families if family["passed"]}
    missing_methods = sorted(set(required_methods) - passed_methods)
    complete = bool(pairs) and not pending and not unmatched and not unreadable
    passed = complete and not rejected and not missing_methods and all(
        family["passed"] for family in families)
    status = ("relative_structure_energy_checks_passed" if passed else
              "needs_attention" if rejected or unreadable or any(f["status"] == "needs_attention" for f in families)
              else "pending")
    return {
        "schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": status, "passed": passed, "coverage_complete": complete,
        "quantum_calculations_started": False, "analyzer_sha256": sha(__file__),
        "allowed_stages": sorted(stages), "required_methods": list(required_methods),
        "relative_energy_change_tolerance_meV_atom": TOL, "normalization_atoms": 12,
        "energy_semantics": "Force-consistent electronic free energy in eV/cell; no vibrational free energy.",
        "scope": "Within-method numerical convergence of two frozen C4H8 structures. "
                 "Do not rank different functionals by their absolute or relative energies. "
                 "This does not establish relaxed-phase thermodynamic stability.",
        "pairing_policy": "Same method, all numerical inputs, native grids/effective density cutoffs, "
                          "PAW multiplicities, runner, software and dispersion. Source SHA selects role. "
                          "Latest successful retry selected; other attempts and duplicate aliases retained.",
        "refinement_policy": "Componentwise nondecreasing cutoff/kmesh/gpts/fine grid/effective density cutoffs, "
                             "with at least one strict increase. Three compatible paired precisions, "
                             "stable last two intervals and candidate agreement with the unique highest reference. "
                             "Incomparable settings are never ordered by energy or by descriptive labels.",
        "completion_policy": "All included provenance families must pass and all planned points must have a "
                             "compatible completed record; no incompatible or unfinished family is hidden "
                             "by a passing family of the same method.",
        "references": {role: {key: val for key, val in ref.items() if key != "atoms"}
                       for role, ref in refs.items()},
        "selected_successful_records": len(records), "complete_pairs": len(pairs),
        "missing_or_unpassed_required_methods": missing_methods,
        "pending_points": pending, "unmatched_pairs": unmatched, "rejected_records": rejected,
        "unreadable_inputs": unreadable, "attempt_history": history,
        "records": records, "pairs": pairs, "families": families,
    }


def output(root, report):
    tables, reports = root / "tables", root / "reports"
    tables.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)
    save(tables / "relative_structure_energies.json", report)
    by_family = {x["family_id"]: x for x in report["families"]}
    by_pair = {x["pair_id"]: x for x in report["pairs"]}
    fields = ["method", "family_id", "pair_id", "ecut_eV", "kmesh", "gpts",
              "Tashiro_energy_eV_cell", "GFN2_energy_eV_cell", "relative_energy_meV_atom",
              "change_from_highest_meV_atom", "candidate_accepted", "family_status",
              "Tashiro_result_path", "Tashiro_result_sha256", "GFN2_result_path",
              "GFN2_result_sha256", "Tashiro_source_sha256", "GFN2_source_sha256"]
    with (tables / "relative_structure_energies.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for pair in sorted(report["pairs"], key=lambda x: (x["method"], refinement_vector(x))):
            family = by_family[pair["family_id"]]
            high = by_pair.get(family.get("highest_pair_id"))
            a, b = pair["records"]["Tashiro"], pair["records"]["GFN2"]
            writer.writerow({
                "method": pair["method"], "family_id": pair["family_id"], "pair_id": pair["pair_id"],
                "ecut_eV": pair["precision"]["ecut_eV"], "kmesh": str(pair["precision"]["kmesh"]),
                "gpts": str(pair["precision"]["gpts"]),
                "Tashiro_energy_eV_cell": a["energy_free_eV_cell"],
                "GFN2_energy_eV_cell": b["energy_free_eV_cell"],
                "relative_energy_meV_atom": pair["relative_energy_meV_atom"],
                "change_from_highest_meV_atom": delta(pair, high)["signed_change_meV_atom"] if high else None,
                "candidate_accepted": pair["pair_id"] in family["accepted_candidate_pair_ids"],
                "family_status": family["status"],
                "Tashiro_result_path": a["result_path"], "Tashiro_result_sha256": a["result_sha256"],
                "GFN2_result_path": b["result_path"], "GFN2_result_sha256": b["result_sha256"],
                "Tashiro_source_sha256": a["source_sha256"], "GFN2_source_sha256": b["source_sha256"]})
    lines = [
        "# 两套冻结PE结构的相对能数值收敛", "",
        f"状态：**{report['status']}**；完整配对数：{report['complete_pairs']}；"
        f"全部所需方法通过：{report['passed']}。", "",
        "定义：ΔE结构=(E_GFN2−E_Tashiro)/12，表中转换为meV/atom。"
        "负号表示在该方法、该精度和这两套冻结坐标下，GFN2结构能量较低。"
        "不能据此跨泛函排名，也不能将冻结坐标的能差当作优化后晶相稳定性。", "",
        "每对严格核对方法、所有数值设置、实际网格与密度截断、PAW、运行器、软件及色散参数。"
        "来源文件和精确几何分别保留哈希。至少三个相容递增精度，末两段变化均≤0.5 meV/atom；"
        "候选及其更高设置还须对共同最高精度满足同一门槛。"
        "最高点本身不会被推荐；缺少数据及不可比较的精度分支保留待完成。", "",
        "| 方法 | 截断(eV) | k网格 | 实空间网格 | ΔE结构(meV/atom) | 对最高点变化 | 状态 |",
        "|---|---:|---|---|---:|---:|---|"]
    for pair in sorted(report["pairs"], key=lambda x: (x["method"], refinement_vector(x))):
        family = by_family[pair["family_id"]]
        high = by_pair.get(family.get("highest_pair_id"))
        change = f"{delta(pair, high)['signed_change_meV_atom']:.8g}" if high else "待共同参照"
        p = pair["precision"]
        lines.append(f"| {pair['method']} | {p['ecut_eV']:g} | {p['kmesh']} | {p['gpts']} | "
                     f"{pair['relative_energy_meV_atom']:.9g} | {change} | {family['status']} |")
    if not report["pairs"]:
        lines.append("| 待计算 | — | — | — | — | — | pending |")
    lines += ["", "## 待完成或需要处理", "",
              f"- 尚未通过的方法：{', '.join(report['missing_or_unpassed_required_methods']) or '无'}。",
              f"- 待完成运行：{len(report['pending_points'])}；缺配对：{len(report['unmatched_pairs'])}；"
              f"拒绝记录：{len(report['rejected_records'])}。"]
    for family in report["families"]:
        for reason in family["reasons"]:
            lines.append(f"- {family['method']}：{reason}")
    for row in report["rejected_records"]:
        lines.append(f"- {row.get('run_id', row.get('pair_id', '记录'))}：{row['reason']}")
    lines += ["", "完整输入/结果/源文件路径、哈希、失败尝试与精度分支见"
              "../tables/relative_structure_energies.json；逐精度能量见同名CSV。",
              "", f"生成时间：{report['created_utc']}", ""]
    (reports / "相对结构能收敛.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--stages", nargs="+", choices=sorted(STAGES), default=sorted(STAGES))
    parser.add_argument("--required-methods", nargs="+", default=list(METHODS))
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    report = analyze(root, set(args.stages), args.required_methods)
    output(root, report)
    print(json.dumps({key: report[key] for key in
                      ("status", "passed", "coverage_complete", "complete_pairs",
                       "selected_successful_records", "missing_or_unpassed_required_methods")}))
    if args.require_complete and not report["passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
