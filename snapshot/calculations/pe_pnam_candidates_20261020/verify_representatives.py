"""Denser Gamma-supercell checks of the exact delivered representative geometry."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import importlib.util
import json
from pathlib import Path
import shutil
import sys

import numpy as np
from ase.io import read

from generate_candidates import ROOT, dump, sha

SCRIPT_SHA256 = sha(__file__)
WORKFLOW = ROOT.parent/"pe_gfn2_first_run_20261020/convergence_relaxation/workflow.py"
spec = importlib.util.spec_from_file_location("verified_workflow", WORKFLOW)
workflow = importlib.util.module_from_spec(spec)
spec.loader.exec_module(workflow)


def one(group, summary, cache_only=False):
    cid, gid = group["representative_id"], group["cluster_id"]
    row = next(r for r in summary["candidates"] if r["candidate_id"] == cid)
    fine = json.loads(Path(row["final_result_path"]).read_text())
    model = (Path(row["final_result_path"]).parent/"final.extxyz" if cache_only else
             ROOT/"representatives"/gid/"representative.extxyz")
    atoms = read(model)
    source_atoms = read(Path(row["final_result_path"]).parent/"final.extxyz")
    if (not np.array_equal(atoms.numbers, source_atoms.numbers)
            or not np.array_equal(atoms.pbc, source_atoms.pbc)
            or not np.allclose(atoms.cell, source_atoms.cell, rtol=0, atol=1e-12)
            or not np.allclose(atoms.positions, source_atoms.positions, rtol=0, atol=1e-10)):
        raise ValueError("Representative differs from the source final geometry or atom order")
    folder = ROOT/("candidate_sampling_checks" if cache_only else "representative_convergence")/(cid if cache_only else gid)
    folder.mkdir(parents=True, exist_ok=True)
    case = folder/"3x4x10"
    cache = ROOT/"candidate_sampling_checks"/cid/"3x4x10"
    if not cache_only and not case.exists() and (cache/"result.json").exists():
        cached_atoms = read(cache/"input.extxyz")
        cached_result = json.loads((cache/"result.json").read_text())
        if (np.array_equal(cached_atoms.numbers, atoms.numbers)
                and np.array_equal(cached_atoms.pbc, atoms.pbc)
                and np.allclose(cached_atoms.cell, atoms.cell, rtol=0, atol=1e-12)
                and np.allclose(cached_atoms.positions, atoms.positions, rtol=0, atol=1e-10)
                and cached_result.get("settings") == dict(workflow.SETTINGS, accuracy=.001)
                and cached_result.get("repeats") == [3, 4, 10]
                and cached_result.get("input_sha256") == sha(cache/"input.extxyz")
                and cached_result.get("script_sha256") == sha(WORKFLOW)):
            shutil.copytree(cache, case)
    if (case/"result.json").exists():
        dense = json.loads((case/"result.json").read_text())
        previous = read(case/"input.extxyz")
        if (not np.array_equal(previous.numbers, atoms.numbers)
                or not np.array_equal(previous.pbc, atoms.pbc)
                or not np.allclose(previous.cell, atoms.cell, atol=1e-12, rtol=0)
                or not np.allclose(previous.positions, atoms.positions, atol=1e-10, rtol=0)
                or dense.get("settings") != dict(workflow.SETTINGS, accuracy=.001)
                or dense.get("repeats") != [3, 4, 10]
                or dense.get("input_sha256") != sha(case/"input.extxyz")
                or dense.get("script_sha256") != sha(WORKFLOW)):
            raise ValueError("Existing representative check belongs to another geometry or configuration")
    else:
        dense = workflow.launch_case(case, atoms, (3, 4, 10), accuracy=.001)
    delta_e = dense["energy_per_CH2_eV"]-fine["energy_eV_per_CH2"]
    delta_f = float(np.linalg.norm(np.asarray(dense["forces_eV_per_A"]) -
                                   np.asarray(fine["raw_forces_eV_A"]), axis=1).max())
    result = dict(group_id=gid, candidate_id=cid, representative_file=str(model),
                  representative_sha256=sha(model), dense_input_sha256=sha(case/"input.extxyz"),
                  source_final_result_sha256=sha(row["final_result_path"]),
                  fine_repeat=[2, 3, 8], dense_repeat=[3, 4, 10], dense_natoms=1440,
                  comparison_accuracy=.001,
                  energy_change_meV_per_CH2=delta_e*1000,
                  max_force_vector_change_eV_A=delta_f,
                  fine_fmax_eV_A=fine["raw_fmax_eV_A"],
                  dense_fmax_eV_A=dense["max_force_norm_eV_per_A"],
                  replica_force_deviation_eV_A=dense["replica_force_max_deviation_eV_per_A"],
                  sampling_change_passed=bool(abs(delta_e)<.001 and delta_f<.005),
                  dense_force_target_passed=bool(dense["max_force_norm_eV_per_A"]<.01),
                  replica_equivalence_passed=bool(dense["replica_force_max_deviation_eV_per_A"]<1e-5),
                  thresholds=dict(delta_E_meV_per_CH2=1, delta_F_eV_A=.005, dense_fmax_eV_A=.01,
                                  replica_force_eV_A=1e-5),
                  script_sha256=SCRIPT_SHA256, workflow_sha256=sha(WORKFLOW))
    result["passed"] = (result["sampling_change_passed"] and result["dense_force_target_passed"]
                        and result["replica_equivalence_passed"])
    dump(folder/"comparison.json", result)
    print(json.dumps(result), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", nargs="+",
                        help="Check completed candidate geometries in a reusable per-candidate cache")
    args = parser.parse_args()
    summary = json.loads((ROOT/"optimization_summary.json").read_text())
    if args.candidates:
        rows = {r["candidate_id"]: r for r in summary["candidates"]}
        if len(set(args.candidates)) != len(args.candidates):
            raise ValueError("Duplicate candidate request")
        for cid in args.candidates:
            if rows[cid]["status"] != "preoptimized" or not rows[cid]["geometry_passed"]:
                raise ValueError(f"Cannot check unfinished or excluded candidate: {cid}")
        groups = [dict(representative_id=cid, cluster_id=f"candidate_{cid}") for cid in args.candidates]
        with ThreadPoolExecutor(max_workers=2) as pool:
            results = list(pool.map(lambda g: one(g, summary, cache_only=True), groups))
        dump(ROOT/"candidate_sampling_checks/summary.json",
             dict(passed=all(r["passed"] for r in results), results=results,
                  scope="Completed candidate geometry checks; final representatives may differ"))
        if not all(r["passed"] for r in results):
            raise SystemExit("One or more candidate sampling checks did not pass")
        return
    if not summary["all_preoptimization_attempts_completed"]:
        raise RuntimeError("Finalize candidate assessment first")
    snapshot = ROOT/"representative_convergence/source_snapshot"
    snapshot.mkdir(parents=True, exist_ok=True)
    for source in (Path(__file__), WORKFLOW, ROOT/"optimization_summary.json"):
        shutil.copy2(source, snapshot/source.name)
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(lambda g: one(g, summary), summary["groups"]))
    dump(ROOT/"representative_convergence/summary.json",
         dict(passed=bool(results) and all(r["passed"] for r in results), results=results,
              scope="Exact delivered representatives only; no transfer of force convergence to every cluster member"))
    if not results or not all(r["passed"] for r in results):
        raise SystemExit("One or more representative convergence checks did not pass")


if __name__ == "__main__":
    main()
