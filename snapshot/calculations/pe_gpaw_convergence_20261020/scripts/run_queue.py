"""Sequential coordinator: immutable point inputs and full child logs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
ENV = Path("/opt/pe_dft/envs/pe-gpaw-20261020")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2) + "\n")


def aggregate():
    rows = []
    iterations = []
    for d in sorted((ROOT / "runs").glob("*")):
        if not d.is_dir() or not (d / "input.json").exists():
            continue
        s = json.loads((d / "input.json").read_text())
        r = json.loads((d / "result.json").read_text()) \
            if (d / "result.json").exists() else {}
        if not r and (d / "interruption_record.json").exists():
            r = json.loads((d / "interruption_record.json").read_text())
        total = r.get("total", {})
        dft = r.get("dft", {})
        disp = r.get("dispersion", {})
        row = {
            "run_id": d.name, "stage": s.get("stage"),
            "method": s["method"], "ecut_eV": s["ecut_eV"],
            "kmesh": "x".join(map(str, s["kmesh"])),
            "gpts": "x".join(map(str, s.get("gpts", []))) or "automatic",
            "xc_stencil": s.get("xc_stencil", 2),
            "geometry_sha256": r.get("source_geometry_sha256"),
            "actual_geometry_sha256": hashlib.sha256(json.dumps({
                "formula": r.get("formula"), "cell": r.get("cell_A"),
                "positions": r.get("positions_A")}, sort_keys=True).encode()).hexdigest()
                if r.get("positions_A") is not None else None,
            "displacement": json.dumps(s.get("displacement")),
            "strain": json.dumps(s.get("strain")),
            "dispersion_parameters": json.dumps(
                disp.get("calculator_parameters", {}), sort_keys=True),
            "status": r.get("status", "no_result"),
            "energy_free_eV_cell": total.get("energy_free_eV"),
            "energy_extrapolated_eV_cell": total.get("energy_extrapolated_eV"),
            "energy_eV_atom": total.get("energy_eV_atom"),
            "energy_eV_CH2": total.get("energy_eV_CH2"),
            "energy_Ry_cell": total.get("energy_Ry"),
            "dft_energy_free_eV": dft.get("energy_free_eV"),
            "dispersion_energy_eV": disp.get("energy_eV"),
            "max_force_norm_eV_A": total.get("max_force_norm_eV_A"),
            "SCF_iterations": r.get("SCF_iterations"),
            "wall_seconds": r.get("wall_seconds"),
            "mpi_ranks": r.get("mpi_world_size"),
            "result_json": str(d / "result.json"),
            "error": r.get("error", r.get("reason"))
        }
        for i, name in enumerate(("xx", "yy", "zz", "yz", "xz", "xy")):
            vals = total.get("stress_kbar")
            row["stress_" + name + "_kbar"] = vals[i] if vals else None
        rows.append(row)
        p = d / "scf_iterations.jsonl"
        if p.exists():
            for line in p.read_text().splitlines():
                try:
                    iterations.append({"run_id": d.name, **json.loads(line)})
                except json.JSONDecodeError:
                    pass
    (ROOT / "tables").mkdir(exist_ok=True)
    dump(ROOT / "tables" / "energies_all_runs.json", rows)
    if rows:
        with (ROOT / "tables" / "energies_all_runs.csv").open(
                "w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    if iterations:
        fields = sorted(set().union(*(d.keys() for d in iterations)))
        with (ROOT / "tables" / "scf_iterations_all.csv").open(
                "w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fields)
            w.writeheader()
            w.writerows(iterations)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("queue", nargs="?")
    ap.add_argument("--aggregate-only", action="store_true")
    args = ap.parse_args()
    if args.aggregate_only:
        aggregate()
        return
    queuepath = Path(args.queue).resolve()
    queue = json.loads(queuepath.read_text())
    # Freeze the complete queue's implementation and protocol at dispatch,
    # even if later independent work edits the source files.
    runner_bytes = (ROOT / "scripts" / "run_point.py").read_bytes()
    runner_hash = hashlib.sha256(runner_bytes).hexdigest()
    protocol_bytes = (ROOT / "protocol.json").read_bytes()
    protocol_hash = hashlib.sha256(protocol_bytes).hexdigest()
    environment = os.environ.copy()
    environment.update({
        "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
        "GPAW_MPI_BACKEND": "cgpaw",
        "OMPI_ALLOW_RUN_AS_ROOT": "1", "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM": "1",
        "PATH": str(ENV / "bin") + ":" + environment.get("PATH", "")
    })
    if queue.get("setup_path"):
        environment["GPAW_SETUP_PATH"] = queue["setup_path"]
    events = queuepath.with_suffix(".events.jsonl")
    for item in queue["points"]:
        ident = item["id"]
        work = ROOT / "runs" / ident
        if work.exists():
            p = work / "result.json"
            if p.exists() and json.loads(p.read_text()).get("status") == "passed_single_point":
                print(f"SKIP passed {ident}", flush=True)
                continue
            # Failed, interrupted and uncertain runs are retained without overwriting.
            n = 2
            while (ROOT / "runs" / (ident + f"__attempt{n}")).exists():
                n += 1
            work = ROOT / "runs" / (ident + f"__attempt{n}")
        work.mkdir(parents=True)
        spec = {**queue.get("defaults", {}), **item}
        spec.setdefault("structure", str(ROOT / "reference" / "POSCAR"))
        spec["protocol_sha256"] = protocol_hash
        spec["queue_sha256"] = sha(queuepath)
        spec["runner_sha256_at_submission"] = runner_hash
        dump(work / "input.json", spec)
        (work / "run_point_snapshot.py").write_bytes(runner_bytes)
        (work / "protocol_snapshot.json").write_bytes(protocol_bytes)
        cmd = [str(ENV / "bin" / "mpiexec"), "--allow-run-as-root",
               "--bind-to", "none", "-n", str(spec.get("mpi_ranks", 4)),
               str(ENV / "bin" / "python"), str(work / "run_point_snapshot.py"),
               str(work / "input.json")]
        started = time.time()
        with (work / "stdout.log").open("wb") as out, \
                (work / "stderr.log").open("wb") as err:
            child = subprocess.Popen(cmd, stdout=out, stderr=err, env=environment)
            record = {"id": work.name, "command": cmd, "pid": child.pid,
                      "start_unix": started, "environment_overrides": {
                          k: environment[k] for k in
                          ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                           "GPAW_SETUP_PATH") if k in environment}}
            dump(work / "process.json", record)
            print(f"START {work.name} pid={child.pid}", flush=True)
            code = child.wait()
        record.update({"exit_code": code, "wall_seconds": time.time() - started})
        dump(work / "process.json", record)
        with events.open("a") as f:
            f.write(json.dumps(record) + "\n")
        aggregate()
        print(f"END {work.name} code={code} seconds={record['wall_seconds']:.1f}",
              flush=True)
        if code and queue.get("stop_on_failure", True):
            sys.exit(code)


if __name__ == "__main__":
    main()
