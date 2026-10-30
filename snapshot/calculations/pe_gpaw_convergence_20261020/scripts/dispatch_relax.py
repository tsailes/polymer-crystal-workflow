"""Launch one explicitly prepared relaxation and preserve its process logs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
ENV = Path("/opt/pe_dft/envs/pe-gpaw-20261020")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    args = parser.parse_args()
    source = Path(args.input).resolve()
    spec = json.loads(source.read_text(encoding="utf-8-sig"))
    work = source.parent
    if not work.is_relative_to(ROOT / "relaxations"):
        raise ValueError("Input must be inside this batch's relaxations directory")
    for name in ("process.json", "result.json", "stdout.log", "stderr.log"):
        if (work / name).exists():
            raise FileExistsError(f"Existing attempt retained: {work / name}")
    environment = os.environ.copy()
    overrides = {
        "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
        "GPAW_MPI_BACKEND": "cgpaw",
        "OMPI_ALLOW_RUN_AS_ROOT": "1", "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM": "1",
    }
    environment.update(overrides)
    environment["PATH"] = str(ENV / "bin") + ":" + environment.get("PATH", "")
    command = [str(ENV / "bin/mpiexec"), "--allow-run-as-root", "--bind-to",
               "none", "-n", str(spec["mpi_ranks"]), str(ENV / "bin/python"),
               str(ROOT / "scripts/run_relax.py"), str(source)]
    started = time.time()
    with (work / "stdout.log").open("xb") as stdout, \
            (work / "stderr.log").open("xb") as stderr:
        process = subprocess.Popen(command, stdout=stdout, stderr=stderr,
                                   env=environment)
        record = {
            "command": command, "pid": process.pid, "start_unix": started,
            "environment_overrides": overrides,
            "dispatcher_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        }
        (work / "process.json").write_text(json.dumps(record, indent=2) + "\n")
        print(f"START {work.name} pid={process.pid}", flush=True)
        code = process.wait()
    record.update(exit_code=code, wall_seconds=time.time() - started)
    (work / "process.json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"END {work.name} code={code} seconds={record['wall_seconds']:.1f}",
          flush=True)
    raise SystemExit(code)


if __name__ == "__main__":
    main()
