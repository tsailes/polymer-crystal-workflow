"""Resume the saved reference geometry after the candidate releases its slot."""
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def write_new(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2)
        stream.write("\n")


def available_mib():
    lines = Path("/proc/meminfo").read_text().splitlines()
    return int(next(s for s in lines if s.startswith("MemAvailable:")).split()[1]) / 1024


def process_is_live(pid):
    path = Path("/proc") / str(pid) / "status"
    if not path.exists():
        return False
    return not path.read_text().split("State:", 1)[1].splitlines()[0].strip().startswith("Z")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    args = parser.parse_args()
    source = args.input.resolve()
    if source.parent.parent != ROOT / "relaxations" or source.parent.name != "cell_reference_resume_01":
        raise ValueError("Only the prepared reference continuation is authorised here.")
    spec = json.loads(source.read_text())
    if spec["method"] != "PBE-D3(BJ)" or spec["ecut_eV"] != 1800 or spec["kmesh"] != [4, 8, 12]:
        raise ValueError("Reference method or precision changed.")
    if any((source.parent / name).exists() for name in ("process.json", "result.json", "scheduling.json")):
        raise FileExistsError("Continuation has already been scheduled or started.")
    candidate_path = ROOT / "relaxations/cell_candidate/process.json"
    candidate = json.loads(candidate_path.read_text())
    candidate_pid = candidate["pid"]
    expected = candidate["command"]
    if expected[-2:] != [
        str(ROOT / "scripts/run_relax.py"),
        str(ROOT / "relaxations/cell_candidate/input.json"),
    ]:
        raise ValueError("Unexpected candidate task identity.")
    if process_is_live(candidate_pid):
        observed = (Path("/proc") / str(candidate_pid) / "cmdline").read_bytes().rstrip(b"\0").decode().split("\0")
        if observed[1:] != expected[1:]:
            raise ValueError("Candidate PID was reused.")
    scheduled = {
        "time_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "state": "waiting_for_candidate_slot",
        "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "candidate_process_path": str(candidate_path),
        "candidate_pid": candidate_pid,
        "minimum_available_mib": 5000,
        "reason": "Restore the reference optimization without interrupting any active quantum job or overcommitting memory.",
        "restart": "Completed saved geometry; fresh electronic solution and fresh LBFGS history.",
    }
    scheduling_path = source.parent / "scheduling.json"
    write_new(scheduling_path, scheduled)
    print(json.dumps(scheduled), flush=True)
    while True:
        candidate = json.loads(candidate_path.read_text())
        memory = available_mib()
        if "exit_code" in candidate and not process_is_live(candidate_pid) and memory >= 5000:
            break
        time.sleep(10)
    if hashlib.sha256(source.read_bytes()).hexdigest() != scheduled["input_sha256"]:
        raise RuntimeError("Prepared input changed while queued.")
    scheduled.update(state="dispatching", dispatch_time_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                     available_mib=memory, candidate_exit_code=candidate["exit_code"])
    scheduling_path.write_text(json.dumps(scheduled, indent=2) + "\n")
    print(json.dumps(scheduled), flush=True)
    code = subprocess.call([sys.executable, str(ROOT / "scripts/dispatch_relax.py"), str(source)])
    scheduled.update(state="dispatcher_finished", exit_code=code,
                     finish_time_utc=dt.datetime.now(dt.timezone.utc).isoformat())
    scheduling_path.write_text(json.dumps(scheduled, indent=2) + "\n")
    raise SystemExit(code)


if __name__ == "__main__":
    main()
