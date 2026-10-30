"""Reproduce the accepted protocol in a new timestamped run (never overwrite a run)."""
from pathlib import Path
import json
import os
import subprocess
import sys

root = Path(__file__).resolve().parent
workflow = root / "workflow.py"
model = root.parent / "model" / "pe_initial.extxyz"
environment = {**os.environ, "OMP_NUM_THREADS": "4", "OPENBLAS_NUM_THREADS": "1",
               "MKL_NUM_THREADS": "1", "PYTHONUNBUFFERED": "1"}


def call(script, *args, log=None):
    command = [sys.executable, "-u", str(script), *map(str, args)]
    if log:
        with Path(log).open("x", encoding="utf-8") as f:
            subprocess.run(command, env=environment, stdout=f, stderr=subprocess.STDOUT, check=True)
    else:
        subprocess.run(command, env=environment, check=True)


call(workflow, "init")
run = Path(json.loads((root / "active_run.json").read_text())["run_directory"])
for name, grids, accuracy in [
    ("initial_sampling", ["1x1x6", "1x1x8", "1x1x10", "1x1x12", "2x1x12", "1x2x12", "2x2x12"], "0.01"),
    ("transverse_sampling", ["2x2x8", "3x2x8", "2x3x8", "3x3x8"], "0.01"),
    ("refined_sampling", ["2x4x8", "2x3x12", "2x3x6"], "0.01"),
    ("initial_accuracy_check", ["2x3x8"], "0.001"),
]:
    call(workflow, "sample", "--run", run, "--input", model, "--name", name,
         "--grids", *grids, "--accuracy", accuracy)
call(root / "check_folded_force.py", "--input", model, "--output", run / "folded_force_check",
     "--repeat", "2", "3", "8", log=run / "folded_force_check_native.log")
call(root / "assess_results.py", "--run", run, "--phase", "initial")
call(workflow, "optimize", "--run", run, "--input", model, "--name", "relaxation_2x3x8",
     "--repeat", "2", "3", "8", "--steps", "300", log=run / "relaxation_2x3x8_native.log")
final = run / "relaxation_2x3x8" / "pe_relaxed.extxyz"
for name, grids, accuracy in [
    ("final_sampling", ["2x3x8", "3x3x8", "2x4x8", "2x3x12"], "0.01"),
    ("final_combined", ["3x4x10"], "0.01"),
    ("final_accuracy_check", ["2x3x8"], "0.001"),
]:
    call(workflow, "sample", "--run", run, "--input", final, "--name", name,
         "--grids", *grids, "--accuracy", accuracy)
call(root / "assess_results.py", "--run", run, "--phase", "final")
call(root / "plot_results.py", "--run", run, "--final", final)
print("COMPLETED", run, flush=True)
