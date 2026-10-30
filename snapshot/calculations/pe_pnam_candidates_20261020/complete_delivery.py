"""Finalize completed candidates in dependency order; stop on any failed check."""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def run(name, *args):
    subprocess.run([sys.executable, str(ROOT / name), *map(str, args)], check=True)


def main():
    run("assess_candidates.py", "--require-complete")
    run("verify_representatives.py")
    summary = json.loads((ROOT / "optimization_summary.json").read_text())
    models = [ROOT / "representatives" / g["cluster_id"] / "representative.extxyz"
              for g in summary["groups"]]
    run("xrd_benchmark.py", "--structures", *models, "--space-group", "Pnam",
        "--output-dir", ROOT / "xrd_benchmark_final")
    run("write_report.py")
    run("plot_candidates.py")
    run("audits/validate_delivery.py")
    print("Numerical and file checks complete. Visually inspect the final plots before packaging.",
          flush=True)


if __name__ == "__main__":
    main()
