"""Verify archived source bytes and selected independently recorded PE invariants.

Uses only the Python standard library. This does not rerun GFN2 or DFT.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Optional verification report outside snapshot/")
    args = parser.parse_args()
    failures = []
    manifest = load(ROOT / "provenance/source_manifest.json")
    paths = set()
    total_bytes = 0
    for entry in manifest["files"]:
        name = entry["path"]
        path = (ROOT / name).resolve()
        if not path.is_relative_to(ROOT / "snapshot"):
            failures.append(f"Outside snapshot: {name}")
            continue
        if name in paths:
            failures.append(f"Duplicate path: {name}")
        paths.add(name)
        if not path.is_file():
            failures.append(f"Missing: {name}")
            continue
        content = path.read_bytes()
        total_bytes += len(content)
        if len(content) != entry["size_bytes"]:
            failures.append(f"Size mismatch: {name}")
        if hashlib.sha256(content).hexdigest() != entry["sha256"]:
            failures.append(f"SHA256 mismatch: {name}")
        if path.suffix == ".json":
            try:
                load(path)
            except Exception as exc:
                failures.append(f"Invalid JSON: {name}: {exc}")
        if path.suffix == ".py":
            try:
                ast.parse(content.decode("utf-8-sig"), filename=name)
            except Exception as exc:
                failures.append(f"Invalid Python: {name}: {exc}")
    actual = {p.relative_to(ROOT).as_posix() for p in (ROOT / "snapshot").rglob("*") if p.is_file()}
    failures.extend(f"Unlisted snapshot file: {p}" for p in sorted(actual - paths))
    calc = ROOT / "snapshot/calculations"
    counts = load(calc / "pe_pnam_candidates_20261020/generation_manifest.json")
    expected = dict(enumerated_count=288, physical_survivor_count=248,
                    unique_initial_count=16, rejected_count=40, exact_duplicates_count=232)
    for key, value in expected.items():
        if counts.get(key) != value:
            failures.append(f"Candidate count {key}: {counts.get(key)} != {value}")
    cold = load(calc / "pe_g4_variable_cell_20261027/final_cold_verification.json")
    if cold.get("passed") is not True or cold.get("independent_density_and_wavefunctions") is not True:
        failures.append("Variable-cell independent cold-start acceptance missing")
    if not 0 <= cold["raw_fmax_eV_A"] <= 0.002:
        failures.append("Variable-cell force threshold exceeded")
    if not 0 <= cold["raw_max_stress_component_kbar"] <= 0.05:
        failures.append("Variable-cell stress threshold exceeded")
    report = dict(passed=not failures, files_verified=len(paths), bytes_verified=total_bytes,
                  candidate_counts={k: counts.get(k) for k in expected},
                  variable_cell_cold_verification=cold,
                  quantum_calculations_rerun=False, failures=failures)
    if args.output:
        out = args.output.resolve()
        if out.is_relative_to(ROOT / "snapshot"):
            parser.error("Write verification reports outside snapshot/")
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    raise SystemExit(0 if report["passed"] else 1)


if __name__ == "__main__":
    main()
