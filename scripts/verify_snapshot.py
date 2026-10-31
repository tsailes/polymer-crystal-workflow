"""Verify archived source bytes and selected PE and PA6 scientific record invariants.

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
    manifests = [
        ("snapshot", load(ROOT / "provenance/source_manifest.json")),
        ("snapshot_pa6", load(ROOT / "provenance/pa6_source_manifest.json")),
    ]
    paths = set()
    total_bytes = 0
    for archive, manifest in manifests:
        manifest_paths = set()
        for entry in manifest["files"]:
            name = entry["path"]
            path = (ROOT / name).resolve()
            if not path.is_relative_to(ROOT / archive):
                failures.append(f"Outside {archive}: {name}")
                continue
            if name in paths:
                failures.append(f"Duplicate path: {name}")
            paths.add(name)
            manifest_paths.add(name)
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
        actual = {p.relative_to(ROOT).as_posix() for p in (ROOT / archive).rglob("*") if p.is_file()}
        failures.extend(f"Unlisted {archive} file: {p}" for p in sorted(actual - manifest_paths))
        if manifest.get("file_count") != len(manifest_paths):
            failures.append(f"Manifest count mismatch: {archive}")
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
    pa6 = ROOT / "snapshot_pa6/calculations/gfn2_single_pilot_20261030"
    selection = load(pa6 / "input/selection_record.json")
    diagnosis = load(pa6 / "diagnosis_20261030/summary.json")
    if selection.get("candidate_count") != 1 or selection.get("atom_count") != 152:
        failures.append("PA6 pilot scope differs from one 152-atom candidate")
    if selection.get("hydrogen_bonds_evaluated") is not False:
        failures.append("PA6 pilot selection unexpectedly uses hydrogen bonds")
    if diagnosis["tight_geometry"]["strict_bond_count"] != 0:
        failures.append("PA6 recorded no-network diagnostic changed")
    if diagnosis["tight_endpoint_b_doubling"]["b_doubled_fmax_eV_A"] <= 0.005:
        failures.append("PA6 recorded periodic-sensitivity failure changed")
    report = dict(passed=not failures, files_verified=len(paths), bytes_verified=total_bytes,
                  candidate_counts={k: counts.get(k) for k in expected},
                  variable_cell_cold_verification=cold,
                  archives={archive: manifest["file_count"] for archive, manifest in manifests},
                  pa6_pilot=dict(candidate_count=selection["candidate_count"],
                                 atom_count=selection["atom_count"],
                                 strict_bond_count=diagnosis["tight_geometry"]["strict_bond_count"],
                                 periodic_convergence_passed=False),
                  quantum_calculations_rerun=False, failures=failures)
    if args.output:
        out = args.output.resolve()
        if any(out.is_relative_to(ROOT / archive) for archive, _ in manifests):
            parser.error("Write verification reports outside source snapshot directories")
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    raise SystemExit(0 if report["passed"] else 1)


if __name__ == "__main__":
    main()
