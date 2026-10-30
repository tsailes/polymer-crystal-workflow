"""Re-run the frozen geometric generator without changing the source snapshot.

This runs no GFN2 or DFT calculation. Install requirements-analysis.txt first.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

REPO = Path(__file__).resolve().parents[1]
SOURCE = REPO / "snapshot/calculations/pe_pnam_candidates_20261020"
EXPECTED_COUNTS = {
    "enumerated_count": 288,
    "physical_survivor_count": 248,
    "unique_initial_count": 16,
    "rejected_count": 40,
    "exact_duplicates_count": 232,
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def identities(manifest: dict) -> list[dict]:
    """Compare all enumerated identities and classifications, not just counts."""
    return [
        {key: entry.get(key) for key in (
            "candidate_id", "parameters", "passed", "status", "representative_id"
        )}
        for entry in manifest["records"]
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("work/candidates"),
                        help="Fresh output directory, relative to the current directory")
    args = parser.parse_args()
    output = args.output.expanduser().resolve()
    if output.is_relative_to((REPO / "snapshot").resolve()):
        parser.error("Output must be outside the immutable snapshot directory")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.error("Output must be absent or empty; choose a new directory")
    generator = SOURCE / "generate_candidates.py"
    manifest_path = SOURCE / "generation_manifest.json"
    archived = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    if any(archived[key] != value for key, value in EXPECTED_COUNTS.items()):
        parser.error("Archived candidate counts differ from this release's baseline")
    versions = {}
    for name in ("numpy", "scipy", "ase", "gemmi", "spglib"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            parser.error(f"Missing {name}; install requirements-analysis.txt")
    source_hash = sha256(generator)
    protocol = json.loads((SOURCE / "protocol.json").read_text(encoding="utf-8-sig"))
    if source_hash != protocol["generator_sha256"]:
        parser.error("Generator hash does not match the frozen protocol")
    output.mkdir(parents=True, exist_ok=True)
    copied = output / generator.name
    shutil.copyfile(generator, copied)
    started = time.perf_counter()
    completed = subprocess.run([sys.executable, str(copied)], cwd=output, check=False)
    report = {
        "purpose": "Geometric candidate reproduction only; no quantum calculations",
        "quantum_calculations_started": 0,
        "source_generator_sha256": source_hash,
        "archived_manifest_sha256": sha256(manifest_path),
        "python": sys.version,
        "platform": platform.platform(),
        "package_versions": versions,
        "exit_code": completed.returncode,
        "wall_seconds": time.perf_counter() - started,
        "expected_counts": EXPECTED_COUNTS,
        "passed": False,
    }
    if completed.returncode == 0:
        fresh = json.loads((output / "generation_manifest.json").read_text(encoding="utf-8"))
        counts = {key: fresh.get(key) for key in EXPECTED_COUNTS}
        report.update({
            "observed_counts": counts,
            "counts_match": counts == EXPECTED_COUNTS,
            "all_candidate_identities_and_classifications_match": identities(fresh) == identities(archived),
            "unique_candidate_ids_match": [r["candidate_id"] for r in fresh["unique_survivors"]]
            == [r["candidate_id"] for r in archived["unique_survivors"]],
            "source_generator_unchanged": sha256(generator) == source_hash,
            "copied_generator_unchanged": sha256(copied) == source_hash,
        })
        report["passed"] = all(report[key] for key in (
            "counts_match", "all_candidate_identities_and_classifications_match",
            "unique_candidate_ids_match", "source_generator_unchanged", "copied_generator_unchanged"
        ))
    report_path = output / "reproduction_verification.json"
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"passed": report["passed"], "report": str(report_path)}, ensure_ascii=False))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
