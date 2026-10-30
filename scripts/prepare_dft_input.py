"""Prepare the archived variable-cell g4 cold single-point input in a fresh folder.

No electronic-structure calculation is launched. Only file paths are relocated;
all scientific parameters and the frozen structure/runner bytes are preserved.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

REPO = Path(__file__).resolve().parents[1]
SOURCE = REPO / "snapshot/calculations/pe_g4_variable_cell_20261027"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("work/g4-cold-point"),
                        help="Fresh output directory, relative to the current directory")
    args = parser.parse_args()
    output = args.output.expanduser().resolve()
    if output.is_relative_to((REPO / "snapshot").resolve()):
        parser.error("Output must be outside the immutable snapshot directory")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.error("Output must be absent or empty; choose a new directory")
    source_input = SOURCE / "cold_verification/input.json"
    source_result = SOURCE / "cold_verification/result.json"
    structure = SOURCE / "runs/attempt_01/final_structure.extxyz"
    runner = SOURCE / "scripts/run_point.py"
    spec = json.loads(source_input.read_text(encoding="utf-8-sig"))
    expected = json.loads(source_result.read_text(encoding="utf-8-sig"))
    runner_hash = sha256(runner)
    structure_hash = sha256(structure)
    if runner_hash != spec["point_runner_sha256"] or runner_hash != expected["runner_sha256"]:
        parser.error("Runner SHA256 differs from the accepted historical input/result")
    if structure_hash != expected["source_geometry_sha256"]:
        parser.error("Structure SHA256 differs from the accepted historical result")
    if expected["status"] != "passed_single_point" or spec != expected["input"]:
        parser.error("Input does not match the archived accepted single-point result")
    new_spec = dict(spec)
    new_spec["structure"] = "structure.extxyz"
    new_spec["point_runner"] = "run_point.py"
    changed = [key for key in spec if spec[key] != new_spec[key]]
    if set(changed) != {"structure", "point_runner"}:
        parser.error("Unexpected changes to the historical input")
    output.mkdir(parents=True, exist_ok=True)
    copies = [(structure, "structure.extxyz"), (runner, "run_point.py"),
              (source_result, "expected_result.json")]
    for source, name in copies:
        shutil.copyfile(source, output / name)
    (output / "input.json").write_text(json.dumps(new_spec, indent=2) + "\n", encoding="utf-8")
    record = {
        "purpose": "Prepare frozen variable-cell g4 cold single-point rerun",
        "quantum_calculations_started": 0,
        "input_source": source_input.relative_to(REPO).as_posix(),
        "input_source_sha256": sha256(source_input),
        "prepared_input_sha256": sha256(output / "input.json"),
        "changed_fields": changed,
        "all_scientific_parameters_preserved": True,
        "files": [{"source": source.relative_to(REPO).as_posix(), "output": name,
                   "sha256": sha256(source)} for source, name in copies],
        "mpi_ranks_in_historical_input": spec["mpi_ranks"],
        "historical_result_is_not_a_new_rerun_result": True,
    }
    (output / "preparation_manifest.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"prepared": str(output), "quantum_calculations_started": 0,
                      "changed_fields": changed}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
