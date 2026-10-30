"""Freeze the declared PE geometries without launching an energy calculation.

Run with the existing GPAW/ASE WSL Python.  Only this calibration directory is
written.  Already frozen scientific files are verified rather than overwritten.
"""
from __future__ import annotations

from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read, write
import spglib


ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parents[1]
CANDIDATES = PROJECT / "calculations/pe_pnam_candidates_20261020"
SUMMARY = PROJECT / "reports/PE计算实验报告_20261021/data/gfn2_compact_summary.json"
CELL = np.diag([7.432, 4.929, 2.545])


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_bytes(value) -> bytes:
    return (json.dumps(value, ensure_ascii=False, indent=2,
                       allow_nan=False) + "\n").encode("utf-8")


def create_or_verify(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != content:
            raise RuntimeError(f"Refusing to overwrite a differing frozen file: {path}")
        return
    path.write_bytes(content)


def pnam_residual(atoms: Atoms) -> float:
    operations = spglib.get_symmetry_from_database(297)
    fractional = atoms.get_scaled_positions(wrap=False)
    numbers = atoms.numbers
    residual = 0.0
    for rotation, translation in zip(operations["rotations"], operations["translations"]):
        transformed = fractional @ np.asarray(rotation).T + translation
        for number in set(numbers):
            source = transformed[numbers == number]
            target = fractional[numbers == number]
            delta = source[:, None, :] - target[None, :, :]
            delta -= np.rint(delta)
            distances = np.linalg.norm(delta @ atoms.cell.array, axis=2)
            residual = max(residual, float(distances.min(axis=1).max()))
    return residual


def freeze(entry: dict) -> dict:
    source = Path(entry["source_path"])
    source_sha = digest(source)
    original = read(source)
    assert len(original) == 12 and Counter(original.get_chemical_symbols()) == {"C": 4, "H": 8}
    assert np.all(original.pbc), (source, original.pbc)
    assert np.allclose(original.cell.array, CELL, rtol=0, atol=1e-10), source
    assert np.isfinite(original.positions).all()
    # Deliberately drop attached SinglePointCalculator, energies, velocities,
    # constraints and info.  The calibration consumes geometry only.
    clean = Atoms(numbers=original.numbers.copy(), positions=original.positions.copy(),
                  cell=original.cell.array.copy(), pbc=original.pbc.copy())
    folder = ROOT / "structures" / entry["id"]
    poscar = folder / "POSCAR"
    stream = io.StringIO()
    write(stream, clean, format="vasp", direct=True, sort=False, vasp5=True)
    create_or_verify(poscar, stream.getvalue().encode("utf-8"))
    frozen = read(poscar)
    assert np.array_equal(frozen.numbers, original.numbers)
    assert np.allclose(frozen.cell.array, original.cell.array, rtol=0, atol=1e-13)
    coordinate_delta = float(np.max(np.abs(frozen.positions - original.positions)))
    assert coordinate_delta < 5e-13, (source, coordinate_delta)
    assert np.all(frozen.pbc) and frozen.calc is None and not frozen.info
    residual = pnam_residual(frozen)
    assert residual < 1e-5, (source, residual)
    dataset = spglib.get_symmetry_dataset((frozen.cell.array,
                                         frozen.get_scaled_positions(), frozen.numbers),
                                        symprec=1e-5)
    assert dataset is not None and int(dataset.number) in {62, 63}
    geometry = {
        "numbers": frozen.numbers.tolist(),
        "symbols": frozen.get_chemical_symbols(),
        "positions_A": frozen.positions.tolist(),
        "cell_A": frozen.cell.array.tolist(),
        "pbc": frozen.pbc.tolist(),
    }
    canonical = json.dumps(geometry, sort_keys=True, separators=(",", ":"),
                           allow_nan=False).encode("utf-8")
    geometry_sha = hashlib.sha256(canonical).hexdigest()
    geometry_path = folder / "geometry.json"
    create_or_verify(geometry_path, json_bytes({
        "id": entry["id"], "geometry": geometry,
        "geometry_sha256": geometry_sha,
        "source_path": str(source), "source_sha256": source_sha,
        "structure_path": str(poscar), "structure_sha256": digest(poscar),
        "source_to_frozen_max_cartesian_difference_A": coordinate_delta,
        "Pnam_Hall297_max_operation_residual_A": residual,
        "detected_spacegroup_number": int(dataset.number),
        "electronic_calculator_removed": True,
        "old_source_modified": False,
    }))
    assert digest(source) == source_sha
    return {
        **entry,
        "structure": str(poscar),
        "structure_sha256": digest(poscar),
        "geometry_path": str(geometry_path),
        "geometry_sha256": geometry_sha,
        "source_sha256": source_sha,
        "formula": frozen.get_chemical_formula(),
        "natoms": len(frozen),
        "N_CH2": 4,
        "cell_lengths_A": frozen.cell.lengths().tolist(),
        "cell_A": frozen.cell.array.tolist(),
        "PBC": frozen.pbc.tolist(),
        "source_to_frozen_max_cartesian_difference_A": coordinate_delta,
        "Pnam_Hall297_max_operation_residual_A": residual,
        "detected_spacegroup_number": int(dataset.number),
        "geometry_only": True,
    }


def main():
    summary = json.loads(SUMMARY.read_text(encoding="utf-8-sig"))
    records = {record["candidate_id"]: record for record in summary["candidate_final_records"]}
    groups = {group["cluster_id"]: group for group in summary["groups"]}
    proposals = []
    for group_number in range(1, 5):
        group_id = f"group_{group_number:03d}"
        group = groups[group_id]
        candidate = group["representative_id"]
        proposals.append({
            "id": f"final_{group_id}",
            "source_path": str(CANDIDATES / "representatives" / group_id / "POSCAR"),
            "role": "preoptimized_representative",
            "stratum": "preoptimized_representatives",
            "source_candidate_id": candidate,
            "initial_angle_deg": records[candidate]["initial_angle_deg"],
            "group_id": group_id,
            "group_member_count": len(group["members"]),
            "split_block": "same_preoptimized_basin",
            "structural_family": "same_preoptimized_basin",
            "source_status": records[candidate]["status"],
            "screening_eligible": True,
            "independent_polymorph_claim": False,
        })
    for index, angle in enumerate(range(0, 80, 5)):
        candidate = f"cand_{1 + index * 8:04d}"
        record = records[candidate]
        assert record["initial_angle_deg"] == angle
        proposals.append({
            "id": f"initial_{candidate}",
            "source_path": str(CANDIDATES / "candidates" / candidate / "POSCAR"),
            "role": "initial_grid", "stratum": "initial_angular_grid",
            "source_candidate_id": candidate, "initial_angle_deg": angle,
            "group_id": record.get("group_id"),
            "group_member_count": None,
            "split_block": "block" + "ABCD"[index // 4],
            "structural_family": "same_rigid_chain_angular_grid",
            "source_status": "initial_geometrically_valid_candidate",
            "screening_eligible": False,
            "independent_polymorph_claim": False,
        })
    unstable = records["cand_0001"]
    assert unstable["status"] == "unstable_stationary" and unstable["excluded_from_retained_pool"]
    proposals.append({
        "id": "unstable_cand_0001", "source_path": unstable["final_model_path"],
        "role": "unstable_control", "stratum": "unstable_control",
        "source_candidate_id": "cand_0001", "initial_angle_deg": 0,
        "group_id": None, "group_member_count": None,
        "split_block": "negative_control", "structural_family": "zero_angle_stationary",
        "source_status": "unstable_stationary", "screening_eligible": False,
        "independent_polymorph_claim": False,
    })
    assert len(proposals) == len({entry["id"] for entry in proposals}) == 21
    hashes_before = {entry["source_path"]: digest(Path(entry["source_path"])) for entry in proposals}
    samples = [freeze(entry) for entry in proposals]
    record = {
        "schema_version": 1,
        "reference_id": "final_group_001",
        "scope": "Frozen geometries only; static fixed Tashiro cell, C4H8, no new quantum calculation",
        "sampling_note": "16 angular-grid initial geometries, 4 near-minimum representative groups, and 1 unstable control are separate strata. The 4 final groups are not 4 independent polymorphs.",
        "primary_authority": "Each structure path points to the full-precision frozen POSCAR; geometry.json is its read-back audit.",
        "origin_summary": str(SUMMARY), "origin_summary_sha256": digest(SUMMARY),
        "samples": samples,
    }
    create_or_verify(ROOT / "samples.json", json_bytes(record))
    scalar_fields = ["id", "role", "stratum", "source_candidate_id", "initial_angle_deg",
                     "group_id", "group_member_count", "split_block", "structural_family",
                     "screening_eligible", "formula", "natoms", "N_CH2", "structure",
                     "structure_sha256", "geometry_path", "geometry_sha256", "source_path",
                     "source_sha256", "source_to_frozen_max_cartesian_difference_A",
                     "Pnam_Hall297_max_operation_residual_A", "detected_spacegroup_number"]
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=scalar_fields, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(samples)
    create_or_verify(ROOT / "samples.csv", stream.getvalue().encode("utf-8-sig"))
    hashes_after = {source: digest(Path(source)) for source in hashes_before}
    assert hashes_after == hashes_before
    audit = {
        "status": "passed_geometry_freeze",
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "sample_count": len(samples),
        "role_counts": dict(Counter(sample["role"] for sample in samples)),
        "unique_frozen_geometry_hashes": len({sample["geometry_sha256"] for sample in samples}),
        "old_sources_checked": len(hashes_before), "old_sources_changed": 0,
        "old_source_sha256": hashes_before,
        "all_frozen_geometry_only": all(sample["geometry_only"] for sample in samples),
        "reference_id": record["reference_id"],
        "samples_json_sha256": digest(ROOT / "samples.json"),
        "samples_csv_sha256": digest(ROOT / "samples.csv"),
        "script_sha256": digest(Path(__file__)),
        "quantum_calculations_started": False,
    }
    (ROOT / "audit").mkdir(exist_ok=True)
    (ROOT / "audit/sample_freeze_audit.json").write_bytes(json_bytes(audit))
    print(json.dumps({key: audit[key] for key in ("status", "sample_count", "role_counts",
                                                "unique_frozen_geometry_hashes", "old_sources_changed")},
                     ensure_ascii=False))


if __name__ == "__main__":
    main()
