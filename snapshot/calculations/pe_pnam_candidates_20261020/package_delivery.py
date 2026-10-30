"""Package the completed, verified batch with a manifest and ZIP integrity check."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import zipfile

ROOT = Path(__file__).resolve().parent


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    summary = json.loads((ROOT / "optimization_summary.json").read_text())
    dense = json.loads((ROOT / "representative_convergence/summary.json").read_text())
    xrd = json.loads((ROOT / "xrd_benchmark_final/summary.json").read_text())
    audit = json.loads((ROOT / "audits/final_delivery_audit.json").read_text())
    if (not summary["all_preoptimization_attempts_completed"] or not dense["passed"]
            or audit.get("passed") is not True
            or xrd["candidate_count"] != len(summary["groups"])):
        raise RuntimeError("Complete all numerical and delivery checks before packaging")
    if any(ROOT.rglob("*.lock")):
        raise RuntimeError("An active lock remains; inspect it before packaging")
    support = ROOT / "supporting_code"
    support.mkdir(exist_ok=True)
    backend = ROOT.parent / "pe_gfn2_first_run_20261020/convergence_relaxation/workflow.py"
    shutil.copy2(backend, support / "validated_gfn2_workflow.py")
    archive = ROOT / "PE_Pnam_candidates_20261020.zip"
    manifest = ROOT / "delivery_manifest.json"
    verification = ROOT / "archive_verification.json"
    checksum = ROOT / "PE_Pnam_candidates_20261020.zip.sha256"
    excluded_files = {archive, manifest, verification, checksum}
    groups = {g["cluster_id"] for g in summary["groups"]}
    files = []
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or path in excluded_files:
            continue
        rel = path.relative_to(ROOT)
        if ("__pycache__" in rel.parts or "representatives_preview" in rel.parts
                or path.suffix in {".pyc", ".lock"}):
            continue
        if rel.parts[0] == "representatives" and len(rel.parts) > 1 and rel.parts[1] not in groups:
            continue
        files.append(dict(path=rel.as_posix(), bytes=path.stat().st_size, sha256=digest(path)))
    payload = dict(
        created_utc=datetime.now(timezone.utc).isoformat(),
        file_count=len(files),
        files=files,
        scope="Current final groups, all enumerated and optimized candidate evidence, protocols, source audits and code",
        exclusions=["Provisional representatives", "Python bytecode", "Archive and its self-referential checksum"],
        source_pdf="Original source PDFs are referenced by local path and SHA-256; extracted observations and audited pages are included",
        manifest_self_hash="Omitted by design; ZIP hash protects this manifest",
    )
    manifest.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for record in files:
            z.write(ROOT / record["path"], record["path"])
        z.write(manifest, manifest.name)
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None:
            raise RuntimeError("ZIP CRC validation failed")
        for record in files:
            if hashlib.sha256(z.read(record["path"])).hexdigest() != record["sha256"]:
                raise RuntimeError(f"Archive hash mismatch: {record['path']}")
    archive_hash = digest(archive)
    checksum.write_text(f"{archive_hash}  {archive.name}\n", encoding="ascii")
    result = dict(passed=True, archive=archive.name, bytes=archive.stat().st_size,
                  sha256=archive_hash, manifest_entries=len(files),
                  archived_entries=len(files) + 1,
                  validation="CRC and SHA-256 of every manifested archived file",
                  completed_utc=datetime.now(timezone.utc).isoformat())
    verification.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
