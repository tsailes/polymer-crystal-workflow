"""Record installed interfaces and provenance without evaluating any structure."""
from __future__ import annotations

import datetime
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tarfile
import warnings

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
HERE = Path(__file__).resolve().parent
PREFIX = Path("/opt/pe_dft/envs/pe-gpaw-20261020")
SOURCE = Path("/opt/pe_dft/src/gpaw-26.7.0-pe-20261020/gpaw-26.7.0")


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path):
    path = Path(path)
    return {"path": str(path), "resolved_path": str(path.resolve()),
            "size": path.stat().st_size, "sha256": sha(path)}


def command(args):
    result = subprocess.run(args, capture_output=True, text=True, timeout=40)
    return {"argv": args, "exit_code": result.returncode,
            "stdout": result.stdout, "stderr": result.stderr}


def main():
    import _gpaw
    import gpaw
    from gpaw import GPAW, PW, cgpaw
    from gpaw.xc import XC
    from gpaw.xc.libvdwxc import libvdwxc_has_stress
    from dftd3.ase import DFTD3
    from dftd4.ase import DFTD4
    from ase.calculators.mixing import SumCalculator

    info = {
        "status": "interfaces_verified_quantum_validation_pending",
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "quantum_calculations_performed_by_this_script": 0,
        "python": sys.executable,
        "platform": platform.platform(),
        "gpaw_source_path": gpaw.__file__,
        "gpaw_installation": "Official 26.7.0 sdist, unmodified scientific sources, "
                             "local build with explicit Libxc/BLAS/MPI/VDWXC, editable install",
        "libxc_version": cgpaw.libxc_version,
        "libvdwxc_stress": libvdwxc_has_stress(),
        "versions": {name: importlib.metadata.version(name) for name in [
            "gpaw", "ase", "numpy", "scipy", "gpaw-data", "dftd3", "dftd4",
            "pybind11", "tomli", "spglib"]},
        "incremental_installs": [
            json.loads(p.read_text()) for p in sorted(HERE.glob("*_incremental_install.json"))
        ],
        "official_urls": [
            "https://gpaw.readthedocs.io/install.html",
            "https://gpaw.readthedocs.io/releasenotes.html#mpichanges",
            "https://gpaw.readthedocs.io/documentation/xc/libvdwxc.html",
            "https://gitlab.com/gpaw/gpaw/-/tree/26.7.0",
            "https://pypi.org/project/gpaw/26.7.0/",
            "https://gitlab.com/libvdwxc/libvdwxc/-/tree/0.5.0",
        ],
        "interfaces": [],
        "important_files": [],
        "setup_files": [],
    }
    for xc, dispersion in [
        ("PBE", DFTD3(method="pbe", damping="d3bj", params_tweaks={"s9": 1.0})),
        ({"name": "optB88-vdW", "backend": "libvdwxc"}, None),
        ("MGGA_X_R2SCAN+MGGA_C_R2SCAN",
         DFTD4(method="r2scan", model="d4", params_tweaks={"s9": 1.0})),
    ]:
        with warnings.catch_warnings(record=True) as caught:
            obj = XC(xc)
            calc = GPAW(mode=PW(300), xc=xc, kpts=(1, 1, 1),
                        setups="PBE", legacy_gpaw=True, txt=None)
            total = SumCalculator([calc, dispersion]) if dispersion else calc
        entry = {
            "xc_input": xc,
            "xc_class": type(obj).__module__ + "." + type(obj).__name__,
            "xc_name": obj.name,
            "xc_type": obj.type,
            "xc_setup_name": obj.get_setup_name(),
            "gpaw_class": type(calc).__module__ + "." + type(calc).__name__,
            "calculation_was_evaluated": False,
            "warnings": [str(w.message) for w in caught],
            "gpaw_properties": list(calc.implemented_properties),
            "combined_properties": list(total.implemented_properties),
        }
        if dispersion:
            entry["dispersion"] = {
                "class": type(dispersion).__module__ + "." + type(dispersion).__name__,
                "parameters": dict(dispersion.parameters),
                "properties": list(dispersion.implemented_properties),
            }
        if obj.type == "MGGA":
            entry["libxc_disable_fhc_requested"] = obj.kernel.disable_fhc
            entry["libxc_disable_fhc_supported"] = bool(obj.kernel.xc.disable_fhc())
        info["interfaces"].append(entry)
    important = [
        PREFIX / "bin/python",
        PREFIX / "bin/mpiexec",
        Path(_gpaw.__file__),
        PREFIX / "lib/libxc.so",
        PREFIX / "lib/libvdwxc.so",
        PREFIX / "lib/libmpi.so",
        PREFIX / "lib/libblas.so",
        PREFIX / "lib/libdftd4.so",
        PREFIX / "lib/libs-dftd3.so",
    ]
    info["important_files"] = [file_record(p) for p in important if p.exists()]
    for root in gpaw.setup_paths:
        for element in ["C", "H"]:
            path = Path(root) / f"{element}.PBE.gz"
            if path.exists():
                info["setup_files"].append(file_record(path))
    info["runtime_info"] = command([str(PREFIX / "bin/gpaw"), "info"])
    info["mpi_version"] = command([str(PREFIX / "bin/mpiexec"), "--version"])
    info["extension_ldd"] = command(["ldd", _gpaw.__file__])
    info["cpu"] = command(["lscpu"])
    info["memory"] = command(["free", "-m"])
    (HERE / "gpaw_info_final.txt").write_text(info["runtime_info"]["stdout"] +
                                            info["runtime_info"]["stderr"])
    (HERE / "extension_ldd_final.txt").write_text(info["extension_ldd"]["stdout"])
    checked, mismatched = 0, []
    with tarfile.open(HERE / "downloads/gpaw-26.7.0.tar.gz") as archive:
        for member in archive:
            if not member.isfile():
                continue
            relative = Path(member.name).relative_to("gpaw-26.7.0")
            if relative.suffix not in {".py", ".c", ".h", ".cpp", ".hpp"}:
                continue
            local = SOURCE / relative
            expected = hashlib.sha256(archive.extractfile(member).read()).hexdigest()
            checked += 1
            if not local.exists() or sha(local) != expected:
                mismatched.append(str(relative))
    info["official_source_comparison"] = {
        "checked_source_files": checked, "mismatches": mismatched,
        "source_archive_sha256": sha(HERE / "downloads/gpaw-26.7.0.tar.gz"),
    }
    records = json.loads((HERE / "resolved_conda_records.json").read_text())
    cached = {p.name: p for p in Path("/opt/pe_dft/micromamba-gpaw-20261020/pkgs").rglob("*.conda")}
    cached.update({p.name: p for p in Path("/opt/pe_dft/micromamba-gpaw-20261020/pkgs").rglob("*.tar.bz2")})
    packages = []
    for rec in records:
        path = cached.get(rec["fn"])
        measured = sha(path) if path else None
        packages.append({
            "name": rec["name"], "version": rec["version"],
            "build": rec["build"], "url": rec["url"],
            "expected_sha256": rec.get("sha256"),
            "cached_archive": str(path) if path else None,
            "measured_sha256": measured,
            "sha256_matches": measured == rec.get("sha256") if path else None,
            "superseded": rec["name"] == "gpaw",
        })
    (HERE / "downloaded_package_hashes.json").write_text(
        json.dumps(packages, indent=2) + "\n")
    info["conda_package_archive_count"] = len(packages)
    info["conda_package_archive_hashes_verified"] = sum(p["sha256_matches"] is True for p in packages)
    info["pip_freeze"] = command([sys.executable, "-m", "pip", "freeze", "--all"])
    (HERE / "pip_freeze_final.txt").write_text(info["pip_freeze"]["stdout"])
    (HERE / "environment_verified.json").write_text(
        json.dumps(info, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({
        "status": info["status"], "versions": info["versions"],
        "libxc": info["libxc_version"], "libvdwxc_stress": info["libvdwxc_stress"],
        "official_source": info["official_source_comparison"],
        "archive_hashes_verified": info["conda_package_archive_hashes_verified"],
        "interfaces": info["interfaces"],
    }, indent=2))


if __name__ == "__main__":
    main()
