"""Rebuild official GPAW 26.7.0 with explicit Libxc and BLAS support."""
from __future__ import annotations

import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import urllib.request

HERE = Path(__file__).resolve().parent
PREFIX = Path("/opt/pe_dft/envs/pe-gpaw-20261020")
SOURCE = Path("/opt/pe_dft/src/gpaw-26.7.0-pe-20261020/gpaw-26.7.0")
PYTHON = str(PREFIX / "bin/python")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    meta = json.load(urllib.request.urlopen("https://pypi.org/pypi/pybind11/3.0.1/json"))
    (HERE / "downloads/pybind11_pypi_metadata.json").write_text(
        json.dumps(meta, indent=2) + "\n")
    wheel = next(x for x in meta["urls"] if x["filename"].endswith("py3-none-any.whl"))
    wheel_path = HERE / "downloads" / wheel["filename"]
    urllib.request.urlretrieve(wheel["url"], wheel_path)
    assert digest(wheel_path) == wheel["digests"]["sha256"]
    config = f"""# Build configuration only; official scientific sources are unmodified.
compiler = '{PREFIX}/bin/mpicxx'
mpi = True
nolibxc = False
noblas = False
fftw = True
libvdwxc = True
scalapack = True
elpa = False
gpu = False
libraries = ['xc', 'blas', 'lapack', 'fftw3', 'scalapack', 'vdwxc']
include_dirs = ['{PREFIX}/include']
library_dirs = ['{PREFIX}/lib']
runtime_library_dirs = ['{PREFIX}/lib']
extra_compile_args += ['-O3', '-fopenmp']
extra_link_args += ['-fopenmp']
"""
    config_path = HERE / "siteconfig.py"
    config_path.write_text(config)
    env = os.environ.copy()
    env.update(
        PATH=str(PREFIX / "bin") + ":" + env["PATH"],
        GPAW_CONFIG=str(config_path),
        GPAW_BUILD_JOBS="4",
        OMPI_CC="/usr/bin/gcc",
        OMPI_CXX="/usr/bin/g++",
        OMP_NUM_THREADS="4",
        OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
    )
    commands = [
        [PYTHON, "-m", "pip", "install", "--no-deps", str(wheel_path)],
        [PYTHON, "-m", "pip", "wheel", "--no-build-isolation", "--no-deps",
         "--wheel-dir", str(HERE / "wheels"), str(SOURCE), "--verbose"],
    ]
    audit = {
        "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "official_source_archive": str(HERE / "downloads/gpaw-26.7.0.tar.gz"),
        "source_sha256": digest(HERE / "downloads/gpaw-26.7.0.tar.gz"),
        "source_directory": str(SOURCE),
        "scientific_source_modifications": [],
        "configuration_sha256": digest(config_path),
        "pybind11_version": meta["info"]["version"],
        "pybind11_url": wheel["url"],
        "pybind11_sha256": digest(wheel_path),
        "commands": commands,
        "environment": {k: env[k] for k in [
            "GPAW_CONFIG", "GPAW_BUILD_JOBS", "OMPI_CC", "OMPI_CXX",
            "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]},
        "note": "Replaces only GPAW inside the newly created isolated environment; "
                "the original conda GPAW package metadata remains historical. "
                "The installed wheel and its hash identify the final GPAW build.",
    }
    audit_path = HERE / "gpaw_build_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2) + "\n")
    with (HERE / "gpaw_build.log").open("w") as log:
        def execute(args):
            log.write(json.dumps(args) + "\n")
            log.flush()
            print("Running:", " ".join(args[:5]), flush=True)
            p = subprocess.run(args, env=env, cwd=SOURCE, stdout=log,
                               stderr=subprocess.STDOUT)
            log.flush()
            if p.returncode:
                audit["exit_code"] = p.returncode
                audit_path.write_text(json.dumps(audit, indent=2) + "\n")
                raise SystemExit(p.returncode)
        for cmd in commands:
            execute(cmd)
        built = next((HERE / "wheels").glob("gpaw-26.7.0-*.whl"))
        install = [PYTHON, "-m", "pip", "install", "--no-deps", "--force-reinstall", str(built)]
        audit["commands"].append(install)
        audit["wheel_sha256"] = digest(built)
        audit["wheel_path"] = str(built)
        execute(install)
    audit["exit_code"] = 0
    audit["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    audit_path.write_text(json.dumps(audit, indent=2) + "\n")
    print("GPAW official source build and installation completed.", flush=True)


if __name__ == "__main__":
    main()
