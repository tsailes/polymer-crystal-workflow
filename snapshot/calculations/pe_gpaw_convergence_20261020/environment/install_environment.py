"""Create an isolated conda-forge GPAW environment; retain all installation evidence."""
from __future__ import annotations

import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tarfile
import urllib.request

HERE = Path(__file__).resolve().parent
PREFIX = Path("/opt/pe_dft/envs/pe-gpaw-20261020")
MAMBA_ROOT = Path("/opt/pe_dft/micromamba-gpaw-20261020")
DOWNLOADS = HERE / "downloads"
BOOTSTRAP = HERE / "bootstrap"
URL = "https://conda.anaconda.org/conda-forge/linux-64/micromamba-2.9.0-0.tar.bz2"
EXPECTED = "8761c382127e6363bd9e0a2451aa3ef90d071a79133f736e2f759a3bf13040dd"
SPECS = [
    "python=3.12",
    "gpaw=26.7.0=*_mpi_openmpi_omp_*",
    "libvdwxc=0.5.0=mpi_openmpi_*",
    "libxc=7.1.2=cpu_*",
    "gpaw-data=1.0.1",
    "dftd3-python=1.6.0",
    "dftd4-python=4.2.0",
    "ase>=3.29",
    "libblas=*=*openblas",
    "pip",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    if PREFIX.exists():
        raise SystemExit(f"Refusing to replace an existing environment: {PREFIX}")
    DOWNLOADS.mkdir(parents=True, exist_ok=True)
    BOOTSTRAP.mkdir(parents=True, exist_ok=True)
    archive = DOWNLOADS / URL.rsplit("/", 1)[1]
    if not archive.exists():
        urllib.request.urlretrieve(URL, archive)
    actual = sha256(archive)
    if actual != EXPECTED:
        raise RuntimeError(f"Bootstrap checksum mismatch: {actual}")
    with tarfile.open(archive) as source:
        member = source.getmember("bin/micromamba")
        if not member.isfile():
            raise RuntimeError("Unexpected micromamba archive entry")
        target = BOOTSTRAP / "micromamba"
        with source.extractfile(member) as incoming, target.open("wb") as outgoing:
            outgoing.write(incoming.read())
    target.chmod(0o755)
    args = [
        str(target), "create", "--yes", "--no-rc", "--override-channels",
        "--strict-channel-priority", "-c", "conda-forge",
        "-p", str(PREFIX), *SPECS,
    ]
    env = os.environ.copy()
    env.update(
        MAMBA_ROOT_PREFIX=str(MAMBA_ROOT),
        MAMBA_DOWNLOAD_THREADS="4",
        OMP_NUM_THREADS="4",
        OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
    )
    metadata = {
        "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "bootstrap_url": URL,
        "bootstrap_sha256": actual,
        "bootstrap_binary_sha256": sha256(target),
        "prefix": str(PREFIX),
        "mamba_root": str(MAMBA_ROOT),
        "command_argv": args,
        "environment_overrides": {
            key: env[key] for key in (
                "MAMBA_ROOT_PREFIX", "MAMBA_DOWNLOAD_THREADS",
                "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            )
        },
        "requested_packages": SPECS,
        "install_script_sha256": sha256(Path(__file__)),
    }
    (HERE / "installation_request.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print("Bootstrap verified; starting isolated conda-forge solve.", flush=True)
    with (HERE / "installation.log").open("w", encoding="utf-8") as log:
        proc = subprocess.Popen(
            args, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1,
        )
        for line in proc.stdout:
            log.write(line)
            log.flush()
            print(line, end="", flush=True)
        code = proc.wait()
    metadata.update(
        exit_code=code,
        finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    )
    (HERE / "installation_result.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    if code:
        raise SystemExit(code)
    records = [
        json.loads(path.read_text())
        for path in sorted((PREFIX / "conda-meta").glob("*.json"))
    ]
    (HERE / "resolved_conda_records.json").write_text(
        json.dumps(records, indent=2) + "\n", encoding="utf-8")
    for name, options in (
        ("conda_explicit_lock.txt", ["list", "--explicit"]),
        ("conda_packages.json", ["list", "--json"]),
    ):
        result = subprocess.run(
            [str(target), *options, "-p", str(PREFIX)],
            env=env, capture_output=True, text=True, check=True,
        )
        (HERE / name).write_text(result.stdout, encoding="utf-8")
    print(f"Installation completed: {PREFIX}", flush=True)


if __name__ == "__main__":
    main()
