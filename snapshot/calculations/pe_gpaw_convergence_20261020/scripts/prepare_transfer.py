"""Prepare explicitly provisional PE precision-transfer queues; no SCF is run."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMMON = {
    "stage": "functional_transfer",
    "scan_parameters": ["ecut_eV", "kmesh", "gpts"],
    "mpi_ranks": 4,
    "point_group": True,
    "checkpoint": True,
    "xc_stencil": 4,
    "scf_energy_eV_electron": 1e-9,
    "scf_density_electron_electron": 1e-9,
    "scf_eigenstates_eV2_electron": 1e-11,
}
LEVELS = {
    "L": {"ecut_eV": 1400, "kmesh": [3, 6, 9], "gpts": [96, 60, 36]},
    "M": {"ecut_eV": 1600, "kmesh": [4, 8, 12], "gpts": [112, 72, 44]},
    "H": {"ecut_eV": 1800, "kmesh": [4, 8, 12], "gpts": [112, 72, 44]},
}
METHODS = {
    "pbe_d3bj": "PBE-D3(BJ)",
    "optb88_vdw": "optB88-vdW",
    "r2scan_d4": "r2SCAN-D4",
}


def write_new(name, queue):
    path = ROOT / name
    with path.open("x", encoding="utf-8") as handle:
        json.dump(queue, handle, indent=2)
        handle.write("\n")
    print(path)


def main():
    for geometry, filename in (("tashiro", "POSCAR"), ("gfn2", "POSCAR_GFN2")):
        points = []
        for slug, method in METHODS.items():
            for level, precision in LEVELS.items():
                points.append({
                    "id": f"transfer_{geometry}_{slug}_{level}",
                    "method": method,
                    "structure": str(ROOT / "reference" / filename),
                    "scan_label": level,
                    **precision,
                })
        write_new(f"transfer_{geometry}_queue.json",
                  {"stop_on_failure": True, "defaults": COMMON, "points": points})
    write_new("parallel4_validation_queue.json", {
        "stop_on_failure": True,
        "defaults": {**COMMON, "stage": "parallel_validation",
                     "scan_parameters": ["mpi_ranks"]},
        "points": [{"id": "parallel4_pbe_E1400_k3x6x9_g96",
                    "method": "PBE", **LEVELS["L"]}],
    })


if __name__ == "__main__":
    main()
