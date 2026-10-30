"""Prepare nonzero PE probe and two-step derivative checks."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parents[1]
METHODS = [("pbe_d3bj", "PBE-D3(BJ)"),
           ("optb88_vdw", "optB88-vdW"), ("r2scan_d4", "r2SCAN-D4")]


def dump(path, data):
    Path(path).write_text(json.dumps(data, indent=2) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--endpoints", action="store_true")
    ap.add_argument("--prefix", default="deriv_smooth")
    ap.add_argument("--reuse-probe", action="store_true")
    ap.add_argument("--gpts", nargs=3, type=int)
    ap.add_argument("--ecut", type=float, default=600)
    ap.add_argument("--xc-stencil", type=int)
    args = ap.parse_args()
    defaults = {
        "stage": "representative_derivative_validation",
        "structure": str(ROOT / "reference" / "derivative_probe.extxyz"),
        "ecut_eV": args.ecut, "kmesh": [1, 2, 3], "mpi_ranks": 4,
        "point_group": False, "checkpoint": False,
        "scf_energy_eV_electron": 1e-10,
        "scf_density_electron_electron": 1e-10,
        "scf_eigenstates_eV2_electron": 1e-12
    }
    if args.gpts:
        defaults["gpts"] = args.gpts
    if args.xc_stencil:
        defaults["xc_stencil"] = args.xc_stencil
    if not args.endpoints:
        a = read(ROOT / "reference" / "POSCAR")
        strain = np.array([[0, .004, .003], [.004, 0, -.002], [.003, -.002, 0]])
        a.set_cell(a.cell.array @ (np.eye(3) + strain).T, scale_atoms=True)
        changes = {0: [.010, .007, .005], 4: [-.006, .009, -.004],
                   8: [.008, -.004, .010]}
        for i, delta in changes.items():
            a.positions[i] += delta
        path = ROOT / "reference" / "derivative_probe.extxyz"
        if path.exists() and not args.reuse_probe:
            raise FileExistsError("Probe is immutable after first preparation.")
        if not args.reuse_probe:
            write(path, a)
            dump(ROOT / "reference" / "derivative_probe_definition.json", {
            "source": "POSCAR", "strain_matrix": strain.tolist(),
            "additional_cartesian_displacements_A": changes,
            "purpose": "Break mirror and orthorhombic symmetries to test nonzero force and shear stress.",
            "not_a_candidate_for_physical_ranking": True,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
        points = [{"id": f"{args.prefix}_{tag}_center_E{args.ecut:g}", "method": method}
                  for tag, method in METHODS]
        dump(ROOT / f"{args.prefix}_centers_queue.json",
             {"defaults": defaults, "points": points, "stop_on_failure": True})
        return
    points = []
    for tag, method in METHODS:
        p = ROOT / "runs" / f"{args.prefix}_{tag}_center_E{args.ecut:g}" / "result.json"
        r = json.loads(p.read_text())
        assert r["status"] == "passed_single_point"
        common = {"method": method, "gpts": r["basis"]["coarse_grid"]["shape"],
                  "center_id": p.parent.name}
        for atom, axis in [(0, 0), (4, 1), (8, 2)]:
            for step in [.001, .002]:
                for sign in [-1, 1]:
                    points.append({
                        **common,
                        "id": f"{args.prefix}_{tag}_F{atom}_{axis}_h{step:g}_{'p' if sign > 0 else 'm'}",
                        "derivative_kind": "force", "component": [atom, axis],
                        "step": step, "sign": sign,
                        "displacement": {"atom": atom, "axis": axis,
                                         "delta_A": sign * step}})
        for voigt, (i, j) in enumerate([(0, 0), (1, 1), (2, 2),
                                        (1, 2), (0, 2), (0, 1)]):
            for step in [1e-5, 2e-5]:
                for sign in [-1, 1]:
                    eps = np.zeros((3, 3))
                    eps[i, j] = sign * step if i == j else sign * step / 2
                    eps[j, i] = eps[i, j]
                    points.append({
                        **common,
                        "id": f"{args.prefix}_{tag}_S{voigt}_h{step:g}_{'p' if sign > 0 else 'm'}",
                        "derivative_kind": "stress", "component": voigt,
                        "step": step, "sign": sign, "strain": eps.tolist()})
    dump(ROOT / f"{args.prefix}_endpoints_queue.json",
         {"defaults": defaults, "points": points, "stop_on_failure": True})
    for tag, method in METHODS:
        dump(ROOT / f"{args.prefix}_endpoints_{tag}_queue.json",
             {"defaults": defaults,
              "points": [p for p in points if p["method"] == method],
              "stop_on_failure": True})


if __name__ == "__main__":
    main()
