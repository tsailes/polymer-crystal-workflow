"""Find a shell-safe technical probe cutoff; no SCF or energy evaluation.

Exact finite FFT grids and the installed GPAW PWDescriptor mask are used.
Forbidden intervals for a vector are [min_state(E_G/factor),
max_state(E_G/factor)]; inside one its membership changes across states.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np
from ase.calculators.calculator import kpts2ndarray
from ase.io import read
from ase.units import Bohr, Hartree

ROOT = Path(__file__).resolve().parents[1]
GPTS = np.array([64, 40, 24])
KMESH = [1, 2, 3]
COMPONENTS = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def integer_grid(shape, real):
    shape = np.array(shape, int)
    used = shape.copy()
    if real:
        used[2] = shape[2] // 2 + 1
    indices = np.indices(used).transpose(1, 2, 3, 0)
    if real:
        indices[..., :2] = (indices[..., :2] + shape[:2] // 2) % shape[:2] - shape[:2] // 2
    else:
        indices = (indices + shape // 2) % shape - shape // 2
    indices = indices.reshape(-1, 3)
    Q = np.arange(len(indices), dtype=np.int32)
    if real:
        mask = ((indices[:, 2] > 0) | (indices[:, 1] > 0)
                | ((indices[:, 0] >= 0) & (indices[:, 1] == 0)))
        return indices[mask], Q[mask]
    return indices, Q


def states(cell, max_step):
    values = [{"name": "center", "cell_A": cell, "strain": np.zeros((3, 3))}]
    for step in [max_step / 2, max_step]:
        for component, (a, b) in enumerate(COMPONENTS):
            for sign in (-1, 1):
                eps = np.zeros((3, 3))
                eps[a, b] = sign * step if a == b else sign * step / 2
                eps[b, a] = eps[a, b]
                values.append({"name": f"S{component}_h{step:g}_{'p' if sign > 0 else 'm'}",
                               "cell_A": cell @ (np.eye(3) + eps).T, "strain": eps})
    return values


def modes(kpoints):
    complex_i, complex_Q = integer_grid(GPTS, False)
    density_i, density_Q = integer_grid(GPTS, True)
    fine_i, fine_Q = integer_grid(2 * GPTS, True)
    out = [
        {"name": f"wfs_q{q}", "factor": 1, "indices": complex_i,
         "Q": complex_Q, "k": np.array(k), "grid": GPTS}
        for q, k in enumerate(kpoints)
    ]
    out.extend([
        {"name": "density_pd2", "factor": 2, "indices": density_i,
         "Q": density_Q, "k": np.zeros(3), "grid": GPTS},
        {"name": "density_pd3", "factor": 8, "indices": fine_i,
         "Q": fine_Q, "k": np.zeros(3), "grid": 2 * GPTS},
    ])
    return out


def boundary_eV(mode, cell):
    reciprocal = 2 * np.pi * np.linalg.inv(cell).T
    Gplusk = (mode["indices"] + mode["k"]) @ reciprocal
    return 0.5 * Hartree * Bohr**2 * np.einsum("ij,ij->i", Gplusk, Gplusk) / mode["factor"]


def safe_intervals(mode_list, geometries, lower, upper):
    forbidden = []
    stats = []
    for mode in mode_list:
        minimum = np.full(len(mode["indices"]), np.inf)
        maximum = np.full(len(mode["indices"]), -np.inf)
        for state in geometries:
            boundary = boundary_eV(mode, state["cell_A"])
            minimum = np.minimum(minimum, boundary)
            maximum = np.maximum(maximum, boundary)
        mask = (maximum >= lower) & (minimum <= upper) & (maximum > minimum)
        intervals = np.column_stack([minimum[mask], maximum[mask]])
        forbidden.extend(intervals.tolist())
        stats.append({"mode": mode["name"], "vectors_examined": len(minimum),
                      "crossing_intervals_in_search_range": len(intervals)})
    merged = []
    for a, b in sorted(forbidden):
        a, b = max(lower, a), min(upper, b)
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    safe, cursor = [], lower
    for a, b in merged:
        if a > cursor:
            safe.append([cursor, a])
        cursor = max(cursor, b)
    if cursor < upper:
        safe.append([cursor, upper])
    return safe, merged, stats


def pick_candidates(intervals, target=600, count=3):
    candidates = []
    for low, high in intervals:
        if high - low < 2e-5:
            continue
        margin = min(0.002, (high - low) / 4)
        nearest = min(max(target, low + margin), high - margin)
        rounded = round(nearest, 3)
        if not low + margin / 2 <= rounded <= high - margin / 2:
            rounded = nearest
        candidates.append({
            "ecut_eV": rounded, "safe_interval_eV": [low, high],
            "distance_from_600_eV": abs(rounded - target),
            "minimum_interval_margin_eV": min(rounded - low, high - rounded),
        })
    return sorted(candidates, key=lambda x: (
        x["distance_from_600_eV"], -x["minimum_interval_margin_eV"]))[:count]


def numpy_membership(mode_list, geometries, ecut):
    reference = None
    rows = []
    for state in geometries:
        current = {}
        for mode in mode_list:
            selected = mode["Q"][boundary_eV(mode, state["cell_A"]) <= ecut]
            current[mode["name"]] = {
                "nG": len(selected), "integer_fft_Q_sha256": hashlib.sha256(
                    selected.astype(np.int64).tobytes()).hexdigest()}
        if reference is None:
            reference = current
        rows.append({"state": state["name"], "same_as_center": current == reference,
                     "modes": current})
    return rows


def actual_descriptors(mode_list, kpoints, geometries, ecut):
    """Construct the official descriptor directly, using FFTW_ESTIMATE; no SCF."""
    from gpaw import fftw
    from gpaw.mpi import serial_comm
    from gpaw.old.grid_descriptor import GridDescriptor
    from gpaw.old.kpt_descriptor import KPointDescriptor
    from gpaw.old.pw.descriptor import PWDescriptor
    kd = KPointDescriptor(kpoints)
    reference = None
    rows = []
    for state in geometries:
        cell = state["cell_A"]
        gd = GridDescriptor(GPTS, cell / Bohr, comm=serial_comm)
        finegd = GridDescriptor(2 * GPTS, cell / Bohr, comm=serial_comm)
        wfs = PWDescriptor(ecut / Hartree, gd, dtype=complex, kd=kd, fftwflags=fftw.ESTIMATE)
        pd2 = PWDescriptor(2 * ecut / Hartree, gd, dtype=float, fftwflags=fftw.ESTIMATE)
        pd3 = PWDescriptor(8 * ecut / Hartree, finegd, dtype=float, fftwflags=fftw.ESTIMATE)
        actual_arrays = [*wfs.Q_qG, pd2.Q_qG[0], pd3.Q_qG[0]]
        current, exact_match = {}, True
        for mode, actual in zip(mode_list, actual_arrays):
            prediction = mode["Q"][boundary_eV(mode, cell) <= ecut]
            exact_match &= np.array_equal(actual, prediction)
            current[mode["name"]] = {
                "nG": len(actual), "integer_fft_Q_sha256": hashlib.sha256(
                    np.asarray(actual, dtype=np.int64).tobytes()).hexdigest()}
        if reference is None:
            reference = current
        rows.append({"state": state["name"], "matches_numpy_integer_sets": bool(exact_match),
                     "same_integer_sets_as_center": current == reference, "modes": current})
        del wfs, pd2, pd3, actual_arrays
    return {
        "status": "passed" if all(x["matches_numpy_integer_sets"]
                                 and x["same_integer_sets_as_center"] for x in rows) else "failed",
        "SCF_evaluations": 0, "energy_evaluations": 0, "states": rows,
    }


def source_evidence():
    import gpaw
    root = Path(gpaw.__file__).parent
    result = []
    for name, ranges in [
        ("old/pw/descriptor.py", [(27, 31), (43, 56), (63, 74), (119, 136)]),
        ("old/pw/density.py", [(19, 29)]),
        ("old/calculator.py", [(1253, 1263)]),
        ("old/domain.py", [(40, 48)]),
    ]:
        path = root / name
        lines = path.read_text().splitlines()
        result.append({
            "path": str(path), "sha256": sha(path),
            "official_url": "https://gitlab.com/gpaw/gpaw/-/blob/26.7.0/gpaw/" + name,
            "excerpts": [{"lines": [a, b], "text": "\n".join(
                f"{i+1}: {lines[i]}" for i in range(a - 1, min(b, len(lines))))}
                for a, b in ranges],
        })
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lower", type=float, default=590)
    parser.add_argument("--upper", type=float, default=650)
    parser.add_argument("--max-strain", type=float, default=2e-5)
    parser.add_argument("--verify-count", type=int, default=3, choices=[1, 2, 3])
    parser.add_argument("--output", type=Path,
                        default=ROOT / "validation_audit/basis_safe_cutoff.json")
    args = parser.parse_args()
    start = time.monotonic()
    probe = ROOT / "reference/derivative_probe.extxyz"
    atoms = read(probe)
    kpoints = kpts2ndarray({"size": KMESH, "gamma": True}, atoms=atoms)
    geometries = states(atoms.cell.array, args.max_strain)
    mode_list = modes(kpoints)
    grid_caps = [
        0.5 * Hartree * np.pi**2 / (
            (s["cell_A"] / Bohr / GPTS[:, None])**2).sum(1).max()
        for s in geometries]
    assert min(grid_caps) > 2 * args.upper, "Density grid cutoff cap is active; model invalid"
    safe, forbidden, counts = safe_intervals(mode_list, geometries, args.lower, args.upper)
    candidates = pick_candidates(safe, count=args.verify_count)
    print(f"Found {len(safe)} safe intervals; candidates "
          f"{[c['ecut_eV'] for c in candidates]}", flush=True)
    for candidate in candidates:
        candidate["descriptor_verification"] = actual_descriptors(
            mode_list, kpoints, geometries, candidate["ecut_eV"])
        print(f"Actual descriptor E={candidate['ecut_eV']}: "
              f"{candidate['descriptor_verification']['status']}", flush=True)
    recommended = next((x["ecut_eV"] for x in candidates
                        if x["descriptor_verification"]["status"] == "passed"), None)
    baseline = numpy_membership(mode_list, geometries, 600.0)
    result = {
        "status": "safe_cutoff_found" if recommended is not None else "no_verified_safe_cutoff",
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "purpose": "Technical finite-difference probe only; not a convergence-selected production cutoff",
        "probe_path": str(probe), "probe_sha256": sha(probe),
        "script_sha256": sha(__file__),
        "gpts": GPTS.tolist(), "fine_gpts": (2 * GPTS).tolist(),
        "gamma_centered_kmesh": KMESH, "full_BZ_kpoints": kpoints.tolist(),
        "point_group": False,
        "strain_convention": "engineering shear gamma uses symmetric eps_ij=eps_ji=gamma/2; cell'=cell@(I+eps).T",
        "strain_steps": [args.max_strain / 2, args.max_strain],
        "states": [{"name": s["name"], "cell_A": s["cell_A"].tolist(),
                    "strain": s["strain"].tolist()} for s in geometries],
        "formula": "G+k=(n+k_fractional) @ (2*pi*inverse(cell_A).T); "
                   "E_vector_eV=Hartree*Bohr**2/2 * |G+k|**2; included if E_vector <= factor*E.",
        "cutoff_factors": {"wfs": 1, "density_pd2": 2, "density_pd3": 8},
        "density_pd2_minimum_grid_cap_eV": min(grid_caps),
        "density_grid_cap_inactive": True,
        "search_range_eV": [args.lower, args.upper],
        "forbidden_cutoff_intervals_eV": forbidden,
        "safe_open_intervals_eV": safe,
        "search_statistics": counts,
        "baseline_E600_numpy_sets": baseline,
        "candidates": candidates, "recommended_technical_ecut_eV": recommended,
        "source_evidence": source_evidence(),
        "SCF_evaluations": 0, "energy_evaluations": 0,
        "wall_seconds": time.monotonic() - start,
        "limitations": [
            "Basis membership stability is a prerequisite, not proof that analytic stress passes finite differences.",
            "Valid only for this exact cell, grid, k mesh, strain scheme, and unscaled cutoff.",
            "No scientific source code was changed; no DFT energy or SCF calculation was run.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "recommended_ecut_eV": recommended,
                      "wall_seconds": result["wall_seconds"], "output": str(args.output)}))


if __name__ == "__main__":
    main()
