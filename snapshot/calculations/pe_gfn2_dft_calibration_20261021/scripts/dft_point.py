"""Isolated, auditable GPAW point. Run all MPI ranks with identical input JSON."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys
import time
import traceback

import numpy as np
from ase.io import read, write
from ase.stress import voigt_6_to_full_3x3_stress
from ase.units import Hartree, GPa, Ry, Bohr
from gpaw import GPAW, PW
from gpaw.mpi import world


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def convert(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(type(obj).__name__)


def save(path, data):
    if world.rank == 0:
        Path(path).write_text(json.dumps(data, indent=2, default=convert) + "\n")


def versions():
    out = {"python": sys.version, "platform": platform.platform()}
    for name in ("gpaw", "ase", "numpy", "scipy", "dftd3", "dftd4"):
        try:
            out[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            out[name] = None
    return out


def dispersion_calculator(method):
    # DFTD3's params_tweaks replaces its default method parameter dictionary.
    if method == "PBE-D3(BJ)":
        from dftd3.ase import DFTD3
        return DFTD3(damping="d3bj",
                     params_tweaks={"method": "PBE", "atm": True},
                     realspace_cutoff={"disp2": 60 * Bohr, "disp3": 40 * Bohr,
                                       "cn": 60 * Bohr, "width2": .05 * Bohr,
                                       "width3": .05 * Bohr})
    if method == "r2SCAN-D4":
        from dftd4.ase import DFTD4
        return DFTD4(method="r2scan", params_tweaks={"atm": True},
                     realspace_cutoff={"disp2": 60 * Bohr, "disp3": 40 * Bohr,
                                       "cn": 30 * Bohr, "width2": .05 * Bohr,
                                       "width3": .05 * Bohr})
    return None


def make_calculator(spec, txt):
    method = spec["method"]
    if method in ("PBE", "PBE-D3(BJ)"):
        xc = "PBE"
    elif method == "optB88-vdW":
        xc = {"name": "optB88-vdW", "backend": "libvdwxc"}
    elif method == "r2SCAN-D4":
        xc = "MGGA_X_R2SCAN+MGGA_C_R2SCAN"
    else:
        raise ValueError(f"Unrecognized method {method}")
    if spec.get("xc_stencil") is not None:
        xc = ({"name": xc} if isinstance(xc, str) else dict(xc))
        xc["stencil"] = int(spec["xc_stencil"])
    mode_kwargs = {}
    if spec.get("reference_cell_A") is not None:
        # GPAW volume-scaled cutoff; not proof of invariant integer G vectors.
        mode_kwargs["cell"] = np.array(spec["reference_cell_A"])
    kwargs = {
        "mode": PW(float(spec["ecut_eV"]), **mode_kwargs),
        "xc": xc,
        # Setup "type" is paw. XC.get_setup_name() selects the PBE dataset
        # for all three requested XC objects; "PBE" here would be a type error.
        "setups": "paw",
        "kpts": {"size": tuple(spec["kmesh"]), "gamma": True},
        "occupations": {"name": "fermi-dirac", "width": 0.0},
        "nbands": int(spec.get("nbands", 16)),
        "convergence": {
            "energy": spec.get("scf_energy_eV_electron", 1e-8),
            "density": spec.get("scf_density_electron_electron", 1e-8),
            "eigenstates": spec.get("scf_eigenstates_eV2_electron", 1e-10)
        },
        "symmetry": {"point_group": spec.get("point_group", True),
                     "time_reversal": True},
        "parallel": {"augment_grids": True},
        "maxiter": 300,
        "txt": str(txt),
        "legacy_gpaw": True
    }
    if spec.get("gpts"):
        kwargs["gpts"] = tuple(spec["gpts"])
    elif spec.get("h_A"):
        kwargs["h"] = float(spec["h_A"])
    return GPAW(**kwargs)


def basis_metadata(calc):
    data = {}
    for name, grid in (("coarse_grid", calc.wfs.gd),
                       ("fine_grid", calc.density.finegd)):
        data[name] = {"shape": grid.N_c.tolist(),
                      "cell_A": (grid.cell_cv * 0.5291772105638411).tolist()}
    pd = getattr(calc.wfs, "pd", None)
    if pd is not None:
        data["ngmin"] = getattr(pd, "ngmin", None)
        data["ngmax"] = getattr(pd, "ngmax", None)
        # Each q array stores integer FFT indices for its actual PW basis.
        qs = getattr(pd, "Q_qG", [])
        data["local_integer_G_index_sets"] = [
            {"local_q": q, "nG": len(vals),
             "sha256": hashlib.sha256(
                 np.asarray(vals, dtype=np.int64).tobytes()).hexdigest()}
            for q, vals in enumerate(qs)]
    data["local_kpoints"] = [
        {"q": int(kpt.q), "k": int(kpt.k), "s": int(kpt.s)}
        for kpt in calc.wfs.kpt_u]
    data["bz_kpoints"] = calc.wfs.kd.bzk_kc.tolist()
    data["ibz_kpoints"] = calc.wfs.kd.ibzk_kc.tolist()
    data["xc_class"] = str(type(calc.hamiltonian.xc))
    data["xc_name"] = getattr(calc.hamiltonian.xc, "name", None)
    descriptors = {}
    for owner_name, owner in (("density", calc.density),
                              ("hamiltonian", calc.hamiltonian)):
        for attr in ("pd2", "pd3"):
            desc = getattr(owner, attr, None)
            if desc is not None:
                descriptors[owner_name + "_" + attr] = {
                    "ecut_Ha": float(desc.ecut),
                    "grid_shape": desc.gd.N_c.tolist(),
                    "Q_qG": [
                        {"q": q, "nG": len(vals), "sha256": hashlib.sha256(
                            np.asarray(vals, dtype=np.int64).tobytes()).hexdigest()}
                        for q, vals in enumerate(desc.Q_qG)]}
    data["reciprocal_descriptors"] = descriptors
    return data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input_json")
    args = ap.parse_args()
    inp = Path(args.input_json).resolve()
    work = inp.parent
    spec = json.loads(inp.read_text())
    os.chdir(work)
    started = time.time()
    result = {
        "status": "running", "input": spec,
        "input_sha256": sha256(inp), "runner_sha256": sha256(__file__),
        "software": versions(), "mpi_world_size": world.size,
        "start_unix": started, "units": {
            "energy": "eV/cell", "forces": "eV/Angstrom",
            "stress": "eV/Angstrom^3; ASE tensile-positive",
            "stress_voigt_order": ["xx", "yy", "zz", "yz", "xz", "xy"],
            "energy_force_consistent": "Free energy (zero occupation width)"
        }}
    save(work / "status.json", result)
    calc = None
    try:
        structure = Path(spec["structure"]).resolve()
        atoms = read(structure)
        atoms.pbc = True
        result["source_geometry_sha256"] = sha256(structure)
        if spec.get("displacement"):
            item = spec["displacement"]
            atoms.positions[item["atom"], item["axis"]] += item["delta_A"]
        if spec.get("strain"):
            # matrix is infinitesimal symmetric Cartesian strain. For an
            # engineering shear gamma, off-diagonals must each be gamma/2.
            eps = np.asarray(spec["strain"], float)
            atoms.set_cell(atoms.cell.array @ (np.eye(3) + eps).T,
                           scale_atoms=True)
        # ASE write performs MPI synchronization itself; call on every rank.
        write(work / "input_structure.extxyz", atoms)
        result["natoms"] = len(atoms)
        result["formula"] = atoms.get_chemical_formula()
        result["cell_A"] = atoms.cell.array
        result["positions_A"] = atoms.positions.copy()
        result["volume_A3"] = atoms.get_volume()
        calc = make_calculator(spec, work / "gpaw.txt")

        def record_iteration():
            if world.rank != 0:
                return
            h = calc.hamiltonian
            data = {"iteration": getattr(calc.scf, "niter", None),
                    "wall_seconds": time.time() - started}
            for name in ("e_total_free", "e_total_extrapolated",
                         "e_kinetic", "e_coulomb", "e_zero",
                         "e_external", "e_xc", "e_entropy"):
                val = getattr(h, name, None)
                if val is not None:
                    data[name + "_eV"] = float(val * Hartree)
            with (work / "scf_iterations.jsonl").open("a") as f:
                f.write(json.dumps(data) + "\n")
                f.flush()

        calc.attach(record_iteration, 1)
        atoms.calc = calc
        e_extrap = atoms.get_potential_energy()
        e_free = atoms.get_potential_energy(force_consistent=True)
        forces = atoms.get_forces()
        stress = atoms.get_stress()
        dft = {"energy_extrapolated_eV": e_extrap,
               "energy_free_eV": e_free, "forces_eV_A": forces.copy(),
               "stress_eV_A3": stress.copy(),
               "stress_tensor_eV_A3": voigt_6_to_full_3x3_stress(stress),
               "stress_kbar": stress / GPa * 10}
        disp = dispersion_calculator(spec["method"])
        contribution = {"energy_eV": 0.0, "forces_eV_A": np.zeros_like(forces),
                        "stress_eV_A3": np.zeros(6), "stress_kbar": np.zeros(6)}
        if disp is not None:
            atoms.calc = disp
            de = atoms.get_potential_energy()
            df = atoms.get_forces()
            ds = atoms.get_stress()
            contribution = {"energy_eV": de, "forces_eV_A": df,
                            "stress_eV_A3": ds, "stress_kbar": ds / GPa * 10,
                            "calculator_class": str(type(disp)),
                            "calculator_parameters": dict(disp.parameters)}
            e_extrap += de
            e_free += de
            forces += df
            stress += ds
        result["dft"] = dft
        result["dispersion"] = contribution
        result["total"] = {
            "energy_extrapolated_eV": e_extrap, "energy_free_eV": e_free,
            "energy_Ry": e_free / Ry, "energy_eV_atom": e_free / len(atoms),
            "energy_eV_CH2": e_free / (len(atoms) / 3),
            "forces_eV_A": forces, "stress_eV_A3": stress,
            "stress_tensor_eV_A3": voigt_6_to_full_3x3_stress(stress),
            "stress_kbar": stress / GPa * 10,
            "max_force_norm_eV_A": np.linalg.norm(forces, axis=1).max(),
            "net_force_eV_A": forces.sum(axis=0)}
        result["SCF_converged"] = bool(calc.scf.converged)
        result["SCF_iterations"] = int(calc.scf.niter)
        result["electrons"] = float(calc.get_number_of_electrons())
        result["fermi_level_eV"] = float(calc.get_fermi_level())
        result["basis"] = basis_metadata(calc)
        # Archive each MPI rank because local k-point partitions may differ.
        (work / f"basis_rank{world.rank:04d}.json").write_text(
            json.dumps(result["basis"], indent=2, default=convert) + "\n")
        result["datasets"] = []
        for s in calc.wfs.setups:
            filename = getattr(s, "filename", None)
            result["datasets"].append({
                "symbol": s.symbol, "type": s.type,
                "filename": filename,
                "setupname": getattr(getattr(s, "data", None), "setupname", None),
                "fingerprint": getattr(s, "fingerprint", None),
                "file_sha256": sha256(filename) if filename and
                Path(filename).is_file() else None})
        isfinite = all(np.isfinite(v).all() for v in
                       (e_free, e_extrap, forces, stress))
        if not isfinite or not result["SCF_converged"]:
            raise RuntimeError("Nonfinite properties or unconverged SCF.")
        if spec.get("checkpoint", True):
            calc.write(str(work / "final.gpw"))
        result["status"] = "passed_single_point"
    except Exception as exc:
        result["status"] = "failed"
        result["error"] = repr(exc)
        result["traceback"] = traceback.format_exc()
        if world.rank == 0:
            print(result["traceback"], file=sys.stderr, flush=True)
    finally:
        result["wall_seconds"] = time.time() - started
        save(work / "result.json", result)
        save(work / "status.json", {"status": result["status"],
                                  "wall_seconds": result["wall_seconds"]})
    if result["status"] == "failed":
        sys.exit(1)
    if world.rank == 0:
        print(json.dumps({"status": result["status"],
                          "energy_free_eV": result["total"]["energy_free_eV"],
                          "wall_seconds": result["wall_seconds"]}), flush=True)


if __name__ == "__main__":
    main()
