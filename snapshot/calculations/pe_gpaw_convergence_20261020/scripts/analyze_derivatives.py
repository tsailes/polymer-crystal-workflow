"""Read-only GPAW derivative analysis; writes reports, never starts calculations.

Run with the existing ASE/NumPy environment:
  python scripts/analyze_derivatives.py
  python scripts/analyze_derivatives.py --require-complete

The default series is deriv_smooth_*_center_E600.  Every selected point must
be the latest successful attempt of its logical id, with matching geometry,
method, numerical settings, PAW data, software and runner provenance.
"""
from __future__ import annotations

import argparse
import copy
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import sys

sys.dont_write_bytecode = True
import numpy as np
from ase.io import read
from ase.units import Bohr, GPa

ROOT = Path(__file__).resolve().parents[1]
METHODS = (("pbe_d3bj", "PBE-D3(BJ)"),
           ("optb88_vdw", "optB88-vdW"), ("r2scan_d4", "r2SCAN-D4"))
VOIGT = ("xx", "yy", "zz", "yz", "xz", "xy")
PAIRS = ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1))
FORCE_COMPONENTS = ((0, 0), (4, 1), (8, 2))
STEPS = {"force": (0.001, 0.002), "stress": (1e-5, 2e-5)}
ABS_TOL = {"force": 0.001, "stress": 0.1}
STEP_TOL = {"force": 0.0005, "stress": 0.05}
LEVELS = ("dft", "dispersion", "total")
KBAR = 10.0 / GPa
META_KEYS = {
    "id", "center_id", "derivative_kind", "component", "step", "sign",
    "displacement", "strain", "stage", "structure", "checkpoint",
    "mpi_ranks", "queue_sha256", "protocol_sha256",
    "runner_sha256_at_submission", "gpts", "h_A",
    "validation_role", "noise_check_parent_id", "noise_check_parent_run",
    "noise_check_group",
}
DEFAULT_SETTINGS = {
    "nbands": 16, "point_group": True,
    "scf_energy_eV_electron": 1e-8,
    "scf_density_electron_electron": 1e-8,
    "scf_eigenstates_eV2_electron": 1e-10,
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical_json(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), allow_nan=False)


def object_digest(obj):
    return hashlib.sha256(canonical_json(obj).encode()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def write_json(path, value):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2,
                               allow_nan=False) + "\n", encoding="utf-8")
    temp.replace(path)


def logical_id(name):
    return re.sub(r"__attempt\d+$", "", str(name))


def attempt_number(name):
    match = re.search(r"__attempt(\d+)$", name)
    return int(match[1]) if match else 1


def local_path(value, base):
    """Resolve the recorded Windows/WSL path on either host, without guessing."""
    value = str(value)
    if os.name == "nt" and re.match(r"^/mnt/[a-zA-Z]/", value):
        value = value[5].upper() + ":/" + value[7:]
    elif os.name != "nt" and re.match(r"^[a-zA-Z]:[/\\]", value):
        value = "/mnt/" + value[0].lower() + "/" + value[3:].replace("\\", "/")
    path = Path(value)
    return path if path.is_absolute() else base / path


def strain_matrix(component, signed_step):
    i, j = PAIRS[int(component)]
    matrix = np.zeros((3, 3))
    matrix[i, j] = signed_step if i == j else signed_step / 2
    matrix[j, i] = matrix[i, j]
    return matrix


def group_key(spec):
    component = spec["component"]
    if spec["derivative_kind"] == "force":
        component = tuple(int(x) for x in component)
    else:
        component = int(component)
    return (logical_id(spec["center_id"]), spec["derivative_kind"], component)


def point_key(spec):
    return (*group_key(spec), float(spec["step"]), int(spec["sign"]))


def point_label(kind, component):
    if kind == "force":
        return f"atom{component[0]}_{'xyz'[component[1]]}"
    return VOIGT[component]


def finite_array(value, shape, label):
    array = np.asarray(value, dtype=float)
    if array.shape != shape or not np.isfinite(array).all():
        raise ValueError(f"Invalid {label}: expected finite shape {shape}")
    return array


def position_error(first, second, cell):
    delta = (np.asarray(first) - np.asarray(second)) @ np.linalg.inv(cell)
    delta -= np.rint(delta)
    return float(np.max(np.abs(delta @ cell)))


def effective_settings(spec):
    out = dict(DEFAULT_SETTINGS)
    out.update({k: v for k, v in spec.items() if k not in META_KEYS})
    return out


def dataset_signature(record):
    fields = ("symbol", "type", "setupname", "fingerprint", "file_sha256")
    rows = [{key: data.get(key) for key in fields}
            for data in record["datasets"]]
    if any(not row["file_sha256"] for row in rows):
        raise ValueError("Dataset file SHA-256 is missing")
    return rows


def energy(record, level):
    key = "energy_eV" if level == "dispersion" else "energy_free_eV"
    return float(record[level][key])


def analytic(record, level, kind, component):
    if kind == "force":
        return float(record[level]["forces_eV_A"][component[0]][component[1]])
    return float(record[level]["stress_eV_A3"][component]) * KBAR


def discover(root, prefix):
    attempts, unreadable = {}, []
    for folder in sorted((root / "runs").glob(prefix + "_*")):
        if not folder.is_dir() or not (folder / "input.json").exists():
            continue
        try:
            spec = load(folder / "input.json")
        except (OSError, ValueError) as exc:
            unreadable.append({"folder": str(folder), "error": str(exc)})
            continue
        ident = logical_id(spec.get("id", folder.name))
        result, reason = {}, None
        try:
            if (folder / "result.json").exists():
                result = load(folder / "result.json")
        except (OSError, ValueError) as exc:
            reason = "Result still being written or unreadable: " + str(exc)
        process = {}
        try:
            if (folder / "process.json").exists():
                process = load(folder / "process.json")
        except (OSError, ValueError):
            pass
        state = "pending"
        if result.get("status") == "passed_single_point":
            # A result can appear just before the MPI child exits.  Do not
            # call an unfinished/nonzero MPI process a successful attempt.
            if process.get("exit_code") == 0:
                state = "success"
            elif process.get("exit_code") is not None:
                state = "failed"
                reason = "MPI process nonzero despite passed result"
            else:
                reason = "Passed result exists; process completion still pending"
        elif result.get("status") == "failed" or process.get("exit_code", 0):
            state = "failed"
            reason = result.get("error", "MPI process failed")
        item = {"id": ident, "folder": folder, "input": spec, "result": result,
                "state": state, "reason": reason,
                "attempt": attempt_number(folder.name),
                "start_unix": result.get("start_unix", process.get("start_unix", 0))}
        attempts.setdefault(ident, []).append(item)
    for items in attempts.values():
        items.sort(key=lambda x: (x["attempt"], x["start_unix"], x["folder"].name))
    return attempts, unreadable


def select(attempts, ident):
    items = attempts.get(logical_id(ident), [])
    successful = [item for item in items if item["state"] == "success"]
    return successful[-1] if successful else None


def cn60_view(root, run):
    """Verified overlay: retain original DFT, use the already recombined total."""
    if run["input"]["method"] != "PBE-D3(BJ)":
        return run
    path = root / "method_audit/derivative_d3_cn60" / (run["folder"].name + ".json")
    if not path.exists():
        return None
    revision = load(path)
    original = run["result"]
    if (revision.get("status") != "passed_dispersion_recombination"
            or revision.get("revision_id") != "D3_CN60_smooth005bohr"
            or revision.get("DFT_recomputed") is not False
            or revision.get("source_result_sha256") != digest(run["folder"] / "result.json")
            or revision.get("same_geometry_sha256") != digest(run["folder"] / "input_structure.extxyz")):
        raise ValueError("CN60 revision identity/source hashes are invalid")
    for key in ("input", "dft", "datasets", "basis", "positions_A", "cell_A"):
        if revision[key] != original[key]:
            raise ValueError("CN60 overlay changed original DFT/geometry/settings: " + key)
    expected_parameters = copy.deepcopy(original["dispersion"]["calculator_parameters"])
    expected_parameters["realspace_cutoff"]["cn"] = 60 * Bohr
    if revision["dispersion"]["calculator_parameters"] != expected_parameters:
        raise ValueError("CN60 changed another dispersion parameter or has the wrong CN cutoff")
    if not re.fullmatch(r"[0-9a-f]{64}", revision["revision_script_sha256"]):
        raise ValueError("Missing CN60 revision-script provenance")
    atoms = read(run["folder"] / "input_structure.extxyz")
    exact_geometry = {
        "atomic_numbers": atoms.numbers.tolist(), "positions_A": original["positions_A"],
        "cell_A": original["cell_A"], "pbc": atoms.pbc.tolist(),
    }
    if object_digest(exact_geometry) != revision["exact_geometry_sha256"]:
        raise ValueError("CN60 exact geometry hash mismatch")
    revised = copy.deepcopy(run)
    revised["result"]["dispersion"] = revision["dispersion"]
    revised["result"]["total"] = revision["total"]
    revised["revision"] = {"path": str(path), "sha256": digest(path),
                           "revision_id": revision["revision_id"],
                           "revision_script_sha256": revision["revision_script_sha256"],
                           "source_result_sha256": revision["source_result_sha256"],
                           "DFT_recomputed": False}
    return revised


def history(attempts, ident):
    return [{"folder": str(x["folder"]), "attempt": x["attempt"],
             "status": x["state"], "reason": x["reason"]}
            for x in attempts.get(logical_id(ident), [])]


def merged_basis(run):
    """Merge rank-local q through the explicit global IBZ k index."""
    record, folder = run["result"], run["folder"]
    size = int(record["mpi_world_size"])
    merged, records, common, descriptors = {}, [], None, None
    descriptor_presence = []
    for rank in range(size):
        path = folder / f"basis_rank{rank:04d}.json"
        if not path.exists():
            raise ValueError(f"Missing rank basis file: {path.name}")
        basis = load(path)
        header = {
            "coarse_shape": basis["coarse_grid"]["shape"],
            "fine_shape": basis["fine_grid"]["shape"],
            "bz_kpoints": basis["bz_kpoints"],
            "ibz_kpoints": basis["ibz_kpoints"],
        }
        if common is None:
            common = header
        elif header != common:
            raise ValueError("Rank-dependent grid or global k-point metadata")
        for grid in ("coarse_grid", "fine_grid"):
            if not np.allclose(basis[grid]["cell_A"], record["cell_A"],
                               atol=1e-8, rtol=0):
                raise ValueError(f"Rank {rank} {grid} cell differs from result")
        rank_descriptors = basis.get("reciprocal_descriptors")
        descriptor_presence.append(rank_descriptors is not None)
        if rank_descriptors is not None:
            required_names = {"density_pd2", "density_pd3",
                              "hamiltonian_pd2", "hamiltonian_pd3"}
            if set(rank_descriptors) != required_names:
                raise ValueError("Incomplete reciprocal density/Hamiltonian descriptor metadata")
            for name, descriptor in rank_descriptors.items():
                grid_name = "coarse_grid" if name.endswith("pd2") else "fine_grid"
                if descriptor["grid_shape"] != basis[grid_name]["shape"]:
                    raise ValueError("Reciprocal descriptor grid differs: " + name)
                if not np.isfinite(descriptor["ecut_Ha"]) or descriptor["ecut_Ha"] <= 0:
                    raise ValueError("Invalid reciprocal descriptor cutoff")
                qsets_d = descriptor["Q_qG"]
                if not qsets_d or len({q["q"] for q in qsets_d}) != len(qsets_d):
                    raise ValueError("Empty/duplicate reciprocal descriptor Q sets")
                for qset in qsets_d:
                    if qset["nG"] <= 0 or not re.fullmatch(r"[0-9a-f]{64}", qset["sha256"]):
                        raise ValueError("Invalid density/Hamiltonian G count/hash")
            if descriptors is None:
                descriptors = rank_descriptors
            elif rank_descriptors != descriptors:
                raise ValueError("Ranks disagree about density/Hamiltonian descriptors")
        qsets = {int(q["local_q"]): q for q in basis["local_integer_G_index_sets"]}
        for kp in basis["local_kpoints"]:
            q, k, spin = (int(kp[x]) for x in ("q", "k", "s"))
            if q not in qsets or not 0 <= k < len(header["ibz_kpoints"]):
                raise ValueError("Unresolvable local-q to global-k mapping")
            value = {"global_k": k, "spin": spin,
                     "fractional_k": header["ibz_kpoints"][k],
                     "nG": int(qsets[q]["nG"]), "sha256": qsets[q]["sha256"]}
            if value["nG"] <= 0 or not re.fullmatch(r"[0-9a-f]{64}", value["sha256"]):
                raise ValueError("Invalid PW count/hash")
            key = (spin, k)
            if key in merged and merged[key] != value:
                raise ValueError("Ranks disagree about a duplicate global k basis")
            merged[key] = value
        records.append({"rank": rank, "path": str(path), "sha256": digest(path)})
    if not merged:
        raise ValueError("Empty merged PW basis")
    spins = sorted({key[0] for key in merged})
    required = {(spin, k) for spin in spins for k in range(len(common["ibz_kpoints"]))}
    if set(merged) != required or spins != [0]:
        raise ValueError("Incomplete global k coverage or unexpected spin state")
    if any(descriptor_presence) and not all(descriptor_presence):
        raise ValueError("Only some ranks have reciprocal descriptor metadata")
    payload = {**common, "global_G_index_sets": [merged[key] for key in sorted(merged)]}
    wavefunction_sha = object_digest(payload)
    # The continuously changing cutoff value is not a change in a discrete
    # basis.  Compare actual Q sets together with their FFT grid shapes.
    if descriptors is not None:
        payload["reciprocal_G_index_sets"] = {
            name: {"grid_shape": val["grid_shape"], "Q_qG": val["Q_qG"]}
            for name, val in sorted(descriptors.items())}
    return {**payload, "merged_sha256": object_digest(payload),
            "wavefunction_sha256": wavefunction_sha,
            "reciprocal_descriptors": descriptors,
            "density_hamiltonian_basis_recorded": descriptors is not None,
            "rank_files": records,
            "interpretation": "Flattened integer FFT-index hashes; compare only together "
                              "with identical grid shapes and global k-point mapping."}


def validate_run(run):
    record, spec, folder = run["result"], run["input"], run["folder"]
    if record.get("input") != spec or record.get("input_sha256") != digest(folder / "input.json"):
        raise ValueError("Stored input or input SHA differs from actual input.json")
    snapshot = folder / "run_point_snapshot.py"
    if (not snapshot.exists() or record.get("runner_sha256") != digest(snapshot)
            or spec.get("runner_sha256_at_submission") != digest(snapshot)):
        raise ValueError("Runner snapshot/submission/result hashes disagree")
    if record.get("SCF_converged") is not True:
        raise ValueError("Native SCF convergence not confirmed")
    if record.get("natoms") != 12 or record.get("formula") != "C4H8":
        raise ValueError("Expected the 12-atom C4H8 probe")
    if record.get("units", {}).get("stress_voigt_order") != list(VOIGT):
        raise ValueError("Unexpected Voigt ordering")
    if record.get("units", {}).get("energy") != "eV/cell":
        raise ValueError("Unexpected energy unit")
    if not np.isclose(float(record["electrons"]), 24, atol=1e-10):
        raise ValueError("Unexpected valence electron count")
    source = local_path(spec["structure"], folder)
    if not source.exists() or digest(source) != record["source_geometry_sha256"]:
        raise ValueError("Original probe file missing or SHA changed")
    atoms = read(source)
    if not all(atoms.pbc) or len(atoms) != 12:
        raise ValueError("Source is not the complete periodic probe")
    if spec.get("displacement"):
        disp = spec["displacement"]
        atoms.positions[disp["atom"], disp["axis"]] += disp["delta_A"]
    if spec.get("strain"):
        eps = finite_array(spec["strain"], (3, 3), "strain")
        atoms.set_cell(atoms.cell.array @ (np.eye(3) + eps).T, scale_atoms=True)
    cell = finite_array(record["cell_A"], (3, 3), "cell")
    positions = finite_array(record["positions_A"], (12, 3), "positions")
    if not np.allclose(cell, atoms.cell.array, atol=1e-8, rtol=0):
        raise ValueError("Calculated cell does not match declared transformation")
    if position_error(positions, atoms.positions, cell) > 2e-8:
        raise ValueError("Calculated positions do not match declared transformation")
    if not np.isclose(record["volume_A3"], abs(np.linalg.det(cell)), atol=1e-8, rtol=0):
        raise ValueError("Recorded volume is inconsistent with the cell")
    saved = read(folder / "input_structure.extxyz")
    if saved.get_chemical_symbols() != atoms.get_chemical_symbols() or not all(saved.pbc):
        raise ValueError("Saved structure species/PBC differ from declared input")
    if (not np.allclose(saved.cell.array, cell, atol=1e-8, rtol=0)
            or position_error(saved.positions, positions, cell) > 2e-8):
        raise ValueError("Saved structure differs from result geometry")
    signature = dataset_signature(record)
    if [x["symbol"] for x in signature] != atoms.get_chemical_symbols():
        raise ValueError("Dataset atom ordering differs from saved structure")
    for level in LEVELS:
        if not np.isfinite(energy(record, level)):
            raise ValueError("Nonfinite energy")
        forces = finite_array(record[level]["forces_eV_A"], (12, 3), "forces")
        stress = finite_array(record[level]["stress_eV_A3"], (6,), "stress")
        kbar = finite_array(record[level]["stress_kbar"], (6,), "stress_kbar")
        if not np.allclose(kbar, stress * KBAR, atol=1e-6, rtol=1e-8):
            raise ValueError("Stress unit conversion inconsistent with ASE")
        if "stress_tensor_eV_A3" in record[level]:
            tensor = finite_array(record[level]["stress_tensor_eV_A3"], (3, 3), "stress tensor")
            vector = tensor.flat[[0, 4, 8, 5, 2, 1]]
            if not np.allclose(tensor, tensor.T, atol=1e-12, rtol=0) or not np.allclose(
                    vector, stress, atol=1e-12, rtol=0):
                raise ValueError("Full tensor/Voigt inconsistency")
    if abs(energy(record, "total") - energy(record, "dft") - energy(record, "dispersion")) > 1e-10:
        raise ValueError("Total energy is not DFT plus dispersion")
    for name in ("forces_eV_A", "stress_eV_A3"):
        difference = (np.asarray(record["total"][name]) -
                      np.asarray(record["dft"][name]) -
                      np.asarray(record["dispersion"][name]))
        if np.max(np.abs(difference)) > 1e-11:
            raise ValueError("Total force/stress is not the sum of contributions")
    return {"basis": merged_basis(run), "datasets": signature,
            "symbols": atoms.get_chemical_symbols(), "settings": effective_settings(spec)}


def compare_endpoint(center, endpoint, checks, kind, component, step, sign):
    a, b = center["result"], endpoint["result"]
    ca, cb = checks[center["folder"].name], checks[endpoint["folder"].name]
    if a["source_geometry_sha256"] != b["source_geometry_sha256"]:
        raise ValueError("Center and endpoint originate from different probe files")
    if ca["settings"] != cb["settings"]:
        differing = sorted(k for k in ca["settings"].keys() | cb["settings"].keys()
                           if ca["settings"].get(k) != cb["settings"].get(k))
        raise ValueError("Numerical/method settings differ: " + ", ".join(differing))
    if ca["datasets"] != cb["datasets"] or ca["symbols"] != cb["symbols"]:
        raise ValueError("PAW datasets or element ordering changed")
    if a["software"] != b["software"] or a["runner_sha256"] != b["runner_sha256"]:
        raise ValueError("Software or runner changed between center and endpoint")
    if a["dispersion"].get("calculator_parameters") != b["dispersion"].get("calculator_parameters"):
        raise ValueError("Dispersion parameters differ between center and endpoint")
    for key in ("xc_name", "xc_class"):
        if a["basis"][key] != b["basis"][key]:
            raise ValueError("Native XC identity differs")
    if center["input"].get("displacement") or center["input"].get("strain"):
        raise ValueError("Center must be the already prepared probe itself")
    spec = endpoint["input"]
    if logical_id(spec["center_id"]) != center["id"] or group_key(spec)[1:] != (kind, component):
        raise ValueError("Endpoint derivative metadata does not name this center/component")
    if float(spec["step"]) != step or int(spec["sign"]) != sign:
        raise ValueError("Endpoint step/sign metadata mismatch")
    if ca["settings"].get("point_group") is not False:
        raise ValueError("Probe derivative checks require point_group=False")
    if kind == "force":
        expected = {"atom": component[0], "axis": component[1], "delta_A": sign * step}
        if spec.get("displacement") != expected or spec.get("strain"):
            raise ValueError("Force endpoint is not the declared single Cartesian displacement")
        xyz = np.array(a["positions_A"])
        xyz[component[0], component[1]] += sign * step
        if (not np.allclose(a["cell_A"], b["cell_A"], atol=1e-8, rtol=0)
                or position_error(xyz, b["positions_A"], np.array(a["cell_A"])) > 2e-8):
            raise ValueError("Force geometry differs from center plus exactly one displacement")
    else:
        eps = strain_matrix(component, sign * step)
        if spec.get("displacement") or not np.allclose(spec.get("strain"), eps, atol=1e-16, rtol=0):
            raise ValueError("Stress endpoint does not use the declared engineering strain")
        expected_cell = np.array(a["cell_A"]) @ (np.eye(3) + eps).T
        frac = np.array(a["positions_A"]) @ np.linalg.inv(a["cell_A"])
        if (not np.allclose(b["cell_A"], expected_cell, atol=1e-8, rtol=0)
                or position_error(frac @ expected_cell, b["positions_A"], expected_cell) > 2e-8):
            raise ValueError("Stress geometry does not preserve center fractional coordinates")
    wave_same = ca["basis"]["wavefunction_sha256"] == cb["basis"]["wavefunction_sha256"]
    density_recorded = (ca["basis"]["density_hamiltonian_basis_recorded"] and
                        cb["basis"]["density_hamiltonian_basis_recorded"])
    density_same = (ca["basis"].get("reciprocal_G_index_sets") ==
                    cb["basis"].get("reciprocal_G_index_sets")) if density_recorded else None
    return {
        # Missing old metadata is unknown, not proof of a shell crossing.
        # The separate scope/recorded fields prevent a full-basis claim.
        "basis_unchanged": wave_same and (density_same if density_recorded else True),
        "wavefunction_basis_unchanged": wave_same,
        "density_hamiltonian_basis_recorded": density_recorded,
        "density_hamiltonian_basis_unchanged": density_same,
        "center_basis_sha256": ca["basis"]["merged_sha256"],
        "endpoint_basis_sha256": cb["basis"]["merged_sha256"],
        "coarse_grid_unchanged": ca["basis"]["coarse_shape"] == cb["basis"]["coarse_shape"],
        "fine_grid_unchanged": ca["basis"]["fine_shape"] == cb["basis"]["fine_shape"],
    }


def expected_points(root, prefix, cutoff, attempts):
    """Generate the agreed coverage and supplement it from actual queue metadata."""
    centers, specs = {}, {}
    available = []
    for ident, items in attempts.items():
        spec = items[-1]["input"]
        if (re.fullmatch(re.escape(prefix) + r"_.*_center_E[0-9.]+", ident)
                and spec.get("validation_role") != "tight_scf_noise_check"
                and not spec.get("derivative_kind")):
            available.append((ident, spec))
    observed_cutoffs = {float(spec["ecut_eV"]) for _, spec in available}
    fallback_cutoff = cutoff if cutoff is not None else (
        next(iter(observed_cutoffs)) if len(observed_cutoffs) == 1 else 600.)
    for tag, method in METHODS:
        matches = [(ident, spec) for ident, spec in available
                   if spec.get("method") == method and
                   (cutoff is None or np.isclose(spec["ecut_eV"], cutoff, atol=1e-9, rtol=0))]
        if len(matches) > 1:
            raise ValueError("Multiple centers for one method; choose --cutoff explicitly: " + method)
        center = matches[0][0] if matches else f"{prefix}_{tag}_center_E{fallback_cutoff:g}"
        centers[center] = method
        for kind, components in (("force", FORCE_COMPONENTS), ("stress", range(6))):
            for component in components:
                for step in STEPS[kind]:
                    for sign in (-1, 1):
                        suffix = (f"F{component[0]}_{component[1]}" if kind == "force"
                                  else f"S{component}")
                        ident = f"{prefix}_{tag}_{suffix}_h{step:g}_{'p' if sign > 0 else 'm'}"
                        spec = {"id": ident, "method": method, "center_id": center,
                                "derivative_kind": kind, "component": component,
                                "step": step, "sign": sign}
                        specs[point_key(spec)] = spec
    for path in sorted(root.glob(prefix + "_endpoints*queue.json")):
        try:
            queue = load(path)
            candidates = [{**queue.get("defaults", {}), **p} for p in queue["points"]]
        except (OSError, ValueError, KeyError):
            continue
        for spec in candidates:
            if logical_id(spec.get("center_id", "")) in centers:
                specs[point_key(spec)] = spec
    for items in attempts.values():
        spec = items[-1]["input"]
        if spec.get("validation_role") == "tight_scf_noise_check":
            continue
        if spec.get("derivative_kind") and logical_id(spec.get("center_id", "")) in centers:
            specs[point_key(spec)] = spec
    return centers, specs


def analyze(root, prefix, cutoff, view="original"):
    attempts, unreadable = discover(root, prefix)
    center_methods, endpoints = expected_points(root, prefix, cutoff, attempts)
    selected, checks, invalid, index, revision_pending = {}, {}, {}, {}, []
    for key, spec in endpoints.items():
        index[key] = logical_id(spec["id"])
    wanted = set(center_methods) | set(index.values())
    for ident in wanted:
        run = select(attempts, ident)
        if run is not None:
            if view == "d3-cn60":
                try:
                    revised = cn60_view(root, run)
                except (OSError, ValueError, KeyError, TypeError) as exc:
                    invalid[run["folder"].name] = "CN60 revision: " + str(exc)
                    revised = run
                if revised is None:
                    revision_pending.append(ident)
                    continue
                run = revised
            selected[ident] = run
            try:
                checks[run["folder"].name] = validate_run(run)
            except (OSError, ValueError, KeyError, TypeError, IndexError) as exc:
                invalid[run["folder"].name] = str(exc)
    groups = []
    keys = sorted({key[:3] for key in endpoints}, key=str)
    for center_id, kind, component in keys:
        center = selected.get(center_id)
        steps = sorted({key[3] for key in endpoints if key[:3] == (center_id, kind, component)})
        group = {
            "center_id": center_id, "method": center_methods[center_id],
            "derivative_kind": kind, "component": list(component) if kind == "force" else component,
            "component_label": point_label(kind, component), "steps": steps,
            "unit": "eV/Angstrom" if kind == "force" else "kbar",
            "absolute_tolerance": ABS_TOL[kind], "two_step_tolerance": STEP_TOL[kind],
            "center_run": center["folder"].name if center else None,
            "center_attempts": history(attempts, center_id),
            "coverage_complete": False, "status": "pending", "reasons": [],
            "points": [], "levels": {},
        }
        missing = center is None
        errors = [invalid[center["folder"].name]] if center and center["folder"].name in invalid else []
        values, bases = {level: {} for level in LEVELS}, []
        if len(steps) != 2 or tuple(steps) != STEPS[kind]:
            errors.append("Unexpected derivative step set; explicit protocol revision required")
        for step in steps:
            pair = {}
            for sign in (-1, 1):
                ident = index[(center_id, kind, component, step, sign)]
                run = selected.get(ident)
                point = {"step": step, "sign": sign, "logical_id": ident,
                         "selected_run": run["folder"].name if run else None,
                         "attempts": history(attempts, ident)}
                group["points"].append(point)
                if run is None:
                    missing = True
                    continue
                pair[sign] = run
                name = run["folder"].name
                if name in invalid:
                    errors.append(name + ": " + invalid[name])
                    continue
                if center and center["folder"].name not in invalid:
                    try:
                        comparison = compare_endpoint(center, run, checks, kind, component, step, sign)
                        point["basis_comparison"] = comparison
                        bases.append(comparison["basis_unchanged"])
                    except (ValueError, KeyError, TypeError, IndexError) as exc:
                        errors.append(name + ": " + str(exc))
            if center and len(pair) == 2 and not errors:
                for level in LEVELS:
                    delta = energy(pair[1]["result"], level) - energy(pair[-1]["result"], level)
                    derivative = (-delta / (2 * step) if kind == "force" else
                                  delta / (2 * step * center["result"]["volume_A3"]) * KBAR)
                    reference = analytic(center["result"], level, kind, component)
                    values[level][step] = {"step": step, "analytic": reference,
                                          "finite_difference": derivative,
                                          "signed_error": derivative - reference,
                                          "absolute_error": abs(derivative - reference),
                                          "energy_plus_eV": energy(pair[1]["result"], level),
                                          "energy_minus_eV": energy(pair[-1]["result"], level),
                                          "energy_difference_eV": delta}
        group["coverage_complete"] = not missing
        group["reasons"] = sorted(set(errors))
        group["basis_unchanged_all_endpoints"] = bool(bases) and all(bases) and len(bases) == 4
        all_descriptor_metadata = bool(center) and all(
            point.get("basis_comparison", {}).get("density_hamiltonian_basis_recorded")
            for point in group["points"])
        group["density_hamiltonian_basis_recorded_all_points"] = all_descriptor_metadata
        group["basis_scope"] = (
            "Wave-function Q_qG, global k mapping, grids and density/Hamiltonian pd2/pd3 Q sets."
            if all_descriptor_metadata else
            "Wave-function Q_qG, global k mapping and grids; complete density/Hamiltonian "
            "pd2/pd3 metadata is unavailable for at least one point.")
        for level in LEVELS:
            entries = [values[level][h] for h in steps if h in values[level]]
            difference = abs(entries[0]["finite_difference"] - entries[1]["finite_difference"]) \
                if len(entries) == 2 else None
            applicable = not (level == "dispersion" and group["method"] == "optB88-vdW")
            numeric_passed = (len(entries) == 2 and
                              all(x["absolute_error"] <= ABS_TOL[kind] for x in entries) and
                              difference <= STEP_TOL[kind])
            group["levels"][level] = {
                "applicable": applicable,
                "note": ("Nonlocal vdW is already inside the DFT contribution; "
                         "the separate additive contribution is identically zero."
                         if not applicable else None),
                "steps": entries, "two_step_difference": difference,
                "numerically_passed": numeric_passed if len(entries) == 2 else None,
            }
        if missing:
            group["status"] = "pending"
            group["reasons"].append("Center or an endpoint lacks a completed successful attempt")
        elif errors:
            group["status"] = "invalid_comparison"
        elif not group["basis_unchanged_all_endpoints"]:
            group["status"] = "inconclusive_basis_change"
            group["reasons"].append("Grid/global-k/PW basis changed; higher-cutoff follow-up is required")
        elif all(x["numerically_passed"] for x in group["levels"].values() if x["applicable"]):
            group["status"] = "passed"
        else:
            group["status"] = "failed"
        groups.append(group)
    counts = {name: sum(g["status"] == name for g in groups)
              for name in ("passed", "pending", "failed", "invalid_comparison", "inconclusive_basis_change")}
    complete = bool(groups) and all(g["coverage_complete"] for g in groups)
    passed = complete and counts["passed"] == len(groups)
    if not complete:
        status = "pending"
    elif passed:
        status = "representative_derivative_checks_passed"
    elif counts["failed"] or counts["invalid_comparison"]:
        status = "needs_attention"
    else:
        status = "inconclusive_basis_change"
    report = {
        "schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
        "data_view": view,
        "view_description": ("Original saved results; dispersion parameters are checked "
                             "within each method and retained in source records."
                             if view == "original" else
                             "PBE-D3 uses separately archived CN60 recombination; DFT is "
                             "unchanged. D4 and optB88 remain their original saved results. "
                             "Use revised total directly; never add dispersion a second time."),
        "status": status, "coverage_complete": complete, "passed": passed,
        "scope": "Three representative C/H force components and six stress components "
                 "at two steps for each of three methods, on the asymmetric probe. "
                 "Not all 36 forces; not a production cutoff/k-point convergence claim.",
        "tight_SCF_noise_validation_included_in_pass_flag": False,
        "follow_up": "A representative tight-SCF endpoint noise check remains separately required.",
        "quantum_calculations_started": False, "series_prefix": prefix,
        "requested_cutoff_eV": cutoff,
        "center_cutoffs_eV": {ident: selected[ident]["input"].get("ecut_eV")
                              if ident in selected else None for ident in center_methods},
        "analyzer_sha256": digest(__file__),
        "thresholds": {"force_absolute_eV_A": .001, "force_two_step_eV_A": .0005,
                       "stress_absolute_kbar": .1, "stress_two_step_kbar": .05},
        "stress_conversion_eV_A3_to_kbar": KBAR,
        "counts": counts, "expected_component_groups": len(groups),
        "expected_unique_endpoint_points": len(endpoints),
        "selected_successful_points": len(selected),
        "unreadable_in_progress_inputs": unreadable,
        "pending_CN60_revisions": revision_pending,
        "selected_runs": {
            ident: {"folder": str(run["folder"]),
                    "input_sha256": digest(run["folder"] / "input.json"),
                    "result_sha256": digest(run["folder"] / "result.json"),
                    "runner_sha256": run["result"].get("runner_sha256"),
                    "dispersion_revision": run.get("revision"),
                    "attempts": history(attempts, ident),
                    "validation_error": invalid.get(run["folder"].name),
                    "merged_basis": checks.get(run["folder"].name, {}).get("basis")}
            for ident, run in sorted(selected.items())
        },
        "groups": groups,
    }
    queue = (noise_queue(center_methods, endpoints, selected, groups) if view == "original" else
             {"status": "not_generated_for_recombined_view", "points": [],
              "reason": "Use the original DFT noise-check proposal with an explicitly "
                        "version-matched runner; the CN60 overlay did not rerun DFT."})
    return report, queue


def noise_queue(centers, endpoints, selected, groups):
    """Recommend center plus one +/- small-strain pair per method; do not submit."""
    points, suggestions = [], []
    for center_id, method in centers.items():
        center = selected.get(center_id)
        if center is None:
            suggestions.append({"method": method, "status": "pending_center"})
            continue
        candidates = [g for g in groups if g["center_id"] == center_id
                      and g["derivative_kind"] == "stress" and int(g["component"]) >= 3]
        # Prefer the largest measured absolute FD error; otherwise the strongest
        # nonzero shear signal.  This choice does not alter acceptance thresholds.
        def score(group):
            entries = group["levels"]["total"]["steps"]
            measured = max((x["absolute_error"] for x in entries), default=-1.)
            signal = abs(analytic(center["result"], "total", "stress", int(group["component"])))
            return measured, signal
        component = int(max(candidates, key=score)["component"])
        ids = [center_id]
        for sign in (-1, 1):
            key = (center_id, "stress", component, STEPS["stress"][0], sign)
            ids.append(logical_id(endpoints[key]["id"]))
        source_specs = [center["input"]]
        for ident in ids[1:]:
            run = selected.get(ident)
            if run is not None:
                source_specs.append(run["input"])
            else:
                source_specs.append(next(spec for spec in endpoints.values()
                                         if logical_id(spec["id"]) == ident))
        ready = all(ident in selected for ident in ids)
        suggestions.append({"method": method, "component": VOIGT[component],
                            "strain_step": STEPS["stress"][0],
                            "parent_ids": ids, "all_parents_completed": ready,
                            "status": "proposed_not_submitted" if ready else "pending_parent_endpoints"})
        if not ready:
            continue
        for spec in source_specs:
            item = {k: v for k, v in spec.items()
                    if k not in {"protocol_sha256", "queue_sha256", "runner_sha256_at_submission"}}
            parent = logical_id(item["id"])
            item.update(id=parent + "_tightscf", stage="scf_noise_validation",
                        validation_role="tight_scf_noise_check",
                        noise_check_parent_id=parent,
                        noise_check_parent_run=selected[parent]["folder"].name,
                        noise_check_group=f"{center_id}_{VOIGT[component]}")
            for key in ("scf_energy_eV_electron", "scf_density_electron_electron",
                        "scf_eigenstates_eV2_electron"):
                item[key] = float(spec.get(key, DEFAULT_SETTINGS[key])) / 10
            points.append(item)
    return {
        "status": "suggested_only_not_submitted",
        "purpose": "SCF-noise control at small strain; these jobs are excluded from ordinary "
                   "derivative pairing by validation_role.",
        "defaults": {}, "points": points, "stop_on_failure": True,
        "suggestions": suggestions,
        "acceptance": {"energy_change_eV_cell_per_endpoint": 1e-8,
                       "ideal_energy_change_eV_cell_per_endpoint": 5e-9,
                       "FD_stress_change_kbar": .02,
                       "rule": "Compare the identical +/- structures at base/tighter SCF; "
                               "a tighter center also checks the analytic stress."},
    }


def output_reports(root, report, queue):
    out = root / "validation_audit"
    out.mkdir(parents=True, exist_ok=True)
    prefix = report["series_prefix"]
    if not re.fullmatch(r"[A-Za-z0-9_-]+", prefix):
        raise ValueError("Unsafe series prefix for report filenames")
    # Preserve the established smooth-series report paths; every other series
    # receives an explicit prefix and cannot overwrite the old analysis.
    suffix = "" if prefix == "deriv_smooth" else "_" + prefix
    if report["data_view"] != "original":
        suffix += "_d3_cn60"
    write_json(out / f"derivative_results{suffix}.json", report)
    write_json(out / f"tight_scf_suggested_queue{suffix}.json", queue)
    rows = []
    for group in report["groups"]:
        for level, values in group["levels"].items():
            entries = values["steps"] or [{"step": step} for step in group["steps"]]
            for entry in entries:
                rows.append({
                    "method": group["method"], "center_id": group["center_id"],
                    "data_view": report["data_view"],
                    "center_run": group["center_run"], "kind": group["derivative_kind"],
                    "component": group["component_label"], "contribution": level,
                    "applicable": values["applicable"], "step": entry["step"],
                    "unit": group["unit"], "analytic": entry.get("analytic"),
                    "finite_difference": entry.get("finite_difference"),
                    "signed_error": entry.get("signed_error"),
                    "absolute_error": entry.get("absolute_error"),
                    "two_step_difference": values["two_step_difference"],
                    "absolute_tolerance": group["absolute_tolerance"],
                    "two_step_tolerance": group["two_step_tolerance"],
                    "energy_plus_eV": entry.get("energy_plus_eV"),
                    "energy_minus_eV": entry.get("energy_minus_eV"),
                    "energy_difference_eV": entry.get("energy_difference_eV"),
                    "basis_unchanged": group["basis_unchanged_all_endpoints"],
                    "coverage_complete": group["coverage_complete"],
                    "status": group["status"],
                    "reasons": " | ".join(group["reasons"]),
                })
    with (out / f"derivative_results{suffix}.csv").open("w", encoding="utf-8-sig", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]) if rows else ["status"])
        writer.writeheader()
        writer.writerows(rows)
    lines = [
        "# GPAW 代表性力与六分量应力差分",
        "",
        f"数据系列：{prefix}；视图：{report['data_view']}。{report['view_description']}",
        "",
        f"状态：**{report['status']}**；数据覆盖完成：{report['coverage_complete']}；"
        f"本轮代表性差分通过：{report['passed']}。",
        "",
        "每种方法检查三个代表性 C/H 力分量和六个独立应力分量，各用两档步长。"
        "这不等于全部36个力分量逐项检查，也不证明生产 cutoff/k 网格收敛。"
        "更严格 SCF 的噪声复核仍需单独完成；本表通过标志不包含该后续检查。",
        "",
        "力步长0.001/0.002 Å，绝对误差≤0.001 eV/Å，两档差≤0.0005 eV/Å。"
        "应变1e−5/2e−5，绝对误差≤0.1 kbar，两档差≤0.05 kbar。"
        "应力采用ASE拉正号及工程剪切，分母2hV使用该探针中心的体积。",
        "",
        "缺少成功中心/端点时列为pending，失败尝试保留在JSON中；"
        "按逻辑ID选择最后成功attempt，绝不把failed最后一帧当成功端点。"
        "G集合比较合并了所有MPI rank，并同时检查全局k映射及粗/细网格形状。"
        "新批次另检查各rank密度/Hamiltonian的pd2/pd3集合；旧批次缺少这些记录时明确注明，"
        "不能依据波函数G集合不变声称所有基组均未变化。",
        "",
        "| 方法 | 类型 | 分量 | 分项 | 最大绝对误差 | 两档差 | 单位 | 状态 |",
        "|---|---|---|---|---:|---:|---|---|",
    ]
    for group in report["groups"]:
        for level, val in group["levels"].items():
            errors = [x["absolute_error"] for x in val["steps"]]
            error = f"{max(errors):.8g}" if errors else "待完成"
            diff = f"{val['two_step_difference']:.8g}" if val["two_step_difference"] is not None else "待完成"
            status = group["status"] if val["applicable"] else "无单独色散项"
            lines.append(f"| {group['method']} | {group['derivative_kind']} | "
                         f"{group['component_label']} | {level} | {error} | {diff} | "
                         f"{group['unit']} | {status} |")
    issues = [g for g in report["groups"] if g["reasons"]]
    if issues:
        lines += ["", "## 待完成或需处理", ""]
        for group in issues:
            lines.append(f"- {group['method']} {group['component_label']}："
                         + "；".join(group["reasons"]))
    lines += ["", "## Tight-SCF 复核建议", "",
              f"tight_scf_suggested_queue{suffix}.json 仅是队列建议，本分析器不会提交。"
              "每种方法选一组敏感非零剪切的小步长±端点及中心，三项SCF阈值再收紧10倍。"
              "只有原始三个点已成功时才写入可运行队列；否则保持待完成。"
              "建议确认每端点能量变化≤1e−8 eV/胞（理想≤5e−9），"
              "并直接比较差分应力变化是否≤0.02 kbar。",
              "", f"生成时间：{report['created_utc']}", ""]
    (out / f"derivative_results{suffix}.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--prefix", default="deriv_smooth")
    parser.add_argument("--cutoff", type=float, default=None,
                        help="Optional exact cutoff selection; by default discover recorded center IDs.")
    parser.add_argument("--view", choices=("original", "d3-cn60"), default="original",
                        help="CN60 creates separate reports; never overwrites original analysis.")
    parser.add_argument("--require-complete", action="store_true",
                        help="Exit nonzero unless all agreed derivative checks pass.")
    args = parser.parse_args()
    report, queue = analyze(args.root.resolve(), args.prefix, args.cutoff, args.view)
    output_reports(args.root.resolve(), report, queue)
    print(json.dumps({key: report[key] for key in
                      ("status", "coverage_complete", "passed", "counts",
                       "expected_unique_endpoint_points", "selected_successful_points")}))
    if args.require_complete and not report["passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
