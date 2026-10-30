"""Independent, finite Pnam PE packing grid; no refined coordinates are inputs."""
from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
from pathlib import Path

import ase
import gemmi
import numpy as np
import spglib
from ase import Atoms
from ase.io import read, write
from ase.neighborlist import neighbor_list
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parent
CELL = np.diag([7.432, 4.929, 2.545])
SG = gemmi.SpaceGroup("P n a m")
OPS = list(SG.operations())
ROTATIONS = np.array([o.rot for o in OPS], float) / 24
TRANSLATIONS = np.array([o.tran for o in OPS], float) / 24
VDW = {"C": 1.70, "H": 1.20}


def dump(path, data):
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2,
                                     allow_nan=False) + "\n", encoding="utf-8")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def assignment_distance(fa, fb, cell, symbols):
    mapping = np.empty(len(fa), int)
    residual = np.empty(len(fa))
    for element in sorted(set(symbols)):
        ids = np.flatnonzero(np.asarray(symbols) == element)
        delta = fa[ids, None, :] - fb[None, ids, :]
        delta -= np.rint(delta)
        distance = np.linalg.norm(delta @ cell, axis=2)
        rows, cols = linear_sum_assignment(distance)
        mapping[ids[rows]] = ids[cols]
        residual[ids[rows]] = distance[rows, cols]
    return mapping, residual


def equivalent(first, second, rms_tol=1e-6, max_tol=2e-6):
    """Exact orthorhombic lattice isometries + origin shift + element bijection.

    All three lengths are distinct; signed diagonal operations are the full
    lattice point group. Test every C anchor; do not merge on a histogram alone.
    """
    if len(first) != len(second) or not np.allclose(first.cell, second.cell, atol=1e-8):
        return None
    a, b = first.get_scaled_positions(), second.get_scaled_positions()
    symbols = first.get_chemical_symbols()
    if symbols != second.get_chemical_symbols():
        return None
    anchor = symbols.index("C")
    best = None
    for signs in itertools.product((-1, 1), repeat=3):
        transformed = a * signs
        for j in np.flatnonzero(np.asarray(symbols) == "C"):
            shift = b[j] - transformed[anchor]
            mapping, residual = assignment_distance(
                transformed + shift, b, first.cell.array, symbols)
            rms = float(np.sqrt(np.mean(residual**2)))
            record = dict(rms_A=rms, max_A=float(residual.max()),
                          signs=list(signs), origin_shift_fractional=shift.tolist(),
                          atom_mapping=mapping.tolist())
            if best is None or rms < best["rms_A"]:
                best = record
    if best["rms_A"] <= rms_tol and best["max_A"] <= max_tol:
        return best
    return None


def symmetry_audit(atoms):
    frac = atoms.get_scaled_positions()
    errors = []
    for rotation, translation in zip(ROTATIONS, TRANSLATIONS):
        _, residual = assignment_distance(frac @ rotation.T + translation, frac,
                                          atoms.cell.array, atoms.get_chemical_symbols())
        errors.append(float(residual.max()))
    dataset = spglib.get_symmetry_dataset(
        (atoms.cell.array, frac, atoms.numbers), symprec=1e-5)
    return dict(target="Pnam", hall_number=297,
                operation_max_displacement_A=errors,
                contains_target_group=bool(max(errors) < 1e-5),
                detected_number=int(dataset.number) if dataset else None,
                detected_symbol=dataset.international if dataset else None,
                detected_hall_number=int(dataset.hall_number) if dataset else None)


def build(theta, center=(0.0, 0.0), phase=0.25, cc=1.54, ch=1.09):
    transverse = math.sqrt(cc**2 - (CELL[2, 2] / 2)**2)
    q = cc / (3 * transverse)
    phi = math.radians(theta)
    e = np.array([math.cos(phi), math.sin(phi), 0])
    n = np.array([-math.sin(phi), math.cos(phi), 0])
    offset = np.array([center[0], center[1], phase]) @ CELL
    carbon = offset + 0.5 * transverse * e
    asu_cart = [carbon, carbon + ch * (q * e + math.sqrt(1-q*q) * n),
                carbon + ch * (q * e - math.sqrt(1-q*q) * n)]
    asu_frac = np.asarray(asu_cart) @ np.linalg.inv(CELL)
    positions, numbers, asu_ids, operation_ids = [], [], [], []
    for site, (number, frac) in enumerate(zip((6, 1, 1), asu_frac)):
        for oi, (rotation, translation) in enumerate(zip(ROTATIONS, TRANSLATIONS)):
            p = (rotation @ frac + translation) % 1.0
            # Normalize near-boundary roundoff before constructing actual Cartesian
            # coordinates; ASE otherwise may wrap 1.0 to 0 only when reading frac.
            p[np.abs(p-1.0) < 1e-12] = 0.0
            p[np.abs(p) < 1e-12] = 0.0
            if any(num == number and np.linalg.norm(
                    (p - old - np.rint(p - old)) @ CELL) < 1e-8
                   for old, num in zip(positions, numbers)):
                continue
            positions.append(p)
            numbers.append(number)
            asu_ids.append(site)
            operation_ids.append(oi)
    order = sorted(range(len(numbers)), key=lambda i: (-numbers[i], asu_ids[i], operation_ids[i]))
    atoms = Atoms(numbers=np.array(numbers)[order],
                  scaled_positions=np.array(positions)[order], cell=CELL, pbc=True)
    centers = np.array([center, [(center[0]+0.5) % 1, (-center[1]+0.5) % 1]])
    frac = atoms.get_scaled_positions()
    chain, images = [], []
    for f in frac:
        deltas = f[:2] - centers
        deltas -= np.rint(deltas)
        ci = int(np.argmin(np.linalg.norm(deltas * np.diag(CELL)[:2], axis=1)))
        unwrapped_xy = centers[ci] + deltas[ci]
        image = [*np.rint(unwrapped_xy-f[:2]).astype(int), 0]
        chain.append(ci)
        images.append(image)
    atoms.set_array("chain_id", np.array(chain, int))
    atoms.set_array("wrap_image", np.array(images, int))
    atoms.set_array("asu_id", np.array(asu_ids, int)[order])
    atoms.info.update(solved_coordinates_used=False, chemical_repeat_count=4,
                      setting_angle_deg=float(theta), target_spacegroup="Pnam")
    return atoms, dict(asu_fractional=asu_frac.tolist(),
                       asu_elements=["C", "H", "H"], chain_centers_fractional=centers.tolist(),
                       cc_A=cc, ch_A=ch,
                       CCC_deg=float(np.degrees(2*np.arcsin(CELL[2, 2]/(2*cc)))))


def same_physical_chain(atoms, i, j, shift):
    wrap = atoms.arrays["wrap_image"]
    return (atoms.arrays["chain_id"][i] == atoms.arrays["chain_id"][j] and
            np.array_equal(np.asarray(shift)[:2], (wrap[j]-wrap[i])[:2]))


def build_topology(atoms):
    ii, jj, shifts, distances = neighbor_list("ijSd", atoms, 1.81)
    symbols = atoms.get_chemical_symbols()
    bonds = {}
    for i, j, s, d in zip(ii, jj, shifts, distances):
        kind = "".join(sorted([symbols[i], symbols[j]]))
        limit = {"CC": 1.8, "CH": 1.3}.get(kind, 0)
        # Classify intended chain, never excuse a close intermolecular collision.
        if d < limit and same_physical_chain(atoms, i, j, s):
            forward, reverse = (int(i), int(j), *map(int, s)), (int(j), int(i), *map(int, -s))
            k = min(forward, reverse)
            bonds[k] = dict(i=k[0], j=k[1], shift=list(k[2:]), kind=kind)
    edges = list(bonds.values())
    if sum(e["kind"] == "CC" for e in edges) != 4 or sum(e["kind"] == "CH" for e in edges) != 8:
        raise ValueError("Incorrect periodic chemical topology")
    return dict(periodic_bonds=edges, chain_id=atoms.arrays["chain_id"].tolist(),
                wrap_image=atoms.arrays["wrap_image"].tolist())


def excluded_pairs(topology, natoms):
    graph = [[] for _ in range(natoms)]
    for edge in topology["periodic_bonds"]:
        s = np.array(edge["shift"], int)
        graph[edge["i"]].append((edge["j"], s))
        graph[edge["j"]].append((edge["i"], -s))
    excluded = set()
    for i in range(natoms):
        for j, s in graph[i]:
            excluded.add((i, j, *map(int, s)))
            for k, t in graph[j]:
                excluded.add((i, k, *map(int, s+t)))
    return excluded


def audit_geometry(atoms, topology, initial=False):
    for name in ("chain_id", "wrap_image"):
        if name not in atoms.arrays:
            atoms.set_array(name, np.array(topology[name], int))
    symbols = atoms.get_chemical_symbols()
    excluded = excluded_pairs(topology, len(atoms))
    ii, jj, shifts, distances = neighbor_list("ijSd", atoms, 5.5)
    clashes, nonbonded, inter_cc = [], [], []
    for i, j, s, d in zip(ii, jj, shifts, distances):
        if (int(i), int(j), *map(int, s)) > (int(j), int(i), *map(int, -s)):
            continue
        record = dict(i=int(i), j=int(j), shift=s.tolist(),
                      pair=symbols[i]+symbols[j], distance_A=float(d))
        if (int(i), int(j), *map(int, s)) not in excluded:
            limit = 0.85 * (VDW[symbols[i]]+VDW[symbols[j]])
            record["vdw_threshold_A"] = limit
            record["ratio_to_vdw_threshold"] = float(d/limit)
            nonbonded.append(record)
            if d < limit-1e-8:
                clashes.append(record)
        if symbols[i] == symbols[j] == "C" and not same_physical_chain(atoms, i, j, s):
            inter_cc.append(record)
    bond_records = []
    degree = {i: [] for i in range(len(atoms))}
    for edge in topology["periodic_bonds"]:
        i, j, s = edge["i"], edge["j"], np.asarray(edge["shift"])
        d = float(np.linalg.norm(atoms.positions[j]+s@atoms.cell.array-atoms.positions[i]))
        bond_records.append({**edge, "distance_A": d})
        degree[i].append(symbols[j])
        degree[j].append(symbols[i])
    valid_degree = all(sorted(degree[i]) == (["C", "C", "H", "H"] if s == "C" else ["C"])
                       for i, s in enumerate(symbols))
    cc = [e["distance_A"] for e in bond_records if e["kind"] == "CC"]
    ch = [e["distance_A"] for e in bond_records if e["kind"] == "CH"]
    cc_approx = all(1.50 <= d <= 1.57 for d in cc)
    chemical_ok = valid_degree and cc_approx and all(1.04 <= d <= 1.15 for d in ch)
    mincc = min(inter_cc, key=lambda x: x["distance_A"])
    density = float(atoms.get_masses().sum()*1.66053906660/atoms.get_volume())
    density_ratio = density/0.997
    symmetry = symmetry_audit(atoms)
    formula_ok = symbols.count("C") == 4 and symbols.count("H") == 8
    passed = (formula_ok and chemical_ok and not clashes and mincc["distance_A"] >= 3.4-1e-8 and
              .95 <= density_ratio <= 1.05 and symmetry["contains_target_group"])
    return dict(passed=bool(passed), formula=atoms.get_chemical_formula(), natoms=len(atoms),
                cell_A=atoms.cell.array.tolist(), volume_A3=float(atoms.get_volume()),
                density_g_cm3=density, density_reference_g_cm3=0.997,
                density_reference_kind="External crystallographic crystalline-PE reference; not Tashiro sample bulk density",
                density_ratio=density_ratio, density_ratio_passed=bool(.95 <= density_ratio <= 1.05),
                chemistry_topology_passed=bool(chemical_ok), bonds=bond_records,
                CC_range_A=[min(cc), max(cc)], CH_range_A=[min(ch), max(ch)],
                CC_strict_1p53_1p54_passed=bool(all(1.53-1e-8 <= d <= 1.54+1e-8 for d in cc)),
                CC_approximate_criterion_A=[1.50, 1.57],
                CC_approximate_criterion_note="Explicit chemical sanity margin; strict interval separately reported",
                minimum_interchain_CC=mincc, interchain_CC_3p4_passed=bool(mincc["distance_A"] >= 3.4-1e-8),
                minimum_nonbonded_ratio=min(nonbonded, key=lambda x: x["ratio_to_vdw_threshold"]),
                vdw_clashes=clashes, exclusion="Periodic chemical graph paths of 1 or 2 bonds only",
                symmetry=symmetry, initial=initial)


def symmetry_fingerprint(atoms):
    """Candidate key only; equivalence requires atom-mapping confirmation."""
    frac = atoms.get_scaled_positions()
    d = frac[:, None, :] - frac[None, :, :]
    d -= np.rint(d)
    dist = np.linalg.norm(d@atoms.cell.array, axis=2)
    symbols = atoms.get_chemical_symbols()
    pairs = {}
    for i in range(len(atoms)):
        for j in range(i):
            kind = "".join(sorted((symbols[i], symbols[j])))
            pairs.setdefault(kind, []).append(float(dist[i, j]))
    audit = symmetry_audit(atoms)
    return dict(composition=atoms.get_chemical_formula(),
                target_group_contained=audit["contains_target_group"],
                detected_spacegroup=audit["detected_number"],
                pair_distributions_A={k: sorted(v) for k, v in pairs.items()})


def export_structure(atoms, folder, stem="input"):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    clean = atoms.copy()
    clean.calc = None
    clean.set_constraint()
    write(folder / f"{stem}.extxyz", clean)
    write(folder / "POSCAR", clean, format="vasp", direct=True, sort=False, vasp5=True)
    # Use precisely the imposed target group, even when a structure has higher symmetry.
    frac = clean.get_scaled_positions()
    symbols = clean.get_chemical_symbols()
    unseen, reps = set(range(len(clean))), []
    while unseen:
        i = min(unseen)
        reps.append(i)
        for rot, tr in zip(ROTATIONS, TRANSLATIONS):
            target = rot@frac[i]+tr
            delta = frac-target
            delta -= np.rint(delta)
            for j in range(len(clean)):
                if symbols[j] == symbols[i] and np.linalg.norm(delta[j]@clean.cell.array) < 1e-5:
                    unseen.discard(j)
    lines = [f"data_{stem}", "_chemical_formula_sum 'C H2'", "_cell_formula_units_Z 4",
             "_space_group_name_H-M_alt 'P n a m'", "_space_group_name_Hall '-P 2c 2n'",
             "_space_group_IT_number 62"]
    for name, val in zip(("a", "b", "c"), clean.cell.lengths()):
        lines.append(f"_cell_length_{name} {val:.12f}")
    for name, val in zip(("alpha", "beta", "gamma"), clean.cell.angles()):
        lines.append(f"_cell_angle_{name} {val:.12f}")
    lines += ["loop_", "_space_group_symop_id", "_space_group_symop_operation_xyz"]
    lines += [f"{i+1} '{o.triplet()}'" for i, o in enumerate(OPS)]
    lines += ["loop_", "_atom_site_label", "_atom_site_type_symbol",
              "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z", "_atom_site_occupancy"]
    for i in reps:
        lines.append(f"{symbols[i]}{i+1} {symbols[i]} " +
                     " ".join(f"{v:.14f}" for v in frac[i]) + " 1")
    cif = folder / f"{stem}.cif"
    cif.write_text("\n".join(lines)+"\n")
    roundtrip = {}
    for path in (cif, folder/"POSCAR", folder/f"{stem}.extxyz"):
        restored = read(path)
        # CIF may change element ordering. Reorder by species for same matcher.
        order = sorted(range(len(restored)), key=lambda k: -restored.numbers[k])
        restored = restored[order]
        if len(restored) != len(clean) or not np.allclose(clean.cell, restored.cell, atol=1e-8):
            raise ValueError(f"Round trip failed: {path}")
        _, residual = assignment_distance(clean.get_scaled_positions(), restored.get_scaled_positions(),
                                          clean.cell.array, symbols)
        if residual.max() > 2e-7:
            raise ValueError(f"Round trip changed coordinates: {path}")
        roundtrip[path.name] = dict(sha256=sha(path), passed=True,
                                   rms_A=float(np.sqrt(np.mean(residual**2))),
                                   max_A=float(residual.max()))
    dump(folder/"file_validation.json", roundtrip)
    return roundtrip


def main():
    if (ROOT/"generation_manifest.json").exists():
        raise FileExistsError("Generation is immutable; use a new directory for changed protocol")
    protocol = dict(
        scope="Finite rigid-chain grid under specified cell, Pnam, two chains, 2-carbon axial period",
        coverage_claim="Exhaustive only over the declared grid; not all continuous polymer conformations",
        cell_A=CELL.tolist(), spacegroup="Pnam", hall_number=297,
        chemistry="C4H8, two infinite all-trans PE chains along c",
        solved_coordinates_used=False,
        angle_grid_deg=list(range(0, 180, 5)),
        angle_180="Same chain plane as 0; endpoint equivalence separately tested",
        transverse_center_grid=[[0, 0], [0, .5], [.5, 0], [.5, .5]],
        carbon_phase_grid=[.25, .75],
        translation_scope="All symmetry-compatible discrete centers/phases for this two-chain generic-orientation model",
        nonbonded_filter=dict(vdw_radii_A=VDW, factor=.85, excluded_graph_distances=[1, 2],
                             source="Bondi 1964, DOI 10.1021/j100785a001; ASE data.vdw"),
        additional_interchain_CC_min_A=3.4,
        density_reference=dict(value_g_cm3=.997,
            kind="Crystallographic crystalline-PE literature reference; not independently measured density of Tashiro sample",
            source="Martin and Passaglia 1966, J Res NBS 70A, 221-225",
            url="https://pmc.ncbi.nlm.nih.gov/articles/PMC6687595/",
            ratio_limits=[.95, 1.05],
            limitation="Fixed cell and composition make density identical across candidate packings"),
        initial_cc_A=1.54, initial_ch_A=1.09,
        target_symmetry=dict(rotations=ROTATIONS.astype(int).tolist(),
                             translations=TRANSLATIONS.tolist(),
                             operations_xyz=[o.triplet() for o in OPS]),
        software=dict(ase=ase.__version__, gemmi=gemmi.__version__, spglib=spglib.__version__,
                      numpy=np.__version__), generator_sha256=sha(__file__),
        xrd_role="Tashiro measured d positions as external benchmark only; no solved atomic coordinates",
        numerical_optimizer_policy="All geometrically distinct survivors receive GFN2; no coarse-energy cutoff")
    dump(ROOT/"protocol.json", protocol)
    representatives, entries, rejection_counts = [], [], {}
    for index, (theta, cx, cy, phase) in enumerate(itertools.product(
            protocol["angle_grid_deg"], (0, .5), (0, .5), (.25, .75)), 1):
        cid = f"cand_{index:04d}"
        atoms, construction = build(theta, (cx, cy), phase)
        topology = build_topology(atoms)
        audit = audit_geometry(atoms, topology, initial=True)
        params = dict(setting_angle_deg=theta, center_fractional=[cx, cy],
                      carbon_phase_fractional=phase, conformation="all-trans")
        entry = dict(candidate_id=cid, parameters=params, passed=audit["passed"],
                     min_interchain_CC_A=audit["minimum_interchain_CC"]["distance_A"],
                     min_nonbonded_ratio=audit["minimum_nonbonded_ratio"]["ratio_to_vdw_threshold"],
                     n_clashes=len(audit["vdw_clashes"]),
                     symmetry=audit["symmetry"], status="rejected")
        if not audit["passed"]:
            entry["reasons"] = []
            if audit["vdw_clashes"]:
                entry["reasons"].append("nonbonded_vdw_overlap")
                entry["limiting_clash"] = min(audit["vdw_clashes"],
                                              key=lambda x: x["ratio_to_vdw_threshold"])
            if not audit["interchain_CC_3p4_passed"]:
                entry["reasons"].append("interchain_CC_below_3p4_A")
            if not audit["chemistry_topology_passed"]:
                entry["reasons"].append("chemical_topology")
            for reason in entry["reasons"]:
                rejection_counts[reason] = rejection_counts.get(reason, 0)+1
        else:
            duplicate = None
            for rep in representatives:
                match = equivalent(atoms, rep["atoms"])
                if match is not None:
                    duplicate = (rep["candidate_id"], match)
                    break
            if duplicate:
                entry.update(status="exact_symmetry_duplicate",
                             representative_id=duplicate[0], equivalence=duplicate[1])
            else:
                folder = ROOT/"candidates"/cid
                export_structure(atoms, folder)
                dump(folder/"reference_topology.json", topology)
                dump(folder/"geometry_validation.json", audit)
                dump(folder/"construction.json", {**construction, **params})
                dump(folder/"symmetry_fingerprint.json", symmetry_fingerprint(atoms))
                entry.update(status="unique_survivor",
                             input_path=f"candidates/{cid}/input.extxyz")
                representatives.append(dict(candidate_id=cid, atoms=atoms))
        entries.append(entry)
    good = [e for e in entries if e["passed"]]
    unique = [e for e in entries if e["status"] == "unique_survivor"]
    manifest = dict(status="geometric_generation_completed", protocol_sha256=sha(ROOT/"protocol.json"),
                    target_symmetry=protocol["target_symmetry"],
                    enumerated_count=len(entries), physical_survivor_count=len(good),
                    unique_initial_count=len(unique), rejected_count=len(entries)-len(good),
                    exact_duplicates_count=len(good)-len(unique),
                    rejection_counts=rejection_counts, unique_survivors=unique,
                    records=entries, quantum_calculations_performed=False)
    dump(ROOT/"generation_manifest.json", manifest)
    with (ROOT/"generation_log.jsonl").open("w", encoding="utf-8") as f:
        for e in entries:
            f.write(json.dumps(e, ensure_ascii=False)+"\n")
    with (ROOT/"candidate_catalog.csv").open("w", newline="", encoding="utf-8-sig") as f:
        fields = ["candidate_id", "status", "setting_angle_deg", "center_fractional",
                  "carbon_phase_fractional", "min_interchain_CC_A", "min_nonbonded_ratio",
                  "n_clashes", "representative_id"]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for e in entries:
            writer.writerow({k: {**e, **e["parameters"]}.get(k, "") for k in fields})
    print(json.dumps({k: manifest[k] for k in ("enumerated_count", "physical_survivor_count",
          "unique_initial_count", "rejected_count", "exact_duplicates_count", "rejection_counts")}))


if __name__ == "__main__":
    main()
