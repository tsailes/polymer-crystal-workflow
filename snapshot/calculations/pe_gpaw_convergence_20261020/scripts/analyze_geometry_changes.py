"""Geometry-only explanation of PE cell relaxation, never a DFT acceptance gate.

The original reference defines periodic covalent edges and transverse chain
images once. Final distances use those same edges and original atom order.
No calculator is constructed and no source geometry is changed.
"""
from __future__ import annotations

import argparse
from collections import deque
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import numpy as np
import spglib
from ase import Atoms
from ase.io import read
from ase.neighborlist import neighbor_list

from analyze_relaxations import inspect_run
from export_pnam_cif import match_periodic, serialize, symmetry_summary

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def record(path):
    path = Path(path)
    return {"path": str(path), "sha256": sha(path)}


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def canonical_edge(i, j, shift):
    forward = (int(i), int(j), *map(int, shift))
    reverse = (int(j), int(i), *map(int, -np.asarray(shift)))
    return min(forward, reverse)


def reference_topology(atoms):
    """Infer connectivity only on the already audited experimental reference."""
    symbols = atoms.get_chemical_symbols()
    bonds = {}
    ii, jj, shifts, distances = neighbor_list("ijSd", atoms, 1.8)
    for i, j, shift, distance in zip(ii, jj, shifts, distances):
        kind = "".join(sorted((symbols[i], symbols[j])))
        if distance < {"CC": 1.8, "CH": 1.3}.get(kind, 0):
            edge = canonical_edge(i, j, shift)
            bonds[edge] = {"i": edge[0], "j": edge[1],
                           "shift": list(edge[2:]), "kind": kind}
    bonds = list(bonds.values())
    if sorted(e["kind"] for e in bonds) != ["CC"] * 4 + ["CH"] * 8:
        raise ValueError("Reference periodic topology is not four CC and eight CH bonds")
    graph = [[] for _ in atoms]
    for edge in bonds:
        i, j, shift = edge["i"], edge["j"], np.asarray(edge["shift"], int)
        graph[i].append((j, shift))
        graph[j].append((i, -shift))
    for i, neighbors in enumerate(graph):
        expected = ["C", "C", "H", "H"] if symbols[i] == "C" else ["C"]
        if sorted(symbols[j] for j, _ in neighbors) != expected:
            raise ValueError("Reference topology does not have PE covalent coordination")
    labels = np.full(len(atoms), -1, int)
    images = np.zeros((len(atoms), 3), int)
    component = -1
    for origin in range(len(atoms)):
        if labels[origin] != -1:
            continue
        component += 1
        labels[origin] = component
        queue = deque([origin])
        while queue:
            i = queue.popleft()
            for j, shift in graph[i]:
                proposal = images[i, :2] + shift[:2]
                if labels[j] == -1:
                    labels[j] = component
                    images[j, :2] = proposal
                    queue.append(j)
                elif labels[j] != component or not np.array_equal(images[j, :2], proposal):
                    raise ValueError("Reference chains are not periodic along c only")
    if component != 1 or sorted(np.bincount(labels).tolist()) != [6, 6]:
        raise ValueError("Reference must contain two six-atom chain repeats")
    return {
        "periodic_bonds": bonds, "chain_id": labels.tolist(),
        "transverse_unwrap_image": images.tolist(),
        "reference_only_neighbor_identification_A": {"CC": 1.8, "CH": 1.3},
        "note": "Distance cutoffs identify only the audited reference graph; "
                "they are not final-geometry acceptance thresholds.",
    }


def same_chain(topology, i, j, shift):
    image = topology["transverse_unwrap_image"]
    return (topology["chain_id"][i] == topology["chain_id"][j]
            and np.array_equal(np.asarray(shift)[:2],
                               np.asarray(image[j])[:2] - np.asarray(image[i])[:2]))


def geometry(atoms, topology):
    if (len(atoms) != 12 or atoms.get_chemical_formula() != "C4H8"
            or not atoms.pbc.all() or not np.isfinite(atoms.positions).all()
            or not np.isfinite(atoms.cell.array).all() or atoms.get_volume() <= 0):
        raise ValueError("Invalid geometry for this C4H8 diagnostic")
    symbols = atoms.get_chemical_symbols()
    bond_records, angles = [], []
    for edge in topology["periodic_bonds"]:
        i, j = edge["i"], edge["j"]
        vector = (atoms.positions[j] + np.asarray(edge["shift"]) @ atoms.cell.array
                  - atoms.positions[i])
        bond_records.append({**edge, "distance_A": float(np.linalg.norm(vector))})
        if edge["kind"] == "CC":
            angle = float(np.degrees(np.arctan2(vector[1], vector[0])) % 180)
            angles.append(min(angle, 180 - angle))
    distances = {kind: [e["distance_A"] for e in bond_records if e["kind"] == kind]
                 for kind in ("CC", "CH")}
    cutoff = float(max(atoms.cell.lengths()) + 1.0)
    contacts = {"CC": {}, "HH": {}}
    ii, jj, shifts, dd = neighbor_list("ijSd", atoms, cutoff)
    for i, j, shift, distance in zip(ii, jj, shifts, dd):
        kind = symbols[i] + symbols[j]
        if kind not in contacts or same_chain(topology, i, j, shift):
            continue
        edge = canonical_edge(i, j, shift)
        contacts[kind][edge] = {"i": edge[0], "j": edge[1], "shift": list(edge[2:]),
                                "distance_A": float(distance)}
    nearest = {}
    for kind, found in contacts.items():
        if not found:
            raise ValueError(f"No interchain {kind} contact within diagnostic search radius")
        ordered = sorted(found.values(), key=lambda v: v["distance_A"])
        minimum = ordered[0]["distance_A"]
        nearest[kind] = {
            "minimum_distance_A": minimum, "search_cutoff_A": cutoff,
            "unique_periodic_pairs_within_cutoff": len(ordered),
            "degenerate_nearest_pairs": [x for x in ordered
                                        if abs(x["distance_A"] - minimum) <= 1e-8],
        }
    operations = spglib.get_symmetry_from_database(297)
    closure = []
    frac = atoms.get_scaled_positions()
    for rotation, translation in zip(operations["rotations"], operations["translations"]):
        transformed = Atoms(numbers=atoms.numbers, cell=atoms.cell, pbc=True,
                            scaled_positions=frac @ rotation.T + translation)
        closure.append(match_periodic(atoms, transformed, np.inf)["maximum_periodic_distance_A"])
    symmetry = {"target_Hall_number": 297,
                "target_operation_residuals_A": closure,
                "target_operation_max_residual_A": max(closure)}
    for key, hall in (("Hall297_input_setting", 297), ("default_standard_setting", 0)):
        try:
            _, summary = symmetry_summary(atoms, 1e-5, hall)
            symmetry[key] = summary
        except ValueError as exc:
            symmetry[key] = {"not_identified": str(exc)}
    return {
        "formula": atoms.get_chemical_formula(), "natoms": len(atoms),
        "cell_A": atoms.cell.array.tolist(), "cell_lengths_A": atoms.cell.lengths().tolist(),
        "volume_A3": float(atoms.get_volume()),
        "density_g_cm3": float(atoms.get_masses().sum() * 1.66053906660 / atoms.get_volume()),
        "periodic_covalent_bonds": bond_records,
        "CC_range_A": [min(distances["CC"]), max(distances["CC"])],
        "CH_range_A": [min(distances["CH"]), max(distances["CH"])],
        "chain_plane_acute_angle_to_a_deg": float(np.mean(angles)),
        "nearest_interchain_contacts": nearest, "symmetry": symmetry,
    }


def accepted_geometry(folder, original):
    run = inspect_run(folder)
    if run["acceptance"] != "verified_relaxation":
        raise ValueError(f"{folder.name} is not an accepted final relaxation: {run['acceptance']}")
    result_path, trajectory = folder / "result.json", folder / "trajectory.traj"
    result = load(result_path)
    if result["source_geometry_sha256"] != sha(original):
        raise ValueError("Relaxation has a different original geometry")
    atoms = read(trajectory, index=-1)
    final = result["final_evaluation"]
    if (atoms.cell.array.tolist() != final["cell_A"]
            or atoms.positions.tolist() != final["positions_A"]
            or atoms.numbers.tolist() != read(original).numbers.tolist()):
        raise ValueError("Final full-precision trajectory or atom order differs from result")
    return atoms, {"run_id": folder.name, "result": record(result_path),
                   "trajectory": record(trajectory),
                   "acceptance": run["acceptance"], "precision": run["precision"]}


def report_markdown(report):
    baseline = report["structures"]["Tashiro"]["geometry"]
    lines = [
        "# PE 优化后几何变化（解释性对照）", "",
        "本表对比原始实验精修模型与两套数值精度的 PBE-D3(BJ) 0 K 优化终态。"
        "这是几何解释，不增加收敛验收门槛，也不把实验精修结构当成无误差真值。", "",
        "| 结构 | a / Å | b / Å | c / Å | 密度 / g·cm⁻³ | C–C / Å | C–H / Å | 最近链间 C···C / Å | 最近链间 H···H / Å |",
        "|---|---:|---:|---:|---:|---|---|---:|---:|",
    ]
    for name, item in report["structures"].items():
        g = item["geometry"]
        a, b, c = g["cell_lengths_A"]
        cc, ch = g["CC_range_A"], g["CH_range_A"]
        inter = g["nearest_interchain_contacts"]
        lines.append(f"| {name} | {a:.7f} | {b:.7f} | {c:.7f} | {g['density_g_cm3']:.7f} "
                     f"| {cc[0]:.7f}–{cc[1]:.7f} | {ch[0]:.7f}–{ch[1]:.7f} "
                     f"| {inter['CC']['minimum_distance_A']:.7f} "
                     f"| {inter['HH']['minimum_distance_A']:.7f} |")
    lines += ["", "周期 C–C/C–H 键由原始模型定义并冻结；优化后按原有原子顺序和周期平移重算，"
              "没有根据终态距离重新猜键。链间接触包含横向周期复制的邻链，排除同一无限链沿 c 的周期像。"
              "密度由晶胞质量/体积计算，不等同于含非晶区、孔隙等因素的样品体密度。", ""]
    for name, item in report["structures"].items():
        g = item["geometry"]
        if name != "Tashiro":
            delta = (np.asarray(g["cell_lengths_A"]) / baseline["cell_lengths_A"] - 1) * 100
            lines.append(f"- {name} 相对原始模型：a、b、c 分别变化 "
                         f"{delta[0]:+.4f}%、{delta[1]:+.4f}%、{delta[2]:+.4f}%；"
                         f"密度变化 {(g['density_g_cm3']/baseline['density_g_cm3']-1)*100:+.4f}%。")
        sym = g["symmetry"]
        hall = sym["Hall297_input_setting"]
        lines.append(f"- {name}：目标 Pnam 操作最大匹配残差 "
                     f"{sym['target_operation_max_residual_A']:.3g} Å；"
                     f"Hall 297 识别结果 {hall.get('international', '未识别')}，"
                     f"操作数 {hall.get('operation_count', '未知')}。")
    lines += ["", "实验模型的氢位置包含精修约束；0 K 静态计算也未包含热膨胀和零点振动。"
              "因此实验与计算差异不能仅归因于泛函误差。这里不使用密度 95%–105% 或"
              "链间 C···C ≥3.4 Å 作为通过条件。", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-only", action="store_true")
    parser.add_argument("--candidate", default="cell_candidate")
    parser.add_argument("--reference", default="cell_reference")
    args = parser.parse_args()
    original = ROOT / "reference/POSCAR"
    atoms = read(original)
    topology = reference_topology(atoms)
    baseline = geometry(atoms, topology)
    audit_path = ROOT / "reference_audit/reference_audit.json"
    audited = load(audit_path)["experimental_reference"]["geometry"]
    for key in ("CC_range_A", "CH_range_A", "density_g_cm3",
                "chain_plane_acute_angle_to_a_deg"):
        if not np.allclose(baseline[key], audited[key], atol=1e-10, rtol=0):
            raise ValueError(f"Reference diagnostic does not reproduce previous audit: {key}")
    report = {
        "schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "reference_selfcheck_only" if args.reference_only else "geometry_comparison_available",
        "quantum_calculations_started": False, "acceptance_gate": False,
        "script": record(__file__), "reused_helpers": [
            record(ROOT / "scripts/export_pnam_cif.py"),
            record(ROOT / "scripts/analyze_relaxations.py")],
        "existing_reference_audit": record(audit_path),
        "reference_values_reproduced_from_existing_audit": True,
        "software": {p: importlib.metadata.version(p) for p in ("ase", "numpy", "spglib", "scipy")},
        "reference_topology": topology,
        "structures": {"Tashiro": {"source": record(original), "geometry": baseline}},
    }
    if not args.reference_only:
        if args.candidate == args.reference:
            raise ValueError("Candidate and reference must be distinct relaxation attempts")
        for label, name in (("candidate", args.candidate), ("reference", args.reference)):
            final, provenance = accepted_geometry(ROOT / "relaxations" / name, original)
            report["structures"][label] = {**provenance, "geometry": geometry(final, topology)}
    output = ROOT / "validation_audit" / (
        "geometry_reference_selfcheck.json" if args.reference_only else "optimized_geometry_changes.json")
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False,
                                 default=serialize) + "\n", encoding="utf-8")
    if not args.reference_only:
        (ROOT / "reports/优化后几何变化.md").write_text(report_markdown(report), encoding="utf-8")
    print(json.dumps({"status": report["status"], "output": str(output),
                      "baseline_density_g_cm3": baseline["density_g_cm3"],
                      "baseline_nearest_interchain_contacts": baseline["nearest_interchain_contacts"]},
                     ensure_ascii=False))


if __name__ == "__main__":
    main()
