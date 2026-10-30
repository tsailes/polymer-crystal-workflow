#!/usr/bin/env python3
"""Make readable, reproducible views of the constructed PE candidates.

Writes only to plots/ (or --output-dir). No calculator, structural modification,
or reference coordinates are used. Periodic bond shifts come from each
candidate's immutable topology; bonds are never joined across the whole box.
Rerunning refreshes the images and provenance for the current saved inputs.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
from ase.io import read


ROOT = Path(__file__).resolve().parent
CARBON = "#303c48"
HYDROGEN = "#6ab0d0"
CELL_EDGE = "#727b82"
CLUSTER_COLORS = ["#28679b", "#dc7d28", "#758c36", "#a95c89", "#856a46"]
MARKERS = ["o", "s", "^", "D", "v", "P", "X"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def setup_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11,
        "axes.titlesize": 12, "axes.labelsize": 11,
        "axes.spines.top": False, "axes.spines.right": False,
        "savefig.facecolor": "white", "figure.facecolor": "white",
        "svg.fonttype": "none",
    })


def load_candidates(root):
    candidates = []
    for folder in sorted((root / "candidates").glob("cand*")):
        if not folder.is_dir():
            continue
        structure = folder / "input.extxyz"
        if not structure.exists():
            structure = folder / "input.cif"
        if not structure.exists():
            continue
        atoms = read(structure)
        construction_path = folder / "construction.json"
        topology_path = folder / "reference_topology.json"
        construction, topology = load_json(construction_path), load_json(topology_path)
        cell = np.asarray(atoms.cell)
        if not np.allclose(cell, np.diag(np.diag(cell)), atol=1e-7):
            raise ValueError("These projected views require the supplied orthorhombic a,b,c setting")
        if not np.all(atoms.pbc):
            raise ValueError("The candidate must have 3D periodic boundaries")
        # Topology is indexed in the archived atom order. Do not sort atoms.
        for bond in topology["periodic_bonds"]:
            i, j = bond["i"], bond["j"]
            vector = atoms.positions[j] + np.asarray(bond["shift"]) @ cell - atoms.positions[i]
            distance = float(np.linalg.norm(vector))
            lower, upper = (1.3, 1.8) if bond["kind"] == "CC" else (0.85, 1.3)
            if not lower <= distance <= upper:
                raise ValueError(f"Topology/order mismatch in {folder.name}: {bond} ({distance:g} A)")
        candidates.append(dict(
            identifier=folder.name, angle=float(construction["setting_angle_deg"]),
            atoms=atoms, topology=topology,
            chain_centers=np.asarray(construction["chain_centers_fractional"]),
            sources=[dict(path=str(p), sha256=sha(p))
                     for p in [structure, construction_path, topology_path]],
        ))
    if not candidates:
        raise ValueError("No archived candidates found")
    candidates.sort(key=lambda c: (c["angle"], c["identifier"]))
    first_cell = np.asarray(candidates[0]["atoms"].cell)
    if any(not np.allclose(c["atoms"].cell, first_cell, atol=1e-7) for c in candidates):
        raise ValueError("A shared-scale candidate sheet requires identical initial cells")
    return candidates


def projected_geometry(candidate, margin):
    """Return short bond segments and displayed atoms, including nearby images."""
    atoms = candidate["atoms"]
    cell = np.asarray(atoms.cell)
    a, b = cell[0, 0], cell[1, 1]
    low, high = np.array([-margin, -margin]), np.array([a + margin, b + margin])
    chains = np.asarray(candidate["topology"]["chain_id"])
    wraps = np.asarray(candidate["topology"]["wrap_image"])
    centers = candidate["chain_centers"] * np.array([a, b])
    # Draw whole chain sections whose axes lie inside/on the primary cell.
    # This includes the four periodic images of a corner chain, while avoiding
    # isolated H dots from distant chain axes just outside the viewing window.
    chain_images = [(chain, np.array([sx, sy, 0]) @ cell)
                    for chain, center in enumerate(centers)
                    for sx, sy in itertools.product((-1, 0, 1), repeat=2)
                    if np.all(center + np.array([sx*a, sy*b]) >= -1e-8)
                    and np.all(center + np.array([sx*a, sy*b]) <= [a+1e-8, b+1e-8])]
    segments, seen_bonds = [], set()
    for chain, translate in chain_images:
        for bond in candidate["topology"]["periodic_bonds"]:
            if chains[bond["i"]] != chain:
                continue
            anchor = wraps[bond["i"]] @ cell + translate
            first = (atoms.positions[bond["i"]] + anchor)[:2]
            second = (atoms.positions[bond["j"]] + anchor
                      + np.asarray(bond["shift"]) @ cell)[:2]
            if np.any(np.maximum(first, second) < low) or np.any(np.minimum(first, second) > high):
                continue
            # Two different z images can project to the same segment.
            key = tuple(sorted((tuple(np.round(first, 7)), tuple(np.round(second, 7)))))
            if key in seen_bonds:
                continue
            seen_bonds.add(key)
            if np.linalg.norm(second - first) > 1.81:
                raise AssertionError("A projected bond is longer than its physical bond")
            midpoint = (first + second) / 2
            inside = 0 <= midpoint[0] <= a and 0 <= midpoint[1] <= b
            segments.append((np.array([first, second]), inside))
    display, seen_atoms = [], set()
    for chain, translate in chain_images:
        for atom in atoms:
            if chains[atom.index] != chain:
                continue
            xy = (atom.position + wraps[atom.index] @ cell + translate)[:2]
            if np.any(xy < low) or np.any(xy > high):
                continue
            key = (atom.symbol, *np.round(xy, 7))
            if key in seen_atoms:
                continue
            seen_atoms.add(key)
            inside = -1e-8 <= xy[0] <= a + 1e-8 and -1e-8 <= xy[1] <= b + 1e-8
            display.append((xy, atom.symbol, inside))
    return segments, display


def save_figure(fig, path, dpi):
    fig.savefig(path.with_suffix(".png"), dpi=dpi, bbox_inches="tight", pad_inches=0.18)
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.18)
    plt.close(fig)
    return [str(path.with_suffix(".png")), str(path.with_suffix(".svg"))]


def plot_initial(candidates, out, dpi):
    ncols = 4
    nrows = int(np.ceil(len(candidates) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(15.8, 3.12 * nrows + 1.55),
                             squeeze=False)
    margin = 1.65
    a, b, c = np.diag(candidates[0]["atoms"].cell)
    records = []
    for index, (ax, candidate) in enumerate(zip(axes.flat, candidates)):
        segments, atoms = projected_geometry(candidate, margin)
        for inside in (False, True):
            lines = [s for s, status in segments if status == inside]
            if lines:
                ax.add_collection(LineCollection(
                    lines, colors=CARBON, linewidths=1.05 if inside else 0.85,
                    alpha=0.67 if inside else 0.22, zorder=1))
            for element, color, size in [("H", HYDROGEN, 23), ("C", CARBON, 47)]:
                positions = np.array([xy for xy, symbol, status in atoms
                                      if status == inside and symbol == element])
                if len(positions):
                    ax.scatter(positions[:, 0], positions[:, 1], s=size,
                               c=color, edgecolors="white", linewidths=0.45,
                               alpha=1.0 if inside else 0.30, zorder=3 if element == "C" else 2)
        ax.add_patch(Rectangle((0, 0), a, b, fill=False, linewidth=1.15,
                               edgecolor=CELL_EDGE, linestyle=(0, (4, 3)), zorder=0))
        ax.set_xlim(-margin, a + margin)
        ax.set_ylim(-margin, b + margin)
        ax.set_aspect("equal")
        ax.set_xticks([0, 2, 4, 6])
        ax.set_yticks([0, 2, 4])
        ax.tick_params(labelsize=9, length=3, colors="#4c555e")
        ax.set_title(f"{candidate['identifier']}  |  $\\phi_0$ = {candidate['angle']:g}°",
                     loc="left", pad=7)
        if index // ncols == nrows - 1:
            ax.set_xlabel("a direction (Å)", labelpad=3)
        if index % ncols == 0:
            ax.set_ylabel("b direction (Å)", labelpad=3)
        for side in ("left", "bottom"):
            ax.spines[side].set_color("#c4c9cd")
            ax.spines[side].set_linewidth(0.7)
        records.append(dict(candidate_id=candidate["identifier"],
                            initial_setting_angle_from_a_deg=candidate["angle"],
                            projected_bonds=len(segments), displayed_projected_atoms=len(atoms)))
    for ax in list(axes.flat)[len(candidates):]:
        ax.set_visible(False)
    fig.suptitle("Constructed PE candidates: view along the chain axis c",
                 fontsize=19, x=0.045, ha="left", y=0.985)
    fig.text(0.045, 0.951,
             f"{len(candidates)} initial candidates · same cell {a:g} × {b:g} × {c:g} Å · "
             "Pnam · C4H8 per cell", fontsize=12, color="#4c555e")
    handles = [
        Line2D([], [], marker="o", color="none", markerfacecolor=CARBON, markersize=7, label="C"),
        Line2D([], [], marker="o", color="none", markerfacecolor=HYDROGEN, markersize=6, label="H"),
        Line2D([], [], color=CELL_EDGE, linestyle="--", linewidth=1.2, label="Periodic cell"),
        Line2D([], [], marker="o", color="none", markerfacecolor=CARBON, alpha=0.30,
               markersize=7, label="Neighbor-cell images"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.04, 0.932),
               frameon=False, ncol=4, fontsize=11, handlelength=1.6, columnspacing=2.3)
    fig.text(0.045, 0.022,
             "$\\phi_0$: constructed chain-plane angle from +a, not the literature angle from b.  "
             "Projection hides z separation; marker sizes are schematic.",
             fontsize=10, color="#4c555e")
    fig.subplots_adjust(left=0.055, right=0.985, bottom=0.065, top=0.88,
                        hspace=0.32, wspace=0.25)
    return save_figure(fig, out / "initial_candidates_ab", dpi), records


def plot_optimization(summary_path, out, dpi):
    if not summary_path.exists():
        return dict(status="not_available", reason="optimization_summary.json does not yet exist")
    summary = load_json(summary_path)
    all_rows = summary.get("candidates", [])
    fields = ["initial_angle_deg", "final_chain_plane_angle_deg",
              "relative_energy_meV_CH2", "cluster_id"]
    rows = [r for r in all_rows if r.get("geometry_passed") is True
            and all(r.get(k) is not None for k in fields)]
    if not rows:
        return dict(status="no_accepted_results_yet", source_path=str(summary_path),
                    source_sha256=sha(summary_path), total_candidate_count=len(all_rows))
    groups = sorted(set(row["cluster_id"] for row in rows))
    styles = {group: dict(color=CLUSTER_COLORS[i % len(CLUSTER_COLORS)],
                          marker=MARKERS[(i // len(CLUSTER_COLORS)) % len(MARKERS)])
              for i, group in enumerate(groups)}
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.7))
    xs = [float(row["initial_angle_deg"]) for row in all_rows]
    x_max = max(xs) if xs else 90
    for group in groups:
        selected = sorted((r for r in rows if r["cluster_id"] == group),
                          key=lambda r: r["initial_angle_deg"])
        x = [r["initial_angle_deg"] for r in selected]
        angle = [r["final_chain_plane_angle_deg"] for r in selected]
        energy = [r["relative_energy_meV_CH2"] for r in selected]
        label = f"{group}  (n={len(selected)})"
        for ax, ys in zip(axes, [angle, energy]):
            ax.scatter(x, ys, s=59, edgecolors="white", linewidths=0.8,
                       zorder=3, label=label, **styles[group])
    axes[0].set_ylabel("Final chain-plane angle from a (°)")
    axes[0].set_ylim(-3, 93)
    axes[0].set_yticks([0, 15, 30, 45, 60, 75, 90])
    axes[0].plot([0, min(x_max, 90)], [0, min(x_max, 90)],
                 linestyle="--", color="#a7afb5", linewidth=0.9, zorder=0)
    axes[0].text(0.97, 0.96, "Dashed: unchanged angle", transform=axes[0].transAxes,
                 ha="right", va="top", fontsize=9, color="#6d757c")
    axes[1].set_ylabel("Relative energy (meV / CH₂)")
    axes[1].axhline(0, color="#a7afb5", linewidth=0.8, zorder=0)
    # A single zero-energy snapshot should not acquire an artificial narrow,
    # symmetric +/- range that suggests measured energy resolution.
    energy_top = max(1.0, max(float(r["relative_energy_meV_CH2"]) for r in rows) * 1.12)
    axes[1].set_ylim(-0.04 * energy_top, energy_top)
    for ax in axes:
        ax.set_xlabel("Initial setting angle from a (°)")
        ax.set_xlim(-3, max(x_max, 15) + 3)
        ax.set_xticks(np.arange(0, max(x_max, 15) + 1, 15))
        ax.grid(axis="y", color="#e4e7e9", linewidth=0.6, zorder=0)
        ax.spines["left"].set_color("#a7afb5")
        ax.spines["bottom"].set_color("#a7afb5")
    fig.suptitle("Saved GFN2 preoptimization results at the same fixed cell",
                 x=0.06, ha="left", fontsize=17, y=0.97)
    fig.text(0.06, 0.91,
             f"Accepted plotted results: {len(rows)}/{len(all_rows)} · "
             f"geometry/energy groups: {len(groups)} · summary status: {summary.get('status', 'unspecified')}",
             color="#4c555e", fontsize=11)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.52, 0.89),
               ncol=min(len(groups), 4), frameon=False, fontsize=10)
    fig.text(0.06, 0.025,
             "Angles are folded into 0–90°; energy zero is the lowest accepted result in this snapshot.\n"
             "Pending or rejected candidates are omitted. Groups are geometric/energy clusters, not experimentally confirmed phases.",
             fontsize=9.5, color="#4c555e", linespacing=1.4)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.19, top=0.76, wspace=0.32)
    outputs = save_figure(fig, out / "optimization_angle_and_energy", dpi)
    return dict(status="plotted", source_path=str(summary_path), source_sha256=sha(summary_path),
                plotted_count=len(rows), total_candidate_count=len(all_rows),
                clusters=groups, outputs=outputs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "plots")
    parser.add_argument("--optimization-summary", type=Path)
    parser.add_argument("--dpi", type=int, default=170)
    args = parser.parse_args()
    root, out = args.root.resolve(), args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    setup_style()
    candidates = load_candidates(root)
    outputs, records = plot_initial(candidates, out, args.dpi)
    optimization = plot_optimization(
        args.optimization_summary or root / "optimization_summary.json", out, args.dpi)
    provenance = dict(
        purpose="User-facing scientific schematic and saved-result diagnostic; no new calculations",
        builder_path=str(Path(__file__).resolve()), builder_sha256=sha(__file__),
        projection="Cartesian a-b plane, looking along c; equal scale for every candidate",
        angle_definition="Initial independent-chain construction angle measured from +a; not literature setting angle from b",
        periodic_drawing="Topology edge shifts plus neighboring xy images; duplicated z projections collapsed",
        image_alpha_policy="Atoms outside the displayed primary cell are translucent",
        marker_radii="Schematic, not van der Waals radii", initial_outputs=outputs,
        candidates=records, source_files=[s for candidate in candidates for s in candidate["sources"]],
        optimization_plot=optimization,
    )
    (out / "plot_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(dict(initial_candidates=len(candidates), outputs=outputs,
                          optimization_plot_status=optimization["status"]), indent=2))


if __name__ == "__main__":
    main()
