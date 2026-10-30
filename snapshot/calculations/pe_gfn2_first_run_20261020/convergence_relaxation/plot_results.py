"""Plot recorded PE calculations; this program never evaluates an energy.

Required inputs: --run RUN_DIRECTORY --final RELAXED_EXTXYZ.
The initial model is ../model/pe_initial.extxyz relative to this script.
Reads native scalar result.json files and optimization history.json files;
summary/decision.json is optional and is not needed to produce the figures.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from ase.io import read
from ase.neighborlist import neighbor_list


HERE = Path(__file__).resolve().parent
INITIAL = HERE.parent / "model" / "pe_initial.extxyz"
COLORS = ["#2875AB", "#DA7C30", "#4D9B70", "#9A6AB0", "#B85D6E"]


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def collect_records(run):
    """Identify records by their schema, not by assumed suite directory names."""
    samples, histories = [], []
    for path in sorted(run.rglob("result.json")):
        data = load_json(path)
        if all(k in data for k in ("repeats", "energy_per_CH2_eV", "forces_eV_per_A")):
            if data.get("status") == "completed":
                samples.append((path, data))
    for path in sorted(run.rglob("history.json")):
        data = load_json(path)
        if not isinstance(data, list) or not data:
            continue
        if not all(all(k in row for k in ("step", "energy_eV_per_CH2", "max_force_eV_A")) for row in data):
            continue
        final_path = path.parent / "final_result.json"
        final = load_json(final_path) if final_path.exists() else None
        histories.append((path, data, final))
    return samples, histories


def save_figure(fig, output, stem, figures):
    for extension in ("png", "svg"):
        destination = output / f"{stem}.{extension}"
        fig.savefig(destination, dpi=220, facecolor="white")
        figures.append(str(destination))
    plt.close(fig)


def plot_histories(histories, run, target, output, figures):
    if not histories:
        return {"generated": False, "reason": "No optimization history is available yet."}
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.6), layout="constrained")
    for ordinal, (path, data, final) in enumerate(histories):
        steps = np.array([row["step"] for row in data])
        energy = np.array([row["energy_eV_per_CH2"] for row in data], float)
        fmax = np.array([row["max_force_eV_A"] for row in data], float)
        if not np.all(np.isfinite(energy)) or not np.all(np.isfinite(fmax)) or np.any(fmax < 0):
            raise ValueError(f"Invalid recorded optimization values: {path}")
        label = str(path.parent.relative_to(run))
        if final is not None:
            label += " (force converged)" if final.get("converged") else " (not force converged)"
        else:
            label += " (in progress)"
        color = COLORS[ordinal % len(COLORS)]
        axes[0].plot(steps, (energy - energy[0]) * 1000, color=color, lw=1.8, label=label)
        axes[1].plot(steps, fmax, color=color, lw=1.8, label=label)
        for axis, values in zip(axes, ((energy - energy[0]) * 1000, fmax)):
            axis.scatter(steps[-1], values[-1], color=color, s=24, zorder=4)
    axes[0].axhline(0, lw=0.8, color="#9AA3AA")
    axes[0].set(ylabel="Energy change from each stage start (meV / CH$_2$)")
    axes[1].set_yscale("symlog", linthresh=max(target / 10, 1e-6))
    axes[1].axhline(target, color="#656565", ls="--", lw=1.1,
                    label=f"Force target: {target:g} eV / Angstrom")
    axes[1].set(ylabel="Maximum atomic force (eV / Angstrom)")
    for axis in axes:
        axis.set_xlabel("BFGS step")
        axis.grid(alpha=0.2)
        axis.legend(fontsize=8, loc="best")
    fig.suptitle("GFN2-xTB local relaxation at a fixed cell", fontsize=13, weight="bold")
    save_figure(fig, output, "fixed_cell_optimization", figures)
    return {"generated": True, "energy_reference": "Each optimization stage's own first recorded energy.",
            "force_target_eV_per_A": target, "history_count": len(histories)}


def draw_xy(axis, atoms, title):
    positions = atoms.get_positions()
    cell = atoms.cell.array
    chain_ids = np.asarray(atoms.arrays["chain_id"])
    symbols = np.array(atoms.get_chemical_symbols())
    a, b = cell[0, 0], cell[1, 1]
    # Draw real periodic bonds, including segments crossing cell boundaries.
    ii, jj, shifts, distances = neighbor_list("ijSd", atoms, 1.85)
    seen = set()
    for i, j, shift, distance in zip(ii, jj, shifts, distances):
        pair = "".join(sorted([symbols[i], symbols[j]]))
        cutoff = {"CC": 1.80, "CH": 1.30}.get(pair, 0)
        if distance >= cutoff:
            continue
        start, end = positions[i, :2], (positions[j] + shift @ cell)[:2]
        projected_key = tuple(sorted((tuple(np.round(start, 7)), tuple(np.round(end, 7)))))
        if projected_key in seen:
            continue
        seen.add(projected_key)
        color = COLORS[int(chain_ids[i]) % len(COLORS)]
        axis.plot([start[0], end[0]], [start[1], end[1]], color=color,
                  alpha=0.7, lw=1.5, zorder=1, clip_on=True)
    for chain in sorted(set(chain_ids)):
        color = COLORS[int(chain) % len(COLORS)]
        for element, size in (("H", 36), ("C", 125)):
            mask = (chain_ids == chain) & (symbols == element)
            axis.scatter(positions[mask, 0], positions[mask, 1], s=size,
                         facecolors=color if element == "C" else "white",
                         edgecolors=color, linewidths=1.2, zorder=3)
    axis.plot([0, a, a, 0, 0], [0, 0, b, b, 0], color="#505965", lw=1.0)
    axis.set(xlim=(-0.08, a + 0.08), ylim=(-0.08, b + 0.08),
             xlabel="x (Angstrom)", ylabel="y (Angstrom)", title=title)
    axis.set_aspect("equal", adjustable="box")
    axis.grid(alpha=0.13)


def plot_structures(initial_path, final_path, output, figures):
    initial, final = read(initial_path), read(final_path)
    if not (np.all(initial.pbc) and np.all(final.pbc)):
        raise ValueError("Both models must have three periodic dimensions")
    if not np.allclose(initial.cell.array, final.cell.array, atol=1e-10, rtol=0):
        raise ValueError("Initial and relaxed cells differ; a same-cell comparison would be misleading")
    if not np.allclose(initial.cell.array, np.diag(initial.cell.lengths()), atol=1e-10, rtol=0):
        raise ValueError("This XY comparison expects the declared orthorhombic axis-aligned cell")
    if initial.get_chemical_symbols() != final.get_chemical_symbols():
        raise ValueError("Initial and relaxed atom compositions/orders differ")
    for atoms in (initial, final):
        if "chain_id" not in atoms.arrays or set(atoms.arrays["chain_id"]) != {0, 1}:
            raise ValueError("Two explicit chain_id groups are required")
        atoms.wrap()
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 5.1), layout="constrained")
    draw_xy(axes[0], initial, "Initial candidate: analytic chains")
    draw_xy(axes[1], final, "Relaxed candidate: GFN2-xTB")
    a, b, c = initial.cell.lengths()
    fig.suptitle(f"Same fixed cell: {a:.3f} x {b:.3f} x {c:.3f} Angstrom | View along c",
                 fontsize=12, weight="bold")
    handles = [Line2D([], [], marker="o", color=COLORS[i], linestyle="none", label=f"Chain {i + 1}")
               for i in (0, 1)]
    handles += [Line2D([], [], marker="o", color="#555555", linestyle="none", markersize=8, label="C"),
                Line2D([], [], marker="o", markerfacecolor="white", markeredgecolor="#555555",
                       color="white", linestyle="none", markersize=5, label="H")]
    axes[1].legend(handles=handles, loc="upper right", fontsize=8, ncol=2)
    fig.supxlabel("Projected model coordinates; overlapping z positions are not resolved. Experimental packing is not shown.",
                  fontsize=9)
    save_figure(fig, output, "initial_relaxed_xy", figures)
    return {"same_cell_verified": True, "cell_lengths_A": [float(a), float(b), float(c)],
            "scope": "XY projection of initial and relaxed candidates, not experimental structure validation."}


def plot_sampling(samples, run, output, figures):
    groups = defaultdict(list)
    for path, record in samples:
        # Compare only identical recorded input bytes and SCF settings.
        key = (str(path.parent.parent.relative_to(run)), record.get("input_sha256"),
               record.get("settings", {}).get("accuracy"))
        groups[key].append((path, record))
    if not groups:
        return []
    fig, axes = plt.subplots(len(groups), 2, figsize=(11.5, max(3.1, 2.8 * len(groups))),
                             squeeze=False, layout="constrained")
    references = []
    for row, (key, records) in enumerate(sorted(groups.items(), key=lambda item: str(item[0]))):
        records.sort(key=lambda item: (np.prod(item[1]["repeats"]), tuple(item[1]["repeats"])))
        ref_path, ref = records[-1]
        ref_forces = np.asarray(ref["forces_eV_per_A"], float)
        labels, energy, force = [], [], []
        for path, record in records:
            labels.append("x".join(map(str, record["repeats"])))
            energy.append((record["energy_per_CH2_eV"] - ref["energy_per_CH2_eV"]) * 1000)
            force.append(float(np.linalg.norm(np.asarray(record["forces_eV_per_A"]) - ref_forces, axis=1).max()))
        x = np.arange(len(labels))
        for col, values in enumerate((energy, force)):
            axes[row, col].plot(x, values, "o-", color=COLORS[row % len(COLORS)], lw=1.5, ms=4)
            axes[row, col].set_xticks(x, labels, rotation=35, ha="right")
            axes[row, col].grid(alpha=0.18)
            axes[row, col].set_xlabel("Exact replica grid (Gamma calculation)")
            axes[row, col].set_title(f"{key[0]} | reference {labels[-1]}", fontsize=10)
        axes[row, 0].set_ylabel("Energy difference (meV / CH$_2$)")
        axes[row, 1].set_ylabel("Max. force-vector difference\n(eV / Angstrom)")
        references.append({"suite": key[0], "reference_path": str(ref_path),
                           "reference_selection": "Largest sampled atom count; tuple order breaks ties. "
                                                  "This choice is not a convergence verdict."})
    fig.suptitle("Recorded replica-sampling diagnostics", fontsize=13, weight="bold")
    save_figure(fig, output, "sampling_diagnostics", figures)
    return references


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--final", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    run, final = args.run.resolve(), args.final.resolve()
    if not run.is_dir() or not final.is_file() or not INITIAL.is_file():
        raise FileNotFoundError("Run directory, --final structure, and initial model must exist")
    output = (args.output_dir or run / "plots").resolve()
    targets = [output / f"{stem}.{ext}" for stem in
               ("fixed_cell_optimization", "initial_relaxed_xy", "sampling_diagnostics") for ext in ("png", "svg")]
    targets.append(output / "plot_manifest.json")
    if any(path.exists() for path in targets):
        raise FileExistsError("Plot artifacts already exist; choose a fresh --output-dir")
    output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.fonttype": "none"})
    samples, histories = collect_records(run)
    protocol_path = run / "protocol.json"
    protocol = load_json(protocol_path) if protocol_path.exists() else {}
    force_target = float(protocol.get("thresholds", {}).get("relaxation_fmax_eV_per_A", 0.01))
    figures = []
    structure_info = plot_structures(INITIAL, final, output, figures)
    optimization_info = plot_histories(histories, run, force_target, output, figures)
    references = plot_sampling(samples, run, output, figures)
    sources = [INITIAL, final, *[path for path, _ in samples], *[path for path, _, _ in histories]]
    sources.extend(path.parent / "final_result.json" for path, _, result in histories if result is not None)
    if protocol_path.exists():
        sources.append(protocol_path)
    manifest = {"script": str(Path(__file__).resolve()), "script_sha256": sha256(__file__),
                "run": str(run), "figures": figures, "structure_comparison": structure_info,
                "optimization": optimization_info, "sampling_references": references,
                "source_sha256": {str(path): sha256(path) for path in sources},
                "scope": "Diagnostic plots from recorded data; no new chemical calculations, "
                         "global search conclusion or experimental-structure claim."}
    (output / "plot_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"output_dir": str(output), "figures": figures,
                      "optimization_history_count": len(histories), "sampling_record_count": len(samples)}, indent=2))


if __name__ == "__main__":
    main()
