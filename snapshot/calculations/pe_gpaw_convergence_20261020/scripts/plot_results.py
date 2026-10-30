"""Generate source-linked scientific figures from archived numerical results."""
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports" / "figures"
OUT.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 220})


def save(fig, name, sources):
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{suffix}", bbox_inches="tight")
    plt.close(fig)
    (OUT / f"{name}.sources.json").write_text(json.dumps([
        {"file": str(p.relative_to(ROOT)),
         "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources
    ], indent=2) + "\n")


def diagnostic_plot():
    sources = [ROOT / "validation_audit" / filename for filename in (
        "derivative_results.csv", "derivative_results_deriv_dense.csv",
        "derivative_results_deriv_stencil4.csv")]
    if not all(p.exists() for p in sources):
        return
    values = []
    for p in sources:
        with p.open(encoding="utf-8-sig", newline="") as f:
            rows = list(csv.DictReader(f))
        matches = [float(r["absolute_error"]) for r in rows
                   if r["method"] == "optB88-vdW" and r["kind"] == "stress"
                   and r["component"] == "yz" and r["contribution"] == "total"
                   and r["absolute_error"]]
        if not matches:
            return
        values.append(max(matches))
    fig, ax = plt.subplots(figsize=(7.1, 4.0))
    ax.bar(range(3), values, color=["#d9834a", "#4778a8", "#33856b"], width=.6)
    ax.axhline(.1, color="#aa3333", linestyle="--", lw=1.5,
               label="Acceptance threshold: 0.1 kbar")
    ax.set_xticks(range(3), ["42 x 28 x 15\nstencil 2",
                            "64 x 40 x 24\nstencil 2",
                            "64 x 40 x 24\nstencil 4"])
    ax.set_ylabel("Maximum analytic - finite-difference error (kbar)")
    ax.set_title("PE optB88-vdW: diagnosing the yz stress")
    for i, v in enumerate(values):
        ax.text(i, v + .009, f"{v:.4f}", ha="center", va="bottom")
    ax.set_ylim(0, max(values) * 1.2)
    ax.legend(frameon=False)
    save(fig, "stress_grid_diagnostic", sources)


def convergence_plots():
    p = ROOT / "tables" / "convergence_summary.json"
    if not p.exists():
        return
    data = json.loads(p.read_text())
    for g in data.get("groups", []):
        if g["method"] != "PBE" or g["stage"] not in (
                "numerical_cutoff", "numerical_kmesh", "numerical_grid"):
            continue
        rows = g["points"]
        if len(rows) < 2:
            continue
        stage = g["stage"].replace("numerical_", "")
        xs = np.arange(len(rows))
        if stage == "cutoff":
            labels = [f"{r['ecut_eV']:g}" for r in rows]
            xlabel = "Plane-wave cutoff (eV)"
        elif stage == "kmesh":
            labels = [" x ".join(map(str, r["kmesh"])) for r in rows]
            xlabel = "Gamma-centered k-point mesh"
        else:
            labels = [" x ".join(map(str, r["coarse_grid"])) for r in rows]
            xlabel = "Coarse real-space grid"
        metrics = [("delta_energy_abs_meV_atom", "Energy (meV/atom)", .5),
                   ("max_atom_force_vector_difference_eV_A",
                    "Max. force-vector\n difference (eV/A)", .005),
                   ("max_stress_component_difference_kbar",
                    "Max. stress-component\n difference (kbar)", .1)]
        fig, axes = plt.subplots(3, 1, figsize=(7.5, 8.4), sharex=True)
        for ax, (key, ylabel, threshold) in zip(axes, metrics):
            vals = [r["against_terminal"][key] for r in rows]
            ax.plot(xs, vals, "o-", color="#326f91", lw=1.5, ms=4)
            ax.axhline(threshold, color="#b2553a", linestyle="--", lw=1)
            ax.set_ylabel(ylabel)
            ax.grid(axis="y", alpha=.2)
        axes[-1].set_xticks(xs, labels, rotation=20 if len(rows) > 5 else 0)
        axes[-1].set_xlabel(xlabel)
        first = rows[0]
        fixed = (
            f"k={'x'.join(map(str, first['kmesh']))}; "
            f"grid={'x'.join(map(str, first['coarse_grid']))}"
            if stage == "cutoff" else
            f"cutoff={first['ecut_eV']:g} eV; "
            f"grid={'x'.join(map(str, first['coarse_grid']))}"
            if stage == "kmesh" else
            f"cutoff={first['ecut_eV']:g} eV; "
            f"k={'x'.join(map(str, first['kmesh']))}"
        )
        axes[0].set_title(
            f"PE PBE: {stage} scan, differences from highest tested setting\n"
            f"{fixed}\nStatus: {g['status']} (fixed geometry)", fontsize=10)
        fig.tight_layout()
        save(fig, f"convergence_{stage}_{g['group_id'][:8]}", [p])


def relaxation_plots():
    """Plot saved evaluations; a low energy alone never implies acceptance."""
    for folder in sorted((ROOT / "relaxations").glob("*")):
        source = folder / "evaluations.jsonl"
        if not source.exists():
            continue
        rows = []
        for line in source.read_text().splitlines():
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                # A running writer may not have completed its final line.
                continue
            if row.get("SCF_converged"):
                rows.append(row)
        if len(rows) < 2:
            continue
        spec_path = folder / "input.json"
        spec = json.loads(spec_path.read_text())
        status_path = folder / "status.json"
        status = json.loads(status_path.read_text())["status"]
        xs = np.array([r["evaluation_id"] for r in rows])
        energies = np.array([r["total"]["energy_free_eV"] for r in rows])
        cell = np.array([r["cell_lengths_A"] for r in rows])
        fig, axes = plt.subplots(4, 1, figsize=(8, 10.3), sharex=True)
        axes[0].plot(xs, (energies - energies[0]) * 1000 / rows[0]["natoms"],
                     "o-", color="#326f91", ms=3)
        axes[0].set_ylabel("Energy change\n(meV/atom)")
        axes[1].semilogy(xs, [r["raw_max_force_norm_eV_A"] for r in rows],
                         "o-", color="#33856b", ms=3)
        axes[1].axhline(spec["fmax_eV_A"], color="#b2553a", ls="--",
                        label=f"Threshold {spec['fmax_eV_A']:g} eV/A")
        axes[1].set_ylabel("Max. raw force\n(eV/A)")
        axes[1].legend(frameon=False, fontsize=9)
        axes[2].semilogy(xs, [r["raw_max_stress_component_kbar"] for r in rows],
                         "o-", color="#77599c", ms=3)
        axes[2].axhline(spec["stress_max_kbar"], color="#b2553a", ls="--",
                        label=f"Threshold {spec['stress_max_kbar']:g} kbar")
        axes[2].set_ylabel("Max. raw stress\n(kbar)")
        axes[2].legend(frameon=False, fontsize=9)
        for i, label in enumerate("abc"):
            axes[3].plot(xs, (cell[:, i] / cell[0, i] - 1) * 100,
                         "o-", ms=3, label=label)
        axes[3].set_ylabel("Cell length change\nfrom initial (%)")
        axes[3].set_xlabel("Saved geometry evaluation")
        axes[3].legend(frameon=False, ncol=3)
        for ax in axes:
            ax.grid(axis="y", alpha=.2)
        axes[0].set_title(
            f"PE {spec['method']}: {folder.name}\n"
            f"{spec['ecut_eV']:g} eV; k={'x'.join(map(str, spec['kmesh']))}; "
            f"grid={'x'.join(map(str, spec['gpts']))}\nStatus: {status}",
            fontsize=10)
        fig.tight_layout()
        save(fig, f"relaxation_{folder.name}", [source, spec_path, status_path])


if __name__ == "__main__":
    diagnostic_plot()
    convergence_plots()
    relaxation_plots()
