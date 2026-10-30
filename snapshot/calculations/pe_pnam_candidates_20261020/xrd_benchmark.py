#!/usr/bin/env python3
"""Compare a cell/target space group with Tashiro's 45 observed d positions.

This is a conditional lattice/absence check, not a powder-profile R factor or
a packing score. No observed/calculated intensities, ADP or reference atom
coordinates are loaded. Candidate coordinates are used only to check whether
the requested space-group operations remain valid.

Example (run in the GFN2 environment; no calculator is invoked):
  python xrd_benchmark.py --candidates-dir candidates \
      --space-group Pnam --output-dir xrd_benchmark_initial
  python xrd_benchmark.py --structures path/to/final.extxyz \
      --space-group Pnam --output-dir xrd_benchmark_final

Dependencies: numpy, scipy, ASE, Gemmi.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path

import ase
from ase.cell import Cell
from ase.geometry import find_mic
from ase.io import read
import gemmi
import numpy as np
from scipy.optimize import linear_sum_assignment


ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parents[1]
SOURCE_DIR = PROJECT / "reports/pe_xrd_seed_search_20260922/source_audit"
SOURCE_PDF = PROJECT / "examples/PE_Tashiro_1997/sources/Tashiro_1997_original_Wiley.pdf"
EXPECTED = {
    "tashiro1997_dobs_only.json": "8bbd18402e44ba9c4bfca364411917a4289a075d161157294177b693333ac4b3",
    "tashiro1997_dobs_only.csv": "1fb6db866d63a81227825b902c3fd4c9f08ede54fe75ff5524cbc816dcbb8086",
    "source_pdf": "4750738ed0c2e29a569d3b402d1a94d46db920f3f5f40f4ba953b0d937ad15f0",
}
LIMITATIONS = [
    "Only published, processed fiber-diffraction positions are compared; this is not a continuous powder profile.",
    "The target space group and cell are supplied priors, not inferred from these observations.",
    "Position matches contain no intensity, ADP, atomic-scattering-factor or refined-reference-coordinate information.",
    "At identical cell and target group, all packings have identical allowed lattice positions and position residuals.",
    "The reported Q tolerance is a descriptive threshold, not an experimental standard deviation or statistical confidence interval.",
    "Matching is observed-to-nearest-allowed-position, permits many-to-one assignments, and does not prove indexing.",
    "No penalty is assigned to unobserved predicted positions: the published list is not known to be complete.",
    "No full-pattern R factor, chi-square, probability, or packing ranking is produced.",
    "Extra structure-specific intensity cancellations and experimental detection thresholds are not modeled.",
]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def archive_reference(destination):
    """Copy the two clean inputs and record original-PDF hash; no PDF is copied."""
    destination.mkdir(parents=True, exist_ok=True)
    sources = []
    for name in ("tashiro1997_dobs_only.json", "tashiro1997_dobs_only.csv"):
        source = SOURCE_DIR / name
        if sha(source) != EXPECTED[name]:
            raise ValueError(f"Audited input hash changed: {source}")
        target = destination / name
        if target.exists() and sha(target) != EXPECTED[name]:
            raise ValueError(f"Refusing to overwrite different archived input: {target}")
        if not target.exists():
            target.write_bytes(source.read_bytes())
        sources.append(dict(source_path=str(source), archive_name=name,
                            sha256=EXPECTED[name], bytes=source.stat().st_size))
    if sha(SOURCE_PDF) != EXPECTED["source_pdf"]:
        raise ValueError("Original PDF hash differs from the independently audited source")
    manifest = dict(
        purpose="Immutable clean position-only input archive; no reference atom coordinates",
        original_pdf=dict(path=str(SOURCE_PDF), sha256=EXPECTED["source_pdf"],
                          table="Table I, PDF page 11 / printed page 1687",
                          wavelength="PDF page 16 / printed page 1692",
                          intensity_limitation="Eq.11-13, PDF page 17: overlapping apparent Fobs use Fcalc-derived weights"),
        clean_files=sources, observed_count=45,
        raw_table_rows=59, deduplication="Exact printed d(obs) values",
        observations_not_included=["Fobs", "Fcalc", "Iobs", "dcalc", "author_hkl",
                                   "refined_atom_coordinates", "ADP"],
        primary_source_doi="10.1002/(SICI)1099-0488(199708)35:11<1677::AID-POLB3>3.0.CO;2-Q",
    )
    path = destination / "source_manifest.json"
    if path.exists() and json.loads(path.read_text(encoding="utf-8")) != manifest:
        raise ValueError("Existing reference source manifest differs; preserve it and choose another archive")
    if not path.exists():
        write_json(path, manifest)
    return destination / "tashiro1997_dobs_only.json"


def load_observations(path):
    # This task is intentionally tied to the independently audited clean table.
    if sha(path) != EXPECTED["tashiro1997_dobs_only.json"]:
        raise ValueError("Observed input is not the audited position-only JSON")
    data = json.loads(path.read_text(encoding="utf-8"))
    rows = data["observations"]
    allowed = {"observation_id", "d_obs_A", "q_A_inv", "two_theta_deg"}
    if len(rows) != 45 or any(set(row) != allowed for row in rows):
        raise ValueError("Observed data violates the 45-row field whitelist")
    wavelength = float(data["wavelength_A"])
    for row in rows:
        d = float(row["d_obs_A"])
        if d <= 0 or not np.isclose(row["q_A_inv"], 2 * math.pi / d, atol=1e-10):
            raise ValueError("Invalid observed d/Q conversion")
        if not np.isclose(row["two_theta_deg"], two_theta(d, wavelength), atol=1e-9):
            raise ValueError("Invalid observed wavelength/angle conversion")
    return data


def two_theta(d, wavelength):
    x = wavelength / (2 * d)
    return None if x > 1 else math.degrees(2 * math.asin(x))


def verify_cell_metric(cell, group):
    metric = cell @ cell.T
    worst = 0.0
    for op in group.operations():
        r = np.asarray(op.rot, dtype=float) / op.DEN
        worst = max(worst, float(np.max(np.abs(r.T @ metric @ r - metric))))
    if worst > 1e-6:
        raise ValueError(f"Cell metric is incompatible with target space group ({worst:g} A^2)")
    return worst


def target_symmetry_check(atoms, group, tolerance):
    """Check the supplied operations in the supplied setting, with no origin fit."""
    fractional = atoms.get_scaled_positions(wrap=True)
    symbols = np.asarray(atoms.get_chemical_symbols())
    cell = np.asarray(atoms.cell)
    errors = []
    for op in group.operations():
        transformed = np.array([op.apply_to_xyz(f.tolist()) for f in fractional])
        worst = 0.0
        for symbol in sorted(set(symbols)):
            indices = np.flatnonzero(symbols == symbol)
            delta = transformed[indices, None, :] - fractional[None, indices, :]
            _, lengths = find_mic((delta @ cell).reshape(-1, 3), cell, pbc=True)
            distances = lengths.reshape(len(indices), len(indices))
            first, second = linear_sum_assignment(distances)
            worst = max(worst, float(distances[first, second].max()))
        errors.append(dict(operation=op.triplet(), max_matching_distance_A=worst))
    maximum = max(row["max_matching_distance_A"] for row in errors)
    return dict(checked=True, contains_supplied_target_operations=maximum <= tolerance,
                tolerance_A=tolerance, maximum_matching_distance_A=maximum,
                origin_policy="Use supplied setting/origin; no origin alignment or symmetry inference",
                operations=errors)


def allowed_positions(cell, group, q_max, wavelength):
    """Enumerate bounded signed hkl, drop absences, group equal d values.

    |h_i| <= |a_i| / d_min bounds every reciprocal vector inside the Q sphere,
    including non-orthogonal cells. One Friedel representative is retained.
    The hkl lists are labels, not powder/fiber intensity multiplicities.
    """
    inverse_cell = np.linalg.inv(cell)
    bounds = np.ceil(np.linalg.norm(cell, axis=1) * q_max / (2 * math.pi)).astype(int)
    ops = group.operations()
    individual = []
    for hkl in itertools.product(*(range(-int(n), int(n) + 1) for n in bounds)):
        nonzero = next((v for v in hkl if v), None)
        if nonzero is None or nonzero < 0 or ops.is_systematically_absent(hkl):
            continue
        inverse_d = float(np.linalg.norm(inverse_cell @ np.asarray(hkl)))
        q = 2 * math.pi * inverse_d
        if q <= q_max + 1e-10:
            individual.append((q, hkl))
    individual.sort()
    positions = []
    for q, hkl in individual:
        if positions and abs(q - positions[-1]["q_A_inv"]) <= 1e-9:
            positions[-1]["hkl_labels"].append(list(hkl))
        else:
            d = 2 * math.pi / q
            positions.append(dict(position_id=f"allowed_{len(positions)+1:04d}",
                                  q_A_inv=q, d_A=d, two_theta_deg=two_theta(d, wavelength),
                                  hkl_labels=[list(hkl)]))
    if not positions:
        raise ValueError("No allowed reciprocal-lattice positions were generated")
    return positions


def compare_positions(observations, predicted, tolerance):
    qs = np.array([p["q_A_inv"] for p in predicted])
    matches = []
    for row in observations:
        near = predicted[int(np.argmin(np.abs(qs - row["q_A_inv"])))]
        dq = near["q_A_inv"] - row["q_A_inv"]
        matches.append(dict(
            **row, allowed_position_id=near["position_id"],
            nearest_allowed_hkl_labels=near["hkl_labels"],
            nearest_allowed_d_A=near["d_A"], nearest_allowed_q_A_inv=near["q_A_inv"],
            nearest_allowed_two_theta_deg=near["two_theta_deg"],
            delta_q_pred_minus_obs_A_inv=dq,
            abs_delta_q_A_inv=abs(dq),
            delta_d_pred_minus_obs_A=near["d_A"] - row["d_obs_A"],
            delta_two_theta_pred_minus_obs_deg=near["two_theta_deg"] - row["two_theta_deg"],
            within_descriptive_q_tolerance=abs(dq) <= tolerance,
        ))
    dq = np.array([row["delta_q_pred_minus_obs_A_inv"] for row in matches])
    used = [row["allowed_position_id"] for row in matches]
    repeated = {p: used.count(p) for p in sorted(set(used)) if used.count(p) > 1}
    stats = dict(
        observed_count=len(matches), q_tolerance_A_inv=tolerance,
        tolerance_interpretation="Descriptive proximity threshold; no uncertainty/confidence claim",
        observations_within_tolerance=sum(row["within_descriptive_q_tolerance"] for row in matches),
        fraction_observations_within_tolerance=float(np.mean(np.abs(dq) <= tolerance)),
        mean_abs_delta_q_A_inv=float(np.mean(np.abs(dq))),
        rms_delta_q_A_inv=float(np.sqrt(np.mean(dq * dq))),
        max_abs_delta_q_A_inv=float(np.max(np.abs(dq))),
        unique_predicted_positions_selected=len(set(used)),
        predicted_positions_reused_by_multiple_observations=repeated,
        assignment_policy="Independent nearest Q; not one-to-one indexing",
        unmatched_predicted_position_penalty=None,
        intensity_or_full_profile_R_factor=None,
        eligible_for_packing_ranking=False,
    )
    return matches, stats


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v) if isinstance(v, (list, dict)) else v
                             for k, v in row.items()})


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--structures", nargs="+", type=Path)
    parser.add_argument("--candidates-dir", type=Path,
                        help="Read direct cand*/input.extxyz children only")
    parser.add_argument("--cell", nargs=6, type=float, metavar=("A", "B", "C", "ALPHA", "BETA", "GAMMA"),
                        help="Optional cell-only input, or assert all structure cells match")
    parser.add_argument("--space-group", default="Pnam")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--reference-dir", type=Path, default=ROOT / "reference_xrd")
    parser.add_argument("--observed-json", type=Path,
                        help="Must hash-match the audited clean input; default archived copy")
    parser.add_argument("--q-tolerance", type=float, default=0.03,
                        help="Descriptive proximity in inverse Angstrom, not experimental sigma")
    parser.add_argument("--q-margin", type=float, default=0.15,
                        help="Enumerate beyond largest observed Q to avoid boundary clipping")
    parser.add_argument("--symmetry-tolerance", type=float, default=1e-3)
    args = parser.parse_args()
    if min(args.q_tolerance, args.q_margin, args.symmetry_tolerance) <= 0:
        parser.error("Tolerances and Q margin must be positive")
    if args.q_margin < args.q_tolerance:
        parser.error("Q margin must be at least the descriptive match tolerance")
    group = gemmi.find_spacegroup_by_name(args.space_group)
    if group is None:
        parser.error("Unknown space-group symbol/setting")
    archived = archive_reference(args.reference_dir.resolve())
    observed_path = (args.observed_json or archived).resolve()
    observed = load_observations(observed_path)
    structures = list(args.structures or [])
    if args.candidates_dir:
        structures += sorted(args.candidates_dir.glob("cand*/input.extxyz"))
    structures = list(dict.fromkeys(p.resolve() for p in structures))
    if not structures and args.cell is None:
        parser.error("Supply --structures, nonempty --candidates-dir, or --cell")
    asserted_cell = np.asarray(Cell.fromcellpar(args.cell)) if args.cell else None
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    # Preserve previous runs; no silent refresh of frozen reports.
    if (output / "summary.json").exists():
        raise FileExistsError("Output already has summary.json; select a new output directory")
    rows, cache, identifiers = [], {}, set()
    q_max = max(row["q_A_inv"] for row in observed["observations"]) + args.q_margin
    for path in structures or [None]:
        if path is None:
            cell, identifier = asserted_cell, "cell_only"
            symmetry = dict(checked=False, reason="No candidate coordinates supplied")
            structure = None
        else:
            atoms = read(path)
            if not np.all(atoms.pbc) or abs(atoms.cell.volume) <= 0:
                raise ValueError(f"Candidate requires a finite 3D periodic cell: {path}")
            cell = np.asarray(atoms.cell)
            if asserted_cell is not None and not np.allclose(cell, asserted_cell, atol=1e-7, rtol=0):
                raise ValueError(f"Candidate cell differs from explicit --cell: {path}")
            identifier = f"{path.parent.name}__{path.stem}"
            symmetry = target_symmetry_check(atoms, group, args.symmetry_tolerance)
            structure = dict(path=str(path), sha256=sha(path), formula=atoms.get_chemical_formula(),
                             atom_count=len(atoms),
                             coordinate_usage="Target-operation closure only; coordinates do not enter peak positions")
        if identifier in identifiers:
            raise ValueError(f"Duplicate output identifier {identifier}; run paths separately")
        identifiers.add(identifier)
        metric_error = verify_cell_metric(cell, group)
        signature_payload = dict(cell_A=np.round(cell, 10).tolist(), hall=group.hall)
        signature = hashlib.sha256(json.dumps(signature_payload, sort_keys=True).encode()).hexdigest()
        if signature not in cache:
            predicted = allowed_positions(cell, group, q_max, observed["wavelength_A"])
            matches, stats = compare_positions(observed["observations"], predicted, args.q_tolerance)
            cache[signature] = (predicted, matches, stats)
        predicted, matches, stats = cache[signature]
        folder = output / identifier
        folder.mkdir(exist_ok=True)
        write_csv(folder / "allowed_positions.csv", predicted)
        write_csv(folder / "observed_nearest_allowed.csv", matches)
        result = dict(
            candidate_id=identifier, candidate=structure, cell_A=cell.tolist(),
            target_space_group=dict(name=group.xhm(), number=group.number, hall=group.hall,
                                    role="Supplied target, not inferred"),
            metric_invariance_max_error_A2=metric_error, target_symmetry_check=symmetry,
            observed_input=dict(path=str(observed_path), sha256=sha(observed_path),
                                count=45, wavelength_A=observed["wavelength_A"]),
            enumeration=dict(q_min_A_inv=0, q_max_A_inv=q_max,
                             distinct_allowed_positions=len(predicted),
                             hkl_policy="Signed integer indices; one Friedel representative; equal-Q labels grouped; no intensity multiplicities",
                             q_grouping_tolerance_A_inv=1e-9),
            cell_group_signature=signature, descriptive_statistics=stats,
            limitations=LIMITATIONS,
        )
        write_json(folder / "result.json", result)
        rows.append(dict(
            candidate_id=identifier, cell_group_signature=signature,
            target_symmetry_pass=symmetry.get("contains_supplied_target_operations"),
            max_target_symmetry_error_A=symmetry.get("maximum_matching_distance_A"),
            **{k: stats[k] for k in ["observations_within_tolerance", "mean_abs_delta_q_A_inv",
                                     "rms_delta_q_A_inv", "max_abs_delta_q_A_inv"]},
            result_path=str(folder / "result.json"),
        ))
    summary = dict(
        status="position_only_diagnostic_completed", candidate_count=len(rows),
        unique_cell_group_position_tables=len(cache),
        all_candidates_share_cell_and_target_group=len(cache) == 1,
        packing_ranking_produced=False, energy_calculations_run=0,
        target_space_group=group.xhm(), observed_count=45, q_tolerance_A_inv=args.q_tolerance,
        builder_sha256=sha(__file__), versions=dict(ase=ase.__version__, gemmi=gemmi.__version__),
        candidates=rows, limitations=LIMITATIONS,
    )
    write_json(output / "summary.json", summary)
    write_csv(output / "summary.csv", rows)
    text = [
        "# Tashiro PE: position-only lattice consistency check", "",
        f"Compared {len(rows)} structures/cells with 45 published observed d spacings.",
        f"Target group: {group.xhm()} (Hall {group.hall}).",
        f"Distinct cell/group tables: {len(cache)}. Q proximity threshold: {args.q_tolerance:g} inverse Angstrom.",
        "", "This output does not rank molecular packings or report a powder-profile R factor.",
        "At fixed cell and target group, orientation/translation changes leave these allowed positions unchanged.",
        "", "Each candidate folder contains allowed_positions.csv, observed_nearest_allowed.csv and result.json.",
        "Signed residuals are predicted minus observed. Nearest-Q matches can reuse a predicted position.",
        "target_symmetry_check verifies the supplied operations in the supplied origin/axis setting.",
        "If that check fails, the allowed-position table is conditional on the requested target group and is not the candidate's verified extinction pattern.",
        "", "Limitations:", *[f"- {item}" for item in LIMITATIONS], "",
    ]
    (output / "README.md").write_text("\n".join(text), encoding="utf-8")
    print(json.dumps({k: summary[k] for k in ["status", "candidate_count",
          "unique_cell_group_position_tables", "packing_ranking_produced",
          "energy_calculations_run"]}, indent=2))


if __name__ == "__main__":
    main()
