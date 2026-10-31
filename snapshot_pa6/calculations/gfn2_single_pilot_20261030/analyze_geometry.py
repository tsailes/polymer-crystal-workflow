"""Independent before/after PA6 geometry diagnostics, without energy calculations.

Atom identity/order and the fixed cell must be unchanged. Relaxation coordinates
remain unwrapped, so the original declared covalent image shifts are retained.
The comparison is a post-run diagnostic; no H-bond property filters candidates.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from ase.io import read
import spglib

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
sys.path.insert(0, str(PROJECT / 'scripts'))
import validate_ordered as validator


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def stats(values):
    values = np.asarray(values, float)
    return {'count': int(values.size),
            'min_A': float(values.min()) if values.size else None,
            'max_A': float(values.max()) if values.size else None,
            'mean_A': float(values.mean()) if values.size else None}


def graph_key(edges):
    return {validator.canonical_edge(e['i'], e['j'], e['shift']) for e in edges}


def symmetry_details(frac, symbols, cell):
    reports = []
    for rotation, translation in validator.OPS:
        transformed = frac @ rotation.T + translation
        assignment = validator.assignment_residual(transformed, frac, symbols, symbols, cell)
        mapping = np.asarray(assignment['mapping'], int)
        residual, images = validator.minimum_images(frac[mapping] - transformed, cell)
        reports.append({'rotation': rotation.tolist(), 'translation': translation.tolist(),
                        **assignment, 'residual_per_atom_A': np.linalg.norm(residual, axis=1).tolist(),
                        'minimum_image_per_atom': images.tolist()})
    return {'operations': reports,
            'max_residual_A': max(op['max_A'] for op in reports),
            'interpretation': 'Species-preserving global bijection per declared operation; no chain-wise alignment'}


def snapshot(atoms, record):
    cell = np.asarray(atoms.cell.array)
    positions = np.asarray(atoms.positions)
    symbols = np.asarray(atoms.get_chemical_symbols())
    frac = positions @ np.linalg.inv(cell)
    edges, adjacency, errors, _ = validator.infer_graph(symbols, frac, cell)
    topology, membership = validator.periodic_components(adjacency)
    chemistry = validator.check_chemistry(symbols, adjacency)
    hb = validator.hydrogen_bonds(symbols, frac, cell, adjacency, membership)
    hb.update(analysis_role='Post-optimization diagnostic only; never used for seed admission, selection, ranking, or convergence',
              used_for_filtering=False)
    roles = {item['index']: item['role'] for item in record['atom_roles']}
    lengths, grouped = [], defaultdict(list)
    for edge in record['bonds']:
        i, j, image = edge['i'], edge['j'], np.asarray(edge['shift'])
        vector = positions[j] + image @ cell - positions[i]
        length = float(np.linalg.norm(vector))
        pair = ''.join(sorted((symbols[i], symbols[j])))
        role_pair = {roles[i], roles[j]}
        category = ('amide_C_N' if role_pair == {'carbonyl_C', 'amide_N'} else
                    'carbonyl_C_O' if role_pair == {'carbonyl_C', 'carbonyl_O'} else
                    'amide_N_H' if role_pair == {'amide_N', 'amide_H'} else pair)
        grouped[category].append(length)
        limits = validator.BOND_LIMITS.get(pair)
        lengths.append({'i': i, 'j': j, 'symbols': [str(symbols[i]), str(symbols[j])],
                        'roles': [roles[i], roles[j]], 'shift': image.tolist(),
                        'category': category, 'length_A': length,
                        'within_broad_covalent_bounds': bool(limits[0] <= length <= limits[1]) if limits else None})
    symmetry_scan = []
    for tolerance in (1e-5, 1e-4, .001, .003, .01, .03, .05):
        data = spglib.get_symmetry_dataset((cell, frac, atoms.numbers), symprec=tolerance)
        symmetry_scan.append({'symprec_A': tolerance,
                              'space_group_number': int(data.number) if data is not None else None,
                              'international': str(data.international) if data is not None else None,
                              'hall_number': int(data.hall_number) if data is not None else None})
    return {'formula': dict(Counter(symbols)), 'atom_count': len(atoms),
            'cell_parameters_A_deg': atoms.cell.cellpar().tolist(),
            'volume_A3': float(atoms.get_volume()), 'declared_bond_lengths': lengths,
            'bond_length_ranges': {name: stats(values) for name, values in grouped.items()},
            'inferred_covalent_edges': edges, 'inferred_bond_errors': errors,
            'chemistry': chemistry, 'periodic_topology': topology,
            'valence_geometry': validator.valence_geometry(symbols, frac, cell, adjacency),
            'declared_symmetry': symmetry_details(frac, symbols, cell),
            'spglib_scan': symmetry_scan, 'hydrogen_bonds': hb,
            'hydrogen_bond_sheet_topology': validator.hydrogen_bond_sheet_topology(adjacency, hb)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--final', type=Path, default=HERE / 'relaxation/final.extxyz')
    parser.add_argument('--output', type=Path, default=HERE / 'geometry_comparison.json')
    args = parser.parse_args()
    record_path, initial_path = HERE / 'input/candidate.json', HERE / 'input/candidate.extxyz'
    if not args.final.is_file():
        raise SystemExit(f'Final geometry is not available yet: {args.final}')
    record = json.loads(record_path.read_text(encoding='utf-8-sig'))
    initial, final = read(initial_path), read(args.final)
    for label, atoms in [('initial', initial), ('final', final)]:
        if atoms.get_chemical_symbols() != record['symbols']:
            raise ValueError(f'{label}: atom order or composition changed')
        if not np.allclose(atoms.cell.array, record['cell'], rtol=0, atol=1e-8):
            raise ValueError(f'{label}: fixed-cell comparison requires unchanged cell')
        if not np.all(atoms.pbc) or not np.isfinite(atoms.positions).all():
            raise ValueError(f'{label}: invalid periodic geometry')
    if not np.allclose(initial.positions, record['positions'], atol=1e-7, rtol=0):
        raise ValueError('Initial extxyz does not match source candidate JSON')
    before, after = snapshot(initial, record), snapshot(final, record)
    cell = initial.cell.array
    raw_delta = final.positions - initial.positions
    mic_delta, mic_images = validator.minimum_images(raw_delta @ np.linalg.inv(cell), cell)
    distances = np.linalg.norm(mic_delta, axis=1)
    declared = graph_key(record['bonds'])
    before_graph = graph_key(before['inferred_covalent_edges'])
    after_graph = graph_key(after['inferred_covalent_edges'])
    graph_preserved = declared == before_graph == after_graph
    topology_preserved = bool(graph_preserved and after['chemistry']['passed'] and
                              after['periodic_topology']['passed'])
    bond_changes = [{k: a[k] for k in ('i', 'j', 'category', 'shift')} |
                    {'initial_A': a['length_A'], 'final_A': b['length_A'],
                     'change_A': b['length_A']-a['length_A']} for a, b in
                    zip(before['declared_bond_lengths'], after['declared_bond_lengths'])]
    hb_keys = lambda snapshot: {(b['N'], b['H'], b['O'], *b['oxygen_image_from_H'])
                               for b in snapshot['hydrogen_bonds']['bonds']}
    hb_before, hb_after = hb_keys(before), hb_keys(after)
    result = {
        'status': 'analysis_completed', 'candidate_id': record['candidate_id'],
        'scope': 'One fixed-cell GFN2 pilot geometry; not energetic, numerical sampling, DFT, experimental, or global-minimum validation',
        'coordinate_convention': 'Fixed atom identity; optimization positions unwrapped; original covalent image shifts retained',
        'sources': {str(path): sha(path) for path in [record_path, initial_path, args.final,
                                                    Path(__file__), Path(validator.__file__)]},
        'software': {'validator_version': validator.VERSION, 'spglib_version': spglib.__version__},
        'hydrogen_bonds_used_for_filtering': False,
        'displacement': {'definition': 'sqrt(mean of squared per-atom minimum-image displacement norms), no rotation, translation, permutation, or chain alignment',
                         'rms_A': float(np.sqrt(np.mean(distances**2))), 'max_A': float(distances.max()),
                         'max_atom_index': int(distances.argmax()),
                         'per_atom_norm_A': distances.tolist(), 'minimum_image_vector_A': mic_delta.tolist(),
                         'minimum_image_shift': mic_images.tolist(),
                         'raw_unwrapped_rms_A': float(np.sqrt(np.mean(np.sum(raw_delta**2, axis=1)))),
                         'mean_unwrapped_translation_A': raw_delta.mean(axis=0).tolist()},
        'covalent_topology_comparison': {
            'preserved': topology_preserved, 'same_declared_and_inferred_periodic_edges': graph_preserved,
            'initial_matches_declared': declared == before_graph,
            'final_matches_declared': declared == after_graph,
            'lost_edges': sorted(before_graph - after_graph), 'new_edges': sorted(after_graph - before_graph),
            'max_absolute_declared_bond_length_change_A': max(abs(b['change_A']) for b in bond_changes)},
        'declared_bond_length_changes': bond_changes,
        'hydrogen_bond_changes': {'initial_bond_count': len(hb_before), 'final_bond_count': len(hb_after),
                                  'lost_geometric_contacts': sorted(hb_before - hb_after),
                                  'new_geometric_contacts': sorted(hb_after - hb_before),
                                  'role': 'Threshold-based post-run geometry diagnostic; not measured hydrogen bonds or an acceptance criterion'},
        'before': before, 'after': after,
    }
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps({'output': str(args.output), 'displacement_rms_A': result['displacement']['rms_A'],
                      'displacement_max_A': result['displacement']['max_A'],
                      'covalent_topology_preserved': topology_preserved,
                      'hydrogen_bond_counts': [len(hb_before), len(hb_after)],
                      'symmetry_max_residual_A': [before['declared_symmetry']['max_residual_A'], after['declared_symmetry']['max_residual_A']]}, indent=2))


if __name__ == '__main__':
    main()
