"""Audit frozen-geometry GFN2/DFT pairs and report scoped calibration results.

This program starts no electronic calculations. It accepts only completed,
hash-consistent points at the frozen coordinates. Initial-angle empirical sigma
and the four related endpoints are separate domains. No energy window is invented.
Run after a batch, or during it to obtain an honest partial summary.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
METHODS = ('PBE-D3(BJ)', 'optB88-vdW', 'r2SCAN-D4')
MODELS = ('identity', 'mean_gap_bias_only', 'affine_gap_model')
POSITION_TOL_A = 1e-8
CELL_TOL_A = 1e-10


def path_on_host(value):
    """Read recorded Windows or WSL paths without changing their target."""
    text = str(value)
    if os.name == 'nt':
        match = re.match(r'^/mnt/([a-zA-Z])/(.*)$', text)
        if match:
            return Path(match[1].upper() + ':/' + match[2])
    else:
        match = re.match(r'^([a-zA-Z]):[\\/](.*)$', text)
        if match:
            return Path('/mnt/' + match[1].lower() + '/' + match[2].replace('\\', '/'))
    return Path(text)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def dump(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f'.{os.getpid()}.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2,
                                    allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def method_name(value):
    return str(value).replace('r²SCAN', 'r2SCAN')


def finite_scalar(value):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError('Non-finite scalar')
    return value


def source(path):
    path = Path(path)
    return {'path': str(path), 'sha256': sha(path), 'bytes': path.stat().st_size}


def geometry_comparison(expected, cell, positions, atomic_numbers=None):
    """Identity atom mapping is required; periodic wrapping alone is allowed."""
    cell, positions = np.asarray(cell, float), np.asarray(positions, float)
    reference_cell = np.asarray(expected['cell_A'], float)
    reference_positions = np.asarray(expected['positions_A'], float)
    if cell.shape != (3, 3) or positions.shape != reference_positions.shape:
        raise ValueError('Wrong cell or coordinate shape')
    if not np.isfinite(cell).all() or not np.isfinite(positions).all():
        raise ValueError('Non-finite geometry')
    # This audit is explicitly restricted to the fixed orthorhombic PE cell.
    if np.max(abs(reference_cell - np.diag(np.diag(reference_cell)))) > CELL_TOL_A:
        raise ValueError('Unexpected nonorthogonal reference cell')
    same_species = (atomic_numbers is None or
                    list(atomic_numbers) == list(expected['numbers']))
    fractional = (positions - reference_positions) @ np.linalg.inv(reference_cell)
    mic = (fractional - np.rint(fractional)) @ reference_cell
    lengths = np.linalg.norm(mic, axis=1)
    cell_error = float(np.max(abs(cell - reference_cell)))
    maximum = float(lengths.max())
    return {'passed': bool(cell_error <= CELL_TOL_A and maximum <= POSITION_TOL_A
                           and same_species),
            'cell_max_abs_difference_A': cell_error,
            'periodic_identity_map_max_displacement_A': maximum,
            'periodic_identity_map_RMS_A': float(np.sqrt(np.mean(lengths ** 2))),
            'same_element_order': bool(same_species),
            'element_order_evidence': ('explicit atomic numbers' if atomic_numbers is not None
                                       else 'source POSCAR hash and unchanged worker atom order'),
            'atom_map': list(range(len(reference_positions))),
            'position_tolerance_A': POSITION_TOL_A, 'cell_tolerance_A': CELL_TOL_A}


def load_samples(root):
    record = load(root / 'samples.json')
    samples, failures = {}, []
    for item in record['samples']:
        sample = dict(item)
        ident = sample['id']
        if ident in samples:
            raise ValueError('Duplicate sample ID: ' + ident)
        try:
            structure_path = path_on_host(sample['structure'])
            geometry_path = path_on_host(sample['geometry_path'])
            frozen = load(geometry_path)
            geometry = frozen['geometry']
            canonical = json.dumps(geometry, sort_keys=True, separators=(',', ':'),
                                   allow_nan=False).encode('utf-8')
            actual_geometry_hash = hashlib.sha256(canonical).hexdigest()
            if sha(structure_path) != sample['structure_sha256']:
                raise ValueError('Frozen POSCAR hash differs')
            if actual_geometry_hash != sample['geometry_sha256']:
                raise ValueError('Frozen geometry canonical hash differs')
            if frozen['geometry_sha256'] != actual_geometry_hash:
                raise ValueError('geometry.json stored canonical hash differs')
            if sample['N_CH2'] != 4 or sample['natoms'] != 12:
                raise ValueError('Expected C4H8 and four CH2 groups')
            if geometry['numbers'] != [6] * 4 + [1] * 8 or not all(geometry['pbc']):
                raise ValueError('Unexpected composition, order or PBC')
            sample['_geometry'] = geometry
            sample['_geometry_file_sha256'] = sha(geometry_path)
            sample['_valid'] = True
        except (OSError, ValueError, KeyError, TypeError) as exc:
            sample['_valid'] = False
            failures.append({'sample_id': ident, 'error': str(exc)})
        samples[ident] = sample
    return record, samples, failures


def discover(root):
    entries = []
    for backend in ('dft', 'gfn2'):
        for input_path in sorted((root / 'runs' / backend).glob('*/input.json')):
            work = input_path.parent
            entries.append({'backend': backend, 'input_path': str(input_path),
                            'result_path': str(work / 'result.json'),
                            'process_path': str(work / 'process.json'), 'external': False})
    reused = root / 'audit/reused_points.json'
    if reused.exists():
        record = load(reused)
        for entry in record.get('points', []):
            entries.append(dict(entry, external=True))
    return entries


def audit_point(entry, samples, protocol):
    """Return a compact point only when its actual result is reusable."""
    audit = {'backend': entry.get('backend'), 'external': entry.get('external', False),
             'result_path': str(entry.get('result_path')), 'accepted': False}
    try:
        backend = entry['backend']
        if backend not in ('dft', 'gfn2'):
            raise ValueError('Unrecognized backend')
        paths = {name: path_on_host(entry[name]) for name in
                 ('input_path', 'result_path', 'process_path')}
        spec = load(paths['input_path'])
        ident = entry.get('sample_id', spec.get('sample_id'))
        label = entry.get('precision_label', spec.get('precision_label', spec.get('scan_label')))
        method = method_name(entry.get('method', spec.get('method', 'GFN2-xTB')))
        audit.update(sample_id=ident, method=method, precision_label=label)
        if not paths['result_path'].exists() or not paths['process_path'].exists():
            audit['reason'] = 'Result or final process record not yet available'
            return None, audit
        result, process = load(paths['result_path']), load(paths['process_path'])
        audit.update(result_status=result.get('status'), exit_code=process.get('exit_code'))
        if ident not in samples or not samples[ident]['_valid']:
            raise ValueError('Unknown or invalid frozen sample')
        sample = samples[ident]
        if process.get('exit_code') != 0 or result.get('status') != 'passed_single_point':
            audit['reason'] = 'Not a completed successful point'
            return None, audit
        if backend == 'dft' and result.get('SCF_converged') is not True:
            raise ValueError('DFT SCF not confirmed converged')
        if backend == 'gfn2' and result.get('scf_converged') is not True:
            raise ValueError('GFN2 SCF not confirmed converged')
        hashes = {name: sha(path) for name, path in paths.items()}
        if result.get('input_sha256') != hashes['input_path']:
            raise ValueError('Result/input hash mismatch')
        if result.get('input') != spec:
            raise ValueError('Result embeds different input settings')
        if spec.get('runner_sha256_at_submission') not in (None, result.get('runner_sha256')):
            raise ValueError('Result runner differs from submitted snapshot hash')
        if entry.get('source_result_sha256') not in (None, hashes['result_path']):
            raise ValueError('External reuse result hash mismatch')
        structure_path = path_on_host(spec['structure'])
        structure_sha = sha(structure_path)
        if structure_sha != sample['structure_sha256']:
            raise ValueError('Calculation source POSCAR differs from frozen sample')
        if result.get('source_geometry_sha256') != structure_sha:
            raise ValueError('Result source geometry hash mismatch')
        if spec.get('source_geometry_sha256_expected') not in (None, structure_sha):
            raise ValueError('Submission expected a different frozen geometry')
        if spec.get('strain') or spec.get('displacement'):
            raise ValueError('Displaced/strained points cannot enter frozen-geometry pairing')
        geom = geometry_comparison(sample['_geometry'], result['cell_A'], result['positions_A'],
                                   result.get('atomic_numbers'))
        if not geom['passed']:
            raise ValueError('Actual result coordinates differ from frozen geometry: ' + json.dumps(geom))
        if backend == 'dft':
            expected = protocol['DFT_production_precision']
            if label == 'H' and any(spec.get(key) != value for key, value in expected.items()):
                raise ValueError('DFT H precision differs from frozen production settings')
            if method not in METHODS:
                raise ValueError('Unrecognized DFT target')
            if spec.get('nbands', 16) != 16 or spec.get('xc_stencil') != 4:
                raise ValueError('DFT basis/method settings differ')
            if spec.get('scf_energy_eV_electron') != 1e-9 or spec.get('scf_density_electron_electron') != 1e-9 or spec.get('scf_eigenstates_eV2_electron') != 1e-11:
                raise ValueError('DFT SCF settings differ')
            if result.get('runner_sha256') != protocol['DFT_science_runner_sha256']:
                raise ValueError('DFT science runner hash differs')
        else:
            expected = protocol.get('GFN2_proposed_sampling', {})
            if label == 'F' and (spec.get('repeat') != expected.get('repeat') or
                                 spec.get('accuracy', .001) != expected.get('accuracy', .001)):
                raise ValueError('GFN2 F sampling/accuracy differs')
            if method != 'GFN2-xTB' or result.get('geometry_optimized') is not False:
                raise ValueError('GFN2 point is not the declared frozen singlepoint')
        total = result['total']
        energy = finite_scalar(total['energy_eV_CH2'])
        cell_energy = finite_scalar(total['energy_free_eV'] if backend == 'dft'
                                    else total.get('energy_free_eV_cell', total['energy_eV_cell']))
        if abs(energy - cell_energy / sample['N_CH2']) > 1e-9:
            raise ValueError('Cell/CH2 energy normalization inconsistent')
        forces, stress = np.asarray(total['forces_eV_A'], float), np.asarray(total['stress_kbar'], float)
        if forces.shape != (12, 3) or stress.shape != (6,) or not np.isfinite(forces).all() or not np.isfinite(stress).all():
            raise ValueError('Missing/non-finite full force or stress arrays')
        sources = {name: {'path': str(paths[name]), 'sha256': hashes[name]} for name in paths}
        sources['structure'] = {'path': str(structure_path), 'sha256': structure_sha}
        sources['geometry'] = {'path': str(path_on_host(sample['geometry_path'])),
                               'sha256': sample['_geometry_file_sha256'],
                               'canonical_geometry_sha256': sample['geometry_sha256']}
        point = {'backend': backend, 'sample_id': ident, 'method': method,
                 'precision_label': label, 'external': bool(entry.get('external')),
                 'energy_eV_CH2': energy, 'energy_free_eV_cell': cell_energy,
                 'forces_eV_A': forces.tolist(), 'stress_kbar': stress.tolist(),
                 'electronic_temperature_K': result.get('electronic_temperature_K', 0 if backend == 'dft' else spec.get('electronic_temperature_K')),
                 'geometry_check': geom, 'sources': sources,
                 'cell_A': result['cell_A'], 'positions_A': result['positions_A'],
                 'atomic_numbers': sample['_geometry']['numbers'],
                 'wall_seconds': result.get('wall_seconds'),
                 'process_wall_seconds': process.get('wall_seconds'),
                 'created_unix': float(process.get('finished_unix', result.get('start_unix', 0))),
                 'runner_sha256': result.get('runner_sha256'),
                 'precision': {key: spec.get(key) for key in ('ecut_eV', 'kmesh', 'gpts', 'repeat', 'accuracy')},
                 'settings_source': str(paths['input_path'])}
        audit.update(accepted=True, source_hashes=hashes, geometry_check=geom,
                     electronic_temperature_K=point['electronic_temperature_K'],
                     energy_eV_CH2=energy)
        return point, audit
    except (OSError, ValueError, KeyError, TypeError, np.linalg.LinAlgError) as exc:
        audit['reason'] = f'{type(exc).__name__}: {exc}'
        return None, audit


def choose_points(entries, samples, protocol):
    chosen, audit, alternatives = {}, [], []
    for entry in entries:
        point, record = audit_point(entry, samples, protocol)
        audit.append(record)
        if point is None:
            continue
        if point['backend']=='gfn2' and point['precision_label']=='F':
            target_temperature=float(protocol.get('GFN2_production_temperature_K',0.0))
            actual_temperature=point['electronic_temperature_K']
            if actual_temperature is None or not math.isclose(float(actual_temperature),target_temperature,abs_tol=1e-10):
                record.update(eligible_for_production_pairing=False,
                              reason='Valid temperature diagnostic, not the explicit GFN2 production pairing temperature',
                              pairing_target_temperature_K=target_temperature)
                continue
        record['eligible_for_production_pairing']=True
        key = (point['backend'], point['sample_id'], point['method'], point['precision_label'])
        old = chosen.get(key)
        priority = (not point['external'], point['created_unix'], point['sources']['result_path']['path'])
        if old is None:
            chosen[key] = point
        else:
            old_priority = (not old['external'], old['created_unix'], old['sources']['result_path']['path'])
            selected, alternate = (point, old) if priority > old_priority else (old, point)
            chosen[key] = selected
            alternatives.append({'key': list(key),
                                 'selected_result': selected['sources']['result_path'],
                                 'alternate_result': alternate['sources']['result_path'],
                                 'selection_rule': 'Successful local replay precedes external reuse; then latest completed successful attempt.'})
    return chosen, audit, alternatives


def vector_errors(first, second):
    force = np.asarray(first['forces_eV_A']) - np.asarray(second['forces_eV_A'])
    stress = np.asarray(first['stress_kbar']) - np.asarray(second['stress_kbar'])
    return {'force_component_RMSE_eV_A': float(np.sqrt(np.mean(force ** 2))),
            'force_atom_vector_RMSE_eV_A': float(np.sqrt(np.mean(np.sum(force ** 2, axis=1)))),
            'force_max_atom_vector_difference_eV_A': float(np.linalg.norm(force, axis=1).max()),
            'stress_component_RMSE_kbar': float(np.sqrt(np.mean(stress ** 2))),
            'stress_max_component_difference_kbar': float(np.max(abs(stress))),
            'force_difference_vectors_eV_A': force.tolist(), 'stress_difference_components_kbar': stress.tolist()}


def paired_rows(samples, points, reference_id):
    rows = []
    for sample in samples.values():
        methods = METHODS if sample['role'] == 'preoptimized_representative' else METHODS[:1]
        for method in methods:
            row = {'sample_id': sample['id'], 'role': sample['role'], 'stratum': sample['stratum'],
                   'structural_family': sample['structural_family'], 'split_block': sample['split_block'],
                   'initial_angle_deg': sample['initial_angle_deg'], 'group_id': sample.get('group_id'),
                   'independent_polymorph_claim': False, 'N_CH2': sample['N_CH2'],
                   'method': method, 'reference_id': reference_id, 'status': 'missing_points',
                   'GFN2_energy_eV_CH2': None, 'DFT_energy_eV_CH2': None,
                   'GFN2_relative_energy_meV_CH2': None, 'DFT_relative_energy_meV_CH2': None,
                   'raw_relative_energy_error_DFT_minus_GFN2_meV_CH2': None,
                   'force_component_RMSE_eV_A': None, 'force_atom_vector_RMSE_eV_A': None,
                   'force_max_atom_vector_difference_eV_A': None,
                   'stress_component_RMSE_kbar': None, 'stress_max_component_difference_kbar': None,
                   'GFN2_precision': 'F', 'DFT_precision': 'H', 'GFN2_temperature_K': None,
                   'DFT_ecut_eV': None, 'DFT_kmesh': None, 'DFT_gpts': None,
                   'pair_geometry_consistent': False, 'source_records': {}, 'missing': []}
            required = {'GFN2': ('gfn2', sample['id'], 'GFN2-xTB', 'F'),
                        'DFT': ('dft', sample['id'], method, 'H'),
                        'GFN2_reference': ('gfn2', reference_id, 'GFN2-xTB', 'F'),
                        'DFT_reference': ('dft', reference_id, method, 'H')}
            located = {}
            for name, key in required.items():
                point = points.get(key)
                if point is None:
                    row['missing'].append(name)
                else:
                    located[name] = point
                    row['source_records'][name] = point['sources']
            for prefix in ('GFN2', 'DFT'):
                if prefix in located:
                    row[prefix + '_energy_eV_CH2'] = located[prefix]['energy_eV_CH2']
            if 'GFN2' in located and 'DFT' in located:
                row.update(vector_errors(located['DFT'], located['GFN2']))
                gfn_geometry = {'cell_A':located['GFN2']['cell_A'],
                                'positions_A':located['GFN2']['positions_A'],
                                'numbers':located['GFN2']['atomic_numbers']}
                pair_geometry = geometry_comparison(gfn_geometry, located['DFT']['cell_A'],
                                                     located['DFT']['positions_A'], located['DFT']['atomic_numbers'])
                row['pair_geometry_consistent'] = pair_geometry['passed']
                row['direct_DFT_GFN2_geometry_check'] = pair_geometry
                row['geometry_checks'] = {name: point['geometry_check'] for name, point in located.items()}
                if not pair_geometry['passed']:
                    row['missing'].append('consistent_DFT_GFN2_geometry')
            if 'GFN2' in located:
                row['GFN2_temperature_K'] = located['GFN2']['electronic_temperature_K']
            if 'DFT' in located:
                for key in ('ecut_eV', 'kmesh', 'gpts'):
                    row['DFT_' + key] = located['DFT']['precision'][key]
            if not row['missing']:
                x = 1000 * (located['GFN2']['energy_eV_CH2'] - located['GFN2_reference']['energy_eV_CH2'])
                y = 1000 * (located['DFT']['energy_eV_CH2'] - located['DFT_reference']['energy_eV_CH2'])
                row.update(status='paired', GFN2_relative_energy_meV_CH2=x,
                           DFT_relative_energy_meV_CH2=y,
                           raw_relative_energy_error_DFT_minus_GFN2_meV_CH2=y-x)
            rows.append(row)
    return rows


def residual_metrics(values):
    values = np.asarray(values, float)
    if not len(values):
        return {'n': 0, 'mean_residual_meV_CH2': None, 'RMSE_meV_CH2': None,
                'sample_sd_meV_CH2': None, 'max_abs_residual_meV_CH2': None}
    return {'n': len(values), 'mean_residual_meV_CH2': float(values.mean()),
            'RMSE_meV_CH2': float(np.sqrt(np.mean(values ** 2))),
            'sample_sd_meV_CH2': float(values.std(ddof=1)) if len(values) > 1 else None,
            'max_abs_residual_meV_CH2': float(np.max(abs(values)))}


def fit_model(model, x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if not len(x):
        raise ValueError('No training samples')
    if model == 'identity':
        return {'intercept_meV_CH2': 0.0, 'slope': 1.0, 'fitted_parameters': 0}
    if model == 'mean_gap_bias_only':
        return {'intercept_meV_CH2': float(np.mean(y-x)), 'slope': 1.0, 'fitted_parameters': 1}
    if model == 'affine_gap_model':
        design = np.column_stack([np.ones(len(x)), x])
        coefficients, _, rank, singular = np.linalg.lstsq(design, y, rcond=None)
        if rank != 2 or float(np.ptp(x)) < 1e-10:
            raise ValueError('Affine training data do not identify slope and intercept')
        return {'intercept_meV_CH2': float(coefficients[0]), 'slope': float(coefficients[1]),
                'fitted_parameters': 2, 'design_condition_number': float(singular[0]/singular[-1])}
    raise ValueError(model)


def empirical_scoped_calibration(rows, samples, reference_id, gates, policy, data_integrity_passed=True):
    """Report the user-accepted finite-domain sigma without a coverage claim."""
    expected = [s['id'] for s in samples.values() if s['role'] == 'initial_grid']
    initial = [r for r in rows if r['role'] == 'initial_grid' and r['method'] == METHODS[0]
               and r['status'] == 'paired']
    ids = [r['sample_id'] for r in initial]
    accepted = policy.get('empirical_scoped_sigma_accepted') is True
    gate_passed = all(gates.get(key, {}).get('passed') is True
                      for key in ('DFT_transfer', 'GFN2_sampling', 'GFN2_temperature'))
    complete = (len(expected) == 16 and len(initial) == 16 and len(set(ids)) == 16
                and set(ids) == set(expected) and data_integrity_passed
                and all(r.get('pair_geometry_consistent') is True for r in initial))
    base = {
        'status': 'waiting_for_pairs', 'n': len(initial), 'expected_n': 16,
        'sample_ids': ids, 'expected_sample_ids': expected,
        'missing_sample_ids': [ident for ident in expected if ident not in ids],
        'method': METHODS[0], 'reference_id': reference_id, 'unit': 'meV/CH2',
        'scope': {
            'polymer': 'PE', 'composition': 'C4H8', 'N_CH2': 4,
            'fixed_Tashiro_cell_A': [7.432, 4.929, 2.545],
            'chain': 'Same all-trans chain construction as the frozen initial samples',
            'setting_angles_deg': list(range(0, 76, 5)),
            'geometry_stage': 'Initial rigid-chain geometries only',
            'applies_to_four_preoptimized_endpoints': False,
        },
        'error_definition': 'e_i = DeltaE_DFT_i - DeltaE_GFN2_i, each relative to R0',
        'bias_definition': 'Arithmetic mean of the 16 initial relative-energy errors',
        'residual_definition': 'r_i = e_i - bias',
        'sigma_definition': 'std(r_i, ddof=1)', 'ddof': 1,
        'bias_meV_CH2': None, 'sigma_meV_CH2': None, 'three_sigma_meV_CH2': None,
        'centered_residuals': [], 'empirical_margin_available': False,
        'empirical_scoped_sigma_accepted': accepted,
        'data_integrity_passed': bool(data_integrity_passed),
        'all_16_initial_pairs_complete': complete,
        'all_required_numerical_gates_passed': gate_passed,
        'numerical_gates': {key: gates.get(key, {}).get('passed') is True
                            for key in ('DFT_transfer', 'GFN2_sampling', 'GFN2_temperature')},
        'independent_validation_required_for_empirical_sigma': False,
        'blocked_CV_is_acceptance_gate': False,
        'demonstrated_DFT_savings_required': False,
        'sample_expansion_for_coverage_required': False,
        'coverage_probability_guarantee': None,
        'independent_population_coverage_validated': False,
        'interpretation': ('The bias and empirical 3-sigma margin describe this fixed-cell initial-angle '
                           'domain. They are usable within that declared scope without a population '
                           'coverage claim. Do not transfer them to the four optimized endpoints.'),
    }
    if not accepted:
        return dict(base, status='not_enabled_by_analysis_policy')
    if not complete:
        return base
    if not gate_passed:
        return dict(base, status='waiting_for_numerical_gates')
    errors = np.asarray([finite_scalar(r['raw_relative_energy_error_DFT_minus_GFN2_meV_CH2'])
                         for r in initial], dtype=float)
    bias = float(errors.mean())
    residuals = errors - bias
    sigma = float(residuals.std(ddof=1))
    return dict(base, status='empirical_sigma_ready', bias_meV_CH2=bias,
                sigma_meV_CH2=sigma, three_sigma_meV_CH2=3*sigma,
                empirical_margin_available=True,
                centered_residuals=[{'sample_id': ident, 'raw_error_meV_CH2': float(error),
                                     'centered_residual_meV_CH2': float(residual)}
                                    for ident, error, residual in zip(ids, errors, residuals)])


def blocked_calibration(rows):
    initial = [r for r in rows if r['role'] == 'initial_grid' and r['method'] == METHODS[0]
               and r['status'] == 'paired']
    blocks = sorted({r['split_block'] for r in initial})
    base = {'domain': 'fixed-cell rigid-chain initial angular grid only',
            'reference': 'Relative to the declared fixed final_group_001 geometry',
            'sample_count_available': len(initial), 'expected_sample_count': 16,
            'initial_structures_are_independent_polymorphs': False,
            'split_strategy': 'Leave one contiguous angular block out',
            'coverage_or_generalization_guarantee': None,
            'production_sigma_validated': False,
            'production_sigma_validated_definition': 'Independent predictive coverage validation; separate from the accepted scoped empirical sigma.',
            'role': 'Supplementary diagnostic; no acceptance gate for empirical_scoped_calibration.',
            'limitations': ['Adjacent angles and all 16 inputs belong to one constrained structural family.',
                           'Blocked residuals are pilot diagnostics, not independent polymorph validation.',
                           'Choosing among models using these same folds does not create an independent validation set.',
                           'Do not apply the initial-domain fit to the four near-minimum endpoints.']}
    if len(initial) != 16 or blocks != ['blockA', 'blockB', 'blockC', 'blockD'] or any(sum(r['split_block']==b for r in initial) != 4 for b in blocks):
        return dict(base, status='insufficient_initial_pair_coverage', models={}, folds_not_run=True)
    models = {}
    for model in MODELS:
        folds, pooled = [], []
        for block in blocks:
            train = [r for r in initial if r['split_block'] != block]
            test = [r for r in initial if r['split_block'] == block]
            parameters = fit_model(model, [r['GFN2_relative_energy_meV_CH2'] for r in train],
                                   [r['DFT_relative_energy_meV_CH2'] for r in train])
            heldout = []
            for row in test:
                prediction = parameters['intercept_meV_CH2'] + parameters['slope'] * row['GFN2_relative_energy_meV_CH2']
                residual = row['DFT_relative_energy_meV_CH2'] - prediction
                record = {'sample_id': row['sample_id'], 'split_block': block,
                          'x_GFN2_meV_CH2': row['GFN2_relative_energy_meV_CH2'],
                          'y_DFT_meV_CH2': row['DFT_relative_energy_meV_CH2'],
                          'prediction_DFT_meV_CH2': prediction, 'heldout_residual_meV_CH2': residual}
                heldout.append(record)
                pooled.append(residual)
            train_residuals = [r['DFT_relative_energy_meV_CH2'] - parameters['intercept_meV_CH2']
                               - parameters['slope']*r['GFN2_relative_energy_meV_CH2'] for r in train]
            folds.append({'heldout_block': block, 'train_ids': [r['sample_id'] for r in train],
                          'test_ids': [r['sample_id'] for r in test], 'train_n': len(train), 'test_n': len(test),
                          'domain': base['domain'], 'parameters_fitted_on_training_only': parameters,
                          'heldout_results': heldout, 'heldout_metrics': residual_metrics([v['heldout_residual_meV_CH2'] for v in heldout]),
                          'training_metrics_for_diagnostic_only': residual_metrics(train_residuals)})
        models[model] = {'folds': folds, 'pooled_heldout_metrics': residual_metrics(pooled),
                         'descriptive_three_times_heldout_sd_meV_CH2': 3*float(np.std(pooled, ddof=1)),
                         'this_is_not_a_validated_production_screening_margin': True}
    selected = min(MODELS, key=lambda model: models[model]['pooled_heldout_metrics']['RMSE_meV_CH2'])
    return dict(base, status='blocked_pilot_diagnostics_complete', models=models,
                provisional_best_model_by_blocked_RMSE=selected,
                model_choice_is_not_independent_validation=True)


def endpoint_diagnostics(samples, points, rows):
    final_ids = [s['id'] for s in samples.values() if s['role']=='preoptimized_representative']
    per_sample, cross_method = [], []
    for ident in final_ids:
        subset = [r for r in rows if r['sample_id']==ident]
        per_sample.append({'sample_id': ident, 'source_candidate_id': samples[ident]['source_candidate_id'],
                           'group_member_count': samples[ident].get('group_member_count'),
                           'relative_energies_by_DFT_method': {r['method']:r['DFT_relative_energy_meV_CH2'] for r in subset},
                           'GFN2_relative_energy_meV_CH2': subset[0]['GFN2_relative_energy_meV_CH2'] if subset else None,
                           'GFN2_pairwise_force_stress_errors': {r['method']:{k:r.get(k) for k in
                               ('force_component_RMSE_eV_A','force_atom_vector_RMSE_eV_A','force_max_atom_vector_difference_eV_A',
                                'stress_component_RMSE_kbar','stress_max_component_difference_kbar')} for r in subset}})
        for first, second in itertools.combinations(METHODS, 2):
            a, b = points.get(('dft',ident,first,'H')), points.get(('dft',ident,second,'H'))
            if a is not None and b is not None:
                first_row = next(r for r in subset if r['method']==first)
                second_row = next(r for r in subset if r['method']==second)
                first_gap, second_gap = first_row['DFT_relative_energy_meV_CH2'], second_row['DFT_relative_energy_meV_CH2']
                cross_method.append({'sample_id':ident,'first_method':first,'second_method':second,
                                     'relative_gap_difference_first_minus_second_meV_CH2': None if first_gap is None or second_gap is None else first_gap-second_gap,
                                     **vector_errors(a,b)})
    rankings = {}
    for method in METHODS:
        data = [r for r in rows if r['sample_id'] in final_ids and r['method']==method and r['status']=='paired']
        rankings[method] = [{'sample_id':r['sample_id'],'relative_energy_meV_CH2':r['DFT_relative_energy_meV_CH2']}
                            for r in sorted(data,key=lambda r:r['DFT_relative_energy_meV_CH2'])]
    gfn_rows = [r for r in rows if r['sample_id'] in final_ids and r['method']==METHODS[0]
                and r['GFN2_relative_energy_meV_CH2'] is not None]
    rankings['GFN2-xTB'] = [{'sample_id':r['sample_id'],'relative_energy_meV_CH2':r['GFN2_relative_energy_meV_CH2']}
                              for r in sorted(gfn_rows,key=lambda r:r['GFN2_relative_energy_meV_CH2'])]
    primary = [r for r in rows if r['sample_id'] in final_ids and r['method']==METHODS[0] and r['status']=='paired']
    errors = [r['raw_relative_energy_error_DFT_minus_GFN2_meV_CH2'] for r in primary]
    return {'domain':'Four closely related fixed-cell near-minimum representative endpoints',
            'independent_polymorph_count':None, 'per_sample':per_sample,
            'descriptive_energy_orderings':rankings, 'orderings_are_not_precision_resolved_rank_claims':True,
            'cross_DFT_method_force_stress_comparisons':cross_method,
            'sigma_final':{'validated':False,'sigma_deltaE_meV_CH2':None,
                           'descriptive_raw_residual_metrics_including_zero_reference':residual_metrics(errors),
                           'reason':'Four related endpoints, one shared reference, no independent endpoint-domain validation; do not infer a coverage guarantee.'},
            'decision':'no_rejection','retained_representatives':final_ids,'rejected_representatives':[]}


def gate_record(root, name):
    path = root/'audit'/name
    if not path.exists():
        return {'path':str(path),'available':False,'passed':False,'record':None}
    try:
        record = load(path)
        return {'path':str(path),'available':True,'passed':record.get('passed') is True,
                'sha256':sha(path),'record':record}
    except (OSError, ValueError, TypeError) as exc:
        return {'path':str(path),'available':True,'passed':False,'error':str(exc),'record':None}


def write_csv(path, rows):
    keys = list(dict.fromkeys(key for row in rows for key in row))
    temporary = path.with_name(path.name + f'.{os.getpid()}.tmp')
    with temporary.open('w',encoding='utf-8-sig',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({key:json.dumps(value,ensure_ascii=False,allow_nan=False) if isinstance(value,(list,dict)) else value
                             for key,value in row.items()})
    temporary.replace(path)


def analyze(root):
    root=Path(root).resolve()
    manifest,samples,sample_errors=load_samples(root)
    protocol=load(root/'protocol.json')
    policy_path=root/'analysis_policy.json'
    policy=load(policy_path) if policy_path.exists() else {}
    if not isinstance(policy,dict):
        raise ValueError('analysis_policy.json must contain an object')
    policy_source=dict(source(policy_path),available=True) if policy_path.exists() else {'path':str(policy_path),'available':False}
    reference_id=manifest['reference_id']
    if reference_id != protocol['reference_id'] or reference_id not in samples:
        raise ValueError('Inconsistent or absent R0 sample')
    manifest_hash=sha(root/'samples.json')
    if protocol.get('sample_manifest_sha256') != manifest_hash:
        raise ValueError('Frozen sample manifest differs from protocol')
    entries=discover(root)
    points,attempts,alternatives=choose_points(entries,samples,protocol)
    rows=paired_rows(samples,points,reference_id)
    expected_dft={(s['id'],m) for s in samples.values() for m in
                  (METHODS if s['role']=='preoptimized_representative' else METHODS[:1])}
    available_dft={(ident,method) for (backend,ident,method,label) in points
                   if backend=='dft' and label=='H' and (ident,method) in expected_dft}
    available_gfn={ident for backend,ident,method,label in points if backend=='gfn2' and method=='GFN2-xTB' and label=='F'}
    successful_pairs=sum(row['status']=='paired' for row in rows)
    gates={key:gate_record(root,name) for key,name in
           [('DFT_transfer','dft_transfer_gate.json'),('GFN2_sampling','gfn2_sampling_gate.json'),
            ('GFN2_temperature','gfn2_temperature_gate.json')]}
    temperatures=sorted({float(point['electronic_temperature_K']) for key,point in points.items()
                         if key[0]=='gfn2' and key[3]=='F' and point['electronic_temperature_K'] is not None})
    data_complete=(len(samples)==21 and len(expected_dft)==29 and available_dft==expected_dft
                   and available_gfn==set(samples) and successful_pairs==29 and not sample_errors)
    final_complete=data_complete and len(temperatures)==1 and all(g['passed'] for g in gates.values())
    empirical=empirical_scoped_calibration(rows,samples,reference_id,gates,policy,not sample_errors)
    window=policy.get('screening_window_meV_CH2',protocol.get('screening_window_meV_CH2'))
    if window is not None and finite_scalar(window)<0:
        raise ValueError('Screening energy window must be nonnegative or null')
    coverage={'expected_frozen_samples':21,'frozen_samples':len(samples),
              'expected_primary_PBE_D3_H_points':21,
              'available_primary_PBE_D3_H_points':sum(method==METHODS[0] for _,method in available_dft),
              'expected_DFT_H_points_including_cross_methods':29,'available_unique_DFT_H_points':len(available_dft),
              'expected_GFN2_F_points':21,'available_unique_GFN2_F_points':len(available_gfn),
              'expected_energy_pairs':29,'complete_energy_pairs':successful_pairs,
              'missing_DFT_H_points':[{'sample_id':i,'method':m} for i,m in sorted(expected_dft-available_dft)],
              'missing_GFN2_F_points':sorted(set(samples)-available_gfn),
              'data_coverage_complete':data_complete,'all_required_numerical_gates_passed':all(g['passed'] for g in gates.values())}
    # Check accepted immutable evidence again before publishing derived tables.
    changed=[]
    for point in points.values():
        for value in point['sources'].values():
            if sha(Path(value['path'])) != value['sha256']:
                changed.append(value['path'])
    if changed:
        raise RuntimeError('Completed evidence changed during analysis: '+repr(sorted(set(changed))))
    if policy_source['available'] and sha(policy_path)!=policy_source['sha256']:
        raise RuntimeError('Analysis policy changed during analysis')
    now=datetime.now(timezone.utc).isoformat()
    summary={'schema_version':1,'generated_utc':now,
             'status':'calibration_complete' if final_complete else 'partial',
             'final_complete':final_complete,'quantum_calculations_started_by_analyzer':False,
             'reference_id':reference_id,'energy_unit':'meV/CH2','relative_energy_definition':'1000*(E(sample)/N_CH2 - E(R0)/N_CH2); N_CH2=4',
             'raw_error_definition':'DFT relative energy minus GFN2 relative energy',
             'primary_DFT_method':METHODS[0], 'coverage':coverage,
             'frozen_sample_errors':sample_errors,'GFN2_production_electronic_temperatures_K':temperatures,
             'GFN2_analysis_target_electronic_temperature_K':float(protocol.get('GFN2_production_temperature_K',0.0)),
             'numerical_uncertainty':{'gates':gates,
                  'handling':'Read and retain reported relative-to-R0 changes; numerical uncertainty is not converted to a fitted sigma or hidden by a near-zero residual.',
                  'screening_margin_validated':False,
                  'screening_margin_validated_definition':'Independent predictive coverage; does not block the scoped empirical margin.',
                  'scoped_empirical_margin_available':empirical['empirical_margin_available']},
             'empirical_scoped_calibration':empirical,
             'initial_domain_calibration':blocked_calibration(rows),
             'preoptimized_endpoint_diagnostics':endpoint_diagnostics(samples,points,rows),
             'negative_control_pairs':[r for r in rows if r['role']=='unstable_control'],
             'screening':{'enabled':False,'automatic_candidate_rejection_enabled':False,'window_W_meV_CH2':window,
                          'decision':'awaiting_energy_window' if window is None else 'margin_reported_no_automatic_rejection',
                          'empirical_margin_available':empirical['empirical_margin_available'],
                          'three_sigma_meV_CH2':empirical['three_sigma_meV_CH2'],
                          'empirical_margin_scope':'The 16 frozen PE initial geometries only; not the four preoptimized endpoints.',
                          'rejected_ids':[],
                          'retained_representatives':[s['id'] for s in samples.values() if s['role']=='preoptimized_representative'],
                          'reason':('No energy window W is specified; report the scoped empirical margin when ready and await a window before removing candidates.'
                                    if window is None else 'The analysis reports the specified window and scoped empirical margin; no automatic deletion is performed.'),
                          'coverage_probability_guarantee':None},
             'point_attempt_audit':attempts,'duplicate_point_selection':alternatives,
             'selected_points':list(points.values()),
             'sources':{'samples':source(root/'samples.json'),'protocol':source(root/'protocol.json'),
                        'analyzer':source(Path(__file__)),'analysis_policy':policy_source}}
    out=root/'tables'
    out.mkdir(parents=True,exist_ok=True)
    dump(out/'paired_energy_results.json',{'schema_version':1,'generated_utc':now,'reference_id':reference_id,
                                          'energy_unit':'meV/CH2','coverage':coverage,'rows':rows})
    write_csv(out/'paired_energy_results.csv',rows)
    dump(out/'calibration_summary.json',summary)
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT)
    parser.add_argument('--allow-incomplete',action='store_true',
                        help='Write an honest partial progress snapshot and allow exit 0; default fails closed until data and gates are complete.')
    args=parser.parse_args()
    summary=analyze(args.root)
    print(json.dumps({'status':summary['status'],'final_complete':summary['final_complete'],'coverage':summary['coverage'],
                      'screening_decision':summary['screening']['decision']},ensure_ascii=False))
    if not summary['final_complete'] and not args.allow_incomplete:
        raise SystemExit(1)


if __name__=='__main__':
    main()
