"""Read-only numerical evaluation of the declared GFN2 calibration gates.

No quantum calculation or geometry mutation is performed. Missing, failed, or
unverifiable points cannot pass. Output scope is the explicitly tested frozen
geometries, not every PE packing. Old failures remain recorded as attempts.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = 'final_group_001'
TEMPERATURE_SAMPLES = [REFERENCE, 'initial_cand_0001', 'initial_cand_0121']
SAMPLING_SAMPLES = [f'final_group_{i:03d}' for i in range(1, 5)] + ['initial_cand_0001', 'initial_cand_0121']
REPEATS = {'F': [2, 3, 8], 'D': [3, 4, 10]}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def resolve_path(value):
    """Accept paths recorded on either shared-workspace host."""
    text = str(value).replace('\\', '/')
    if os.name == 'nt' and text.startswith('/mnt/') and len(text) > 7:
        text = text[5].upper() + ':' + text[6:]
    elif os.name != 'nt' and re.match(r'^[A-Za-z]:/', text):
        text = '/mnt/' + text[0].lower() + text[2:]
    path = Path(text)
    return path if path.is_absolute() else ROOT / path


def dump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f'.{os.getpid()}.tmp')
    temporary.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def finite(values):
    if isinstance(values, (int, float)):
        return math.isfinite(values)
    if isinstance(values, list):
        return all(finite(value) for value in values)
    return False


def flatten(values):
    if isinstance(values, list):
        for item in values:
            yield from flatten(item)
    else:
        yield float(values)


def max_component_difference(first, second):
    a, b = list(flatten(first)), list(flatten(second))
    if len(a) != len(b):
        raise ValueError('Array shape mismatch')
    return max(abs(x - y) for x, y in zip(a, b))


def max_vector_difference(first, second):
    if len(first) != len(second) or any(len(row) != 3 for row in first + second):
        raise ValueError('Force vector shape mismatch')
    return max(math.sqrt(sum((x-y)**2 for x, y in zip(a, b))) for a, b in zip(first, second))


def geometry_difference(first, second):
    r1, r2 = first['result'], second['result']
    return {
        'atomic_numbers_equal': r1['atomic_numbers'] == r2['atomic_numbers'],
        'max_cell_component_difference_A': max_component_difference(r1['cell_A'], r2['cell_A']),
        'max_position_component_difference_A': max_component_difference(r1['positions_A'], r2['positions_A']),
        'same_frozen_structure_hash': first['frozen_structure_sha256'] == second['frozen_structure_sha256'],
    }


def same_geometry(first, second):
    data = geometry_difference(first, second)
    return data, (data['atomic_numbers_equal'] and data['same_frozen_structure_hash']
                  and data['max_cell_component_difference_A'] <= 1e-12
                  and data['max_position_component_difference_A'] <= 1e-12)


def evidence(point):
    """Compact provenance; preserve the full result at its hashed source path."""
    return {key: point[key] for key in (
        'origin', 'sample_id', 'precision_label', 'temperature_K',
        'result_path', 'result_sha256', 'input_path', 'input_sha256',
        'process_path', 'process_sha256', 'frozen_structure_path', 'frozen_structure_sha256',
        'runner_sha256', 'energy_eV_CH2', 'replica_force_max_deviation_eV_A')}


def verify_point(result_path, input_path, process_path, sample_id, precision, samples, origin,
                 expected_result_hash=None):
    result_path, input_path, process_path = map(resolve_path, (result_path, input_path, process_path))
    result, spec, process = load(result_path), load(input_path), load(process_path)
    if result.get('status') != 'passed_single_point' or process.get('exit_code') != 0:
        raise ValueError('Point does not have passed_single_point plus exit code 0')
    if sample_id not in samples or precision not in REPEATS:
        raise ValueError('Unknown sample or sampling precision')
    if expected_result_hash and sha(result_path) != expected_result_hash:
        raise ValueError('Reused result source hash changed')
    if result.get('input_sha256') != sha(input_path) or result.get('input') != spec:
        raise ValueError('Input fingerprint or embedded input mismatch')
    if spec.get('sample_id', sample_id) != sample_id or spec.get('precision_label', precision) != precision:
        raise ValueError('Input sample/precision label mismatch')
    if result.get('scf_converged') is not True or result.get('cold_start') is not True or result.get('geometry_optimized') is not False:
        raise ValueError('Cold singlepoint/SCC metadata not confirmed')
    temperature = float(spec['electronic_temperature_K'])
    if not math.isfinite(temperature) or float(result['electronic_temperature_K']) != temperature:
        raise ValueError('Temperature metadata mismatch')
    if spec.get('repeat') != REPEATS[precision] or result.get('repeat') != REPEATS[precision]:
        raise ValueError('Repeated supercell does not match F/D precision')
    settings = result['settings']
    if (settings.get('method') != 'GFN2-xTB' or settings.get('charge') != 0
            or settings.get('multiplicity') != 1 or settings.get('accuracy') != .001
            or settings.get('electronic_temperature') != temperature):
        raise ValueError('GFN2 physics/accuracy settings differ from the declared gate')
    frozen = resolve_path(samples[sample_id]['structure'])
    expected_geometry_hash = samples[sample_id]['structure_sha256']
    if sha(frozen) != expected_geometry_hash:
        raise ValueError('Frozen sample changed after the manifest')
    source = resolve_path(spec['structure'])
    if sha(source) != expected_geometry_hash or result.get('source_geometry_sha256') != expected_geometry_hash:
        raise ValueError('Point source is not the frozen sample, including reuse')
    if spec.get('source_geometry_sha256_expected', expected_geometry_hash) != expected_geometry_hash:
        raise ValueError('Submission geometry hash mismatch')
    if spec.get('runner_sha256_at_submission', result.get('runner_sha256')) != result.get('runner_sha256'):
        raise ValueError('Runner fingerprint mismatch')
    snapshot = result_path.parent / 'runner_snapshot.py'
    if snapshot.exists() and sha(snapshot) != result['runner_sha256']:
        raise ValueError('Saved runner snapshot hash mismatch')
    total = result['total']
    if result.get('primitive_natoms') != 12 or result.get('n_CH2_per_primitive_cell') != 4:
        raise ValueError('Primitive-cell normalization mismatch')
    if result.get('replica_count') != math.prod(REPEATS[precision]) or result.get('supercell_natoms') != 12*math.prod(REPEATS[precision]):
        raise ValueError('Replica count/atom count mismatch')
    if (len(result['atomic_numbers']) != 12 or sorted(result['atomic_numbers']) != [1]*8+[6]*4
            or len(result['positions_A']) != 12 or len(result['cell_A']) != 3
            or any(len(row) != 3 for row in result['positions_A'] + result['cell_A'])
            or len(total['forces_eV_A']) != 12 or any(len(row) != 3 for row in total['forces_eV_A'])
            or len(total['stress_kbar']) != 6):
        raise ValueError('Atomic/geometry/force/stress shape mismatch')
    arrays = [result['positions_A'], result['cell_A'], total['energy_eV_CH2'], total['energy_eV_per_cell'],
              total['forces_eV_A'], total['stress_kbar'], result['replica_force_max_deviation_eV_A']]
    if not all(finite(array) for array in arrays):
        raise ValueError('Non-finite point values')
    if abs(total['energy_eV_per_cell']/4-total['energy_eV_CH2']) > 1e-10:
        raise ValueError('Energy per CH2 and per-cell disagree')
    replicas = result['replica_count']
    if abs(result['supercell']['energy_eV']/replicas-total['energy_eV_per_cell']) > 1e-9:
        raise ValueError('Supercell normalization disagrees')
    return dict(origin=origin, sample_id=sample_id, precision_label=precision,
                temperature_K=temperature, result=result, input=spec,
                result_path=str(result_path), result_sha256=sha(result_path),
                input_path=str(input_path), input_sha256=sha(input_path),
                process_path=str(process_path), process_sha256=sha(process_path),
                frozen_structure_path=str(frozen), frozen_structure_sha256=expected_geometry_hash,
                runner_sha256=result['runner_sha256'], energy_eV_CH2=float(total['energy_eV_CH2']),
                replica_force_max_deviation_eV_A=float(result['replica_force_max_deviation_eV_A']),
                finished_unix=float(process.get('finished_unix', result.get('start_unix', 0))))


def collect_points(samples):
    accepted, rejected = {}, []
    for work in sorted((ROOT/'runs/gfn2').glob('*')):
        if not (work/'input.json').is_file():
            continue
        try:
            spec = load(work/'input.json')
            sample_id = spec.get('sample_id')
            precision = spec.get('precision_label')
            point = verify_point(work/'result.json', work/'input.json', work/'process.json',
                                 sample_id, precision, samples, 'current_run')
            accepted.setdefault((sample_id, precision, point['temperature_K']), []).append(point)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            rejected.append({'path': str(work), 'error': f'{type(exc).__name__}: {exc}'})
    reused_path = ROOT/'audit/reused_points.json'
    if reused_path.is_file():
        reuse = load(reused_path)
        for item in reuse.get('points', []):
            if item.get('backend') != 'gfn2':
                continue
            try:
                point = verify_point(item['result_path'], item['input_path'], item['process_path'],
                                     item['sample_id'], item['precision_label'], samples, 'explicit_reused_point',
                                     expected_result_hash=item['source_result_sha256'])
                accepted.setdefault((point['sample_id'], point['precision_label'], point['temperature_K']), []).append(point)
            except (OSError, ValueError, KeyError, TypeError) as exc:
                rejected.append({'reused_item': item, 'error': f'{type(exc).__name__}: {exc}'})
    selected, duplicate_records = {}, []
    for key, values in accepted.items():
        values.sort(key=lambda p: (p['finished_unix'], p['result_path']))
        selected[key] = values[-1]
        if len(values) > 1:
            duplicate_records.append({'key': list(key), 'selected_result': values[-1]['result_path'],
                                      'accepted_attempts': [evidence(point) for point in values]})
    return selected, rejected, duplicate_records


def pair_or_missing(points, sample_id, first_key, second_key):
    first = points.get((sample_id, *first_key))
    second = points.get((sample_id, *second_key))
    missing = []
    if first is None:
        missing.append({'precision': first_key[0], 'temperature_K': first_key[1]})
    if second is None:
        missing.append({'precision': second_key[0], 'temperature_K': second_key[1]})
    if missing:
        return first, second, {'sample_id': sample_id, 'status': 'missing_or_unverified', 'passed': False,
                               'missing_points': missing}
    return first, second, None


def make_gate(kind, checks, thresholds, context, reference_relative=None):
    complete = all(check.get('status') != 'missing_or_unverified' for check in checks)
    passed = complete and all(check['passed'] for check in checks) and not context['errors']
    gate = dict(schema_version=1, created_utc=datetime.now(timezone.utc).isoformat(),
                kind=kind, passed=passed,
                status='passed' if passed else ('incomplete' if not complete else 'failed'),
                required_sample_count=len(checks), verified_pair_count=sum(check.get('status') != 'missing_or_unverified' for check in checks),
                thresholds=thresholds, checks=checks, provenance=context,
                quantum_calculations_started=False,
                scope_note='Only the listed frozen PE geometries were tested. This gate does not prove convergence for every candidate or a different cell/phase.')
    if reference_relative is not None:
        gate['relative_energy_sampling'] = reference_relative
    return gate


def evaluate_temperature(points, context):
    thresholds = dict(energy_abs_meV_CH2=.1, max_force_vector_eV_A=.001,
                      max_stress_component_kbar=.1, comparison='each numerical difference <= threshold')
    checks = []
    for sample_id in TEMPERATURE_SAMPLES:
        zero, warm, missing = pair_or_missing(points, sample_id, ('F', 0.0), ('F', 300.0))
        if missing:
            checks.append(missing)
            continue
        geo, geometry_ok = same_geometry(zero, warm)
        metrics = dict(energy_300_minus_0_meV_CH2=1000*(warm['energy_eV_CH2']-zero['energy_eV_CH2']),
                       max_force_vector_difference_eV_A=max_vector_difference(zero['result']['total']['forces_eV_A'], warm['result']['total']['forces_eV_A']),
                       max_stress_component_difference_kbar=max_component_difference(zero['result']['total']['stress_kbar'], warm['result']['total']['stress_kbar']))
        passes = dict(same_geometry=geometry_ok,
                      energy=abs(metrics['energy_300_minus_0_meV_CH2']) <= .1,
                      forces=metrics['max_force_vector_difference_eV_A'] <= .001,
                      stress=metrics['max_stress_component_difference_kbar'] <= .1)
        checks.append(dict(sample_id=sample_id, status='verified_pair', passed=all(passes.values()),
                           threshold_checks=passes, geometry=geo, metrics=metrics,
                           zero_K_source=evidence(zero), warm_300K_source=evidence(warm)))
    return make_gate('GFN2_electronic_temperature_0_vs_300K', checks, thresholds, context)


def evaluate_sampling(points, context):
    thresholds = dict(energy_abs_meV_CH2=1., max_force_vector_eV_A=.005,
                      replica_force_max_deviation_eV_A=1e-5,
                      comparison='strictly smaller than each threshold',
                      residual_force_requirement='not imposed; initial geometries are not optimized',
                      stress_requirement='recorded only; no new stress threshold invented')
    checks = []
    for sample_id in SAMPLING_SAMPLES:
        fine, dense, missing = pair_or_missing(points, sample_id, ('F', 0.0), ('D', 0.0))
        if missing:
            checks.append(missing)
            continue
        geo, geometry_ok = same_geometry(fine, dense)
        metrics = dict(energy_D_minus_F_meV_CH2=1000*(dense['energy_eV_CH2']-fine['energy_eV_CH2']),
                       max_force_vector_difference_eV_A=max_vector_difference(fine['result']['total']['forces_eV_A'], dense['result']['total']['forces_eV_A']),
                       max_stress_component_difference_kbar=max_component_difference(fine['result']['total']['stress_kbar'], dense['result']['total']['stress_kbar']),
                       F_replica_force_max_deviation_eV_A=fine['replica_force_max_deviation_eV_A'],
                       D_replica_force_max_deviation_eV_A=dense['replica_force_max_deviation_eV_A'])
        passes = dict(same_geometry=geometry_ok,
                      energy=abs(metrics['energy_D_minus_F_meV_CH2']) < 1.,
                      forces=metrics['max_force_vector_difference_eV_A'] < .005,
                      F_replicas=metrics['F_replica_force_max_deviation_eV_A'] < 1e-5,
                      D_replicas=metrics['D_replica_force_max_deviation_eV_A'] < 1e-5)
        checks.append(dict(sample_id=sample_id, status='verified_pair', passed=all(passes.values()),
                           threshold_checks=passes, geometry=geo, metrics=metrics,
                           fine_F_source=evidence(fine), dense_D_source=evidence(dense)))
    reference_relative = {'reference_id': REFERENCE, 'normalization': 'meV/CH2',
                          'status': 'incomplete', 'records': [], 'max_absolute_change_meV_CH2': None,
                          'acceptance_role': 'diagnostic; no additional tolerance is silently imposed'}
    ref_f = points.get((REFERENCE, 'F', 0.0))
    ref_d = points.get((REFERENCE, 'D', 0.0))
    if ref_f and ref_d:
        for sample_id in SAMPLING_SAMPLES:
            fine, dense = points.get((sample_id, 'F', 0.0)), points.get((sample_id, 'D', 0.0))
            if fine and dense:
                gap_f = 1000*(fine['energy_eV_CH2']-ref_f['energy_eV_CH2'])
                gap_d = 1000*(dense['energy_eV_CH2']-ref_d['energy_eV_CH2'])
                reference_relative['records'].append(dict(sample_id=sample_id, is_reference=sample_id == REFERENCE,
                        relative_energy_F_meV_CH2=gap_f, relative_energy_D_meV_CH2=gap_d,
                        D_minus_F_relative_energy_change_meV_CH2=gap_d-gap_f,
                        F_result_path=fine['result_path'], F_result_sha256=fine['result_sha256'],
                        D_result_path=dense['result_path'], D_result_sha256=dense['result_sha256']))
        reference_relative['reference_F_source'] = evidence(ref_f)
        reference_relative['reference_D_source'] = evidence(ref_d)
        if len(reference_relative['records']) == len(SAMPLING_SAMPLES):
            reference_relative['status'] = 'complete'
        if reference_relative['records']:
            reference_relative['max_absolute_change_meV_CH2'] = max(abs(row['D_minus_F_relative_energy_change_meV_CH2']) for row in reference_relative['records'])
    return make_gate('GFN2_Gamma_repeated_cell_F_vs_D', checks, thresholds, context, reference_relative)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--allow-incomplete', action='store_true', help='Return exit code 0 for an incomplete progress snapshot; passed remains false')
    args = parser.parse_args()
    manifest_path = ROOT/'samples.json'
    manifest = load(manifest_path)
    samples = {row['id']: row for row in manifest['samples']}
    protocol_path = ROOT/'protocol.json'
    protocol = load(protocol_path)
    errors = []
    if protocol.get('sample_manifest_sha256') != sha(manifest_path):
        errors.append('Sample manifest no longer matches frozen protocol')
    points, rejected, duplicates = collect_points(samples)
    context = dict(evaluator_path=str(Path(__file__).resolve()), evaluator_sha256=sha(__file__),
                   sample_manifest_path=str(manifest_path), sample_manifest_sha256=sha(manifest_path),
                   protocol_path=str(protocol_path), protocol_sha256=sha(protocol_path),
                   rejected_or_incomplete_attempts=rejected, accepted_duplicates=duplicates, errors=errors)
    reused = ROOT/'audit/reused_points.json'
    if reused.exists():
        context.update(reused_points_manifest_path=str(reused), reused_points_manifest_sha256=sha(reused))
    temperature = evaluate_temperature(points, context)
    sampling = evaluate_sampling(points, context)
    dump(ROOT/'audit/gfn2_temperature_gate.json', temperature)
    dump(ROOT/'audit/gfn2_sampling_gate.json', sampling)
    print(json.dumps({'temperature_gate': temperature['status'], 'temperature_passed': temperature['passed'],
                      'sampling_gate': sampling['status'], 'sampling_passed': sampling['passed'],
                      'verified_temperature_pairs': temperature['verified_pair_count'],
                      'verified_sampling_pairs': sampling['verified_pair_count']}, ensure_ascii=False))
    if not (temperature['passed'] and sampling['passed']):
        statuses = (temperature['status'], sampling['status'])
        if not (args.allow_incomplete and 'incomplete' in statuses and 'failed' not in statuses):
            raise SystemExit(1)


if __name__ == '__main__':
    main()
