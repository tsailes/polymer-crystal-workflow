"""Select one first-passing v1.2 seed; no energy or H-bond screening."""
from pathlib import Path
import hashlib
import itertools
import json
import sys
import time
import numpy as np
from ase import Atoms
from ase.geometry import cellpar_to_cell
from ase.io import read, write

ROOT = Path('/mnt/e/PA6_DFT')
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'scripts'))
import build_ordered_candidates as builder
import validate_ordered as validator

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))

def dump(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

def main():
    if (OUT / 'candidate.json').exists():
        raise RuntimeError('Existing selected input; refusing to overwrite')
    meta_path = ROOT / 'inputs/allowed_literature_inputs.json'
    manifest_path = ROOT / 'inputs/chemistry_template_manifest.json'
    selection_path = ROOT / 'branches/O/phase24_first_pass/selection.json'
    policy_path = ROOT / 'inputs/candidate_physical_screen_policy.json'
    meta, manifest, selection = map(load, (meta_path, manifest_path, selection_path))
    for item in manifest['source_hashes']:
        assert sha(ROOT / item['path']) == item['sha256'], item['path']
    assert selection['phase_count'] == 24
    assert selection['phase_fractional_c'] == (np.arange(24)/48).tolist()
    cell = cellpar_to_cell(meta['cell']['lengths_A'] + meta['cell']['angles_deg'])
    templates = [load(ROOT / f'inputs/independent_template/chain_{j:02d}.json') for j in range(2)]
    hashes = {str(path.relative_to(ROOT)): sha(path) for path in
              [meta_path, manifest_path, selection_path, policy_path,
               ROOT / 'scripts/build_ordered_candidates.py', ROOT / 'scripts/validate_ordered.py',
               ROOT / 'inputs/independent_template/chain_00.json',
               ROOT / 'inputs/independent_template/chain_01.json']}
    protocol = {
        'status': 'declared_before_generation', 'branch': 'O',
        'generator_version': builder.VERSION, 'source_hashes': hashes,
        'selector_sha256': sha(__file__),
        'selection_rule': 'First builder- and independent-validator-passing candidate in v1.2 deterministic enumeration; stop at one',
        'ordering': 'x then rotation then phase, round robin template 0 then 1',
        'chain_x_fractional': [.25, .225, .275, .20, .30],
        'chain_rotation_about_c_deg': [0., -10., 10.],
        'chain_c_phase_fractional': selection['phase_fractional_c'],
        'template_indices': [0, 1],
        'hydrogen_bonds_evaluated': False,
        'hydrogen_bonds_used_for_admission_ranking_or_sampling': False,
        'energy_evaluations': 0, 'source_coordinates_read': False,
        'coordinate_derived_geometry_read': False,
        'coverage_scope': 'one locally generated pilot input, not full candidate-library coverage',
        'existing_candidate_libraries_modified': False,
        'user_authorization': '2026-10-30 user requested one PA6 candidate GFN coarse optimization'
    }
    dump(OUT / 'selection_protocol.json', protocol)
    checked = []
    start = time.perf_counter()
    for x, phi, phase in itertools.product(protocol['chain_x_fractional'],
                                          protocol['chain_rotation_about_c_deg'],
                                          protocol['chain_c_phase_fractional']):
        for tj, template in enumerate(templates):
            assert np.allclose(template['period_vector_A'], cell[2], rtol=0, atol=1e-10)
            record = builder.expand(np.asarray(template['positions_unwrapped_A']),
                                    template['atom_roles'], template['bonds'], template['symbols'],
                                    cell, x, phi, phase)
            build_check = builder.validate_geometry(record)
            item = {'enumeration_index_1based': len(checked)+1,
                    'parameters': dict(record['generation_parameters'], template_index=tj,
                                       template_id=template['template_id']),
                    'builder_passed': build_check['passed'],
                    'builder_clashes': build_check['clashes']}
            checked.append(item)
            if not build_check['passed']:
                continue
            check = validator.validate(record, analyze_hbonds=False)
            item['independent_validator_passed'] = check['passed']
            if not check['passed']:
                item['independent_errors'] = check
                continue
            cid = 'O-v12-pilot-20261030-firstpass'
            original_k = int(round(phase*96))
            record.update(candidate_id=cid, branch_id='O', generator_version=builder.VERSION,
                          chemistry_template_id=template['template_id'], template_index=tj,
                          input_hashes=hashes,
                          initial_symmetry_hypothesis='P1121/n No.14 c_unique',
                          geometry_status='builder_and_independent_validator_passed',
                          numerical_status='not_started', experimental_status='not_compared',
                          geometry_report=build_check,
                          selection_protocol_sha256=sha(OUT / 'selection_protocol.json'),
                          coverage_status='single_local_pilot_input_not_full_library',
                          enumeration_index_1based=len(checked), original_48_phase_index=original_k)
            dump(OUT / 'candidate.json', record)
            atoms = Atoms(record['symbols'], positions=record['positions'], cell=cell, pbc=True)
            atoms.info.update(candidate_id=cid, branch_id='O', template_id=template['template_id'])
            write(OUT / 'candidate.extxyz', atoms)
            reread = read(OUT / 'candidate.extxyz')
            assert reread.get_chemical_symbols() == record['symbols']
            assert np.allclose(reread.positions, record['positions'], rtol=0, atol=1e-7)
            assert np.allclose(reread.cell.array, cell, rtol=0, atol=1e-7)
            assert all(reread.pbc)
            check.update(input_file=str(OUT / 'candidate.json'), input_sha256=sha(OUT / 'candidate.json'),
                         structure_file=str(OUT / 'candidate.extxyz'),
                         structure_sha256=sha(OUT / 'candidate.extxyz'), extxyz_matches_validated_JSON=True)
            dump(OUT / 'candidate.validation.json', check)
            dump(OUT / 'selection_record.json', {
                'status': 'passed', 'candidate_id': cid, 'candidate_count': 1,
                'parameters': item['parameters'], 'atom_count': len(atoms),
                'formula': atoms.get_chemical_formula(), 'cell_parameters': atoms.cell.cellpar().tolist(),
                'checks_in_enumeration_order': checked, 'elapsed_seconds': time.perf_counter()-start,
                'energy_evaluations': 0, 'hydrogen_bonds_evaluated': False,
                'coverage_scope': 'single local pilot input; enumeration stopped at first passing candidate',
                'input_sha256': check['input_sha256'], 'structure_sha256': check['structure_sha256']})
            print(json.dumps({'candidate_id': cid, 'parameters': item['parameters'],
                              'enumerated': len(checked), 'atom_count': len(atoms),
                              'validation_passed': check['passed'], 'output': str(OUT)}, indent=2), flush=True)
            return
    raise RuntimeError('No passing candidate in declared grid')

if __name__ == '__main__':
    main()
