"""Freeze a tighter-force continuation from four accepted DFT endpoints; no quantum work."""
from pathlib import Path
import copy
import hashlib
import json
import shutil
import numpy as np
import spglib
from ase.io import read

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT.parent / 'pe_dft_representative_refinement_20261022'
TARGET = 0.0002

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def put(p, text):
    p.parent.mkdir(parents=True, exist_ok=True)
    if p.exists():
        assert p.read_text(encoding='utf-8') == text, f'Refuse overwrite: {p}'
    else:
        p.write_text(text, encoding='utf-8')

def save(p, value):
    put(p, json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')

def freeze(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        assert sha(src) == sha(dst), f'Refuse overwrite: {dst}'
    else:
        shutil.copy2(src, dst)

def load(p):
    return json.loads(p.read_text(encoding='utf-8-sig'))

def main():
    old_protocol = load(OLD / 'protocol.json')
    old_summary = load(OLD / 'tables/refinement_summary.json')
    assert old_summary['complete'] is True and old_summary['accepted_count'] == 4
    assert sha(OLD / 'protocol.json') == '430393ee12839e0a749a13809ae9df17bb9827ed1c87a13127ddfba2fbaa824a'
    assert old_protocol['optimizer']['fmax_eV_A'] == 0.002
    for name in ['run_fixedcell_relax.py', 'run_point.py', 'run_refinement.py',
                 'launch_workflow.py', 'structure_match.py', 'keepawake_process_probe.py']:
        freeze(OLD / 'scripts' / name, ROOT / 'scripts' / name)
    analyzer = (OLD / 'scripts/analyze_refinement.py').read_text(encoding='utf-8')
    assert analyzer.count('!= 0.002:') == 1
    analyzer = analyzer.replace('!= 0.002:', '!= 0.0002:')
    analyzer = analyzer.replace('Expected raw force convergence of 0.002 eV/Angstrom',
                                'Expected raw force convergence of 0.0002 eV/Angstrom')
    analyzer = analyzer.replace('**0.002 eV/Å**', '**0.0002 eV/Å**')
    put(ROOT / 'scripts/analyze_refinement.py', analyzer)
    helper = (OLD / 'scripts/keepawake_windows.ps1').read_text(encoding='utf-8-sig')
    helper = helper.replace(OLD.name, ROOT.name)
    helper = helper.replace('PE_DFT_Refinement_KeepAwake_20261022',
                            'PE_DFT_TightForce_KeepAwake_20261022')
    put(ROOT / 'scripts/keepawake_windows.ps1', helper)
    for folder in ['reference', 'environment']:
        for src in (OLD / folder).iterdir():
            if src.is_file():
                freeze(src, ROOT / folder / src.name)
    lock = load(ROOT / 'reference/stage2_parameters_lock.json')
    for paw in lock['environment']['PAW']:
        assert sha(paw['path']) == paw['sha256'], 'PAW changed'
    for src, name in [('protocol.json', 'previous_refinement_protocol.json'),
                      ('sample_manifest.json', 'previous_sample_manifest.json'),
                      ('tables/refinement_summary.json', 'previous_refinement_summary.json')]:
        freeze(OLD / src, ROOT / 'reference' / name)
    protocol = copy.deepcopy(old_protocol)
    protocol.update(batch_id=ROOT.name,
        scope='Continue all four accepted fixed-cell Pnam PBE-D3(BJ) endpoints with a tenfold tighter raw force criterion.',
        authorization='User requested tighter-force continuation of all four DFT endpoints, all other parameters unchanged, with full energy and log retention.',
        continuation={'previous_batch': str(OLD), 'previous_protocol_sha256': sha(OLD / 'protocol.json'),
            'previous_force_target_eV_A': 0.002,
            'only_scientific_parameter_changed': 'optimizer.fmax_eV_A',
            'initialization': 'Previous independently verified final_POSCAR coordinates; fresh electronic state and LBFGS history, same optimizer settings.'})
    protocol['optimizer']['fmax_eV_A'] = TARGET
    save(ROOT / 'protocol.json', protocol)
    samples, queue, checks = [], [], []
    old_samples = load(OLD / 'sample_manifest.json')['samples']
    old_endpoints = {r['sample_id']: r for r in old_summary['endpoints']}
    for old_sample in old_samples:
        sid = old_sample['id']
        endpoint = old_endpoints[sid]
        assert endpoint['accepted'] is True
        work = Path(endpoint['attempt_path'])
        result = load(work / 'result.json')
        process = load(work / 'process.json')
        cold = load(work / 'final_cold_verification.json')
        assert result['status'] == 'passed_relaxation' and process['exit_code'] == 0
        assert cold['passed'] is True and cold['independent_density_and_wavefunctions'] is True
        assert cold['raw_fmax_eV_A'] <= 0.002
        assert sha(work / 'result.json') == endpoint['result_sha256']
        src, dst = work / 'final_POSCAR', ROOT / 'structures' / sid / 'POSCAR'
        atoms = read(src, format='vasp')
        ev = result['final_evaluation']
        delta = (atoms.positions - np.asarray(ev['positions_A'])) @ np.linalg.inv(atoms.cell.array)
        delta -= np.rint(delta)
        position_delta = float(np.abs(delta @ atoms.cell.array).max())
        assert position_delta < 1e-12, 'Previous POSCAR differs from accepted terminal geometry'
        assert atoms.get_chemical_formula() == 'C4H8' and len(atoms) == 12 and atoms.pbc.all()
        assert np.max(abs(atoms.cell.array - np.diag(protocol['cell_lengths_A']))) < 1e-12
        ds = spglib.get_symmetry_dataset((atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers),
                                        symprec=1e-5, hall_number=297)
        assert ds.number == 62 and len(ds.rotations) == 8
        freeze(src, dst)
        sources = {}
        for name in ['result.json', 'process.json', 'final_cold_verification.json', 'final_POSCAR']:
            freeze(work / name, ROOT / 'reference/previous_endpoints' / sid / name)
            sources[name] = {'path': str(work / name), 'sha256': sha(work / name)}
        sample = copy.deepcopy(old_sample)
        sample.update(source_sample_id=sid, source_path=str(src), source_sha256=sha(src),
                      structure=str(dst), structure_sha256=sha(dst),
                      previous_DFT_endpoint=sources,
                      previous_DFT_energy_eV_cell=endpoint['energy_eV_cell'],
                      previous_DFT_fmax_eV_A=cold['raw_fmax_eV_A'])
        samples.append(sample)
        previous_spec = load(OLD / 'inputs' / f'{sid}.json')
        spec = copy.deepcopy(previous_spec)
        spec.update(structure=str(dst), source_geometry_sha256_expected=sha(dst),
                    protocol_sha256_expected=sha(ROOT / 'protocol.json'),
                    point_runner=str(ROOT / 'scripts/run_point.py'), fmax_eV_A=TARGET)
        changed = [k for k in spec if spec[k] != previous_spec[k]]
        allowed = {'structure', 'source_geometry_sha256_expected', 'protocol_sha256_expected',
                   'point_runner', 'fmax_eV_A'}
        assert set(changed) == allowed, changed
        save(ROOT / 'inputs' / f'{sid}.json', spec)
        queue.append({'id': sid, 'input': str(ROOT / 'inputs' / f'{sid}.json'),
                      'input_sha256': sha(ROOT / 'inputs' / f'{sid}.json')})
        checks.append({'id': sid, 'previous_result_sha256': endpoint['result_sha256'],
                       'position_delta_max_A': position_delta,
                       'changed_input_keys': changed, 'source_POSCAR_sha256': sha(src)})
    assert len(samples) == 4
    save(ROOT / 'sample_manifest.json', {'schema_version': 1, 'samples': samples})
    save(ROOT / 'queue.json', {'max_workers': 2, 'protocol_sha256': sha(ROOT / 'protocol.json'), 'points': queue})
    save(ROOT / 'audit/preparation.json', {'passed': True, 'PAW_hashes_verified': 2,
         'source_terminal_checks': checks, 'force_target_before_eV_A': 0.002,
         'force_target_after_eV_A': TARGET, 'merge_criteria_unchanged': True,
         'runner_byte_identical': sha(OLD / 'scripts/run_fixedcell_relax.py') == sha(ROOT / 'scripts/run_fixedcell_relax.py'),
         'quantum_calculations_started': False})
    for name in ['runs', 'tables', 'control']:
        (ROOT / name).mkdir(exist_ok=True)
    put(ROOT / 'README.md', '''# 四代表严格受力继续优化

从上一轮四个已独立验收的 DFT 终态继续。最大原始受力由 0.002 收紧到 **0.0002 eV/Å**；PBE-D3(BJ)、1800 eV、Γ 中心 k=4×8×12、gpts=112×72×44、SCF 阈值、PAW、固定晶胞 7.432×4.929×2.545 Å、Pnam、LBFGS maxstep=0.05 Å、最多200步均沿用。每项4 MPI，最多2项并行。

起始坐标取上一轮 final_POSCAR，已逐原子核对与接受的终态结果一致。电子态和 LBFGS 历史重新初始化；不是重用旧优化器 checkpoint 的逐步续跑。旧批次文件保持不变，新批次保留全部原始能量、原子力、应力、每步结构、SCF日志和独立冷启动复核。

验收后对全部六对使用同一原子映射下 RMS≤0.02 Å、最大位移≤0.05 Å、能差≤1 meV/CH₂ 的联合条件，按完全链接分组。不会为得到单组而放宽条件。最终是否成为一组以计算结果为准；不据此声称全局最低能或独立晶型。
''')
    print(json.dumps({'status': 'prepared', 'root': str(ROOT), 'fmax_eV_A': TARGET,
                      'samples': len(samples), 'protocol_sha256': sha(ROOT / 'protocol.json')}))

if __name__ == '__main__':
    main()
