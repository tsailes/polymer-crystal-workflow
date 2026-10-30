"""Build the PE route evidence index; standard library only, no quantum jobs.

Usage: python freeze_route.py
Only writes route_lock.json, source_manifest.json and validation.json beside this file.
"""
import csv
import hashlib
import json
import re
import statistics
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

def read_json(path):
    return json.loads((ROOT / path).read_text(encoding='utf-8-sig'))

def rows(path):
    with (ROOT / path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def write(name, data):
    (HERE / name).write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

sources = {
    'candidate_protocol': 'calculations/pe_pnam_candidates_20261020/protocol.json',
    'input_provenance': 'calculations/pe_pnam_candidates_20261020/input_provenance.json',
    'candidate_catalog': 'calculations/pe_pnam_candidates_20261020/candidate_catalog.csv',
    'optimization_summary': 'calculations/pe_pnam_candidates_20261020/optimization_summary.json',
    'optimization_catalog': 'calculations/pe_pnam_candidates_20261020/optimization_catalog.csv',
    'dft_original_release': 'calculations/pe_gpaw_convergence_20261020/releases/20261021T043247_514339Z/parameters_lock.json',
    'dft_release_acceptance': 'calculations/pe_gpaw_convergence_20261020/releases/20261021T043247_514339Z/release_check.json',
    'calibration_protocol': 'calculations/pe_gfn2_dft_calibration_20261021/protocol.json',
    'calibration_policy': 'calculations/pe_gfn2_dft_calibration_20261021/analysis_policy.json',
    'initial_calibration': 'reports/PE_GFN2_3sigma_20261021/data/initial_16_energy_calibration.csv',
    'paired_production': 'calculations/pe_gfn2_dft_calibration_20261021/tables/paired_energy_results.csv',
    'normal_fixed_protocol': 'calculations/pe_dft_representative_refinement_20261022/protocol.json',
    'normal_fixed_results': 'calculations/pe_dft_representative_refinement_20261022/tables/final_energy_results.csv',
    'tight_fixed_protocol': 'calculations/pe_dft_tight_force_refinement_20261022/protocol.json',
    'tight_fixed_results': 'calculations/pe_dft_tight_force_refinement_20261022/tables/final_energy_results.csv',
    'tight_pairwise': 'calculations/pe_dft_tight_force_refinement_20261022/tables/pairwise_geometry_energy.csv',
    'fixed_g4': 'calculations/pe_dft_tight_force_refinement_20261022/runs/group_004/attempt_01/final_structure.extxyz',
    'variable_protocol': 'calculations/pe_g4_variable_cell_20261027/protocol.json',
    'variable_result': 'calculations/pe_g4_variable_cell_20261027/runs/attempt_01/result.json',
    'variable_g4': 'calculations/pe_g4_variable_cell_20261027/runs/attempt_01/final_structure.extxyz',
    'variable_cold': 'calculations/pe_g4_variable_cell_20261027/final_cold_verification.json',
    'variable_convergence_protocol': 'calculations/pe_g4_variable_cell_20261027/convergence_checks/protocol.json',
    'variable_convergence_result': 'calculations/pe_g4_variable_cell_20261027/convergence_checks/result.json',
    'variable_convergence_review': 'calculations/pe_g4_variable_cell_20261027/convergence_checks/independent_review.json',
    'angle_coverage': 'reports/pe_questions_audit_20261028/全取向角覆盖与阈值敏感性.md',
    'supplement_not_run': 'reports/pe_questions_audit_20261028/supplemental_angles_080_085_090/README_补充输入说明.md',
    'pbe0_status': 'calculations/pe_pbe0_singlepoint_review_20261023/tables/singlepoint_summary.json',
    'cost_comparison_status': 'calculations/pe_energy_models_benchmark_20261020/state.json',
    'latest_fixed_supplement': 'reports/PE技术文档_正文与SI_20261028/data/fixed_cell/supplement_20261029/summary.json',
    'latest_technical_report': 'reports/PE技术文档_正文与SI_20261028/outputs/PE技术文档_正文与SI_固定胞结果补充版_20261029.md',
    'methods_audit': 'reports/PE技术文档_正文与SI_20261028/audit/methods_results.md',
    'search_audit': 'reports/PE技术文档_正文与SI_20261028/audit/search_calibration.md',
}

manifest = []
for role, path in sources.items():
    p = ROOT / path
    if not p.is_file():
        raise FileNotFoundError(p)
    manifest.append({'role': role, 'path_relative_to_project': path,
                     'bytes': p.stat().st_size, 'sha256': sha(p)})

opt = read_json(sources['optimization_summary'])
original = read_json(sources['dft_original_release'])
tight = read_json(sources['tight_fixed_protocol'])
variable = read_json(sources['variable_protocol'])
cold = read_json(sources['variable_cold'])
review = read_json(sources['variable_convergence_review'])
policy = read_json(sources['calibration_policy'])
cal = rows(sources['initial_calibration'])
fixed = rows(sources['tight_fixed_results'])
pairs = rows(sources['tight_pairwise'])
errors = [float(r['PBE_D3BJ_relative_energy_meV_CH2']) - float(r['GFN2_relative_energy_meV_CH2']) for r in cal]
bias, sigma = statistics.mean(errors), statistics.stdev(errors)
fixed_map = {r['sample_id']: r for r in fixed}
pbe0 = read_json(sources['pbe0_status'])
cost = read_json(sources['cost_comparison_status'])

checks = []
def check(name, condition, value=None):
    checks.append({'name': name, 'passed': bool(condition), 'observed': value})

check('candidate_catalog_288', len(rows(sources['candidate_catalog'])) == 288)
check('optimization_catalog_16', len(rows(sources['optimization_catalog'])) == 16)
check('generation_counts', opt['generation_counts'] == {
    'enumerated_count': 288, 'physical_survivor_count': 248, 'unique_initial_count': 16,
    'rejected_count': 40, 'exact_duplicates_count': 232}, opt['generation_counts'])
check('retained_15_representatives_4', opt['retained_count'] == 15 and opt['geometrically_distinct_final_count'] == 4)
check('paired_DFT_29', len(rows(sources['paired_production'])) == 29)
check('sigma_16_initials', len(cal) == 16 and all(r['sample_id'].startswith('initial_') for r in cal))
check('bias_recomputed', abs(bias - 4.638333875888634) < 1e-9, bias)
check('sigma_recomputed_ddof1', abs(sigma - 1.920697507166551) < 1e-9, sigma)
check('four_tight_fixed_accepted', len(fixed) == 4 and all(r['accepted'] == 'True' for r in fixed))
check('four_tight_raw_forces', all(float(r['raw_fmax_eV_A']) <= 0.0002 for r in fixed))
check('four_cold_verified', all(json.loads(r['cold_verification'])['passed'] for r in fixed))
check('six_pairs_complete', len(pairs) == 6 and len({tuple(sorted((r['first'],r['second']))) for r in pairs}) == 6)
check('six_pairs_merge_pass', all(r['merge_pass'] == 'True' and float(r['rms_A']) <= .02 and float(r['max_A']) <= .05 and float(r['energy_difference_meV_CH2']) <= 1 for r in pairs))
check('fixed_g4_identity', sha(ROOT / sources['fixed_g4']) == fixed_map['group_004']['final_geometry_sha256'] == variable['source_sha256'])
check('variable_g4_identity', sha(ROOT / sources['variable_g4']) == cold['geometry_sha256'] == review['geometry_sha256'])
check('variable_cold_accepted', cold['passed'] and cold['independent_density_and_wavefunctions'] and cold['raw_fmax_eV_A'] <= .002 and cold['raw_max_stress_component_kbar'] <= .05)
check('variable_numerical_review', review['status'] == 'passed')
check('no_sigma_window', policy['screening_window_meV_CH2'] is None)
check('pbe0_cancelled_without_energy', pbe0['status'] == 'cancelled_by_user' and pbe0['accepted_count'] == 0 and not pbe0['valid_final_energy_available'])
check('cost_production_not_started', cost['production_attempts_started'] == 0)
check('production_H', tight['ecut_eV'] == 1800 and tight['kmesh'] == [4,8,12] and tight['gpts'] == [112,72,44])
check('variable_same_H', all(variable['parameters'][k] == tight[k] for k in ('ecut_eV','kmesh','gpts')))

lock = {
    'schema_version': 1, 'route_id': 'PE_CRYSTAL_ROUTE_V1_0', 'version': '1.0',
    'frozen_date_Asia_Shanghai': '2026-10-29',
    'status': 'documented_baseline_from_completed_PE_mainline',
    'canonical_document': 'reports/PE技术路线_v1_0_20261029/技术路线.md',
    'purpose_zh': '以PE为典型对象，固定聚合物晶体结构重建、分层精修、数值验收及实验对照流程',
    'scope': 'PE C4H8; experimental fixed cell and Pnam prior; finite all-trans short-period candidate space; constrained local optimization.',
    'not_an_executable_queue': True, 'new_quantum_calculations_started_by_freeze': 0,
    'historical_protocols_unchanged': True,
    'stages': ['S0_input_and_prior_freeze','S1_candidate_generation','S2_geometry_filter_and_equivalence',
               'S3_GFN2_preoptimization','S4_same_geometry_DFT_calibration','S5_fixed_cell_DFT_and_complete_link_merge',
               'S6_representative_variable_cell','S7_independent_checks_and_experiment','S8_archive'],
    'dependencies': {'S4': ['S3', 'N_DFT_numerical_validation'], 'S5': ['S4', 'N_DFT_numerical_validation']},
    'PE_model': {'formula':'C4H8','atoms':12,'N_CH2':4,'chains':2,'chain_axis':'c',
                 'cell_A':[7.432,4.929,2.545],'angles_deg':[90,90,90],
                 'spacegroup':'Pnam','number':62,'hall_number':297,
                 'initial_CC_A':1.54,'initial_CH_A':1.09},
    'candidate_counts': {**opt['generation_counts'], 'GFN2_completed':16, 'directionally_unstable':1,
                         'retained_endpoints':15,'GFN2_representatives':4,'tight_DFT_groups':1,'energy_rejections':0},
    'GFN2': {'backend':'tblite 0.7.0','model':'GFN2-xTB','sampling':'Gamma repeated supercell, not explicit k mesh',
             'coarse_repeat':[2,2,6],'fine_repeat':[2,3,8],'check_repeat':[3,4,10],
             'coarse_fmax_eV_A':.03,'fine_fmax_eV_A':.01,'optimization_accuracy':.01,
             'cold_accuracy':.001,'optimization_electronic_K':300,'calibration_electronic_K':0,
             'check_abs_delta_energy_meV_CH2_lt':1,'check_max_delta_force_eV_A_lt':.005},
    'DFT': {'engine':'GPAW 26.7.0 legacy CPU PW','PAW':'gpaw-data 1.0.1 native PBE C/H',
            'method_primary':tight['method'],'ATM':True,
            'same_geometry_comparison_methods':['PBE-D3(BJ)','optB88-vdW','r2SCAN-D4'],
            'historical_L':original['candidate_precision'],
            'production_H':{k:tight[k] for k in ('ecut_eV','kmesh','gpts','xc_stencil','nbands','mpi_ranks','scf_energy_eV_electron','scf_density_electron_electron','scf_eigenstates_eV2_electron')},
            'electronic_temperature_K':0, 'PAW_sha256':original['physics']['paw_sha256'],
            'method_details_source_role':'dft_original_release'},
    'fixed_cell': {'normal_fmax_eV_A':.002,'tight_fmax_eV_A':.0002,'stress_is_acceptance_target':False,
                   'optimizer':tight['optimizer'],'merge':tight['merge'],'accepted_count':4,'representative':'group_004',
                   'accepted_attempts':{r['sample_id']:r['attempt'] for r in fixed},
                   'max_pairwise_rms_A':max(float(r['rms_A']) for r in pairs),
                   'source_role':'fixed_g4','sha256':sha(ROOT/sources['fixed_g4'])},
    'variable_cell': {'source_role':'fixed_g4','cell_mask':variable['cell_mask'],'external_pressure_eV_A3':0,
                      'angles_degrees':variable['angles_degrees'],'fmax_eV_A':.002,'stress_max_kbar':.05,
                      'terminal_cold':cold,'numerical_check_status':review['status'],
                      'ecut_comparison_rule':'compare E2000/g128 with E1800/g128 grid_only; same geometry',
                      'stricter_settings_reoptimization_completed':False},
    'calibration': {'frozen_geometries':21,'production_DFT_pairs':29,'initial_n':16,'reference_id':'final_group_001',
                    'error':'DeltaE_DFT - DeltaE_GFN2','unit':'meV/CH2','bias':bias,'sigma':sigma,
                    'ddof':1,'three_sigma':3*sigma,'screening_window':None,'energy_rejections':0,
                    'scope':policy['domain'],'coverage_probability_guarantee':None,
                    'independent_coverage_validation_required':False,'DFT_savings_demonstration_required':False},
    'open_or_cancelled_branches': {'PBE0':pbe0['status'],'PBE0_accepted':0,'formal_cost_comparison_started':0,
                                  'supplemental_80_85_90_deg':'inputs_only_not_run'},
    'claim_limits': ['not a global-minimum proof','not a phonon or Hessian stability proof',
                     'no ZPE or finite-temperature nuclear free energy','fixed and variable cells reported separately',
                     'no final-g4 cross-functional geometry optimization','indexed XRD positions only; no full intensity fit',
                     '3sigma does not imply probabilistic coverage or executed pruning'],
    'reuse_rule': 'Reuse stages and evidence rules; re-establish material-specific search domain, numerical convergence and calibration.',
    'sources': sources,
}
write('route_lock.json', lock)
write('source_manifest.json', {'route_id':lock['route_id'],'project_root':str(ROOT),
                              'indexed_at':datetime.now(timezone.utc).isoformat(),
                              'source_count':len(manifest),'sources':manifest})

# Validate authored artifact links and table shape, rather than launching calculations.
link_failures = []
table_failures = []
for p in [ROOT/'PE技术路线.md', HERE/'技术路线.md', HERE/'验收清单.md']:
    body = p.read_text(encoding='utf-8')
    for target in re.findall(r'\]\(([^)]+)\)', body):
        if target == 'validation.json':
            continue  # emitted after these checks
        if not target.startswith(('http:', 'https:', '#')) and not (p.parent / target).exists():
            link_failures.append({'file':p.name,'target':target})
    width = None
    for n,line in enumerate(body.splitlines(),1):
        if line.startswith('|'):
            count = len(re.findall(r'(?<!\\)\|', line))
            if width is not None and count != width:
                table_failures.append({'file':p.name,'line':n,'pipes':count,'expected':width})
            width = count
        else:
            width = None
check('local_links_resolve', not link_failures, link_failures)
check('markdown_table_shapes', not table_failures, table_failures)
ET.parse(HERE/'技术路线图.svg')
check('diagram_valid_XML', True)
outputs = ['技术路线.md','验收清单.md','技术路线图.svg','route_lock.json','source_manifest.json','freeze_route.py']
write('validation.json', {
    'scope':'source identity, local table arithmetic, terminal status, authored links and table structure; no quantum recalculation',
    'checked_at':datetime.now(timezone.utc).isoformat(),
    'passed':all(c['passed'] for c in checks),'checks':checks,
    'artifacts':[{'path':name,'sha256':sha(HERE/name)} for name in outputs],
})
failed = [c for c in checks if not c['passed']]
print(json.dumps({'checks':len(checks),'passed':not failed,'failed':failed,'sources':len(manifest)},ensure_ascii=False))
if failed:
    raise SystemExit(1)
