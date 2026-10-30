"""Freeze a pilot protocol and queues; never change the stage-2 release."""
from pathlib import Path
import hashlib, json, shutil

ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT.parent/'pe_gpaw_convergence_20261020'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def save_new(path,obj):
    payload=json.dumps(obj,ensure_ascii=False,indent=2)+'\n'
    if path.exists():
        assert json.loads(path.read_text())==obj, 'Existing frozen object differs: '+str(path)
    else: path.write_text(payload,encoding='utf-8')
def main():
    for name in ['queues','audit','runs/dft','runs/gfn2','tables','control']:
        (ROOT/name).mkdir(parents=True,exist_ok=True)
    samples=json.loads((ROOT/'samples.json').read_text())
    ids={x['id']:x for x in samples['samples']}
    assert len(ids)==21 and samples['reference_id']=='final_group_001'
    src=OLD/'scripts/run_point.py'
    dest=ROOT/'scripts/dft_point.py'
    if dest.exists(): assert sha(dest)==sha(src)
    else: shutil.copy2(src,dest)
    assert sha(dest)=='cc4cd9d563a1d5eb69bc72346b096473ea785832487a26d984c441c72e2f57ae'
    release=OLD/'releases/20261021T043247_514339Z'
    protocol={
      'schema_version':1,'status':'pilot_protocol_frozen_before_dispatch',
      'scope':'Fixed Tashiro cell C4H8 PE; 16 rigid-angle inputs, four preoptimized representatives, and one excluded stationary control.',
      'reference_id':'final_group_001','normalization':'meV per CH2; 4 CH2 per primitive cell',
      'same_geometry_pairing_required':True,'same_cell_required':True,
      'main_DFT_method':'PBE-D3(BJ)','cross_methods':['optB88-vdW','r2SCAN-D4'],
      'DFT_science_runner_sha256':sha(dest),
      'stage2_lock':str(release/'parameters_lock.json'),'stage2_lock_sha256':sha(release/'parameters_lock.json'),
      'stage2_release_check_sha256':sha(release/'release_check.json'),
      'DFT_production_precision':{'ecut_eV':1800,'kmesh':[4,8,12],'gpts':[112,72,44]},
      'DFT_previous_candidate_precision':{'ecut_eV':1400,'kmesh':[3,6,9],'gpts':[96,60,36]},
      'DFT_verification_upper_point':{'ecut_eV':2000,'kmesh':[5,10,15],'gpts':[128,80,48]},
      'upper_point_role':'New-geometry transfer verification only; does not edit the prior locked candidate template.',
      'GFN2_proposed_sampling':{'repeat':[2,3,8],'accuracy':0.001,'threads':4,'max_iterations':1000},
      'temperature_gate':'Explicit T=0 and 300 K electronic singlepoints; no silent substitution. Native energy/free_energy alias is not evidence of zero entropy.',
      'temperature_change_thresholds':{'energy_meV_CH2':0.1,'max_force_vector_eV_A':0.001,'max_stress_component_kbar':0.1},
      'numerical_transfer_thresholds':{'energy_abs_meV_atom':0.5,'force_vector_eV_A':0.005,'stress_component_kbar':0.1},
      'split_strategy':'Four contiguous initial-angle blocks; leave one block out. Preoptimized endpoints and unstable control are separate diagnostic domains, never pooled to inflate independent n.',
      'calibration_models_to_audit':['identity','mean_gap_bias_only','affine_gap_model'],
      'model_selection_rule':'Report blocked held-out residuals and bias; no training-residual sigma presented as prediction error. Do not apply a fitted high-energy input model to near-minimum endpoints without separate evidence.',
      'sigma_definition':'Standard deviation of held-out DFT-minus-predicted relative-energy residuals, for a declared domain and fixed R0. Not the spread of candidate energies.',
      'screening_rule':'Only consider rejection if predicted relative energy minus 3*sigma_deltaE exceeds a predeclared low-energy window; requires independent validation and numerical uncertainty below the relevant signal.',
      'production_rejection_enabled':False,'screening_window_meV_CH2':None,
      'small_sample_policy':'Four closely related endpoint groups cannot establish a robust 3sigma coverage claim. Keep all four unless later evidence validates a conservative decision.',
      'cross_functional_policy':'Differences between DFT methods are reported separately, never mixed into a single GFN2 sigma.',
      'raw_data_retention':'All inputs, source/frozen structure hashes, stdout, stderr, SCF records, complete energy/force/stress JSON and failed attempts retained.',
      'sample_manifest_sha256':sha(ROOT/'samples.json')}
    save_new(ROOT/'protocol.json',protocol)
    common={'mpi_ranks':4,'point_group':True,'xc_stencil':4,'nbands':16,'scf_energy_eV_electron':1e-9,'scf_density_electron_electron':1e-9,'scf_eigenstates_eV2_electron':1e-11,'checkpoint':False}
    def point(ident,precision,method='PBE-D3(BJ)'):
        tag={'PBE-D3(BJ)':'pbe_d3bj','optB88-vdW':'optb88_vdw','r2SCAN-D4':'r2scan_d4'}[method]
        return {'id':f'{ident}_{tag}_{precision}','sample_id':ident,'backend':'dft','spec':dict(common,**protocol[{'H':'DFT_production_precision','L':'DFT_previous_candidate_precision','U':'DFT_verification_upper_point'}[precision]],method=method,structure=ids[ident]['structure'],stage='stage3_calibration_'+precision,precision_label=precision)}
    endpoints=['initial_cand_0001','initial_cand_0121','final_group_004']
    pilot=[point('final_group_001','H')]
    for ident in endpoints:
        pilot.extend([point(ident,'H'),point(ident,'U'),point(ident,'L')])
    pilot.append(point('final_group_001','L'))
    pilot_ids={p['id'] for p in pilot}
    production=[point(s['id'],'H') for s in samples['samples'] if point(s['id'],'H')['id'] not in pilot_ids]
    for method in protocol['cross_methods']:
        production.extend(point(f'final_group_{n:03d}','H',method) for n in range(1,5))
    save_new(ROOT/'queues/01_dft_precision_pilot.json',{'phase':'01_dft_precision_pilot','max_workers':2,'stop_on_failure':True,'points':pilot})
    save_new(ROOT/'queues/03_dft_paired_production.json',{'phase':'03_dft_paired_production','max_workers':2,'stop_on_failure':True,'requires':['audit/dft_transfer_gate.json','audit/gfn2_temperature_gate.json'],'points':production})
    print(json.dumps({'samples':len(ids),'pilot_DFT_points':len(pilot),'later_DFT_points':len(production),'protocol_sha256':sha(ROOT/'protocol.json')}))
if __name__=='__main__': main()
