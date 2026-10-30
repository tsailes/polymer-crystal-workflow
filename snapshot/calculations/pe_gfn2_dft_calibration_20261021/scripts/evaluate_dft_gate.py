"""Evaluate new-geometry numerical transfer at fixed, exactly matched geometry."""
from pathlib import Path
import hashlib,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def point(ident,label):
    matches=sorted((ROOT/'runs/dft').glob(f'{ident}_pbe_d3bj_{label}*'),key=lambda p:p.stat().st_mtime,reverse=True)
    for folder in matches:
        if not(folder/'result.json').exists() or not(folder/'process.json').exists():continue
        r=load(folder/'result.json');p=load(folder/'process.json')
        if r.get('status')=='passed_single_point' and p.get('exit_code')==0:
            return r,{'result_path':str(folder/'result.json'),'sha256':sha(folder/'result.json')}
    raise ValueError(f'Missing passed {ident}/{label}')
def compare(ident,a,b):
    first,sa=point(ident,a);second,sb=point(ident,b)
    assert np.allclose(first['cell_A'],second['cell_A'],atol=1e-12,rtol=0)
    assert np.allclose(first['positions_A'],second['positions_A'],atol=1e-12,rtol=0)
    e=abs(first['total']['energy_free_eV']-second['total']['energy_free_eV'])*1000/12
    f=float(np.linalg.norm(np.array(first['total']['forces_eV_A'])-second['total']['forces_eV_A'],axis=1).max())
    s=float(np.max(np.abs(np.array(first['total']['stress_kbar'])-second['total']['stress_kbar'])))
    return {'sample_id':ident,'precisions':[a,b],'energy_abs_meV_atom':e,'max_force_vector_eV_A':f,'max_stress_component_kbar':s,
            'passed':bool(e<=0.5 and f<=0.005 and s<=0.1),'sources':[sa,sb]}
def main():
    ids=['final_group_001','initial_cand_0001','initial_cand_0121','final_group_004']
    result={'status':'incomplete','passed':False,'checks':[],'low_to_high_diagnostics':[]}
    try:
        for ident in ids:
            result['checks'].append(compare(ident,'H','U'))
            result['low_to_high_diagnostics'].append(compare(ident,'L','H'))
        refH,_=point('final_group_001','H');refU,_=point('final_group_001','U')
        relative=[]
        for ident in ids:
            h,_=point(ident,'H');u,_=point(ident,'U')
            dh=1000*(h['total']['energy_free_eV']-refH['total']['energy_free_eV'])/4
            du=1000*(u['total']['energy_free_eV']-refU['total']['energy_free_eV'])/4
            relative.append({'sample_id':ident,'relative_H_meV_CH2':dh,'relative_U_meV_CH2':du,'difference_meV_CH2':dh-du,
                             'threshold_meV_CH2':1.5,'passed':bool(abs(dh-du)<=1.5)})
        relative_passed=all(r['passed'] for r in relative)
        result.update(status='evaluated',passed=all(c['passed'] for c in result['checks']) and relative_passed,relative_energy_checks=relative,
            all_relative_energy_checks_passed=relative_passed,
            relative_energy_threshold_basis='Existing 0.5 meV/atom criterion converted to 1.5 meV/CH2 (3 atoms per CH2), applied separately to relative-energy differences. This alone does not resolve sub-meV endpoint ordering.',
            observed_relative_energy_change_max_meV_CH2=max(abs(r['difference_meV_CH2']) for r in relative),
            numerical_uncertainty_note='Observed finite-precision difference on four declared points, not a statistical confidence bound; retain separately from model sigma.',
            prior_locked_parameters_changed=False,production_precision='H:1800eV,k4x8x12,g112x72x44')
    except Exception as exc:result['error']=repr(exc)
    (ROOT/'audit/dft_transfer_gate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'passed':result['passed'],'status':result['status'],'checks':len(result['checks'])}))
    if not result['passed']:raise SystemExit(1)
if __name__=='__main__':main()
