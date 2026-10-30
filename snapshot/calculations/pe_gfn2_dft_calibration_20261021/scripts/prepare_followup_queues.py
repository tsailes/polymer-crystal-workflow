"""Declare the remaining same-geometry points and explicit reused evidence."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT.parent/'pe_gpaw_convergence_20261020'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def save(path,obj):
    body=json.dumps(obj,ensure_ascii=False,indent=2)+'\n'
    if path.exists(): assert json.loads(path.read_text())==obj, 'Frozen queue differs: '+str(path)
    else:path.write_text(body,encoding='utf-8')
def main():
    entries={s['id']:s for s in json.loads((ROOT/'samples.json').read_text())['samples']}
    refs=[]
    for temperature in [0,300]:
        folder=ROOT/f'pilot/gfn2_temperature_group001/T{temperature:03d}'
        result=json.loads((folder/'result.json').read_text())
        assert result['status']=='passed_single_point' and result['repeat']==[2,3,8]
        refs.append({'backend':'gfn2','sample_id':'final_group_001','method':'GFN2-xTB','precision_label':'F',
             'electronic_temperature_K':temperature,'result_path':str(folder/'result.json'),'input_path':str(folder/'input.json'),
             'process_path':str(folder/'process.json'),'source_result_sha256':sha(folder/'result.json'),
             'note':'Independent temperature pilot; exactly the frozen R0 geometry, validated by source audit.'})
    design=json.loads((ROOT/'audit/design_review.json').read_text())
    for tag,method in [('optb88_vdw','optB88-vdW'),('r2scan_d4','r2SCAN-D4')]:
        folder=OLD/f'runs/transfer_gfn2_{tag}_H'
        verified=next(r for r in design['existing_H_results'] if r['method']==method)
        assert verified['reusable_as_final_group001_H'] and sha(folder/'result.json')==verified['result_sha256']
        refs.append({'backend':'dft','sample_id':'final_group_001','method':method,'precision_label':'H',
              'result_path':str(folder/'result.json'),'input_path':str(folder/'input.json'),
              'process_path':str(folder/'process.json'),'source_result_sha256':sha(folder/'result.json'),
              'note':'Same actual atomic coordinates/cell, PAW, scientific runner, SCF and H precision. Source audit confirms identity mapping.'})
    save(ROOT/'audit/reused_points.json',{'points':refs,'no_old_source_modified':True})
    def gfn(ident,kind='F',temperature=0):
        return {'id':f'{ident}_gfn2_{kind}_T{temperature:03d}','backend':'gfn2','sample_id':ident,
                'spec':{'id':f'{ident}_gfn2_{kind}_T{temperature:03d}','structure':entries[ident]['structure'],
                        'repeat':[2,3,8] if kind=='F' else [3,4,10],'accuracy':0.001,'max_iterations':1000,
                        'threads':4,'electronic_temperature_K':temperature,'charge':0,'multiplicity':1,'precision_label':kind}}
    ends=['initial_cand_0001','initial_cand_0121']
    points=[gfn(i,temperature=t) for i in ends for t in [0,300]]
    points += [gfn(f'final_group_{i:03d}') for i in range(2,5)]
    points += [gfn(i,'D') for i in ['final_group_001','final_group_002','final_group_003','final_group_004',*ends]]
    save(ROOT/'queues/02_gfn2_validation.json',{'phase':'02_gfn2_validation','max_workers':1,'stop_on_failure':True,'points':points})
    provided={p['sample_id'] for p in points if p['spec']['precision_label']=='F' and p['spec']['electronic_temperature_K']==0}|{'final_group_001'}
    save(ROOT/'queues/04_gfn2_paired_production.json',{'phase':'04_gfn2_paired_production','max_workers':1,'stop_on_failure':True,
         'requires':['audit/gfn2_temperature_gate.json','audit/gfn2_sampling_gate.json'],
         'points':[gfn(i) for i in entries if i not in provided]})
    production=json.loads((ROOT/'queues/03_dft_paired_production.json').read_text())
    production['points']=[p for p in production['points'] if not(p['sample_id']=='final_group_001' and p['spec']['method'] in ['optB88-vdW','r2SCAN-D4'])]
    production['requires'].append('audit/gfn2_sampling_gate.json')
    production['supersedes_queue']='03_dft_paired_production.json; removed only two verified reusable R0 points before dispatch'
    save(ROOT/'queues/03_dft_paired_production_reuse_v2.json',production)
    pilot=json.loads((ROOT/'queues/01_dft_precision_pilot.json').read_text())
    reference=next(p for p in pilot['points'] if p['id']=='final_group_001_pbe_d3bj_H')
    reference=json.loads(json.dumps(reference));reference['id']='final_group_001_pbe_d3bj_U'
    reference['spec'].update(ecut_eV=2000,kmesh=[5,10,15],gpts=[128,80,48],precision_label='U',stage='stage3_calibration_U')
    save(ROOT/'queues/01b_dft_reference_upper.json',{'phase':'01b_dft_reference_upper','max_workers':1,'stop_on_failure':True,
       'reason':'Matched R0 upper-precision point needed to resolve relative-energy numerical differences, not just absolute energies.', 'points':[reference]})
    print(json.dumps({'reused_points':len(refs),'gfn2_validation':len(points),'gfn2_later':len(entries)-len(provided),'DFT_later':len(production['points'])}))
if __name__=='__main__':main()
