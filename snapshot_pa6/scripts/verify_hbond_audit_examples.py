"""Reproduce the historical HB audit with archived validator 0.2; no energy.

This is not the current physical admission entry point. Do not use its
historical HB-containing overall pass label for the new candidate pool.
"""
import importlib.util
import json
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import write
import audit_hbond_filter as a

_source=Path(__file__).resolve().parents[1]/'branches/O/validation/validate_ordered_v0_2_source.py'
_spec=importlib.util.spec_from_file_location('historical_validator_v0_2',_source)
v=importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(v)

a.init_worker();protocol,cell,templates=a.DATA
records=a.read(a.OUT/'grid_records.json')
selected=[]
for predicate in (
    lambda r:r['audit_id']=='grid-0329',
    lambda r:r['original_full_HB_coverage'],
    lambda r:r['donors_covered']==0 and r['original_builder_geometry_without_HB_passed'],
    lambda r:r['independent_geometry_without_HB_passed'] and not r['original_builder_geometry_without_HB_passed'],
    lambda r:not r['independent_geometry_without_HB_passed'],
):
    selected.append(next(r for r in records if predicate(r)))
out=[]
for rec in selected:
    p=rec['parameters'];t=templates[p['template']]
    data=a.b.expand(np.array(t['positions_unwrapped_A']),t['atom_roles'],t['bonds'],t['symbols'],cell,p['x'],p['rotation_deg'],p['phase_c'])
    check=v.validate(data)
    good=check['geometry_topology_symmetry_passed']==rec['independent_geometry_without_HB_passed']
    assert good,rec['audit_id']
    out.append({'audit_id':rec['audit_id'],'full_independent_validator_geometry_passed':check['geometry_topology_symmetry_passed'],
                'audit_agrees_with_full_validator':good,'independent_hb_donors':check['hydrogen_bonds']['covered_donor_count'],
                'original_gate_donors':rec['donors_covered']})
    assert check['hydrogen_bonds']['covered_donor_count']==rec['donors_covered']
    if rec['audit_id']=='grid-0329':
        example_dir=a.OUT/'examples';example_dir.mkdir(exist_ok=True)
        a.save(example_dir/'grid-0329.json',data)
        write(example_dir/'grid-0329.extxyz',Atoms(data['symbols'],positions=data['positions'],cell=cell,pbc=True))
        cell_,symbols,_,frac=v.structure_arrays(data)
        _,adj,_,_=v.infer_graph(symbols,frac,cell)
        contacts=rec['closest_8_donor_contacts']
        modified=[]
        for distance_cutoff in (3.5,3.51,3.52):
            picked=[c for c in contacts if c['NO_A']<=distance_cutoff and c['HO_A']<=2.6 and c['NHO_deg']>=120]
            hb=[]
            for c in picked:
                nh_shift=next(s for j,s in adj[c['N']] if j==c['H'])
                hb.append({'H':c['H'],'O':c['O'],'oxygen_image_from_H':(np.array(c['O_shift'])-nh_shift).tolist()})
            sheets=v.hydrogen_bond_sheet_topology(adj,{'bonds':hb})
            modified.append({'NO_cutoff_A':distance_cutoff,'HO_cutoff_A':2.6,'NHO_min_deg':120,
                'donors_covered':len({c['N'] for c in picked}),'acceptors_covered':len({c['O'] for c in picked}),
                'two_ac_infinite_sheets':sheets['passed'],
                'note':'Identical coordinates; only diagnostic cutoff changes. No optimization or prediction of relaxation.'})
        a.save(example_dir/'grid-0329_full_independent_validation.json',check)
        a.save(example_dir/'grid-0329_threshold_sensitivity.json',{'parameters':p,'contacts':contacts,'cutoff_sweep':modified})
a.save(a.OUT/'crosscheck_results.json',{'all_passed':all(r['audit_agrees_with_full_validator'] for r in out),'examples':out})
print(json.dumps({'crosschecks':out,'example_threshold_sweep':modified},indent=2))
