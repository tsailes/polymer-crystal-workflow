"""Additional read-only evidence for retained-item x reflection and cell size."""
import json
from pathlib import Path
from collections import defaultdict
import numpy as np
import spglib
from audit_phase_symmetry import make, match, ROOT, OUT, CELL, sha

groups=defaultdict(list)
for p in sorted((ROOT/'branches/O/generated').glob('O-c*.json')):
    d=json.loads(p.read_text());v=d['generation_parameters']
    key=(d['chemistry_template_id'],round(min(v['chain_x_fractional'],.5-v['chain_x_fractional']),3),
         v['chain_rotation_about_c_deg'],round(v['chain_c_phase_fractional']*96))
    groups[key].append(d)
records=[]
for key,members in groups.items():
    primary=members[0]; links=[]
    for other in members[1:]:
        def construct(d):
            v=d['generation_parameters']; tid=int(d['chemistry_template_id'].split('-')[-1])
            return make(tid,v['chain_x_fractional'],v['chain_rotation_about_c_deg'],v['chain_c_phase_fractional'])
        result=match(construct(primary),construct(other),origins=[np.array([0,.5,0])],tol=1e-7)
        assert result['equivalent']
        links.append({'candidate_id':other['candidate_id'],**result})
    records.append({'representative':primary['candidate_id'],'members':[v['candidate_id'] for v in members],
                    'relation':'x_to_half_minus_x_with_template_screw_numerical_tolerance','links':links})
result={'original_retained_entries':sum(map(len,groups.values())),'groups_under_verified_relation':len(groups),
        'two_member_groups':sum(len(x)==2 for x in groups.values()),'singleton_groups':sum(len(x)==1 for x in groups.values()),
        'tolerance_A':1e-7,'full_minimal_crystal_deduplication_proven':False,'groups':records}
(OUT/'retained_x_reflection_groups.json').write_text(json.dumps(result,indent=2)+'\n')
cellcheck=[]
for cid in ['O-c0010','O-c0011','O-c0057']:
    d=json.loads((ROOT/f'branches/O/generated/{cid}.json').read_text())
    nums=[{'H':1,'C':6,'N':7,'O':8}[s] for s in d['symbols']]
    for tol in [1e-5,1e-7,1e-8]:
        dataset=spglib.get_symmetry_dataset((d['cell'],d['fractional_positions'],nums),symprec=tol)
        cellcheck.append({'candidate':cid,'x':d['generation_parameters']['chain_x_fractional'],'tolerance_A':tol,
            'number':int(dataset.number),'operations':len(dataset.rotations),
            'primitive_cell_volume_ratio':float(abs(np.linalg.det(dataset.primitive_lattice)/np.linalg.det(d['cell']))),
            'transformation_matrix':dataset.transformation_matrix.tolist()})
(OUT/'primitive_cell_audit.json').write_text(json.dumps(cellcheck,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='groups'},indent=2));print(json.dumps(cellcheck,indent=2))
