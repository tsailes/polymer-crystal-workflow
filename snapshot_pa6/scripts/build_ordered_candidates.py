"""Independent ordered PA6 seeds, without importing any reference coordinates.

Geometry uses MMFF94 equilibrium parameters of a SMILES-only cyclic dimer.
The cyclic graph is unwrapped into a translational polymer; no molecular
embedding, crystal databases, paper coordinates, or force-field energy ranking.
Bond lengths stay fixed while internal bond angles satisfy chain closure.
The active entry point reuses the frozen independent templates. Hydrogen-bond
geometry and connectivity never participate in physical admission or ordering.
"""
from __future__ import annotations
import argparse, hashlib, itertools, json, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, linear_sum_assignment

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'inputs/independent_template/dependencies'))
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from ase import Atoms
from ase.io import write
from ase.geometry import cellpar_to_cell, find_mic
import spglib

VERSION = 'O-independent-mmff-angles-v1.2-no-hbond-screen'
OPS_R = np.array([np.eye(3), np.diag([-1,-1,1]), -np.eye(3), np.diag([1,1,-1])], int)
OPS_T = np.array([[0,0,0],[.5,.5,.5],[0,0,0],[.5,.5,.5]])
SHIFTS = np.array(list(itertools.product([-1,0,1], repeat=3)))

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def save(path, data):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')

def molecule_parameters():
    smiles = 'N1C(=O)CCCCCNC(=O)CCCCC1'
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94')
    if props is None:
        raise RuntimeError('MMFF parameters unavailable')
    backbone = [0,1,3,4,5,6,7,8,9,11,12,13,14,15]
    params = {'smiles': smiles, 'rdkit_version': rdBase.rdkitVersion,
              'method':'MMFF94 equilibrium parameters from chemical connectivity only',
              'documentation':'https://rdkit.org/docs/source/rdkit.ForceField.rdForceField.html',
              'coordinate_source':None, 'bond_parameters':[], 'angle_parameters':[]}
    for bond in mol.GetBonds():
        i,j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        p = props.GetMMFFBondStretchParams(mol,i,j)
        if p is None:
            raise RuntimeError(f'Missing bond parameter {i},{j}')
        params['bond_parameters'].append({'i':i,'j':j,'order':bond.GetBondTypeAsDouble(),
                                         'bond_type':p[0],'kb':p[1],'equilibrium_A':p[2]})
    for a in mol.GetAtoms():
        j=a.GetIdx()
        for i,k in itertools.combinations([v.GetIdx() for v in a.GetNeighbors()],2):
            p=props.GetMMFFAngleBendParams(mol,i,j,k)
            if p is None:
                raise RuntimeError(f'Missing angle parameter {i},{j},{k}')
            params['angle_parameters'].append({'i':i,'j':j,'k':k,'angle_type':p[0],
                                              'ka':p[1],'equilibrium_deg':p[2]})
    return mol,props,backbone,params

def build_chain(c, bias_deg, mol, props, backbone):
    n=len(backbone)
    lengths=np.array([props.GetMMFFBondStretchParams(mol,backbone[k],backbone[(k+1)%n])[2]
                      for k in range(n)])
    angle_params=np.array([props.GetMMFFAngleBendParams(mol,backbone[(k-1)%n],backbone[k],backbone[(k+1)%n])[1:]
                           for k in range(n)])
    target=angle_params[:,1].copy()
    for k in [0,1,7,8]: target[k]+=bias_deg
    signs=(-1.)**np.arange(n)
    initial=signs*np.arccos(c/sum(lengths))
    def bondvec(t): return lengths[:,None]*np.column_stack([np.sin(t),np.zeros(n),np.cos(t)])
    def angles(t):
        return 180-np.degrees(np.abs(t-np.roll(t,1)))
    def objective(t):
        return float(np.sum(angle_params[:,0]*(angles(t)-target)**2))
    def closure(t): return np.array([sum(lengths*np.sin(t)),sum(lengths*np.cos(t))-c])
    result=minimize(objective,initial,method='SLSQP',constraints=[{'type':'eq','fun':closure}],
                    bounds=[(s*.2,s*.9) if s>0 else (s*.9,s*.2) for s in signs],
                    options={'ftol':1e-11,'maxiter':1000})
    if not result.success or np.linalg.norm(closure(result.x))>1e-7:
        raise RuntimeError(f'Periodic closure failed: {result.message}; {closure(result.x)}')
    vec=bondvec(result.x)
    p=np.zeros((mol.GetNumAtoms(),3))
    chain=np.vstack([np.zeros(3),np.cumsum(vec,axis=0)])[:-1]
    for k,i in enumerate(backbone): p[i]=chain[k]
    roles={}
    for k,i in enumerate(backbone):
        roles[i]={'index':i,'repeat_index':k//7,'role':('amide_N' if k%7==0 else 'carbonyl_C' if k%7==1 else f'methylene_C{k%7-1}')}
        u=-vec[(k-1)%n]/lengths[(k-1)%n]; v=vec[k]/lengths[k]
        bis=-(u+v); bis/=np.linalg.norm(bis)
        attached=[a.GetIdx() for a in mol.GetAtomWithIdx(i).GetNeighbors() if a.GetIdx() not in backbone]
        if len(attached)==1:
            j=attached[0]; r=props.GetMMFFBondStretchParams(mol,i,j)[2]
            p[j]=p[i]+r*bis
            roles[j]={'index':j,'repeat_index':k//7,'role':'amide_H' if mol.GetAtomWithIdx(j).GetSymbol()=='H' else 'carbonyl_O'}
        elif len(attached)==2:
            h=attached[0]
            target_cos=np.mean([np.cos(np.deg2rad(props.GetMMFFAngleBendParams(mol,h,i,b)[2]))
                               for b in [backbone[(k-1)%n],backbone[(k+1)%n]]])
            cos_alpha=np.clip(target_cos/np.dot(bis,u),0.01,.99)
            for sign,j in zip([-1,1],attached):
                direction=cos_alpha*bis+sign*np.sqrt(1-cos_alpha**2)*np.array([0.,1.,0.])
                r=props.GetMMFFBondStretchParams(mol,i,j)[2]
                p[j]=p[i]+r*direction
                roles[j]={'index':j,'repeat_index':k//7,'role':f'methylene_H{k%7-1}_{"a" if sign<0 else "b"}'}
        else: raise RuntimeError('Unexpected sidegroup degree')
    p[:,0]-=np.mean(p[backbone,0]); p[:,2]-=np.mean(p[backbone,2])
    bonds=[]
    for b in mol.GetBonds():
        i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx(); shift=[0,0,0]
        if {i,j}=={backbone[-1],backbone[0]}:
            shift[2]=1 if i==backbone[-1] else -1
        bonds.append({'i':i,'j':j,'shift':shift,'order':b.GetBondTypeAsDouble(),
                      'equilibrium_A':props.GetMMFFBondStretchParams(mol,i,j)[2]})
    # This is a geometric constrained seed construction, not a full MMFF minimum.
    report={'success':bool(result.success),'closure_residual_A':closure(result.x).tolist(),
            'bond_lengths_fixed':True,'uniform_bond_scaling':False,'backbone_angles_deg':angles(result.x).tolist(),
            'MMFF_backbone_equilibrium_angles_deg':angle_params[:,1].tolist(),
            'amide_angle_bias_deg':bias_deg,'full_force_field_optimization':False,
            'backbone_max_angle_departure_deg':float(np.max(np.abs(angles(result.x)-angle_params[:,1])))}
    return p,[roles[k] for k in range(len(p))],bonds,report

def expand(p, roles, bonds, symbols, cell, x, phi, phase):
    angle=np.deg2rad(phi)
    rot=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1.]])
    f=(p@rot.T)@np.linalg.inv(cell)+[x,0,phase]
    allf=[]; allroles=[]; allbonds=[]
    for g,(R,t) in enumerate(zip(OPS_R,OPS_T)):
        fu=f@R.T+t
        wraps=np.floor(fu+1e-12).astype(int)
        fw=fu-wraps
        allf.extend(fw)
        for role in roles:
            allroles.append(dict(role,index=role['index']+g*len(p),chain_id=f'O-chain-{g}',symmetry_operation_index=g))
        for b in bonds:
            shift=np.array(b['shift'])@R.T+wraps[b['j']]-wraps[b['i']]
            allbonds.append(dict(b,i=b['i']+g*len(p),j=b['j']+g*len(p),shift=shift.tolist()))
    allf=np.array(allf)
    return {'cell':cell.tolist(),'symbols':symbols*4,'fractional_positions':allf.tolist(),
            'positions':(allf@cell).tolist(),'atom_roles':allroles,'bonds':allbonds,
            'generation_parameters':{'chain_x_fractional':x,'chain_rotation_about_c_deg':phi,'chain_c_phase_fractional':phase}}

def hydrogen_bonds(data):
    """Legacy audit helper only; not called by generation or physical screening."""
    p=np.array(data['positions']); cell=np.array(data['cell']); roles=data['atom_roles']
    ns=[a['index'] for a in roles if a['role']=='amide_N']; os=[a['index'] for a in roles if a['role']=='carbonyl_O']
    hs={}
    for b in data['bonds']:
        i,j=b['i'],b['j']
        if i in ns and roles[j]['role']=='amide_H': hs[i]=(j,np.array(b['shift']))
        if j in ns and roles[i]['role']=='amide_H': hs[j]=(i,-np.array(b['shift']))
    hits=[]
    for n in ns:
        h,hshift=hs[n]; hp=p[h]+hshift@cell
        for o in os:
            if roles[n]['chain_id']==roles[o]['chain_id']: continue
            disp=p[o]+SHIFTS@cell-p[n]
            no=np.linalg.norm(disp,axis=1)
            for k in np.where(no<=3.5)[0]:
                op=p[o]+SHIFTS[k]@cell; ho=np.linalg.norm(op-hp)
                angle=np.degrees(np.arccos(np.clip(np.dot(p[n]-hp,op-hp)/(np.linalg.norm(p[n]-hp)*ho),-1,1)))
                if ho<=2.6 and angle>=120:
                    hits.append({'N':n,'H':h,'O':o,'O_shift':SHIFTS[k].tolist(),'NO_A':float(no[k]),'HO_A':float(ho),'NHO_deg':float(angle)})
    return {'bonds':hits,'donors_covered':len({x['N'] for x in hits}),'acceptors_covered':len({x['O'] for x in hits})}

def validate_geometry(data):
    atoms=Atoms(data['symbols'],positions=data['positions'],cell=data['cell'],pbc=True)
    p=atoms.positions; cell=atoms.cell.array
    adjacency={i:set() for i in range(len(p))}
    deviations=[]
    for b in data['bonds']:
        adjacency[b['i']].add(b['j']); adjacency[b['j']].add(b['i'])
        d=np.linalg.norm(p[b['j']]+np.array(b['shift'])@cell-p[b['i']])
        deviations.append(abs(d-b['equilibrium_A']))
    d=atoms.get_all_distances(mic=True)
    clashes=[]; min_ratio=1e9
    lower={('C','C'):2.45,('C','H'):1.6,('C','N'):2.3,('C','O'):2.25,('H','H'):1.25,
           ('H','N'):1.3,('H','O'):1.2,('N','N'):2.3,('N','O'):2.2,('O','O'):2.3}
    for i in range(len(p)):
        excluded=adjacency[i].copy()
        for j in adjacency[i]: excluded.update(adjacency[j])
        for j in range(i+1,len(p)):
            if j in excluded: continue
            bound=lower[tuple(sorted([data['symbols'][i],data['symbols'][j]]))]
            min_ratio=min(min_ratio,d[i,j]/bound)
            if d[i,j]<bound: clashes.append({'i':i,'j':j,'distance_A':float(d[i,j]),'lower_A':bound})
    residuals=[]; f=np.array(data['fractional_positions'])
    for R,t in zip(OPS_R,OPS_T):
        transformed=(f@R.T+t)@cell
        max_res=0.
        for sym in set(data['symbols']):
            ids=np.where(np.array(data['symbols'])==sym)[0]
            v=transformed[ids,None,:]-p[None,ids,:]
            _,dist=find_mic(v.reshape(-1,3),cell,pbc=True)
            cost=dist.reshape(len(ids),len(ids)); ii,jj=linear_sum_assignment(cost)
            max_res=max(max_res,float(cost[ii,jj].max()))
        residuals.append(max_res)
    dataset=spglib.get_symmetry_dataset((cell,f,atoms.numbers),symprec=1e-5)
    report={'atom_count':len(p),'formula':atoms.get_chemical_formula(),
            'max_bond_length_error_A':float(max(deviations)),'clashes':clashes,
            'minimum_nonbonded_lower_bound_ratio':float(min_ratio),'operation_max_residuals_A':residuals,
            'spglib_number':None if dataset is None else int(dataset.number),
            'spglib_international':None if dataset is None else dataset.international,
            'spglib_pointgroup':None if dataset is None else dataset.pointgroup,
            'hydrogen_bonds_evaluated':False,
            'hydrogen_bonds_used_for_admission':False}
    report['passed']=len(p)==152 and atoms.get_chemical_formula()=='C48H88N8O8' and max(deviations)<1e-7 and not clashes and max(residuals)<1e-7
    return report

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--max-candidates',type=int,default=2000)
    args=parser.parse_args()
    if args.max_candidates < 2:
        parser.error('--max-candidates must allow both independent templates (at least 2)')
    inp=ROOT/'inputs/allowed_literature_inputs.json'
    meta=json.loads(inp.read_text(encoding='utf-8-sig'))
    assert meta['chemistry']['repeat_formula']=='C6H11NO'
    assert meta['symmetry_hypothesis']['operations']==['x,y,z','-x+1/2,-y+1/2,z+1/2','-x,-y,-z','x+1/2,y+1/2,-z+1/2']
    assert meta['cell']['angles_deg'][:2]==[90,90]
    template_dir=ROOT/'inputs/independent_template'; out=ROOT/'branches/O/generated_v1_2_no_hbond'
    if out.exists() and any(out.iterdir()):
        raise RuntimeError('Existing versioned output: refusing to overwrite a generation record')
    selection_path=ROOT/'branches/O/phase24_first_pass/selection.json'
    selection=json.loads(selection_path.read_text(encoding='utf-8-sig'))
    if selection['phase_count'] != 24 or selection['phase_fractional_c'] != (np.arange(24)/48).tolist():
        raise RuntimeError('First-pass phase selection differs from the approved 24-point grid')
    manifest_path=ROOT/'inputs/chemistry_template_manifest.json'
    manifest=json.loads(manifest_path.read_text(encoding='utf-8-sig'))
    for source in manifest['source_hashes']:
        if sha(ROOT/source['path']) != source['sha256']:
            raise RuntimeError(f'Frozen independent template provenance changed: {source["path"]}')
    cell=cellpar_to_cell(meta['cell']['lengths_A']+meta['cell']['angles_deg'])
    protocol={'generator_version':VERSION,'seed':0,'stochastic':False,'branch':'O',
              'literature_input_sha256':sha(inp),'generator_sha256':sha(__file__),
              'source_coordinates_read':False,'coordinate_derived_geometry_read':False,
              'template_policy':'reuse frozen independent templates without reconstructing or overwriting them',
              'template_manifest_sha256':sha(manifest_path),
              'phase_selection_sha256':sha(selection_path),
              'physical_screen_policy':'inputs/candidate_physical_screen_policy.json',
              'physical_screen_policy_sha256':sha(ROOT/'inputs/candidate_physical_screen_policy.json'),
              'hydrogen_bonds_evaluated':False,'hydrogen_bonds_used_for_admission':False,
              'hydrogen_bonds_used_for_ranking_or_sampling':False,
              'chain_x_fractional':[.25,.225,.275,.20,.30],
              'chain_position_center_rationale':'x=.25 gives uniform same-sheet chain-axis spacing a/2 under the declared group; derived from cell/group, not atomic coordinates',
              'chain_rotation_about_c_deg':[0.,-10.,10.],
              'chain_c_phase_fractional':selection['phase_fractional_c'],
              'amide_angle_bias_deg':[0.,1.5],
              'max_retained_candidates':args.max_candidates,
              'ordering':'round robin between internal templates; x then rotation then phase',
              'selection':'first geometry-valid candidates by declared deterministic enumeration, not energy pruning',
              'coverage_status':'finite_declared_grid_only_not_global_structure_search',
              'stop_rule':'stop after requested retained candidates; retain at least one per internal template',
              'geometry_thresholds':{'nonbonded_screen':'unchanged v1.1 element-pair overlap lower bounds; excludes 1-2 and 1-3; no H-bond classification',
                                     'interpretation':'operational construction filters, not energetic stability or DFT validation'}}
    save(out/'generation_protocol.json',protocol) # written before all geometry generation
    (out/'build_ordered_candidates_v1_2_source.py').write_bytes(Path(__file__).read_bytes())
    templates=[]
    for j in range(len(protocol['amide_angle_bias_deg'])):
        template=json.loads((template_dir/f'chain_{j:02d}.json').read_text(encoding='utf-8-sig'))
        if not np.allclose(template['period_vector_A'],cell[2],rtol=0,atol=1e-10):
            raise RuntimeError('Frozen template period does not match the declared cell')
        templates.append((np.asarray(template['positions_unwrapped_A']),template['atom_roles'],template['bonds'],template))
    counts={'generated':0,'detailed_geometry_failed':0,'retained':0}
    candidate_ids=[]; per_template=[0]*len(templates); failures=[]; started=time.time()
    for x,phi,phase in itertools.product(protocol['chain_x_fractional'],protocol['chain_rotation_about_c_deg'],protocol['chain_c_phase_fractional']):
        for tj,(p,roles,bonds,template) in enumerate(templates):
            if per_template[tj]>=int(np.ceil(args.max_candidates/len(templates))):continue
            counts['generated']+=1
            data=expand(p,roles,bonds,template['symbols'],cell,x,phi,phase)
            param=dict(data['generation_parameters'],template_id=template['template_id'])
            report=validate_geometry(data)
            if not report['passed']:
                counts['detailed_geometry_failed']+=1
                failures.append({'parameters':param,'reason':'detailed_geometry','clashes':report['clashes']})
                continue
            counts['retained']+=1; per_template[tj]+=1
            cid=f'O-v12-c{counts["retained"]:04d}'
            data.update(candidate_id=cid,branch_id='O',generator_version=VERSION,
                        random_seed=0,chemistry_template_id=template['template_id'],
                        input_hashes={'allowed_inputs':sha(inp),'template':sha(template_dir/f'chain_{tj:02d}.json'),
                                      'protocol':sha(out/'generation_protocol.json'),'generator':sha(__file__)},
                        initial_symmetry_hypothesis='P1121/n No.14 c_unique',
                        geometry_status='builder_checks_passed_independent_validation_pending',
                        numerical_status='not_started',coverage_status='finite_declared_grid_only',experimental_status='not_compared',
                        geometry_report=report)
            save(out/f'{cid}.json',data)
            atoms=Atoms(data['symbols'],positions=data['positions'],cell=cell,pbc=True)
            atoms.info.update(candidate_id=cid,branch_id='O',template_id=template['template_id'])
            write(out/f'{cid}.extxyz',atoms)
            candidate_ids.append(cid)
            print(json.dumps({'candidate_id':cid,'parameters':param,'spglib':report['spglib_international']},ensure_ascii=False),flush=True)
        if counts['retained']>=args.max_candidates:break
    save(out/'generation_failures.json',failures)
    manifest={'counts':counts,'per_template_retained':per_template,'candidates':candidate_ids,
              'elapsed_s':time.time()-started,'protocol_sha256':sha(out/'generation_protocol.json'),
              'full_declared_grid_exhausted':counts['generated']==len(templates)*len(protocol['chain_x_fractional'])*len(protocol['chain_rotation_about_c_deg'])*len(protocol['chain_c_phase_fractional']),
              'status':'geometry_seeds_ready' if counts['retained'] else 'no_valid_seed',
              'independent_validator_status':'pending','energy_evaluations':0,'DFT_runs':0}
    save(out/'generation_manifest.json',manifest)
    print(json.dumps(manifest,indent=2),flush=True)

if __name__=='__main__':main()
