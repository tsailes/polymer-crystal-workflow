"""Audit the frozen 1440-point grid without optimizing or evaluating energy.

Original files are read-only. Full original builder thresholds and independent
chemical/topological screens are evaluated before considering hydrogen bonds.
"""
from __future__ import annotations
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import argparse, hashlib, itertools, json, os, time
from pathlib import Path
import numpy as np
from ase.geometry import cellpar_to_cell
import build_ordered_candidates as b
import validate_ordered as v

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'branches/O/audit_20261029/hbond'
OLD_BOUNDS={'CC':2.45,'CH':1.6,'CN':2.3,'CO':2.25,'HH':1.25,
            'HN':1.3,'HO':1.2,'NN':2.3,'NO':2.2,'OO':2.3}
DATA=None

def read(p): return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d): Path(p).write_text(json.dumps(d,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')

def init_worker():
    global DATA
    protocol=read(ROOT/'branches/O/generated/generation_protocol.json')
    meta=read(ROOT/'inputs/allowed_literature_inputs.json')
    cell=cellpar_to_cell(meta['cell']['lengths_A']+meta['cell']['angles_deg'])
    templates=[read(ROOT/f'inputs/independent_template/chain_{i:02d}.json') for i in range(2)]
    DATA=(protocol,cell,templates)

def contact_geometry(record):
    """Record all relevant N-H...O images within finite exact search bounds.

    Here the period cell is fixed and well-conditioned; use the same explicit
    [-1,1]^3 images as frozen builder so its selection is reproduced verbatim.
    Near-miss labels are observational, never admission rules.
    """
    pos=np.asarray(record['positions']); cell=np.asarray(record['cell'])
    roles=record['atom_roles']; bonds=record['bonds']
    ns=[x['index'] for x in roles if x['role']=='amide_N']
    os_=[x['index'] for x in roles if x['role']=='carbonyl_O']
    contacts=[]
    for n in ns:
        nh=[]
        for edge in bonds:
            i,j=edge['i'],edge['j']; s=np.array(edge['shift'])
            if i==n and roles[j]['role']=='amide_H': nh=[j,s]; break
            if j==n and roles[i]['role']=='amide_H': nh=[i,-s]; break
        h,s=nh; hp=pos[h]+s@cell; nv=pos[n]-hp
        possible=[]
        for o in os_:
            if roles[o]['chain_id']==roles[n]['chain_id']: continue
            op=pos[o]+b.SHIFTS@cell; hv=op-hp; no=np.linalg.norm(op-pos[n],axis=1)
            ho=np.linalg.norm(hv,axis=1)
            angle=np.degrees(np.arccos(np.clip((hv@nv)/(np.maximum(ho,1e-12)*np.linalg.norm(nv)),-1,1)))
            # Keep nearest useful candidate; positive excess expresses proximity
            # to original cutoffs in A and degrees, never an energy/probability.
            for j in np.flatnonzero(no<5.0):
                merit=max(0,no[j]-3.5)+max(0,ho[j]-2.6)+max(0,120-angle[j])/60
                possible.append({'N':n,'H':h,'O':o,'O_shift':b.SHIFTS[j].tolist(),
                    'NO_A':float(no[j]),'HO_A':float(ho[j]),'NHO_deg':float(angle[j]),
                    'geometric_cutoff_excess_score':float(merit)})
        if possible: contacts.append(min(possible,key=lambda x:x['geometric_cutoff_excess_score']))
    return contacts

def evaluate(item):
    ident,x,rot,phase,tidx=item
    _,cell,templates=DATA; template=templates[tidx]
    record=b.expand(np.array(template['positions_unwrapped_A']),template['atom_roles'],template['bonds'],
        template['symbols'],cell,x,rot,phase)
    record.update(candidate_id=ident,branch_id='O',chemistry_template_id=template['template_id'])
    cell,symbols,pos,frac=v.structure_arrays(record)
    edges,adj,errors,pairs=v.infer_graph(symbols,frac,cell)
    chem=v.check_chemistry(symbols,adj)
    val=v.valence_geometry(symbols,frac,cell,adj)
    topo,membership=v.periodic_components(adj)
    roles=v.validate_roles(record,symbols,adj,topo)
    declared={v.canonical_edge(e['i'],e['j'],e['shift']) for e in record['bonds']}
    inferred={v.canonical_edge(e['i'],e['j'],e['shift']) for e in edges}
    graph_ok=declared==inferred and len(declared)==len(record['bonds'])
    d_adj=[set() for _ in symbols]
    for e in record['bonds']: d_adj[e['i']].add(e['j']);d_adj[e['j']].add(e['i'])
    excluded=[]; inferred_excluded=[]
    for i in range(len(symbols)):
        near=set(d_adj[i]);near.update(k for j in d_adj[i] for k in d_adj[j]);excluded.append(near)
        near={j for j,_ in adj[i]};near.update(k for j,_ in adj[i] for k,_ in adj[j]);inferred_excluded.append(near)
    old_clashes=[]; new_clashes=[]
    for i,j,d in zip(*pairs):
        pair=''.join(sorted((symbols[i],symbols[j])))
        if j not in excluded[i] and d<OLD_BOUNDS[pair]: old_clashes.append({'i':int(i),'j':int(j),'pair':pair,'distance_A':float(d),'floor_A':OLD_BOUNDS[pair]})
        if j not in inferred_excluded[i] and d<v.COLLISION_FLOORS[pair]: new_clashes.append({'i':int(i),'j':int(j),'pair':pair,'distance_A':float(d),'floor_A':v.COLLISION_FLOORS[pair]})
    bond_error=max(abs(np.linalg.norm(pos[e['j']]+np.array(e['shift'])@cell-pos[e['i']])-e['equilibrium_A']) for e in record['bonds'])
    # Fast exact orbit check: component rounding is NOT used for distances;
    # it merely gives candidate images, which have negligible Cartesian residual.
    symmax=0.
    for R,t in zip(b.OPS_R,b.OPS_T):
        tf=frac@R.T+t
        for element in ('C','H','N','O'):
            ids=np.flatnonzero(symbols==element)
            delta=tf[ids,None,:]-frac[None,ids,:]
            delta-=np.rint(delta)
            distance=np.linalg.norm(delta@cell,axis=2)
            rows,cols=b.linear_sum_assignment(distance)
            symmax=max(symmax,float(distance[rows,cols].max()))
    symok=symmax<1e-7
    hb=b.hydrogen_bonds(record)
    full=hb['donors_covered']==8 and hb['acceptors_covered']==8
    # HB graph is always defined on DECLARED chemistry; malformed overlapping
    # geometries are still diagnostic records, never accepted through HB alone.
    declared_adj=[[] for _ in symbols]
    for e in record['bonds']:
        s=np.array(e['shift']);declared_adj[e['i']].append((e['j'],s));declared_adj[e['j']].append((e['i'],-s))
    hp=[]
    for hit in hb['bonds']:
        n,h,o=hit['N'],hit['H'],hit['O']
        nhshift=next(s for j,s in declared_adj[n] if j==h)
        hp.append({'H':h,'O':o,'oxygen_image_from_H':(np.array(hit['O_shift'])-nhshift).tolist()})
    sheets=v.hydrogen_bond_sheet_topology(declared_adj,{'bonds':hp})
    chemistry_ok=chem['passed'] and val['passed'] and topo['passed'] and roles['passed'] and graph_ok and not errors
    indep_ok=chemistry_ok and not new_clashes and symok
    original_geometry_ok=bond_error<1e-7 and not old_clashes and symok
    closest=contact_geometry(record)
    near_01=all(c['NO_A']<=3.6 and c['HO_A']<=2.7 and c['NHO_deg']>=110 for c in closest) and len(closest)==8
    near_03=all(c['NO_A']<=3.8 and c['HO_A']<=2.9 and c['NHO_deg']>=100 for c in closest) and len(closest)==8
    return {'audit_id':ident,'parameters':{'x':x,'rotation_deg':rot,'phase_c':phase,'template':tidx},
        'original_full_HB_coverage':full,'donors_covered':hb['donors_covered'],'acceptors_covered':hb['acceptors_covered'],
        'HB_edge_count':len(hb['bonds']),'infinite_ac_sheet':sheets['passed'],
        'sheet_component_ranks':[c['winding_rank'] for c in sheets['components']],
        'original_builder_geometry_without_HB_passed':original_geometry_ok,
        'independent_geometry_without_HB_passed':indep_ok,
        'independent_chemistry_topology_valence_passed':chemistry_ok,
        'graph_agreement':graph_ok,'bond_error_max_A':bond_error,'symmetry_max_A':symmax,
        'original_collision_count':len(old_clashes),'independent_collision_count':len(new_clashes),
        'original_worst_collision':min(old_clashes,key=lambda k:k['distance_A']-k['floor_A']) if old_clashes else None,
        'independent_worst_collision':min(new_clashes,key=lambda k:k['distance_A']-k['floor_A']) if new_clashes else None,
        'closest_8_donor_contacts':closest,'all_donors_near_cutoffs_plus_0p1A_minus10deg':near_01,
        'all_donors_near_cutoffs_plus_0p3A_minus20deg':near_03,
        'chemistry_errors':chem['errors'][:5],'inferred_extra_edges':len(inferred-declared)}

def summarize(records):
    total=len(records);oldfull=[r for r in records if r['original_full_HB_coverage']]
    rejected=[r for r in records if not r['original_full_HB_coverage']]
    def counts(rs):
        return {'count':len(rs),'original_builder_geometry_without_HB_passed':sum(r['original_builder_geometry_without_HB_passed'] for r in rs),
        'independent_geometry_without_HB_passed':sum(r['independent_geometry_without_HB_passed'] for r in rs),
        'both_geometry_screens_passed':sum(r['original_builder_geometry_without_HB_passed'] and r['independent_geometry_without_HB_passed'] for r in rs),
        'independent_chemistry_topology_valence_passed':sum(r['independent_chemistry_topology_valence_passed'] for r in rs),
        'donor_coverage_counts':dict(sorted(Counter(r['donors_covered'] for r in rs).items())),
        'infinite_ac_sheet_count':sum(r['infinite_ac_sheet'] for r in rs)}
    byx={str(x):counts([r for r in records if r['parameters']['x']==x]) for x in (.2,.225,.25,.275,.3)}
    legal_missed=[r for r in rejected if r['original_builder_geometry_without_HB_passed'] and r['independent_geometry_without_HB_passed']]
    examples=sorted(legal_missed,key=lambda r:(-r['donors_covered'],sum(c['geometric_cutoff_excess_score'] for c in r['closest_8_donor_contacts'])))[:12]
    return {'schema_version':1,'scope':'Frozen original grid only; no structural optimization, energy, GFN2 or DFT',
        'all':counts(records),'old_complete_HB_subset':counts(oldfull),'old_rejected_HB_subset':counts(rejected),
        'by_x':byx,
        'original_geometry_passed_HB_classes':dict(sorted(Counter(r['donors_covered'] for r in records if r['original_builder_geometry_without_HB_passed']).items())),
        'independent_geometry_passed_HB_classes':dict(sorted(Counter(r['donors_covered'] for r in records if r['independent_geometry_without_HB_passed']).items())),
        'missed_both_geometry_passed_near_0p1A_minus10deg':sum(r['all_donors_near_cutoffs_plus_0p1A_minus10deg'] for r in legal_missed),
        'missed_both_geometry_passed_near_0p3A_minus20deg':sum(r['all_donors_near_cutoffs_plus_0p3A_minus20deg'] for r in legal_missed),
        'examples':examples,'limitations':['No prediction of relaxation basin or proof of eventual sheet formation',
        'Both collision screens are operational starting-geometry bounds; differing thresholds are reported separately',
        'Near-cutoff score is geometric only and is not energy or optimization probability',
        'All structures reconstructed directly from frozen template geometry; no fresh template solving'],
        'input_hashes':{str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'branches/O/generated/generation_protocol.json',ROOT/'inputs/independent_template/chain_00.json',ROOT/'inputs/independent_template/chain_01.json',Path(b.__file__),Path(v.__file__),Path(__file__)]}}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=4);ap.add_argument('--limit',type=int,default=0);args=ap.parse_args()
    OUT.mkdir(parents=True,exist_ok=True);init_worker();protocol,_,_=DATA
    items=[]
    for x,rot,phase in itertools.product(protocol['chain_x_fractional'],protocol['chain_rotation_about_c_deg'],protocol['chain_c_phase_fractional']):
        for ti in range(2): items.append((f'grid-{len(items)+1:04d}',x,rot,phase,ti))
    if args.limit:items=items[:args.limit]
    started=time.time();records=[]
    with ProcessPoolExecutor(max_workers=args.workers,initializer=init_worker) as pool:
        for result in pool.map(evaluate,items,chunksize=8):
            records.append(result)
            if len(records)%100==0:print(json.dumps({'done':len(records),'total':len(items),'elapsed_s':time.time()-started}),flush=True)
    save(OUT/'grid_records.json',records);summary=summarize(records);summary['elapsed_s']=time.time()-started
    save(OUT/'summary.json',summary)
    print(json.dumps({k:summary[k] for k in ('all','old_complete_HB_subset','old_rejected_HB_subset','missed_both_geometry_passed_near_0p1A_minus10deg','missed_both_geometry_passed_near_0p3A_minus20deg')},indent=2),flush=True)

if __name__=='__main__':main()
