"""Read-only PA6 post-run H-bond, periodic-image and chain-registry diagnostics.

No energy calculation or candidate filtering. No reference crystal coordinates.
All output artifacts are confined to geometry_* files in this directory.
"""
from __future__ import annotations
import argparse
from collections import defaultdict, deque
import hashlib
import itertools
import json
from pathlib import Path
import sys
import numpy as np
from ase.io import read
from ase.io.trajectory import Trajectory

OUT = Path(__file__).resolve().parent
RUN = OUT.parent
PROJECT = RUN.parents[1]
sys.path.insert(0, str(PROJECT / 'scripts'))
import validate_ordered as v


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def jdefault(obj):
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(str(type(obj)))


def graph_edges(record):
    adjacency = [[] for _ in record['symbols']]
    for b in record['bonds']:
        i,j,s = b['i'],b['j'],np.asarray(b['shift'], int)
        adjacency[i].append((j,s)); adjacency[j].append((i,-s))
    return adjacency


def unwrapping_offsets(record, adjacency):
    """BFS branch cut retained for every frame, only used for chain centroids."""
    offsets = np.zeros((len(record['symbols']),3),int)
    groups = defaultdict(list)
    for role in record['atom_roles']:
        groups[role['chain_id']].append(role['index'])
    for indices in groups.values():
        root=min(indices); seen={root}; queue=deque([root])
        while queue:
            i=queue.popleft()
            for j,shift in adjacency[i]:
                if j not in seen:
                    offsets[j]=offsets[i]+shift; seen.add(j); queue.append(j)
    return offsets,groups


def angle_deg(nh,ho):
    return float(np.degrees(np.arccos(np.clip(np.dot(-nh,ho)/(np.linalg.norm(nh)*np.linalg.norm(ho)),-1,1))))


def contact(n,h,o,nh_shift,oi,p,cell,roles):
    nh = p[h]+nh_shift@cell-p[n]
    no = p[o]+oi@cell-p[n]
    ho = no-nh
    angle=angle_deg(nh,ho)
    return {'N':n,'H':h,'O':o,'donor_chain':roles[n]['chain_id'],
            'acceptor_chain':roles[o]['chain_id'],
            'donor_NH_image_shift':nh_shift.tolist(), 'oxygen_image_from_N':oi.tolist(),
            'N_O_A':float(np.linalg.norm(no)),'H_O_A':float(np.linalg.norm(ho)),
            'N_H_O_deg':angle,'donor_NH_unit':(nh/np.linalg.norm(nh)).tolist(),
            'H_to_O_unit':(ho/np.linalg.norm(ho)).tolist(),
            'N_to_O_vector_A':no.tolist(),'H_to_O_vector_A':ho.tolist(),
            'direction_pass_120deg':angle>=120,
            'distance_pass_3p5_2p6':bool(np.linalg.norm(no)<=3.5 and np.linalg.norm(ho)<=2.6),
            'strict_geometric_contact':bool(angle>=120 and np.linalg.norm(no)<=3.5 and np.linalg.norm(ho)<=2.6)}


def all_contacts(p,cell,record,adjacency,radius=6.0,extra=0):
    symbols=record['symbols']; roles={a['index']:a for a in record['atom_roles']}
    frac=p@np.linalg.inv(cell)
    # Include the union of spheres NO <= 6 and HO <= 6; NH < 1.3 A.
    search_r=radius+1.3
    reach=np.ceil(search_r*np.linalg.norm(np.linalg.inv(cell),axis=0)+.5).astype(int)+extra
    offsets=np.asarray(list(itertools.product(*(range(-int(k),int(k)+1) for k in reach))),int)
    rows=[]
    for n,s in enumerate(symbols):
        if s!='N': continue
        for h,nh_shift in adjacency[n]:
            if symbols[h]!='H': continue
            contacts=[]
            for o,os in enumerate(symbols):
                if os!='O' or roles[o]['chain_id']==roles[n]['chain_id']:continue
                base=-np.rint(frac[o]-frac[n]).astype(int)
                for oi in base+offsets:
                    no=p[o]+oi@cell-p[n]
                    ho=p[o]+oi@cell-(p[h]+nh_shift@cell)
                    if np.linalg.norm(no)>radius and np.linalg.norm(ho)>radius:continue
                    contacts.append(contact(n,h,o,nh_shift,oi,p,cell,roles))
            contacts.sort(key=lambda q:(q['N_O_A'],q['O'],q['oxygen_image_from_N']))
            distance_pass=[q for q in contacts if q['distance_pass_3p5_2p6']]
            direction_pass=[q for q in contacts if q['direction_pass_120deg']]
            rows.append({'N':n,'H':h,'donor_chain':roles[n]['chain_id'],
                         'nearest_by_NO':min(contacts,key=lambda q:q['N_O_A']),
                         'nearest_by_HO':min(contacts,key=lambda q:q['H_O_A']),
                         'nearest_angle_qualified':min(direction_pass,key=lambda q:q['H_O_A']) if direction_pass else None,
                         'best_angle_distance_qualified':max(distance_pass,key=lambda q:q['N_H_O_deg']) if distance_pass else None,
                         'all_cross_chain_O_images_within_6A_NO_or_HO':contacts,
                         'strict_contact_count':sum(q['strict_geometric_contact'] for q in contacts),
                         'angle_pass_distance_fail_count':sum(q['direction_pass_120deg'] and not q['distance_pass_3p5_2p6'] for q in contacts),
                         'distance_pass_angle_fail_count':sum(q['distance_pass_3p5_2p6'] and not q['direction_pass_120deg'] for q in contacts),
                         'neither_pass_count':sum(not q['distance_pass_3p5_2p6'] and not q['direction_pass_120deg'] for q in contacts)})
    return rows


def chain_geometry(atoms,record,offsets,groups,initial):
    p=atoms.positions+offsets@atoms.cell.array
    p0=initial.positions+offsets@initial.cell.array
    inv=np.linalg.inv(atoms.cell.array)
    roles={r['index']:r for r in record['atom_roles']}
    normals=np.cross(atoms.cell.array[2],atoms.cell.array[0]);normals/=np.linalg.norm(normals)
    chain_rows=[]
    for label,indices in groups.items():
        ids=np.asarray(indices); masses=atoms.get_masses()[ids]
        center=np.average(p[ids],axis=0,weights=masses)
        center0=np.average(p0[ids],axis=0,weights=masses)
        backbone=[i for i in ids if roles[i]['role']=='amide_N' or roles[i]['role']=='carbonyl_C' or roles[i]['role'].startswith('methylene_C')]
        bc=p[backbone].mean(axis=0); bc0=p0[backbone].mean(axis=0)
        amides=[]
        adj=graph_edges(record)
        for repeat in (0,1):
            aids=[i for i in ids if roles[i]['repeat_index']==repeat and roles[i]['role'] in ('amide_N','carbonyl_C','carbonyl_O')]
            ni=next(i for i in aids if roles[i]['role']=='amide_N')
            ci=next(i for i in aids if roles[i]['role']=='carbonyl_C')
            oi=next(i for i in aids if roles[i]['role']=='carbonyl_O')
            sc=next(s for j,s in adj[ni] if j==ci)
            so=sc+next(s for j,s in adj[ci] if j==oi)
            # Keep each amide coherent even if the polymer BFS cut falls nearby.
            local=np.array([atoms.positions[ni],atoms.positions[ci]+sc@atoms.cell.array,
                            atoms.positions[oi]+so@atoms.cell.array])
            ac=local.mean(axis=0)
            amides.append({'repeat':repeat,'atom_indices':aids,'center_unwrapped_A':ac.tolist(),
                           'c_phase_fractional_mod1':float((ac@inv)[2]%1),
                           'amide_N_index':int(ni),
                           'amide_N_c_phase_fractional_mod1':float((atoms.positions[ni]@inv)[2]%1),
                           'amide_N_delta_c_A':float(atoms.positions[ni,2]-initial.positions[ni,2]),
                           'mean_delta_c_A':float((p[aids]-p0[aids]).mean(axis=0)[2])})
        chain_rows.append({'chain':label,'atom_indices':indices,
                           'center_of_mass_unwrapped_A':center.tolist(),
                           'center_of_mass_delta_A':(center-center0).tolist(),
                           'backbone_centroid_unwrapped_A':bc.tolist(),
                           'backbone_centroid_fractional_mod1':((bc@inv)%1).tolist(),
                           'backbone_mean_delta_A':(bc-bc0).tolist(),
                           'collective_c_phase_change_fraction':float((bc-bc0)[2]/atoms.cell.lengths()[2]),
                           'normal_to_ac_centroid_A':float(bc@normals),
                           'backbone_normal_to_ac_range_A':[float((p[backbone]@normals).min()),float((p[backbone]@normals).max())],
                           'amide_groups':amides})
    pairs=[]
    for a,b in itertools.combinations(chain_rows,2):
        delta=np.array(b['backbone_centroid_unwrapped_A'])-a['backbone_centroid_unwrapped_A']
        mic,shift=v.minimum_images(delta@inv,atoms.cell.array)
        pairs.append({'chain_a':a['chain'],'chain_b':b['chain'],
                      'centroid_minimum_image_vector_A':mic.tolist(),
                      'transverse_axis_distance_A':float(np.linalg.norm(mic[:2])),
                      'normal_to_ac_spacing_A':abs(float(mic@normals)),
                      'relative_c_centroid_phase_A':float(mic[2]),'image':shift.tolist()})
    return {'centroid_convention':'Fixed initial covalent-graph BFS branch cut; temporal changes are independent of that fixed cut. Absolute axial centroids depend on chosen chain segment; use locally connected amide groups and N anchors for registry.',
            'chain_axis':'c from periodic winding, not a fitted finite-segment direction',
            'normal_to_ac_unit':normals.tolist(),'chains':chain_rows,'chain_pairs':pairs}


def summary(rows):
    contacts=[q for r in rows for q in r['all_cross_chain_O_images_within_6A_NO_or_HO']]
    return {'donor_count':len(rows),'O_image_contact_count_within_6A_union':len(contacts),
            'strict_bond_count':sum(q['strict_geometric_contact'] for q in contacts),
            'strict_covered_donors':sum(r['strict_contact_count']>0 for r in rows),
            'angle_pass_distance_fail':sum(r['angle_pass_distance_fail_count'] for r in rows),
            'distance_pass_angle_fail':sum(r['distance_pass_angle_fail_count'] for r in rows),
            'neither_pass':sum(r['neither_pass_count'] for r in rows),
            'minimum_crosschain_NO_A':min(q['N_O_A'] for q in contacts),
            'minimum_crosschain_HO_A':min(q['H_O_A'] for q in contacts),
            'nearest_NO_angle_range_deg':[min(r['nearest_by_NO']['N_H_O_deg'] for r in rows),max(r['nearest_by_NO']['N_H_O_deg'] for r in rows)],
            'nearest_angle_qualified_HO_range_A':[min(r['nearest_angle_qualified']['H_O_A'] for r in rows if r['nearest_angle_qualified']), max(r['nearest_angle_qualified']['H_O_A'] for r in rows if r['nearest_angle_qualified'])] if any(r['nearest_angle_qualified'] for r in rows) else None}


def sensitivity(rows):
    contacts=[q for r in rows for q in r['all_cross_chain_O_images_within_6A_NO_or_HO']]
    table=[]
    for no,ho,angle in itertools.product((3.5,3.7,4.0,4.5),(2.6,2.8,3.0,3.5),(90,100,110,120,130)):
        hits=[q for q in contacts if q['N_O_A']<=no and q['H_O_A']<=ho and q['N_H_O_deg']>=angle]
        table.append({'max_NO_A':no,'max_HO_A':ho,'min_NHO_deg':angle,
                      'contact_count':len(hits),'covered_donors':len({q['N'] for q in hits}),
                      'covered_acceptors':len({q['O'] for q in hits})})
    return {'role':'Post-run sensitivity only, not candidate filtering or chemical proof', 'table':table}


def validate_detection(atoms,record,rows,adjacency):
    cell=atoms.cell.array;p=atoms.positions;frac=p@np.linalg.inv(cell);sym=np.array(record['symbols'])
    inferred,adj2,errors,_=v.infer_graph(sym,frac,cell)
    graph=lambda edges:{v.canonical_edge(q['i'],q['j'],q['shift']) for q in edges}
    topology,membership=v.periodic_components(adj2)
    legacy=v.hydrogen_bonds(sym,frac,cell,adj2,membership)
    enumerated_hits=[q for r in rows for q in r['all_cross_chain_O_images_within_6A_NO_or_HO'] if q['strict_geometric_contact']]
    keys=lambda rr:{(q['N'],q['H'],q['O'],*q['oxygen_image_from_N']) for q in rr}
    extended=all_contacts(p,cell,record,adjacency,extra=2)
    extended_contacts=[q for r in extended for q in r['all_cross_chain_O_images_within_6A_NO_or_HO']]
    regular_contacts=[q for r in rows for q in r['all_cross_chain_O_images_within_6A_NO_or_HO']]
    # A separate ASE MIC route verifies each closest H-O image without using v.minimum_images.
    from ase.geometry import find_mic
    nearest_errors=[]
    for row in rows:
        q=row['nearest_by_HO'];n,h,o=q['N'],q['H'],q['O']
        hp=p[h]+np.array(q['donor_NH_image_shift'])@cell
        _,d=find_mic(p[o]-hp,cell,pbc=True)
        nearest_errors.append(abs(float(d)-q['H_O_A']))
    return {'atom_count_and_order_matches_source':atoms.get_chemical_symbols()==record['symbols'],
            'inferred_periodic_graph_matches_declared':graph(inferred)==graph(record['bonds']),
            'inferred_bond_errors':errors,'four_primitive_c_chains':topology['passed'],
            'legacy_validator_strict_contact_count':len(legacy['bonds']),
            'all_image_enumeration_strict_contact_count':len(enumerated_hits),
            'legacy_and_all_image_strict_counts_agree':len(legacy['bonds'])==len(enumerated_hits),
            'expanded_image_window_same_all_6A_contacts':keys(regular_contacts)==keys(extended_contacts),
            'independent_ASE_MIC_max_HO_difference_A':max(nearest_errors),
            'angle_definition':'angle N-H-O at H = acos((-N_to_H) dot H_to_O / lengths); collinear N-H...O gives 180 degrees',
            'linear_angle_control_deg':angle_deg(np.array([1.,0,0]),np.array([1.,0,0])),
            'right_angle_control_deg':angle_deg(np.array([1.,0,0]),np.array([0.,1,0])),
            'backward_angle_control_deg':angle_deg(np.array([1.,0,0]),np.array([-1.,0,0]))}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--initial',type=Path,default=RUN/'input/candidate.extxyz')
    parser.add_argument('--final',type=Path,default=RUN/'relaxation/final.extxyz')
    parser.add_argument('--trajectory',type=Path,default=RUN/'relaxation/trajectory.traj')
    parser.add_argument('--output',type=Path,default=OUT/'geometry_diagnosis.json')
    args=parser.parse_args()
    record=load(RUN/'input/candidate.json');initial=read(args.initial)
    final=read(args.final);trajectory=Trajectory(args.trajectory)
    adjacency=graph_edges(record);offsets,groups=unwrapping_offsets(record,adjacency)
    selected=sorted({min(i,len(trajectory)-1) for i in (0,1,5,10,20,40,62)}|{len(trajectory)-1})
    frames=[('initial',initial)]+[(f'step_{i}',trajectory[i]) for i in selected]+[('final',final)]
    results={}
    for label,atoms in frames:
        rows=all_contacts(atoms.positions,atoms.cell.array,record,adjacency)
        results[label]={'donors':rows,'summary':summary(rows),
                        'chain_geometry':chain_geometry(atoms,record,offsets,groups,initial)}
        if label in ('initial','final'):
            results[label]['sensitivity']=sensitivity(rows)
            results[label]['detector_audit']=validate_detection(atoms,record,rows,adjacency)
    sources=[Path(__file__),RUN/'input/candidate.json',args.initial,args.final,args.trajectory,Path(v.__file__)]
    result={'scope':'Post-run geometry diagnosis; no quantum calculation, candidate filtering, literature atomic coordinates, or reference-derived geometry',
            'source_sha256':{str(p):sha(p) for p in sources},'atom_index_base':0,
            'image_convention':'N at stored coordinate; H = stored H + donor_NH_image_shift @ cell; O = stored O + oxygen_image_from_N @ cell',
            'neighborhood':'Every distinct cross-chain O image with N...O <= 6 A OR H...O <= 6 A; same-chain O excluded by source chain identity',
            'strict_geometry_thresholds':{'max_NO_A':3.5,'max_HO_A':2.6,'min_NHO_deg':120},
            'threshold_role':'Only diagnoses this already-selected trajectory; does not alter filtering, ranking or eligibility',
            'input_paths':{'initial':str(args.initial),'final':str(args.final),'trajectory':str(args.trajectory)},
            'trajectory_frame_count':len(trajectory),'frames':results}
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2,default=jdefault)+'\n',encoding='utf-8')
    print(json.dumps({k:r['summary'] for k,r in results.items()},indent=2,default=jdefault))


if __name__=='__main__': main()
