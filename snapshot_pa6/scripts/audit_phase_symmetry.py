"""Read-only geometry audit of phase equivalence, no energy/relaxation.

Rebuilds the frozen finite grid from frozen template positions and group ops.
Every accepted equivalence uses ONE global translation, a role-preserving
bijective atomic map, and the complete periodic covalent graph.
"""
from pathlib import Path
import json, hashlib, itertools, time, sys
import numpy as np
from scipy.optimize import linear_sum_assignment

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from validate_ordered import minimum_images, canonical_edge
R=np.array([np.eye(3),np.diag([-1,-1,1]),-np.eye(3),np.diag([1,1,-1])],int)
T=np.array([[0,0,0],[.5,.5,.5],[0,0,0],[.5,.5,.5]])
OUT=ROOT/'branches/O/audit_20261029/phase'
OUT.mkdir(parents=True,exist_ok=True)

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(name,obj):
    (OUT/name).write_text(json.dumps(obj,indent=2)+'\n')
def load(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))

base=load(ROOT/'branches/O/generated/O-c0010.json')
CELL=np.array(base['cell']); INV=np.linalg.inv(CELL)
SVMIN=np.linalg.svd(CELL,compute_uv=False).min()
templates=[load(ROOT/f'inputs/independent_template/chain_{j:02d}.json') for j in range(2)]

def make(tid,x,phi,phase):
    tp=templates[tid]; p=np.array(tp['positions_unwrapped_A'])
    a=np.deg2rad(phi); rot=np.array([[np.cos(a),-np.sin(a),0],[np.sin(a),np.cos(a),0],[0,0,1.]])
    f=p@rot.T@INV+[x,0,phase]
    ff=[]; rr=[]; edges=[]
    for g,(r,t) in enumerate(zip(R,T)):
        z=f@r.T+t; w=np.floor(z+1e-12).astype(int); ff.extend(z-w)
        rr.extend([z['role'].removesuffix('_a').removesuffix('_b') for z in tp['atom_roles']])
        for b in tp['bonds']:
            sh=np.array(b['shift'])@r.T+w[b['j']]-w[b['i']]
            edges.append((*canonical_edge(b['i']+38*g,b['j']+38*g,sh),float(b['order'])))
    return {'f':np.array(ff),'roles':np.array(rr),'edges':set(edges),'params':(tid,x,phi,phase)}

def match(a,b,tol=1e-5,origins=None):
    fa,fb=a['f'],b['f']; ra,rb=a['roles'],b['roles']
    # To reject efficiently, periodic fractional minimum bounds Cartesian MIC.
    if origins is None:
        ia=np.flatnonzero(ra=='amide_N')[0]
        origins=[fb[j]-fa[ia] for j in np.flatnonzero(rb=='amide_N')]
    for delta in origins:
        mapping=np.empty(len(fa),int); possible=True
        for role in ['amide_N','carbonyl_O']+[r for r in sorted(set(ra)) if r not in ('amide_N','carbonyl_O')]:
            ii=np.flatnonzero(ra==role); jj=np.flatnonzero(rb==role)
            df=fa[ii,None,:]+delta-fb[None,jj,:]; df-=np.rint(df)
            d=np.linalg.norm(df,axis=2)*SVMIN
            if np.max(d.min(axis=1))>tol or np.max(d.min(axis=0))>tol:
                possible=False;break
            # For pairs that survive this very small tolerance, fractional
            # image is unique; full exact MIC is separately verified below.
            rows,cols=linear_sum_assignment(np.linalg.norm(df@CELL,axis=2))
            mapping[ii[rows]]=jj[cols]
        if not possible:continue
        dif=fa+delta-fb[mapping]
        vec,shifts=minimum_images(dif,CELL)
        residual=np.linalg.norm(vec,axis=1)
        if residual.max()>tol:continue
        images=-shifts
        mapped=set()
        for i,j,sx,sy,sz,order in a['edges']:
            image=np.array([sx,sy,sz])+images[j]-images[i]
            mapped.add((*canonical_edge(mapping[i],mapping[j],image),order))
        if mapped!=b['edges']:continue
        return {'equivalent':True,'global_translation_fractional':np.asarray(delta).tolist(),
                'atom_mapping':mapping.tolist(),'atomic_lattice_images':images.tolist(),
                'max_residual_A':float(residual.max()),'rms_residual_A':float(np.sqrt(np.mean(residual**2))),
                'role_map_passed':True,'periodic_covalent_graph_passed':True}
    return {'equivalent':False}

def main():
    start=time.time(); xs=[.2,.225,.25,.275,.3]; angles=[-10.,0.,10.]
    within=[]; extra=[]; allrecords={}; controls=[]
    for tid,x,angle in itertools.product(range(2),xs,angles):
        records=[make(tid,x,angle,k/96) for k in range(48)]
        for k,rec in enumerate(records):allrecords[(tid,x,angle,k)]=rec
        pairs=[]
        for k,l in itertools.combinations(range(48),2):
            result=match(records[k],records[l])
            if result['equivalent']:pairs.append({'k':k,'l':l,**result})
        within.append({'template':tid,'x':x,'angle':angle,'phase_count':48,
                       'equivalent_pairs_count':len(pairs),'equivalent_pairs':pairs})
        print(json.dumps({'same_parameter':(tid,x,angle),'equivalent_phase_pairs':len(pairs)}),flush=True)
    save('within_phase_pairs.json',within)
    for tid,x,angle,k in itertools.product(range(2),[.2,.225,.25],angles,range(48)):
        otherx=round(.5-x,3)
        a=allrecords[(tid,x,angle,k)];b=allrecords[(tid,otherx,angle,k)]
        result=match(a,b,origins=[np.array([0,.5,0])])
        extra.append({'relation':'x_to_half_minus_x','template':tid,'x':x,'other_x':otherx,'angle':angle,'phase_index':k,**result})
    save('x_mirror_equivalences.json',extra)
    # Cross-angle comparisons include every phase pair, not only same phase.
    anglepairs=[]
    for tid,x in itertools.product(range(2),xs):
        found=[]
        for k,l in itertools.product(range(48),repeat=2):
            result=match(allrecords[(tid,x,-10.,k)],allrecords[(tid,x,10.,l)])
            if result['equivalent']:found.append({'k':k,'l':l,**result})
        anglepairs.append({'template':tid,'x':x,'equivalent_pairs':found,'count':len(found)})
    save('opposite_angle_pairs.json',anglepairs)
    for tid,x,angle,k in [(0,.25,0.,17),(0,.225,-10.,11),(1,.2,10.,23)]:
        a=allrecords[(tid,x,angle,k)]
        b=make(tid,x,angle,k/96+.5)
        result=match(a,b,origins=[np.array([0,0,.5])])
        controls.append({'relation':'phase_plus_half','parameters':(tid,x,angle,k),**result})
    # Tolerance sensitivity of approximate template screw-induced matches.
    sensitivity=[]
    for tid,angle in itertools.product(range(2),angles):
        a=allrecords[(tid,.225,angle,17)];b=allrecords[(tid,.275,angle,17)]
        sensitivity.append({'template':tid,'angle':angle,'tests':{str(t):match(a,b,tol=t,origins=[np.array([0,.5,0])]) for t in [1e-5,1e-6,1e-7,1e-8,1e-9]}})
    save('positive_controls.json',controls);save('tolerance_sensitivity.json',sensitivity)
    summary={'scope':'same-cell global translation + chemical-role bijection + periodic covalent graph; no optimization',
             'input_grid':1440,'within_parameter_phase_tests':30*1128,
             'within_parameter_equivalent_phase_pairs':sum(z['equivalent_pairs_count'] for z in within),
             'x_mirror_tests':len(extra),'x_mirror_passed':sum(z['equivalent'] for z in extra),
             'cross_opposite_angle_phase_tests':10*48*48,'cross_opposite_angle_equivalences':sum(z['count'] for z in anglepairs),
             'phase_plus_half_controls_passed':all(z['equivalent'] for z in controls),
             'mapping_tolerance_A':1e-5,'elapsed_seconds':time.time()-start,
             'hashes':{'audit_script':sha(__file__),'template0':sha(ROOT/'inputs/independent_template/chain_00.json'),'template1':sha(ROOT/'inputs/independent_template/chain_01.json')},
             'limits':['No arbitrary change-of-cell or supercell matching','No energy or optimization-basin equivalence','x-mirror reduction uses the current template near-screw symmetry; examine tighter tolerance results','Same declared point rotations cannot add equivalences beyond translations because both structures already possess these operations']}
    save('phase_symmetry_summary.json',summary);print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
