"""Small, predeclared rigid-chain direction probes; not new candidate screening."""
from pathlib import Path
import hashlib, json, sys, time
import numpy as np
from ase.io import read, write
from ase.constraints import FixSymmetry

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT.parents[1]/'scripts'))
from run_gfn2_single_pilot import evaluate, save, resources

atoms=read(ROOT/'relaxation/final.extxyz')
meta=json.loads((ROOT/'input/candidate.json').read_text())
roles=meta['atom_roles']
chain0=[r['index'] for r in roles if r['chain_id']=='O-chain-0']
constraint=FixSymmetry(atoms.copy(),symprec=1e-5,adjust_cell=False)
modes=[]
for axis,label in enumerate('xyz'):
    vector=np.zeros((len(atoms),3)); vector[chain0,axis]=1
    vector-=vector.mean(axis=0)
    allowed=vector.copy(); constraint.adjust_forces(atoms,allowed)
    forbidden=vector-allowed
    for kind,v in [('allowed',allowed),('forbidden',forbidden)]:
        norm=np.linalg.norm(v,axis=1).max()
        if norm<1e-8: continue
        v=v/norm
        again=v.copy(); constraint.adjust_forces(atoms,again)
        modes.append(dict(name=f'chain_{label}_{kind}',kind=kind,direction=v,
                          max_direction_norm=float(np.linalg.norm(v,axis=1).max()),
                          projected_norm_ratio=float(np.linalg.norm(again)/np.linalg.norm(v))))
work=HERE/'directional_probes'
work.mkdir(exist_ok=False)
save(work/'protocol.json',dict(input_sha256=hashlib.sha256((ROOT/'relaxation/final.extxyz').read_bytes()).hexdigest(),
     method='GFN2-xTB',accuracy=.01,repeat=[1,1,1],amplitudes_A=[.02,.04],
     mode_definition='single chain Cartesian translations, global translation removed, decomposed by symmetry projector; max atom norm normalized to one',
     selection='xyz axes and symmetry split declared before evaluations; no HB selection',modes=modes,
     limitations=['directional finite differences only, not full Hessian','base geometry has residual force','Gamma sampling not converged']))
base=evaluate(atoms,.01); save(work/'base.json',base)
results=[]; started=time.perf_counter()
for mode in modes:
    v=mode['direction']; pair=[]
    for h in [.02,.04]:
        evals=[]
        for sign in [-1,1]:
            moved=atoms.copy(); moved.positions+=sign*h*v
            tag=f"{mode['name']}_{h:.2f}_{'minus' if sign<0 else 'plus'}"
            write(work/f'{tag}.extxyz',moved)
            r=evaluate(moved,.01); save(work/f'{tag}.json',r)
            evals.append(r)
        curvature=(evals[1]['energy']+evals[0]['energy']-2*base['energy'])/h**2
        curvature_force=-float(np.sum((np.asarray(evals[1]['forces'])-np.asarray(evals[0]['forces']))*v))/(2*h)
        pair.append(dict(amplitude_A=h,minus_delta_eV=evals[0]['energy']-base['energy'],
            plus_delta_eV=evals[1]['energy']-base['energy'],energy_curvature_eV_A2=curvature,
            force_curvature_eV_A2=curvature_force,
            secant_gradient_eV_A=(evals[1]['energy']-evals[0]['energy'])/(2*h)))
    row=dict(name=mode['name'],kind=mode['kind'],projected_norm_ratio=mode['projected_norm_ratio'],
             base_gradient_eV_A=-float(np.sum(np.asarray(base['forces'])*v)),probes=pair)
    results.append(row)
    save(work/'summary.json',dict(status='in_progress',results=results,wall_seconds=time.perf_counter()-started,**resources()))
    print(json.dumps(row),flush=True)
save(work/'summary.json',dict(status='completed',results=results,wall_seconds=time.perf_counter()-started,**resources()))
