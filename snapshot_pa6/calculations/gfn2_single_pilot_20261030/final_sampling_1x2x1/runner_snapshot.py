"""One-candidate, fixed-cell PA6 GFN2 pilot; no batch selection or DFT."""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import resource
import shutil
import sys
import time
import traceback

import numpy as np
from ase.constraints import FixSymmetry
from ase.io import read, write
from ase.optimize import BFGS
from tblite.ase import TBLite

SETTINGS = dict(method='GFN2-xTB', charge=0, multiplicity=1, accuracy=0.1,
                max_iterations=500, electronic_temperature=300., verbosity=1)

def encode(x):
    if isinstance(x, np.ndarray): return x.tolist()
    if isinstance(x, np.generic): return x.item()
    if isinstance(x, Path): return str(x)
    raise TypeError(type(x).__name__)

def save(path, data):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(data, default=encode, indent=2, allow_nan=False)+'\n')
    tmp.replace(path)

def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def resources():
    r = resource.getrusage(resource.RUSAGE_SELF)
    return dict(peak_rss_MiB=r.ru_maxrss/1024, cpu_seconds=r.ru_utime+r.ru_stime)

def evaluate(atoms, accuracy=0.1):
    atoms = atoms.copy()
    atoms.set_constraint()
    calc = TBLite(**dict(SETTINGS, accuracy=accuracy))
    t = time.perf_counter()
    cpu_start = resources()['cpu_seconds']
    calc.calculate(atoms, properties=['energy', 'forces', 'stress', 'charges'])
    res = calc.results
    out = {k:res[k] for k in ['energy','forces','stress','charges']}
    out.update(wall_seconds=time.perf_counter()-t,
               fmax_eV_A=float(np.linalg.norm(res['forces'],axis=1).max()),
               evaluation_cpu_seconds=resources()['cpu_seconds']-cpu_start, **resources())
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('mode',choices=['point','preflight','relax','audit'])
    ap.add_argument('--input',required=True,type=Path)
    ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--repeat',nargs=3,type=int,default=[1,1,1])
    ap.add_argument('--accuracy',type=float,default=.1)
    ap.add_argument('--max-steps',type=int,default=200)
    ap.add_argument('--fmax',type=float,default=.05)
    args=ap.parse_args()
    if args.mode != 'point' and args.repeat != [1,1,1]:
        raise ValueError('Only point implements repeat; relaxation and audit use the base cell')
    out=args.output
    out.mkdir(parents=True,exist_ok=False)
    atoms=read(args.input)
    if len(atoms)!=152 or atoms.get_chemical_formula()!='C48H88N8O8' or not atoms.pbc.all():
        raise ValueError('Expected full-occupancy periodic PA6 C48H88N8O8')
    shutil.copy2(args.input,out/'input.extxyz')
    shutil.copy2(__file__,out/'runner_snapshot.py')
    save(out/'protocol.json',dict(mode=args.mode,input_sha256=digest(args.input),
         script_sha256=digest(__file__),settings=dict(SETTINGS,accuracy=args.accuracy),
         threads={k:os.getenv(k) for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']},
         versions={p:importlib.metadata.version(p) for p in ['ase','tblite','spglib','numpy']},
         cell=atoms.cell.array,atom_count=len(atoms),repeat=args.repeat,
         fmax_eV_A=args.fmax,max_steps=args.max_steps,
         sampling='Gamma only; numerical convergence unvalidated',
         scope='One ordered candidate, fixed cell, ordinary PA6; no DFT, no energy pruning',
         symmetry='FixSymmetry at 1e-5 A only in relax; raw and projected forces both recorded'))
    started=time.perf_counter()
    if args.mode=='point':
        result=evaluate(atoms.repeat(args.repeat),args.accuracy)
        result['energy_per_base_cell_eV']=result['energy']/np.prod(args.repeat)
        result['energy_per_repeat_eV']=result['energy_per_base_cell_eV']/8
        replica_forces=np.asarray(result['forces']).reshape(-1,len(atoms),3)
        result['folded_forces_eV_A']=replica_forces.mean(axis=0)
        result['replica_force_max_deviation_eV_A']=float(np.linalg.norm(replica_forces-replica_forces.mean(axis=0),axis=2).max())
        save(out/'result.json',result)
    elif args.mode=='preflight':
        base=evaluate(atoms,args.accuracy)
        checks=[]
        for element in ['N','O','H']:
            ids=np.flatnonzero(atoms.numbers=={'N':7,'O':8,'H':1}[element])
            ix,jx=np.unravel_index(np.abs(base['forces'][ids]).argmax(),(len(ids),3))
            atom=int(ids[ix]); step=1e-4; es=[]
            for direction in [-1,1]:
                shifted=atoms.copy(); shifted.positions[atom,jx]+=direction*step
                es.append(evaluate(shifted,args.accuracy)['energy'])
            fd=-(es[1]-es[0])/(2*step)
            checks.append(dict(element=element,atom=atom,axis=int(jx),finite_difference=fd,
                               analytic=base['forces'][atom,jx],error_eV_A=abs(fd-base['forces'][atom,jx])))
        shifted=atoms.copy(); shifted.positions += [.271,-.163,.329]; shifted.wrap()
        translated=evaluate(shifted,args.accuracy)
        de=abs(translated['energy']-base['energy'])
        df=float(np.linalg.norm(translated['forces']-base['forces'],axis=1).max())
        passed=all(c['error_eV_A']<.003 for c in checks) and de<1e-5 and df<.001
        save(out/'result.json',dict(passed=passed,base=base,force_checks=checks,
                translation_energy_error_eV=de,translation_force_error_eV_A=df,
                wall_seconds=time.perf_counter()-started,**resources()))
        if not passed: raise RuntimeError('PA6 selected force/translation preflight failed')
    elif args.mode=='relax':
        before_constraint=atoms.copy()
        atoms.set_constraint(FixSymmetry(atoms,symprec=1e-5,adjust_cell=False))
        constraint_change=dict(max_position_change_A=float(np.abs(atoms.positions-before_constraint.positions).max()),
                               max_cell_change_A=float(np.abs(atoms.cell.array-before_constraint.cell.array).max()))
        save(out/'symmetry_initialization_change.json',constraint_change)
        if max(constraint_change.values())>1e-6:
            raise ValueError('Symmetry refinement changed input beyond rounding precision')
        atoms.calc=TBLite(**dict(SETTINGS,accuracy=args.accuracy))
        opt=BFGS(atoms,trajectory=str(out/'trajectory.traj'),logfile=str(out/'optimizer.log'),
                 restart=str(out/'bfgs_restart.json'),maxstep=.04)
        history=[]
        def snapshot():
            raw=atoms.get_forces(apply_constraint=False)
            projected=atoms.get_forces()
            energy=atoms.get_potential_energy()
            row=dict(step=opt.nsteps,wall_seconds=time.perf_counter()-started,energy_eV=energy,
                     energy_eV_per_repeat=energy/8,raw_fmax_eV_A=float(np.linalg.norm(raw,axis=1).max()),
                     projected_fmax_eV_A=float(np.linalg.norm(projected,axis=1).max()),**resources())
            history.append(row)
            with (out/'history.jsonl').open('a') as f: f.write(json.dumps(row)+'\n')
            current=atoms.copy(); current.set_constraint()
            write(out/'latest.extxyz',current)
            save(out/'status.json',dict(status='running',latest=row))
            print('PROGRESS '+json.dumps(row),flush=True)
            if row['wall_seconds']>1800: raise TimeoutError('30 minute pilot limit reached')
            if row['peak_rss_MiB']>6000: raise MemoryError('Pilot 6 GiB process limit reached')
        opt.attach(snapshot,interval=1)
        converged=bool(opt.run(fmax=args.fmax,steps=args.max_steps))
        final=atoms.copy(); final.set_constraint()
        write(out/'final.extxyz',final); write(out/'final.cif',final)
        raw=atoms.get_forces(apply_constraint=False)
        save(out/'final_raw_results.json',dict(energy_eV=atoms.get_potential_energy(),forces_eV_A=raw,
                      stress_eV_A3=atoms.get_stress(apply_constraint=False),charges=atoms.get_charges()))
        save(out/'result.json',dict(status='coarse_converged' if converged and history[-1]['raw_fmax_eV_A']<=args.fmax else 'step_limit_or_raw_force_not_converged',
              projected_converged=converged,raw_converged=history[-1]['raw_fmax_eV_A']<=args.fmax,
              initial=history[0],final=history[-1],steps=opt.nsteps,
              energy_change_eV=history[-1]['energy_eV']-history[0]['energy_eV'],
              wall_seconds=time.perf_counter()-started,**resources()))
        save(out/'status.json',json.loads((out/'result.json').read_text()))
    elif args.mode=='audit':
        # Independent cold start at tighter electronic convergence; no geometry change.
        result=evaluate(atoms,args.accuracy)
        save(out/'result.json',result)
    print('FINISHED '+str(out),flush=True)

if __name__=='__main__':
    try: main()
    except Exception:
        failure=traceback.format_exc()
        if '--output' in sys.argv:
            failure_dir=Path(sys.argv[sys.argv.index('--output')+1])
            if failure_dir.is_dir():
                save(failure_dir/'failure.json',dict(status='failed',traceback=failure,**resources()))
                if (failure_dir/'status.json').is_file():
                    previous=json.loads((failure_dir/'status.json').read_text())
                    save(failure_dir/'status.json',dict(status='failed',previous=previous,traceback=failure))
        sys.stderr.write(failure)
        sys.exit(1)
