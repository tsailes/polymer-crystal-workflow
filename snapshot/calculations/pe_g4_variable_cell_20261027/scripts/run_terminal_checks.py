"""Independent, one-factor terminal numerical checks of the accepted g4 geometry."""
from pathlib import Path
import fcntl,hashlib,json,os,subprocess,time,traceback
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'convergence_checks'
ENVROOT=Path('/opt/pe_dft/envs/pe-gpaw-20261020')
GEOM=ROOT/'runs/attempt_01/final_structure.extxyz'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):p.write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
OUT.mkdir(exist_ok=True)
lock=(OUT/'workflow.lock').open('a+')
try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
except BlockingIOError:raise SystemExit('Terminal checks already running')
if (OUT/'status.json').exists() or (OUT/'result.json').exists():
    raise SystemExit('Existing terminal-check evidence found; refusing duplicate launch')
cold_ver=json.loads((ROOT/'final_cold_verification.json').read_text())
assert cold_ver['passed'] and cold_ver['geometry_sha256']==sha(GEOM)
base_spec=json.loads((ROOT/'cold_verification/input.json').read_text())
base=json.loads((ROOT/'cold_verification/result.json').read_text())
assert base['SCF_converged'] and base['status']=='passed_single_point'
env=dict(os.environ)
env.update(PATH=str(ENVROOT/'bin')+':'+env.get('PATH',''),GPAW_SETUP_PATH=str(ENVROOT/'lib/python3.12/site-packages/gpaw_data/setups'),GPAW_MPI_BACKEND='cgpaw',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',OMPI_ALLOW_RUN_AS_ROOT='1',OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1')
plan=[('k_only',1800,[5,10,15],[112,72,44],'baseline'),('grid_only',1800,[4,8,12],[128,80,48],'baseline'),('ecut_only',2000,[4,8,12],[128,80,48],'grid_only')]
protocol={'purpose':'Separate terminal k-point, grid, and Ecut checks; no cross-functional calculations or geometry changes','baseline':'../cold_verification/result.json','geometry_sha256':sha(GEOM),'point_runner_sha256':sha(ROOT/'scripts/run_point.py'),'plan':[{'name':n,'ecut_eV':e,'kmesh':k,'gpts':g,'comparison':r}for n,e,k,g,r in plan],'strict_absolute_endpoint_thresholds':{'fmax_eV_A':.002,'stress_max_kbar':.05},'historical_difference_diagnostics_only':{'energy_meV_atom':.5,'force_vector_eV_A':.005,'stress_component_kbar':.1},'note':'Historical difference criteria do not relax the strict absolute terminal force/stress thresholds. Ecut comparison uses the same denser grid at both cutoffs.'}
save(OUT/'protocol.json',protocol)
state={'status':'starting','coordinator_pid':os.getpid(),'started_unix':time.time(),'geometry_sha256':sha(GEOM),'stages':[]}
save(OUT/'status.json',state)
results={'baseline':base}
try:
    for name,ecut,kmesh,gpts,refname in plan:
        work=OUT/name;work.mkdir(exist_ok=False)
        spec=dict(base_spec,structure=str(GEOM),ecut_eV=ecut,kmesh=kmesh,gpts=gpts,checkpoint=False)
        save(work/'input.json',spec)
        cmd=[str(ENVROOT/'bin/mpiexec'),'--allow-run-as-root','--bind-to','none','-n','4',str(ENVROOT/'bin/python'),str(ROOT/'scripts/run_point.py'),str(work/'input.json')]
        with (work/'stdout.log').open('x') as stdout,(work/'stderr.log').open('x') as stderr:
            proc=subprocess.Popen(cmd,cwd=work,env=env,stdin=subprocess.DEVNULL,stdout=stdout,stderr=stderr)
            receipt={'pid':proc.pid,'command':cmd,'started_unix':time.time(),'input_sha256':sha(work/'input.json'),'geometry_sha256':sha(GEOM)}
            save(work/'process.json',receipt)
            state.update(status='running',current_stage=name,child_pid=proc.pid);save(OUT/'status.json',state)
            code=proc.wait()
            receipt.update(exit_code=code,finished_unix=time.time());save(work/'process.json',receipt)
        if code:raise RuntimeError(f'{name}: exit {code}')
        result=json.loads((work/'result.json').read_text())
        if result.get('status')!='passed_single_point' or not result.get('SCF_converged'):raise RuntimeError(f'{name}: failed SCF or single point')
        assert sha(GEOM)==protocol['geometry_sha256'],'Geometry changed during checks'
        ref=results[refname];total=result['total'];rt=ref['total']
        f=np.asarray(total['forces_eV_A']);rf=np.asarray(rt['forces_eV_A'])
        s=np.asarray(total['stress_kbar']);rs=np.asarray(rt['stress_kbar'])
        n=len(f)
        de=(total['energy_extrapolated_eV']-rt['energy_extrapolated_eV'])*1000/n
        df=float(np.linalg.norm(f-rf,axis=1).max());ds=float(np.abs(s-rs).max())
        fmax=float(np.linalg.norm(f,axis=1).max());smax=float(np.abs(s).max())
        row={'name':name,'reference':refname,'delta_energy_meV_atom':de,'max_force_vector_difference_eV_A':df,'max_stress_component_difference_kbar':ds,'raw_fmax_eV_A':fmax,'raw_stress_max_kbar':smax,'all_six_stress_components_kbar':s.tolist(),'strict_absolute_endpoint_pass':fmax<=.002 and smax<=.05,'historical_difference_diagnostics_pass':abs(de)<=.5 and df<=.005 and ds<=.1,'wall_seconds':result.get('wall_seconds'),'input':spec}
        save(work/'comparison.json',row);state['stages'].append(row);results[name]=result
        save(OUT/'status.json',state)
    state.update(status='completed',all_strict_absolute_endpoint_checks_pass=all(x['strict_absolute_endpoint_pass']for x in state['stages']),all_historical_difference_diagnostics_pass=all(x['historical_difference_diagnostics_pass']for x in state['stages']))
except Exception as e:
    state.update(status='needs_review',error=repr(e),traceback=traceback.format_exc())
finally:
    state.update(finished_unix=time.time(),wall_seconds=time.time()-state['started_unix'])
    save(OUT/'result.json',state);save(OUT/'status.json',state)
