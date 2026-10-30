"""Run one authorized g4 cell relaxation, then an independent cold electronic check."""
from pathlib import Path
import fcntl, hashlib, json, os, subprocess, sys, time, traceback
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CONTROL = ROOT / 'control'
ENVROOT = Path('/opt/pe_dft/envs/pe-gpaw-20261020')
def save(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

lock = (CONTROL / 'workflow.lock').open('a+')
try:
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
except BlockingIOError:
    raise SystemExit('A g4 cell workflow already holds the exclusive lock')
run = ROOT / 'runs' / 'attempt_01'
if (run / 'status.json').exists() or (CONTROL / 'workflow_result.json').exists():
    raise SystemExit('Existing run evidence found; explicit continuation in a fresh attempt is required')
env = dict(os.environ)
env.update(PATH=str(ENVROOT / 'bin') + ':' + env.get('PATH', ''),
           GPAW_SETUP_PATH=str(ENVROOT / 'lib/python3.12/site-packages/gpaw_data/setups'),
           GPAW_MPI_BACKEND='cgpaw', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
           MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1',
           OMPI_ALLOW_RUN_AS_ROOT='1', OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1')
started = time.time()
record = {'status':'starting', 'coordinator_pid':os.getpid(), 'started_unix':started,
          'root':str(ROOT), 'source_sha256':sha(ROOT/'source/final_structure.extxyz')}
save(CONTROL/'workflow_status.json', record)

def stage(name, script, spec):
    work=spec.parent
    cmd=[str(ENVROOT/'bin/mpiexec'),'--allow-run-as-root','--bind-to','none','-n','4',
         str(ENVROOT/'bin/python'),str(script),str(spec)]
    with (work/'stdout.log').open('x') as out, (work/'stderr.log').open('x') as err:
        proc=subprocess.Popen(cmd, cwd=work, env=env, stdin=subprocess.DEVNULL,
                              stdout=out, stderr=err)
        receipt={'stage':name,'pid':proc.pid,'coordinator_pid':os.getpid(),
                 'command':cmd,'started_unix':time.time(),
                 'input_sha256':sha(spec), 'source_sha256':record['source_sha256']}
        save(work/'process.json',receipt)
        record.update(status=name, child_pid=proc.pid, current_work=str(work))
        save(CONTROL/'workflow_status.json',record)
        code=proc.wait()
        receipt.update(exit_code=code,finished_unix=time.time())
        save(work/'process.json',receipt)
        return code

try:
    code=stage('optimizing_cell',ROOT/'scripts/run_relax.py',run/'input.json')
    result=json.loads((run/'result.json').read_text()) if (run/'result.json').exists() else {}
    if code or result.get('status')!='passed_relaxation':
        raise RuntimeError(f'Relaxation not accepted: exit={code}, status={result.get("status")}')
    cold=ROOT/'cold_verification'
    cold.mkdir(exist_ok=False)
    spec=json.loads((run/'input.json').read_text())
    spec['structure']=str(run/'final_structure.extxyz')
    spec['checkpoint']=False
    save(cold/'input.json',spec)
    endpoint_hash=sha(run/'final_structure.extxyz')
    code=stage('independent_cold_verification',ROOT/'scripts/run_point.py',cold/'input.json')
    check=json.loads((cold/'result.json').read_text())
    maximum_stress=float(np.abs(check.get('total',{}).get('stress_kbar',[float('inf')])).max())
    maximum_force=float(check.get('total',{}).get('max_force_norm_eV_A',float('inf')))
    passed=(code==0 and check.get('status')=='passed_single_point'
            and check.get('SCF_converged') is True and maximum_force<=spec['fmax_eV_A']
            and maximum_stress<=spec['stress_max_kbar']
            and endpoint_hash==sha(run/'final_structure.extxyz'))
    verification={'passed':passed,'independent_density_and_wavefunctions':True,
                  'geometry_sha256':endpoint_hash,'raw_fmax_eV_A':maximum_force,
                  'raw_max_stress_component_kbar':maximum_stress,
                  'energy_delta_cold_minus_relax_eV':check['total']['energy_extrapolated_eV']-
                         result['final_evaluation']['total']['energy_extrapolated_eV'],
                  'source_relaxation_status':result['status']}
    save(ROOT/'final_cold_verification.json',verification)
    record.update(status='accepted_after_independent_cold_verification' if passed else 'cold_verification_failed',
                  verification=verification)
except Exception as exc:
    record.update(status='needs_review',error=repr(exc),traceback=traceback.format_exc())
    print(record['traceback'],file=sys.stderr,flush=True)
finally:
    record.update(finished_unix=time.time(),wall_seconds=time.time()-started)
    save(CONTROL/'workflow_result.json',record)
    save(CONTROL/'workflow_status.json',record)
