"""Resumable coordinator with immutable inputs and complete per-attempt logs."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import argparse, csv, datetime, fcntl, hashlib, json, os, subprocess, sys, time

ROOT=Path(__file__).resolve().parents[1]
DFT_ENV=Path('/opt/pe_dft/envs/pe-gpaw-20261020')
GFN_ENV=Path('/opt/pe_dft/envs/gfn2-pe-20261020')
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p): return json.loads(Path(p).read_text())
def dump(p,obj):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    temp=p.with_name(p.name+f'.{os.getpid()}.{time.time_ns()}.tmp')
    temp.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    temp.replace(p)
def aggregate():
    rows=[]
    for backend in ['dft','gfn2']:
        for work in sorted((ROOT/'runs'/backend).glob('*')):
            if not (work/'input.json').exists(): continue
            spec=load(work/'input.json'); proc=load(work/'process.json') if (work/'process.json').exists() else {}
            result=load(work/'result.json') if (work/'result.json').exists() and proc.get('exit_code') is not None else {}
            total=result.get('total',{})
            row={'run_id':work.name,'backend':backend,'sample_id':spec.get('sample_id'),'method':spec.get('method','GFN2-xTB'),
                 'status':result.get('status','running_or_interrupted'),'exit_code':proc.get('exit_code'),
                 'precision':spec.get('precision_label'),'ecut_eV':spec.get('ecut_eV'),'kmesh':json.dumps(spec.get('kmesh')),
                 'gpts':json.dumps(spec.get('gpts')),'repeat':json.dumps(spec.get('repeat')),
                 'electronic_temperature_K':spec.get('electronic_temperature_K',0 if backend=='dft' else None),
                 'energy_eV_cell':total.get('energy_free_eV',total.get('energy_eV_per_cell')),
                 'energy_eV_CH2':total.get('energy_eV_CH2'),'energy_extrapolated_eV_cell':total.get('energy_extrapolated_eV'),
                 'DFT_part_eV_cell':result.get('dft',{}).get('energy_free_eV'),'dispersion_eV_cell':result.get('dispersion',{}).get('energy_eV'),
                 'max_force_eV_A':total.get('max_force_norm_eV_A'),'wall_seconds':result.get('wall_seconds'),
                 'input_sha256':sha(work/'input.json'),'result_sha256':sha(work/'result.json') if result else None,
                 'result_path':str(work/'result.json'),'error':result.get('error')}
            rows.append(row)
    dump(ROOT/'tables/energies_all_attempts.json',rows)
    if rows:
        with (ROOT/'tables/energies_all_attempts.csv').open('w',encoding='utf-8-sig',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    return rows
def run_one(item,queue_hash,protocol_hash,runner_bytes):
    backend=item['backend']; base=ROOT/'runs'/backend/item['id']; spec=dict(item['spec'],sample_id=item['sample_id'],protocol_sha256=protocol_hash,queue_sha256=queue_hash)
    spec['runner_sha256_at_submission']=hashlib.sha256(runner_bytes[backend]).hexdigest()
    spec['source_geometry_sha256_expected']=sha(spec['structure'])
    index=1; work=base
    while work.exists():
        if (work/'input.json').exists() and load(work/'input.json')==spec and (work/'result.json').exists():
            result=load(work/'result.json'); proc=load(work/'process.json') if (work/'process.json').exists() else {}
            if result.get('status')=='passed_single_point' and proc.get('exit_code')==0:
                print('REUSE '+work.name,flush=True); return {'id':item['id'],'work':str(work),'status':'reused'}
        index+=1;work=base.with_name(base.name+f'__attempt{index}')
    work.mkdir(parents=True)
    dump(work/'input.json',spec)
    dump(work/'protocol_snapshot.json',load(ROOT/'protocol.json'))
    runner=work/'runner_snapshot.py';runner.write_bytes(runner_bytes[backend])
    environment=os.environ.copy()
    environment.update({'OMP_NUM_THREADS':str(spec.get('threads',4) if backend=='gfn2' else 1),'OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1','OMP_MAX_ACTIVE_LEVELS':'1'})
    if backend=='dft':
        environment.update({'GPAW_MPI_BACKEND':'cgpaw','OMPI_ALLOW_RUN_AS_ROOT':'1','OMPI_ALLOW_RUN_AS_ROOT_CONFIRM':'1',
                            'GPAW_SETUP_PATH':str(DFT_ENV/'lib/python3.12/site-packages/gpaw_data/setups'),
                            'PATH':str(DFT_ENV/'bin')+':'+environment.get('PATH','')})
        command=[str(DFT_ENV/'bin/mpiexec'),'--allow-run-as-root','--bind-to','none','-n',str(spec.get('mpi_ranks',4)),str(DFT_ENV/'bin/python'),str(runner),str(work/'input.json')]
    else:
        command=[str(GFN_ENV/'bin/python'),str(runner),str(work/'input.json')]
        # A dense 1440-atom point must not compete with two large PW jobs
        # when the available WSL memory is low. This changes scheduling only.
        minimum_kib=(5 if spec.get('repeat')==[3,4,10] else 2)*1024*1024
        while True:
            memory={line.split(':')[0]:int(line.split()[1]) for line in Path('/proc/meminfo').read_text().splitlines() if len(line.split())>=2}
            if memory['MemAvailable']>=minimum_kib: break
            time.sleep(3)
    started=time.time()
    print('START '+work.name,flush=True)
    with (work/'stdout.log').open('wb') as out,(work/'stderr.log').open('wb') as err:
        child=subprocess.Popen(command,cwd=work,env=environment,stdout=out,stderr=err,stdin=subprocess.DEVNULL,start_new_session=True)
        process={'command':command,'pid':child.pid,'started_unix':started,'backend':backend,'sample_id':item['sample_id'],'environment_overrides':{k:environment[k] for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']}}
        dump(work/'process.json',process)
        code=child.wait()
    process.update(exit_code=code,finished_unix=time.time(),wall_seconds=time.time()-started)
    dump(work/'process.json',process)
    result=load(work/'result.json') if (work/'result.json').exists() else {}
    assert sha(spec['structure'])==spec['source_geometry_sha256_expected'],'Frozen geometry changed during calculation'
    accepted=code==0 and result.get('status')=='passed_single_point'
    print(('DONE ' if accepted else 'FAILED ')+work.name+f' {process["wall_seconds"]:.1f}s',flush=True)
    return {'id':item['id'],'work':str(work),'status':'passed' if accepted else 'failed','exit_code':code}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('queue',nargs='?');ap.add_argument('--aggregate-only',action='store_true');args=ap.parse_args()
    if args.aggregate_only: aggregate();return
    path=Path(args.queue).resolve(); queue=load(path);phase=queue['phase']; protocol=ROOT/'protocol.json'
    lock=(ROOT/'control'/f'{phase}.lock').open('a+')
    fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
    for requirement in queue.get('requires',[]):
        assert load(ROOT/requirement).get('passed') is True, 'Gate not passed: '+requirement
    runner_bytes={}
    for kind in {p['backend'] for p in queue['points']}:
        runner_bytes[kind]=(ROOT/'scripts'/('dft_point.py' if kind=='dft' else 'gfn2_point.py')).read_bytes()
    qhash,phash=sha(path),sha(protocol)
    state={'phase':phase,'status':'running','pid':os.getpid(),'started_unix':time.time(),'queue_sha256':qhash,'protocol_sha256':phash,'total':len(queue['points']),'completed':[]}
    state_path=ROOT/'control'/f'{phase}_state.json';dump(state_path,state)
    pending=iter(queue['points']);active={};failed=False
    with ThreadPoolExecutor(max_workers=queue['max_workers']) as pool:
        while True:
            while not failed and len(active)<queue['max_workers']:
                item=next(pending,None)
                if item is None:break
                future=pool.submit(run_one,item,qhash,phash,runner_bytes);active[future]=item
            if not active:break
            done,_=wait(active,return_when=FIRST_COMPLETED)
            for future in done:
                item=active.pop(future)
                try:record=future.result()
                except Exception as exc:record={'id':item['id'],'status':'failed','coordinator_error':repr(exc)}
                state['completed'].append(record)
                if record['status']=='failed' and queue.get('stop_on_failure',True):failed=True
            state['updated_unix']=time.time();dump(state_path,state);aggregate()
            assert sha(path)==qhash and sha(protocol)==phash,'Protocol/queue mutated while running'
    state.update(status='failed' if failed else 'completed',finished_unix=time.time())
    dump(state_path,state);aggregate()
    if failed:sys.exit(1)
if __name__=='__main__':main()
