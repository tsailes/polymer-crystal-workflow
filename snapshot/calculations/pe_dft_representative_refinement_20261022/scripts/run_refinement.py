"""Run only the frozen four-point fixed-cell queue; preserve every attempt."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import datetime, fcntl, hashlib, json, os, shutil, subprocess, sys, time, traceback

ROOT = Path(__file__).resolve().parents[1]
ENV = Path('/opt/pe_dft/envs/pe-gpaw-20261020')
CONTROL = ROOT / 'control'

def load(p):
    return json.loads(Path(p).read_text())

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def dump(p, obj):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    temp = p.with_name(p.name + f'.{os.getpid()}.{time.time_ns()}.tmp')
    temp.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
    temp.replace(p)

def live_workers():
    running = []
    for p in Path('/proc').iterdir():
        if not p.name.isdigit() or int(p.name) == os.getpid():
            continue
        try:
            cmd = (p / 'cmdline').read_bytes().replace(b'\0', b' ').decode(errors='replace')
            state = (p / 'stat').read_text().rsplit(')', 1)[1].split()[0]
            if str(ROOT) in cmd and any(s in cmd for s in
                    ('runner_snapshot.py', 'run_fixedcell_relax.py')) and state != 'Z':
                running.append({'pid': int(p.name), 'command': cmd})
        except (FileNotFoundError, ProcessLookupError):
            continue
    return running

def environment():
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               NUMEXPR_NUM_THREADS='1', OMP_MAX_ACTIVE_LEVELS='1', GPAW_MPI_BACKEND='cgpaw',
               OMPI_ALLOW_RUN_AS_ROOT='1', OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1',
               GPAW_SETUP_PATH=str(ENV / 'lib/python3.12/site-packages/gpaw_data/setups'),
               PATH=str(ENV / 'bin') + ':' + env.get('PATH', ''))
    return env

def complete(work):
    if not (work / 'result.json').exists() or not (work / 'process.json').exists():
        return False
    r, p = load(work / 'result.json'), load(work / 'process.json')
    return (p.get('exit_code') == 0 and r.get('status') == 'passed_relaxation'
            and r.get('final_cold_verification', {}).get('passed') is True
            and all((work / x).exists() for x in ('final_POSCAR', 'final.cif', 'final.gpw',
                'evaluations.jsonl', 'optimizer.log', 'gpaw_final_cold.txt')))

def run_one(item, frozen, runner_bytes):
    base = ROOT / 'runs' / item['id']
    old = sorted(base.glob('attempt_*')) if base.exists() else []
    for work in old:
        if not (work / 'input.json').exists() or sha(work / 'input.json') != item['input_sha256']:
            raise RuntimeError(f'Input mismatch in previous attempt: {work}')
        if complete(work):
            print('REUSE ' + item['id'], flush=True)
            return {'id': item['id'], 'status': 'passed', 'reused': True, 'work': str(work)}
        if (work / 'result.json').exists() or (work / 'process.json').exists() and load(work / 'process.json').get('exit_code') is not None:
            raise RuntimeError(f'Previous failed calculation preserved; no automatic scientific retry: {work}')
    work = base / f'attempt_{len(old) + 1:02d}'
    work.mkdir(parents=True, exist_ok=False)
    source = Path(item['input'])
    if sha(source) != item['input_sha256']:
        raise RuntimeError('Frozen input changed: ' + str(source))
    shutil.copyfile(source, work / 'input.json')
    shutil.copyfile(ROOT / 'protocol.json', work / 'protocol_snapshot.json')
    runner = work / 'runner_snapshot.py'
    runner.write_bytes(runner_bytes)
    spec = load(source)
    cmd = [str(ENV / 'bin/mpiexec'), '--allow-run-as-root', '--bind-to', 'none', '-n',
           str(spec['mpi_ranks']), str(ENV / 'bin/python'), str(runner), str(work / 'input.json')]
    process = {'sample_id': item['id'], 'command': cmd, 'started_unix': time.time(),
               'input_sha256': sha(source), 'runner_sha256': sha(runner), 'frozen_files': frozen,
               'scheduler_recovery_of_interrupted_attempt': bool(old),
               'environment': {k: v for k, v in environment().items() if k in
                   ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS',
                    'GPAW_SETUP_PATH', 'GPAW_MPI_BACKEND']}}
    dump(work / 'process.json', process)
    print('START ' + item['id'] + ' ' + str(work), flush=True)
    with (work / 'stdout.log').open('xb') as out, (work / 'stderr.log').open('xb') as err:
        child = subprocess.Popen(cmd, cwd=work, env=environment(), stdout=out, stderr=err,
                                 stdin=subprocess.DEVNULL, start_new_session=True)
        process['pid'] = child.pid
        dump(work / 'process.json', process)
        code = child.wait()
    process.update(exit_code=code, finished_unix=time.time(), wall_seconds=time.time() - process['started_unix'])
    dump(work / 'process.json', process)
    ok = complete(work)
    if sha(spec['structure']) != spec['source_geometry_sha256_expected']:
        ok = False
        process['frozen_geometry_changed'] = True
        dump(work / 'process.json', process)
    print(('DONE ' if ok else 'FAILED ') + item['id'] + f' {process["wall_seconds"]:.1f}s', flush=True)
    return {'id': item['id'], 'status': 'passed' if ok else 'failed', 'work': str(work),
            'exit_code': code, 'wall_seconds': process['wall_seconds']}

def progress(item):
    attempts = sorted((ROOT / 'runs' / item['id']).glob('attempt_*'))
    if not attempts:
        return {'id': item['id'], 'status': 'queued'}
    work = attempts[-1]
    record = {'id': item['id'], 'work': str(work)}
    for filename, key in [('status.json', 'calculation'), ('process.json', 'process')]:
        try:
            obj = load(work / filename)
            record[key] = ({k: obj.get(k) for k in ('pid', 'started_unix', 'exit_code', 'wall_seconds')}
                           if key == 'process' else obj)
        except (FileNotFoundError, json.JSONDecodeError):
            pass
    return record

def main():
    CONTROL.mkdir(exist_ok=True)
    lock = (CONTROL / 'coordinator.lock').open('a+')
    fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    orphans = live_workers()
    if orphans:
        raise RuntimeError('Existing quantum workers found; no duplicate start: ' + repr(orphans))
    queue = load(ROOT / 'queue.json')
    if sha(ROOT / 'protocol.json') != queue['protocol_sha256']:
        raise RuntimeError('Protocol hash differs from frozen queue')
    paths = [ROOT / 'queue.json', ROOT / 'protocol.json', ROOT / 'scripts/run_point.py',
             ROOT / 'scripts/run_fixedcell_relax.py'] + [Path(p['input']) for p in queue['points']]
    frozen = {str(p): sha(p) for p in paths}
    runner_bytes = (ROOT / 'scripts/run_fixedcell_relax.py').read_bytes()
    state = {'status': 'running', 'stage': 'fixed_cell_PBE_D3BJ_relaxation', 'pid': os.getpid(),
             'started_unix': time.time(), 'total': len(queue['points']), 'completed': [],
             'frozen_files': frozen, 'max_workers': queue['max_workers']}
    dump(CONTROL / f'coordinator_start_{time.time_ns()}.json', state)
    def save():
        state['updated_unix'] = time.time()
        state['points'] = [progress(p) for p in queue['points']]
        dump(CONTROL / 'workflow_state.json', state)
    save()
    try:
        pending, active, failed = iter(queue['points']), {}, False
        with ThreadPoolExecutor(max_workers=queue['max_workers']) as pool:
            while True:
                if any(sha(p) != h for p, h in frozen.items()):
                    failed = True
                    state['integrity_error'] = 'Frozen input or runner changed; no new submissions'
                while not failed and len(active) < queue['max_workers']:
                    item = next(pending, None)
                    if item is None:
                        break
                    active[pool.submit(run_one, item, frozen, runner_bytes)] = item
                if not active:
                    break
                done, _ = wait(active, timeout=15, return_when=FIRST_COMPLETED)
                for future in done:
                    item = active.pop(future)
                    try:
                        record = future.result()
                    except Exception:
                        record = {'id': item['id'], 'status': 'failed', 'traceback': traceback.format_exc()}
                    state['completed'].append(record)
                    failed |= record['status'] != 'passed'
                state['stop_new_submissions'] = failed
                save()
        state['stage'] = 'final_analysis'
        save()
        with (CONTROL / f'analysis_{time.time_ns()}.log').open('w') as out:
            analysis = subprocess.run([str(ENV / 'bin/python'), str(ROOT / 'scripts/analyze_refinement.py')],
                                      cwd=ROOT, env=environment(), stdout=out, stderr=subprocess.STDOUT)
        state['analysis_exit_code'] = analysis.returncode
        summary = load(ROOT / 'tables/refinement_summary.json') if (ROOT / 'tables/refinement_summary.json').exists() else {}
        state['analysis_summary'] = summary
        all_passed = len(state['completed']) == state['total'] and all(x['status'] == 'passed' for x in state['completed'])
        state['status'] = 'completed' if all_passed and not failed and analysis.returncode == 0 and summary.get('complete') is True else 'needs_attention'
        state['stage'] = 'complete' if state['status'] == 'completed' else 'preserved_for_review'
    except Exception:
        state['status'] = 'needs_attention'
        state['error'] = traceback.format_exc()
        print(state['error'], flush=True)
    state['finished_unix'] = time.time()
    save()
    print(json.dumps({'status': state['status'], 'completed': state['completed']}), flush=True)
    if state['status'] != 'completed':
        sys.exit(1)

if __name__ == '__main__':
    main()
