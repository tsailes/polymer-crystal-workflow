"""Advance this authorized batch through explicit numerical gates, resumably."""
from pathlib import Path
import fcntl,json,os,subprocess,sys,time,traceback
ROOT=Path(__file__).resolve().parents[1]
PY='/opt/pe_dft/envs/pe-gpaw-20261020/bin/python'
def load(p):return json.loads(p.read_text())
def save(obj):
    p=ROOT/'control/workflow_state.json';tmp=p.with_suffix('.tmp')
    tmp.write_text(json.dumps(obj,indent=2)+'\n');tmp.replace(p)
def start(queue):
    q=ROOT/'queues'/queue;phase=load(q)['phase']; state=ROOT/'control'/f'{phase}_state.json'
    if state.exists():
        previous=load(state)
        if previous.get('status')=='completed':return None,state
        pid=previous.get('pid')
        if pid and Path(f'/proc/{pid}/cmdline').exists() and str(ROOT).encode() in Path(f'/proc/{pid}/cmdline').read_bytes():
            return None,state
    out=(ROOT/'control'/f'{phase}.stdout.log').open('ab');err=(ROOT/'control'/f'{phase}.stderr.log').open('ab')
    child=subprocess.Popen([PY,str(ROOT/'scripts/run_batch.py'),str(q)],stdout=out,stderr=err,stdin=subprocess.DEVNULL,start_new_session=True)
    out.close();err.close();return child,state
def wait_phase(item):
    child,path=item
    while True:
        if path.exists():
            state=load(path)
            if state.get('status')=='completed':return
            if state.get('status')=='failed':raise RuntimeError('Calculation phase failed: '+str(path))
        if child is not None and child.poll() is not None:
            if child.returncode:raise RuntimeError(f'Coordinator exit {child.returncode}: {path}')
        time.sleep(5)
def command(script):
    with (ROOT/'control/workflow_analysis.stdout.log').open('ab') as out,(ROOT/'control/workflow_analysis.stderr.log').open('ab') as err:
        subprocess.run([PY,str(ROOT/'scripts'/script)],stdout=out,stderr=err,check=True)
def main():
    lock=(ROOT/'control/supervisor.lock').open('a+')
    fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
    status={'status':'running','pid':os.getpid(),'started_unix':time.time(),'stage':'temperature_sampling_and_DFT_precision_pilots','hard_rejection_enabled':False}
    save(status)
    try:
        dft=start('01_dft_precision_pilot.json')
        gfn=start('02_gfn2_validation.json')
        wait_phase(dft)
        status.update(stage='matched_reference_upper_precision',updated_unix=time.time());save(status)
        wait_phase(start('01b_dft_reference_upper.json'))
        command('evaluate_dft_gate.py')
        wait_phase(gfn)
        command('evaluate_gfn2_gates.py')
        status.update(stage='same_geometry_paired_production',updated_unix=time.time());save(status)
        dft=start('03_dft_paired_production_reuse_v2.json');gfn=start('04_gfn2_paired_production.json')
        wait_phase(dft);wait_phase(gfn)
        status.update(stage='calibration_analysis',updated_unix=time.time());save(status)
        command('analyze_calibration.py')
        status.update(status='completed_pilot_calibration',stage='finished',finished_unix=time.time(),
           summary_path=str(ROOT/'tables/calibration_summary.json'),hard_rejection_enabled=False,
           note='Pilot errors and blocked validation reported. No unsupported final-basin 3sigma coverage claim or forced rejection.')
        save(status)
    except Exception as exc:
        status.update(status='needs_attention',error=repr(exc),traceback=traceback.format_exc(),updated_unix=time.time())
        save(status);raise
if __name__=='__main__':main()
