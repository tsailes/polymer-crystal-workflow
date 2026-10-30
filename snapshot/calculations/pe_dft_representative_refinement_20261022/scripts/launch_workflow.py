"""Detach the existing coordinator only after exclusive-lock/proc checks inside it."""
from pathlib import Path
import json, subprocess, time
root = Path(__file__).resolve().parents[1]
stamp = time.time_ns()
with (root / 'control' / f'coordinator_{stamp}.stdout.log').open('xb') as out, \
     (root / 'control' / f'coordinator_{stamp}.stderr.log').open('xb') as err:
    child = subprocess.Popen(['/opt/pe_dft/envs/pe-gpaw-20261020/bin/python',
        str(root / 'scripts/run_refinement.py')], cwd=root, stdin=subprocess.DEVNULL,
        stdout=out, stderr=err, start_new_session=True)
record = {'pid': child.pid, 'launched_unix': time.time(), 'log_stamp': stamp}
(root / 'control' / f'launch_{stamp}.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record))
