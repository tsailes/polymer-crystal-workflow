"""Summarize recorded pilot observations, without launching calculations."""
from pathlib import Path
from statistics import median
import json, numpy as np
root=Path(__file__).resolve().parent
def read(p): return json.loads((root/p).read_text(encoding='utf-8-sig'))
bench=[]
for n in (1,2,4,8):
    rows=[json.loads(p.read_text()) for p in root.glob(f'benchmark_t{n}_*/result.json')]
    bench.append(dict(threads=n,runs=len(rows),wall_seconds_median=median(r['wall_seconds'] for r in rows),
        process_cpu_seconds_median=median(r['cpu_seconds'] for r in rows),peak_rss_MiB=max(r['peak_rss_MiB'] for r in rows)))
for row in bench:
    row['speedup_vs_1']=bench[0]['wall_seconds_median']/row['wall_seconds_median']
    row['parallel_efficiency']=row['speedup_vs_1']/row['threads']
relax=read('relaxation/result.json')
raw=read('relaxation/final_raw_results.json')
tight=read('final_cold_tight/result.json')
sampling=[]
for grid in ('2x1x1','1x2x1','1x1x2'):
    r=read(f'final_sampling_{grid}/result.json')
    sampling.append(dict(repeat=grid,energy_change_meV_per_repeat=(r['energy_per_repeat_eV']-tight['energy']/8)*1000,
        max_force_vector_change_eV_A=float(np.linalg.norm(np.array(r['folded_forces_eV_A'])-tight['forces'],axis=1).max()),
        replica_force_max_deviation_eV_A=r['replica_force_max_deviation_eV_A'],wall_seconds=r['wall_seconds'],peak_rss_MiB=r['peak_rss_MiB']))
geometry=read('geometry_comparison.json')
result=dict(benchmark=bench,relaxation=relax,
    cold_tight_check=dict(energy_change_eV=tight['energy']-raw['energy_eV'],fmax_eV_A=tight['fmax_eV_A'],
      max_force_vector_change_eV_A=float(np.linalg.norm(np.array(tight['forces'])-raw['forces_eV_A'],axis=1).max())),
    sampling_sensitivity=sampling,geometry_summary=dict(displacement_rms_A=geometry['displacement']['rms_A'],
      displacement_max_A=geometry['displacement']['max_A'],covalent_topology_preserved=geometry['covalent_topology_comparison']['preserved'],
      hydrogen_bond_counts=[geometry['hydrogen_bond_changes']['initial_bond_count'],geometry['hydrogen_bond_changes']['final_bond_count']]),
    batch_singlepoint=read('batch_4x4_summary.json'),
    sampling_convergence_status='axis-doubling sensitivity only; no combined or higher-grid convergence claim',
    production_parameters_validated=False,local_DFT_started=False)
top_results=[p for p in root.glob('*/result.json')]
result['resources']=dict(recorded_process_CPU_seconds=sum(json.loads(p.read_text())['cpu_seconds'] for p in top_results),
    peak_individual_process_RSS_MiB=max(json.loads(p.read_text())['peak_rss_MiB'] for p in top_results),
    run_files_bytes=sum(p.stat().st_size for p in root.rglob('*') if p.is_file()),
    project_files_bytes=sum(p.stat().st_size for p in root.parents[1].rglob('*') if p.is_file()),
    energy_kWh=None,electricity_cost=None,cloud_calculation_cost=0,
    notes='CPU seconds sum terminal process accounting only; excludes geometry-only setup. Disk measured before summary writes. No power measurement.')
(root/'summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,ensure_ascii=False,indent=2))
