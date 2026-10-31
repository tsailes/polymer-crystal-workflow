"""Publish completed single-pilot status, preserving previous project records."""
from pathlib import Path
import json, shutil
root=Path(__file__).resolve().parent
project=root.parents[1]
s=json.loads((root/'summary.json').read_text())
r=s['relaxation']; g=s['geometry_summary']; res=s['resources']
text=f'''# PA6 单候选 GFN2 粗优化试跑（2026-10-30）

已完成一个有序候选的固定胞、保对称、基础胞 Gamma 点粗优化。**流程跑通，但尚未获得目标氢键片层，也未通过周期采样收敛验收。** 本次不启动更多候选或本机 DFT。

## 输入与范围

候选 `O-v12-pilot-20261030-firstpass`：152 原子 C48H88N8O8、4 条周期链、8 个 PA6 重复单元。按 v1.2 已声明的24相位网格顺序取首个几何合格项：模板0、x=0.25、转角0°、相位0。没有按氢键或能量挑选，没有读取文献原子坐标。仅局部生成该候选，未生成整个新版库。

固定白名单晶胞 a=9.43、b=7.62、c=17.46 Å，γ=67°。使用 tblite 0.7.0 / ASE 3.29.0 的周期 GFN2-xTB、accuracy=0.1、闭壳层中性、电子占据温度300 K。该温度不表示核温度或分子动力学。采用 FixSymmetry，初始位置修正约3e-8 Å，胞修正为机器精度；因此保持的 No.14 是约束结果，不是自由稳定性证明。粗优化停止门槛为原始最大力≤0.05 eV/Å，与正式精修门槛不同。

## 实测结果

- 8线程粗优化：{r['steps']}步，{r['wall_seconds']:.2f}秒，原始最大力 {r['initial']['raw_fmax_eV_A']:.4f}→{r['final']['raw_fmax_eV_A']:.5f} eV/Å。
- 能量下降 {-r['energy_change_eV']:.5f} eV/晶胞，即 {-r['energy_change_eV']/8:.5f} eV/重复单元；仅表示同一起点的弛豫降能。
- 152条周期共价边保持，4条周期链完整。原子位移 RMS {g['displacement_rms_A']:.5f} Å、最大 {g['displacement_max_A']:.5f} Å。
- 氢键事后诊断：初、终态均0条满足 N···O≤3.5 Å、H···O≤2.6 Å、N—H···O≥120°的跨链接触；独立邻胞枚举复核一致。不能称已恢复目标氢键片层，也不据此删除候选。
- accuracy=0.01 独立冷启动最大力 {s['cold_tight_check']['fmax_eV_A']:.5f} eV/Å，最大力向量变化 {s['cold_tight_check']['max_force_vector_change_eV_A']:.2e} eV/Å。选定N/O/H分量有限差分及整体平移检查通过，非全部导数验证。

## 并行计时

同一初始结构、独立冷启动，各线程数做两次，顺序1/2/4/8/8/4/2/1，BLAS固定1线程。

|线程|单点计算中位时间/s|相对1线程加速|
|---:|---:|---:|
'''
for b in s['benchmark']:
    text+=f"|{b['threads']}|{b['wall_seconds_median']:.3f}|{b['speedup_vs_1']:.2f}×|\n"
text+=f'''
同一候选另做4个独立单点、每个4线程：合计墙钟 {s['batch_singlepoint']['total_wall_seconds']:.3f} 秒，全部成功。该结果证明并发单点可运行，不等于已测试4个完整优化。单候选可先用8线程；批量优先测试4任务×4线程，给32逻辑线程的桌面保留余量。优化步骤自身有前后依赖，不能把同一次优化的不同步独立并发。

## 周期采样限制

同一终态在accuracy=0.01下分别沿三个方向加倍，能量按8个重复单元的基础胞折算，力取副本平均。

|复制胞|能量变化/meV每重复单元|最大力向量变化/eV Å⁻¹|单点/s|
|---|---:|---:|---:|
'''
for a in s['sampling_sensitivity']:
    text+=f"|{a['repeat']}|{a['energy_change_meV_per_repeat']:.3f}|{a['max_force_vector_change_eV_A']:.5f}|{a['wall_seconds']:.2f}|\n"
text+=f'''
横向采样影响明显，基础胞粗收敛不能作为精细能量排序或生产参数锁。未做组合加密、更密网格、约束释放、DFT或实验符合性验证。后续若推进，应先用同一候选确定可接受周期采样，并保持不同起点的独立性；不把这一个结果推广到整个候选库。

## 消耗与文件

- 粗优化峰值内存 {r['peak_rss_MiB']:.1f} MiB，累计CPU时间 {r['cpu_seconds']/60:.2f} CPU分钟。
- 全部量化计时/复核/优化进程累计 {res['recorded_process_CPU_seconds']/60:.2f} CPU分钟；最大单进程内存 {res['peak_individual_process_RSS_MiB']:.1f} MiB。
- 本轮文件约 {res['run_files_bytes']/1024**2:.2f} MiB；PA6项目含既有依赖约 {res['project_files_bytes']/1024**2:.1f} MiB。统计时点在最终摘要写入前。
- 全部在本机CPU执行，未调用GPU或付费云计算。没有测量整机功率，不能提供实际kWh/电费；也未获取聊天额度/Token账单。

结构文件：[起点](input/candidate.extxyz)、[终态CIF](relaxation/final.cif)、[终态extxyz](relaxation/final.extxyz)。
主要证据：[汇总](summary.json)、[逐步记录](relaxation/history.jsonl)、[几何变化](geometry_comparison.json)、[选择与来源](input/selection_record.json)。
运行入口及代码快照均随结果保留。所有输出目录均不可覆盖复跑；重跑需使用新的输出目录。本次只完成一个候选的资源与粗优化试验。
'''
(root/'README.md').write_text(text,encoding='utf-8')
backup=root/'previous_project_state'
backup.mkdir(exist_ok=True)
for name in ('project.json','README.md','LOCAL_ORDERED_RUN.md','runtime/STATUS.md'):
    source=project/name
    target=backup/name
    target.parent.mkdir(parents=True,exist_ok=True)
    if not target.exists(): shutil.copy2(source,target)
spec=json.loads((project/'project.json').read_text(encoding='utf-8-sig'))
spec['date_local']='2026-10-30'
spec['status']='single_gfn2_pilot_completed_with_limitations'
spec['branches']['O']['status']='one_v1_2_candidate_coarse_optimized_full_new_library_not_generated'
policy=spec['calculation_policy']
policy.update(quantum_runs_started=1,gfn2_runs_started=1,gfn2_processes_completed=len(list(root.glob('*/result.json'))),
              run_count_definition='One pilot campaign; standalone process count includes repeated timing and validation points')
policy['first_pilot'].update(selected_candidates=['O-v12-pilot-20261030-firstpass'],allocation_status='single_candidate_pilot_completed')
spec['gates']['G4_numerics_resources'].update(status='inconclusive',
    evidence=['calculations/gfn2_single_pilot_20261030/summary.json'],
    scope='Single-candidate Gamma coarse optimization and resource pilot completed; sampling convergence and target H-bond sheet not established')
active=spec['active_execution']
active.update(authorized_date_local='2026-10-30',status='single_gfn2_pilot_completed_with_limitations',
    latest_user_instruction='先挑1个候选结构GFN粗优化，实测耗时、并行与消耗；本次仅一个候选',
    next_compute_launch_allowed=False,
    completed_scope='Single independently generated v1.2 candidate, 62-step Gamma fixed-cell coarse optimization, timing and sensitivity checks',
    latest_result='calculations/gfn2_single_pilot_20261030/README.md')
(project/'project.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
note='''\n> **2026-10-30 最新状态：** 用户已授权并完成一个候选的GFN2粗优化（62步、135.4秒）。周期链保持；尚未形成满足诊断阈值的目标氢键片层，周期采样未收敛。只完成这一项试跑，新版完整候选库及DFT未运行。详见 [单候选试跑结果](E:/PA6_DFT/calculations/gfn2_single_pilot_20261030/README.md)。下文2026-10-29的“暂停/尚未启动”表述为历史状态，由本条及project.json更新。\n'''
for name in ('README.md','LOCAL_ORDERED_RUN.md','runtime/STATUS.md'):
    path=project/name
    contents=path.read_text(encoding='utf-8-sig')
    first,rest=contents.split('\n',1)
    path.write_text(first+'\n'+note+'\n'+rest,encoding='utf-8')
print(root/'README.md')
