"""Assess deliberately paired refinements; never infer convergence from row order."""
import argparse
from pathlib import Path
import json
import numpy as np
from ase.io import read
from workflow import compare, dump, digest, THRESHOLDS, SETTINGS

parser = argparse.ArgumentParser()
parser.add_argument("--run", type=Path, required=True)
parser.add_argument("--phase", choices=["initial", "final"], required=True)
parser.add_argument("--relaxation", default="relaxation_2x3x8")
args = parser.parse_args()
run = args.run


def load(path):
    return json.loads((run / path).read_text(encoding="utf-8"))


def replicas_pass(data):
    return bool(data["replica_force_max_deviation_eV_per_A"] <= THRESHOLDS["replica_force_eV_per_A"])


base_path = "transverse_sampling/2x3x8/result.json"
base = load(base_path)
checks = []
for label, other in [
    ("a: 2x3x8 -> 3x3x8", "transverse_sampling/3x3x8/result.json"),
    ("b: 2x3x8 -> 2x4x8", "refined_sampling/2x4x8/result.json"),
    ("c: 2x3x8 -> 2x3x12", "refined_sampling/2x3x12/result.json"),
    ("accuracy: 0.01 -> 0.001", "initial_accuracy_check/2x3x8/result.json"),
]:
    other_data = load(other)
    if base["input_sha256"] != other_data["input_sha256"]:
        raise ValueError("Refinements must use identical primitive geometry")
    comparison = compare(base, other_data)
    comparison["replica_forces_passed"] = replicas_pass(base) and replicas_pass(other_data)
    comparison["passed"] = comparison["passed"] and comparison["replica_forces_passed"]
    checks.append(dict(label=label, first=base_path, second=other, **comparison))
fd = load("folded_force_check/result.json")
passed = all(r["passed"] for r in checks) and fd["passed"]
decision = dict(
    status="accepted_for_fixed_cell_local_relaxation" if passed else "not_accepted",
    selected_repeat=[2, 3, 8], expanded_atoms=576, primitive_atoms=12,
    settings=SETTINGS, thresholds=THRESHOLDS,
    checks=checks, folded_force_check=fd,
    rationale="a, b and c refinements and tighter electronic accuracy pass declared energy and atomic-force-vector criteria. c=8 retained for margin over c=6.",
    scope="Practical convergence at this geometry and fixed electronic smearing; repeat after relaxation.",
)
dump(run / "initial_decision.json", decision)
print(json.dumps(decision, ensure_ascii=False, indent=2), flush=True)
if not passed:
    raise RuntimeError("Initial parameter selection did not pass")
if args.phase == "initial":
    raise SystemExit(0)

relax = load(f"{args.relaxation}/final_result.json")
final_checks = []
fresh = load("final_sampling/2x3x8/result.json")
for name, path in [
    ("a", "final_sampling/3x3x8/result.json"),
    ("b", "final_sampling/2x4x8/result.json"),
    ("c", "final_sampling/2x3x12/result.json"),
    ("combined a,b,c", "final_combined/3x4x10/result.json"),
    ("accuracy", "final_accuracy_check/2x3x8/result.json"),
]:
    data = load(path)
    if fresh["input_sha256"] != data["input_sha256"]:
        raise ValueError("Final comparisons must use identical primitive geometry")
    final_checks.append(dict(label=name, path=path, **compare(fresh, data),
        replica_forces_passed=replicas_pass(data),
        refined_max_force_eV_A=data["max_force_norm_eV_per_A"],
        refined_force_passed=data["max_force_norm_eV_per_A"] < THRESHOLDS["relaxation_fmax_eV_per_A"]))
history = load(f"{args.relaxation}/history.json")
audit = load(f"{args.relaxation}/geometry_validation.json")
final_model = read(run / args.relaxation / "pe_relaxed.extxyz")
fresh_model = read(run / "final_sampling" / "2x3x8" / "input.extxyz")
initial_model = read(run.parent.parent.parent / "model" / "pe_initial.extxyz")
geometry_bound = bool(np.array_equal(final_model.numbers, fresh_model.numbers)
    and np.allclose(final_model.positions, fresh_model.positions, rtol=0, atol=1e-10)
    and np.array_equal(final_model.cell.array, fresh_model.cell.array)
    and np.array_equal(final_model.pbc, fresh_model.pbc))
cell_preserved = bool(np.array_equal(final_model.cell.array, initial_model.cell.array))
final_pass = all(r["passed"] and r["refined_force_passed"] and r["replica_forces_passed"] for r in final_checks)
final_pass = bool(final_pass and relax["converged"] and audit["passed"] and relax["cell_unchanged"]
                  and replicas_pass(fresh) and geometry_bound and cell_preserved
                  and fresh["max_force_norm_eV_per_A"] < THRESHOLDS["relaxation_fmax_eV_per_A"])
summary = dict(
    status="accepted" if final_pass else "needs_further_work", initial_decision=decision,
    relaxation_folder=args.relaxation, optimization=relax, final_checks=final_checks,
    fresh_final=fresh, geometry_validation=audit,
    fresh_geometry_matches_final_file=geometry_bound, final_cell_matches_initial_file=cell_preserved,
    energy_drop_eV_per_cell=fresh["energy_per_primitive_cell_eV"] - base["energy_per_primitive_cell_eV"],
    energy_drop_eV_per_CH2=fresh["energy_per_CH2_eV"] - base["energy_per_CH2_eV"],
    final_model_sha256=digest(run / args.relaxation / "pe_relaxed.extxyz"),
    limitations=["Local optimization from one constructed packing", "Experimental-prior cell held fixed",
                 "Primitive translation period retained", "No lattice equilibrium or global minimum established",
                 "No DFT error calibration, phonon stability or temperature extrapolation performed"],
)
dump(run / "summary.json", summary)
dump(run.parent.parent / "latest_completed_run.json", {"run_directory": str(run), "status": summary["status"]})
rows = ["| 检查 | ΔE（meV/CH₂） | 最大原子力变化（eV/Å） | 通过 |",
        "|---|---:|---:|---|"]
for r in checks:
    rows.append(f'| {r["label"]} | {r["energy_change_meV_per_CH2"]:.6f} | {r["max_force_vector_change_eV_per_A"]:.6f} | {"是" if r["passed"] else "否"} |')
final_rows = ["| 终态复核 | ΔE（meV/CH₂） | 最大原子力变化（eV/Å） | 复核最大受力（eV/Å） |",
              "|---|---:|---:|---:|"]
for r in final_checks:
    final_rows.append(f'| {r["label"]} | {r["energy_change_meV_per_CH2"]:.6f} | {r["max_force_vector_change_eV_per_A"]:.6f} | {r["refined_max_force_eV_A"]:.6f} |')
stress_gpa = np.asarray(fresh["stress_eV_per_A3_voigt"]) * 160.21766208
report = f"""# PE 周期采样与固定晶胞优化

本轮状态：**{summary["status"]}**。这是从化学几何构造初猜得到的、保持原胞周期和固定晶胞的 GFN2-xTB 局部优化候选；尚未验证为实验正交相，也未证明全局最低能量。

## 每一步做了什么

1. 保持初始原子坐标不动，改变沿 a、b、c 的重复数，比较同一结构的能量和力。这样采样误差不会与结构变化混在一起。
2. 选用 **2×3×8 的 Γ 点重复超胞**，计算时共576个原子，优化变量仍为12原子的原胞坐标。优化器使用超胞能量除48，原胞力采用48个等价副本力的平均。eV/CH₂仅作报告归一化。
3. 用 ±0.0001 Å 的原胞坐标扰动验证受力。每个等价副本一起移动；差分受力误差为 **{fd["absolute_error_eV_A"]:.3e} eV/Å**，通过0.001 eV/Å门槛。
4. 固定 a=7.432 Å、b=4.929 Å、c=2.545 Å 和90°夹角，仅移动原子。ASE BFGS 最大步长0.06 Å，受力目标0.01 eV/Å；没有施加实验原子坐标或额外空间群约束。
5. 优化后重新独立计算，并加密三个方向和组合采样，复核受力与能量；检查周期链连接、键长、原子接触以及晶胞是否保持。

## 初态采样选择依据

门槛预先写入 protocol.json：|ΔE| <1 meV/CH₂，maxᵢ‖ΔFᵢ‖ <0.005 eV/Å。它们是本阶段的实用数值门槛，不代表方法对实验的误差。

{chr(10).join(rows)}

较粗2×2×8到2×3×8的能量变化约1.113 meV/CH₂、力变化0.01731 eV/Å，未通过，因此没有用2×2×8优化。完整逐次数据见各sampling文件夹。表中每次只比较对应的细化方向，不能把所有行的相邻差都解释成单方向收敛。

## 优化结果

| 量 | 初态 | 终态（重新独立计算） |
|---|---:|---:|
| 能量（eV/C₄H₈原胞） | {base["energy_per_primitive_cell_eV"]:.9f} | {fresh["energy_per_primitive_cell_eV"]:.9f} |
| 能量（eV/CH₂） | {base["energy_per_CH2_eV"]:.9f} | {fresh["energy_per_CH2_eV"]:.9f} |
| 最大原子受力（eV/Å） | {base["max_force_norm_eV_per_A"]:.6f} | {fresh["max_force_norm_eV_per_A"]:.6f} |

优化共 **{relax["steps"]}步**。同一采样下能量降低 **{-summary["energy_drop_eV_per_cell"]:.6f} eV/原胞**（{-summary["energy_drop_eV_per_CH2"]:.6f} eV/CH₂）。

初末晶胞逐值相同。周期键连接检查：{"通过" if audit["passed"] else "未通过"}；C—C范围 {audit["C_C_range_A"]} Å，C—H范围 {audit["C_H_range_A"]} Å，最近链间原子距离 {audit["minimum_interchain_distance_A"]:.6f} Å。

## 终态采样复核

{chr(10).join(final_rows)}

每项还要求复核后的最大受力仍低于0.01 eV/Å；只看力变化小不够。所有检查均针对同一终态结构。更密采样仍受原胞平移周期限制，未探索长周期畸变。

终态应力（ASE拉伸为正，xx,yy,zz,yz,xz,xy，GPa）：{stress_gpa.tolist()}。本次固定晶胞，不要求应力为零；残余应力会在后续晶胞优化/方法校准时研究。

## 方法、软件与温度

- 计算方法：GFN2-xTB，属于半经验紧束缚量子方法。
- 能量/力引擎：tblite 0.7.0；Python/ASE 3.29.0 负责建模、调用、优化和保存数据。没有用Python重写xTB能量公式。
- 参数：中性闭壳层，accuracy=0.01、最大自洽迭代500次、电子占据展宽300 K；accuracy=0.001作精度复核。
- 这里的300 K只设置电子占据，不是分子动力学温度，也未包含晶格振动或构象热涨落。报告能量是当前tblite电子结果（包含该占据方案的熵处理），不是完整晶体吉布斯自由能。
- `xtb` 是另一个独立程序入口，也可在编译启用后调用 `tblite`。本机发现的Windows独立xtb未在本轮用于这些结果，不能将不同版本的输出混合标记。
- Γ点重复超胞与原胞k点采样存在能带折叠关系，但本轮没有调用显式k点接口。采用“Γ点重复超胞采样检查”名称，以保留实现和截断的实际边界。

## 下一步与结果边界

此结果完成了数值采样和一次固定胞局部优化。它尚未完成多初猜结构搜索、晶胞优化、稳定性验证或GFN2-xTB误差标定。后续应在多个链取向/相对位移初猜上重复搜索，再用已选PBE-D3(BJ)、optB88-vdW、r²SCAN-D4交叉比较。当前实验来源晶胞属于先验，不能据此声称已经完全从零预测晶格常数。

## 师兄路线中前三个级别的含义

1. **候选结构生成**：提出许多不同的周期排列初猜。空间群约束对称性，链构象描述链内扭转，堆积方式描述链之间的位置、取向和错位。三者并非可以任意相乘：生成后要检查无限链能否接续、是否原子碰撞、是否重复。10¹–10³是计划规模；本轮只有一个通过筛查的初猜进入局部优化，未实施该规模的候选生成。
2. **DFT收敛和参考方法**：先验证截断能、k点等数值参数，再比较选定三种泛函的结果。数值收敛解决参数太粗的问题，不能消除泛函本身的近似误差。“一次性”只适用于经过验证的研究对象和范围；更换聚合物、晶胞或结构范围要复核。本轮只做GFN2-xTB采样检查，没有完成DFT参数收敛。
3. **GFN2-xTB粗筛和3σ标定**：先用便宜方法评估候选，再选有代表性的候选做DFT，比较同一几何/相同归一化下的相对能量，估计并校正偏差。用未参与拟合的结构检验剩余误差；接近筛选边界的候选应保留。3σ意指某种误差标准差的三倍，但必须先明确误差定义、样本覆盖和统计假设；当前尚未确定或执行这个规则。

“±3σ覆盖约99.73%”只对应正态分布等相应假设，不能直接宣传为误删概率保证。对不同晶型、链构象的误差可能不同，不能用一两个算例估计的σ可靠淘汰整个候选池。若DFT泛函之间排序冲突，也应保留争议候选并用另一个经验证的参考层级研究，而非把任一结果自动当真值。

三步实际有回路：先产生初步候选和xTB结果，才能抽取有代表性的DFT标定集；标定后再回到整个候选池筛选。本轮的单算例验收是这套循环的准备工作。

10¹–10³个初始候选、10¹–10²个存活候选是计算规模的预估。实际保留数量应由去重、结构多样性、能量差和标定误差共同决定，并非必须凑到某个固定数量。

## 文件

- `{args.relaxation}/pe_relaxed.cif`：终态晶体模型（P1写出，未强行指定空间群）。
- `{args.relaxation}/pe_relaxed.extxyz`：带周期标记、链标签的终态坐标。
- `{args.relaxation}/atomic_forces.csv`：优化末步原子力；独立终态/细采样力见对应result.json。
- `{args.relaxation}/optimization.traj`、`history.json`、`optimizer.log`：完整优化轨迹和收敛记录。
- `summary.json`：最终数值与判定；`protocol.json`：预设阈值；`initial_decision.json`：选参依据。
- 各采样子目录的 `native_output.log`、`input.extxyz`、`result.json`：原始输出与可核查输入。

官方依据：[tblite API](https://tblite.readthedocs.io/en/latest/api/python.html)、[xtb的tblite集成说明](https://github.com/grimme-lab/xtb/issues/514)、[xtb v6.7.1发布页](https://github.com/grimme-lab/xtb/releases/tag/v6.7.1)、[VASP关于截断能收敛](https://vasp.at/wiki/ENCUT)、[NIST关于标准差区间](https://itl.nist.gov/div898/handbook/prc/section2/prc261.htm)。
"""
(run / "result_report_zh.md").write_text(report, encoding="utf-8")
print("FINAL_ACCEPTANCE", final_pass, flush=True)
if not final_pass:
    raise RuntimeError("Final acceptance failed; refine or continue optimization before declaring completion")
