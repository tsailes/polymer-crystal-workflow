"""Collect diagnosis evidence without launching any electronic calculations."""
from pathlib import Path
import json
import numpy as np
root=Path(__file__).resolve().parent
def load(p): return json.loads((root/p).read_text(encoding='utf-8-sig'))
g=load('geometry_diagnosis.json')['frames']
t=load('tight_continuation/result.json')
gt=load('geometry_tight_diagnosis.json')['frames']
sy=load('symmetry_diagnosis.json')
scan=load('geometry_phase24_summary.json')
dense=load('tight_sampling_1x2x1/result.json')
traw=load('tight_continuation/final_raw_results.json')
tight_sampling=dict(base_fmax_eV_A=t['final']['raw_fmax_eV_A'],
    b_doubled_fmax_eV_A=dense['fmax_eV_A'],
    max_force_change_eV_A=float(np.linalg.norm(np.asarray(dense['folded_forces_eV_A'])-traw['forces_eV_A'],axis=1).max()),
    energy_change_meV_per_repeat=(dense['energy_per_base_cell_eV']-traw['energy_eV'])/8*1000)
base=root.parent
paths=lambda p:'E:/PA6_DFT/calculations/gfn2_single_pilot_20261030/diagnosis_20261030/'+p
text='''# PA6 单候选未形成氢键网络：原因诊断

2026-10-30。结论分为已观察到的几何事实与尚待因果验证的解释。最直接的证据指向当前链注册与局部优化范围；周期超胞敏感性明确存在；空间群本身不禁止目标网络；暂不能判断为GFN2方法失效。

## 1 先排除检测问题与粗停止门槛

原子角色、152条周期共价边、四条周期链均复核。全邻胞枚举、扩大图像窗口、独立ASE最小镜像距离与180°/90°/0°角定义控制一致，没有发现漏检或图像/角度错误。八个NH供体等价。

|状态|最近跨链N···O / Å|最近跨链H···O / Å|最近配对N—H···O / °|诊断接触数|
|---|---:|---:|---:|---:|
'''
for name,row in [('独立生成初态',g['initial']),('62步粗优化终态',g['final']),('收紧门槛续算终态',gt['final'])]:
    q=row['summary']
    text+=f"|{name}|{q['minimum_crosschain_NO_A']:.4f}|{q['minimum_crosschain_HO_A']:.4f}|{sum(q['nearest_NO_angle_range_deg'])/2:.2f}|{q['strict_bond_count']}|\n"
text+=f'''
诊断条件是N···O≤3.5 Å、H···O≤2.6 Å、角≥120°。原粗终态的近邻O来自跨层链，角度只有约97.7°；同层方向合适的O仍很远（H···O约5.73 Å）。单独放宽距离或角度并不能把原粗终态变成常规氢键片层，不能把“0”简单归为单一阈值边缘。

收紧试验保持同一固定胞、GFN2模型、基础胞Gamma与实际8操作，将accuracy设为0.01、原始Fmax门槛设为0.005 eV/Å。从原粗终态继续{t['steps']}步，耗时{t['wall_seconds']:.1f}秒，最终原始Fmax={t['final']['raw_fmax_eV_A']:.6f} eV/Å，状态为`{t['status']}`；进一步降能{-t['energy_change_eV']:.6f} eV/晶胞。以上终态几何为本次实际结果，是否达到门槛与是否形成网络分别判断。这一试验不解除对称约束，也不解决周期超胞敏感性。

## 2 初始注册有直接几何证据，但势垒尚未测得

固定同一独立模板、x=0.25、转角0°，对批准的24个相位做了完整几何对照。全部条目保留，不按氢键淘汰、排序或挑选量化代表；19项同时通过生成器及独立几何检查。phase=0.125至0.229167的6个相位各有8个诊断接触、8/8供受体覆盖，并形成两张沿a/c无限延伸的网络。例如phase=1/6为N···O≈2.876 Å、H···O≈1.876 Å、角≈167.79°。

这些结果证明同一化学模板及对称条件下，另一些轴向注册在几何上允许网络。原粗优化的局部酰胺中心沿c只改变约±0.022 Å，未完成该类大幅注册切换。但没有测定完整滑移能量面或路径，不能宣称已经证明能垒、能量排序、稳定局部极小或全局最低能。

## 3 实际约束更强，仍不是网络存在性的禁令

方案声明4个操作，自动FixSymmetry实际保留8个操作，多出初猜的纯(0,1/2,0)平移及其组合。两者都可标成No.14，所以只看空间群号不够。456个位移自由度中，仅声明4操作允许114维，实际8操作允许57维；扣除整体平移后的9维相对刚性链平移分别允许3维与1维。

24个相位全部具有相同8操作，包括有网络的6个相位。因此，“No.14或额外b/2平移直接禁止氢键网络”不成立。更强约束仍可能影响到达网络的路径、局部自由度或分岔，需显式4操作和无约束微扰对照。

原轨迹被投影去掉的最大原子力仅7.34e-13 eV/Å，说明该对称点的一阶力没有被明显削去，不能据此排除被禁方向的不稳定性。另对4个独立刚性链位移方向进行了±0.02/±0.04 Å单点：所试点均升能；部分有限差分曲率对步长敏感。该小范围结果不构成完整Hessian、能垒或稳定性验证。

## 4 周期超胞敏感性确实影响受力判断

以下为原62步粗终态的同几何加倍检查，非收紧续算终态。

|复制胞|原始Fmax / eV Å⁻¹|相对基础胞最大逐原子ΔF / eV Å⁻¹|能量变化 / meV每重复单元|
|---|---:|---:|---:|
|1×1×1|0.038286|—|—|
|2×1×1|0.039983|0.022862|−10.533|
|1×2×1|0.059770|0.029858|−22.206|
|1×1×2|0.039081|0.002101|−1.654|

b加倍已使最大力超过原0.05粗门槛；N/O/NH的受力也有明显变化，不能视为与氢键无关的常数能量平移。最近H/O配对接近坐标的微小残余总力投影还发生符号变化，但它不是单独氢键力，也不能据此预测优化路径。

对本次收紧续算终态又做了一次b加倍复核：基础胞Fmax={tight_sampling['base_fmax_eV_A']:.6f}→{tight_sampling['b_doubled_fmax_eV_A']:.6f} eV/Å，最大逐原子力变化{tight_sampling['max_force_change_eV_A']:.6f} eV/Å，能量变化{tight_sampling['energy_change_meV_per_repeat']:.3f} meV/重复单元。该敏感性没有随基础胞几何加严而消失，也仍超过本次0.005门槛。这不是已完成组合及更密网格收敛。

应称“周期超胞/采样敏感性”，不把全部变化都解释成纯k点误差。当前独立tblite接口的显式k点路径有局限，且v0.7周期二阶静电采用与晶胞有关的Wigner–Seitz处理；需要组合及更密超胞检验，而非把PE的采样参数直接移用。[tblite官方k点讨论](https://github.com/tblite/tblite/issues/4)、[周期静电实现说明](https://github.com/tblite/tblite/pull/350)。

## 5 GFN2目前不能被定为原因

电子accuracy加严后原粗终态最大力变化仅1.23e-5 eV/Å，已有N/O/H抽样导数检查也通过，因此电子迭代精度或这些已抽查的单位/符号不是主要证据指向的问题。不过，少量导数一致不证明物理模型正确。

GFN2本来就以包括氢键的非共价作用为设计对象，含各向异性静电；分子基准不能替代周期PA6验证。[GFN2原始论文](https://doi.org/10.1021/acs.jctc.8b01176)。要判断方法误差，应在同一几何、数值收敛与可比较晶胞条件下对照GFN2和独立DFT的N/O/NH力及相关集体位移能量变化。

## 6 下一步按单因素验证

1. **先统一周期处理。** 在冻结终态比较组合2×2×1、2×2×2与必要的更密横向复制胞；关键力变化应低于计划的优化停止尺度后，再统一几何优化设置。这些是待验证网格，不是已收敛生产参数。
2. **再分离约束。** 同一固定胞、同一已检查的采样与停止门槛，比较实际8操作、仅声明4操作、完全无空间群约束。用预先固定的小扰动和正负成对重启；4操作分支须显式打破额外b/2平移，自由分支还须包含破坏4操作的分量。直接删除FixSymmetry但保留精确对称坐标不足以验证。[ASE官方示例](https://docs.ase-lib.org/ase/constraints.html#the-fixsymmetry-class)。
3. **单独检验初始注册。** 按既有网格均匀预选少量相位，使用同一方法与约束协议；仍只按原有化学/防重叠规则放行，保留所有终态。网络仅作事后分析。不能只挑已有氢键的起点来证明流程成功。
4. **最后隔离电子模型。** 将相同冻结几何与小幅集体位移交给云端PBE-D3(BJ)+ATM，对关键分歧再用r2SCAN-D4；先比较力和能量变化，再决定是否精修。单一DFT方法也不是未经验证的真值。

固定实验晶胞亦可能影响结果，但应在以上比较后独立测试变胞，避免同时改变初猜、约束、周期处理和晶胞导致无法归因。

详细可执行条件保存在[下一步验证协议]('''+paths('next_validation_protocol.json')+''')；它是提案，本次未启动这些批量或DFT对照。所有几何对照都未使用论文原子坐标，未改变氢键不参与候选筛选的规则。
'''
(root/'README.md').write_text(text,encoding='utf-8')
summary=dict(status='diagnosis_completed',tight_continuation=t,
             tight_geometry=gt['final']['summary'],detector_audit=gt['final']['detector_audit'],
             tight_endpoint_b_doubling=tight_sampling,
             actual_operation_count=sy['operation_count'],declared_operation_count=sy['declared_operation_count'],
             extra_translation_does_not_geometrically_forbid_network=True,
             next_validation='next_validation_protocol.json',method_error_established=False,
             conclusions='Initial registry/local search is supported as explanation; periodic sensitivity is confirmed; constraint-path and model causality remain open')
(root/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
project=base.parents[1]
spec=json.loads((project/'project.json').read_text(encoding='utf-8-sig'))
spec['latest_diagnosis']={'date_local':'2026-10-30','report':'calculations/gfn2_single_pilot_20261030/diagnosis_20261030/README.md',
    'status':'completed','production_parameters_validated':False,'DFT_started':False,
    'scope':'Geometry/symmetry/force audit, small directional points and same-candidate tighter Gamma continuation; no quantum screening of additional phase candidates'}
(project/'project.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2))
