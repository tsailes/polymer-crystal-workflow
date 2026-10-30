"""Write a user-facing report using recorded evidence, including pending work."""
from __future__ import annotations
import json
from pathlib import Path
import importlib.metadata
import sys

ROOT = Path(__file__).resolve().parent


def load(name):
    return json.loads((ROOT/name).read_text(encoding="utf-8"))


def main():
    generation = load("generation_manifest.json")
    summary = load("optimization_summary.json") if (ROOT/"optimization_summary.json").exists() else {}
    xrd = load("xrd_benchmark_initial/summary.json")
    dense = load("representative_convergence/summary.json") if (ROOT/"representative_convergence/summary.json").exists() else None
    final_xrd = load("xrd_benchmark_final/summary.json") if (ROOT/"xrd_benchmark_final/summary.json").exists() else None
    audit = load("candidates/cand_0001/geometry_validation.json")
    paths = []
    for entry in generation["unique_survivors"]:
        p = ROOT/"candidates"/entry["candidate_id"]/"POSCAR"
        paths.append("E:"+str(p).removeprefix("/mnt/e"))
    (ROOT/"POSCAR_list_initial.txt").write_text("\n".join(paths)+"\n", encoding="utf-8")
    if summary.get("all_preoptimization_attempts_completed"):
        final_paths = ["E:"+str(ROOT/"representatives"/g["cluster_id"]/"POSCAR").removeprefix("/mnt/e")
                       for g in summary["groups"]]
        (ROOT/"POSCAR_list_final.txt").write_text("\n".join(final_paths)+"\n", encoding="utf-8")
    (ROOT/"runtime_manifest.json").write_text(json.dumps(
        dict(python=sys.executable,
             packages={p:importlib.metadata.version(p)
                       for p in ("ase","numpy","scipy","gemmi","spglib","tblite","matplotlib")}),
        indent=2)+"\n")
    lines = [
        "# Tashiro PE：Pnam 约束候选生成与 GFN2-xTB 预优化",
        "",
        f"候选预优化状态：**{summary.get('status', '几何生成完成，预优化进行中')}**；"
        f"代表结构密采样复核：**{'通过' if dense and dense['passed'] else '尚未通过或待完成'}**。",
        "本批次目录：`E:/PE_DFT/calculations/pe_pnam_candidates_20261020`。",
        "",
        "## 先回答：以前做过了吗",
        "",
        "此前只有单个独立构造的 P1 候选完成了 GFN2-xTB 数值验证与固定胞优化。历史216候选是30°粗网格，"
        "未严格施加本次Pnam、0.85倍范德华接触过滤和批量GFN2流程；不能计作本次四步已经完成。"
        "本目录是新的独立批次，所有几何和能量均有各自来源，没有套用旧P1结果。",
        "",
        "## 输入约束与用途",
        "",
        "- 晶胞：a=7.432、b=4.929、c=2.545 Å，α=β=γ=90°；固定Tashiro实验晶胞。",
        "- 空间群：Pnam，国际编号62，Hall编号297，Hall符号 −P 2c 2n；链轴沿c。",
        "- 化学模型：两条无限PE链，共C₄H₈；每条链在c方向含两个碳。",
        "- 只用化学键长/键角与对称操作生成坐标；不读取或复制文献精修原子坐标。",
        "- setting角定义：全反式锯齿链平面在ab平面的投影相对于+a轴的角度；它不是链内扭转角。",
        "- 这是给定轴向周期、两链、指定空间群下的有限网格覆盖，不宣称穷尽连续空间、其他晶胞、缺陷或全部PE多晶型。",
        "",
        "## 第1步：链构象枚举",
        "",
        "按t=180°、g⁺=+60°、g⁻=−60°列出两扭转的9种状态组合，并先施加周期闭合约束。"
        "在两碳纯平移周期中，四个相邻碳必为A、B、A+c、B+c，四点共面；"
        "合理、不断键的链对应全反式。含gauche的8种组合因不满足这个指定周期而解析排除，"
        "不是先算能量再淘汰。更长周期的gauche缺陷需要另一个搜索模型。",
        "",
        "初始C—C=1.54 Å，C—H=1.09 Å。c与C—C共同确定C—C—C角，氢的位置按相容四面体几何建立。"
        "详见`audits/conformation_enumeration.json`及独立几何审计。",
        "",
        "## 第2步：空间群约束放置",
        "",
        "取向角取0°、5°……175°，180°与0°等价；36个角。"
        "枚举4个允许的半格横向链轴位置和2个碳镜面相位(¼、¾)，共36×4×2=288个参数组合。"
        "链轴的另外一些位置只是整体平移等价，仍保留枚举日志再去重。",
        "",
        "从一个C和两个H的独立位点应用Pnam的8个操作，归入晶胞并去除特殊位置上的重合像，"
        "得到4C+8H。第二条链的角度与第一条关联，不能独立转动。"
        "每个候选均逐操作核对，另用spglib识别空间群；0°候选具有更高的#63对称性，仍包含Pnam。",
        "",
        "## 第3步：物理预过滤",
        "",
        "使用Bondi范德华半径C=1.70 Å、H=1.20 Å，阈值为0.85×半径和："
        "C···C=2.89 Å、C···H=2.465 Å、H···H=2.04 Å。"
        "仅排除周期化学图中相隔1或2条共价键的原子对；所有其他近邻包括跨晶胞像均检查。"
        "另外保留你此前提出的更严格链间C···C≥3.4 Å条件。"
        "原子包裹到晶胞时保留其晶格位移，避免误把另一条周期像链当作同一链。",
        "",
        f"- 枚举：{generation['enumerated_count']}个；物理通过：{generation['physical_survivor_count']}个；"
        f"剔除：{generation['rejected_count']}个。",
        "- 被剔除的是80°–100°区间的取向组合，原因是H···H过近；所有这些组合的链间C···C仍≥3.4 Å。",
        f"- 等价去重：{generation['exact_duplicates_count']}个重复，保留"
        f"**{generation['unique_initial_count']}个不同起点：0°、5°……75°**。",
        "- 初始去重同时允许晶格点群操作、整体原点平移和同元素原子置换，要求逐原子坐标匹配，"
        "不是只看空间群名或成对距离直方图。",
        "",
        f"全部初态密度为{audit['density_g_cm3']:.9f} g/cm³。与外部文献给出的结晶相晶体学参考密度"
        f"0.997 g/cm³相比，比值{audit['density_ratio']:.6f}，落在0.95–1.05。"
        "这个参考不是Tashiro同一样品的独立实测体密度；后者目前未取得。固定晶胞和组成下密度不随堆积改变，"
        "因此它是输入合理性检查，不是堆积优劣评分。",
        "",
        "密度来源：Martin & Passaglia (1966), *Density of Polyethylene Crystals Grown From Solution*, "
        "[J. Res. NBS 70A, 221–224](https://pmc.ncbi.nlm.nih.gov/articles/PMC6687595/)。"
        "该文明确区分样品实测密度与0.997 g/cm³晶体学密度。"
        "范德华半径来源：[Bondi 1964](https://doi.org/10.1021/j100785a001)。",
        "",
        "## 第4步：GFN2-xTB预优化与再次去重",
        "",
        "计算后端为tblite 0.7.0的GFN2-xTB，ASE只负责构建、优化和记录。"
        "这是实际的紧束缚量子化学能量/受力计算；本批调用tblite库接口，未调用独立xtb.exe。"
        "GFN2-xTB是[tblite官方支持的内置方法](https://tblite.readthedocs.io/en/latest/spec/methods.html)。"
        "先用2×2×6精确复制胞作较粗预优化（力阈值0.03 eV/Å），随后每个起点都用2×3×8（576原子）"
        "继续至最大受力<0.01 eV/Å；粗阶段不按能量删除候选。"
        "晶胞固定，优化位移和力只按目标Pnam的8个操作投影；不额外锁住偶然出现的更高对称群。",
        "",
        "2×3×8是Γ点重复超胞计算，不是显式k点设置。超胞能量除以复制数，等价原子力取平均；"
        "每次记录原始力、投影力与副本偏差。末态再冷启动、把accuracy从0.01收紧至0.001检查能量和力。"
        "电子展宽300 K只属于电子占据设置，不代表晶体做了300 K分子动力学。",
        "",
        "加密检查曾发现15°起点的第一版终态虽满足采样变化阈值，"
        "但3×4×10下最大力为0.01152542 eV/Å，超过0.01目标。"
        "因此对该近重复组的15°与40°两个起点，从各自保存终态续优化，"
        "把2×3×8力目标收紧至0.003 eV/Å；随后按新终态重新分组和复核。"
        "旧结果和未通过记录保留在历史attempt及candidate_sampling_checks中，"
        "决策依据见audits/refinement_decision.json。",
        "",
        f"当前已预优化{summary.get('preoptimized_count', 0)}/16个；"
        f"终态仍通过几何过滤{summary.get('accepted_geometry_count', 0)}个；"
        f"检出方向不稳定{summary.get('directionally_unstable_count', 0)}个，"
        f"保留{summary.get('retained_count', 0)}个进入几何/能量合并；"
        f"当前合并后{summary.get('geometrically_distinct_final_count', 0)}组。"
        + ("全部完成后才把分组视为本批次最终结果。" if summary.get("pending_or_failed") else ""),
        "",
        "终态去重需要同时满足：每CH₂能量差≤1 meV、周期逐原子RMS≤0.02 Å、"
        "最大位移≤0.05 Å，并确认化学组成、目标空间群与原子映射。"
        "对8种正交晶格轴符号变换和每个碳锚点，交替优化整体平移与同元素原子指派，"
        "比较最小像距离。能量接近本身不能证明同一结构；所有成员与代表的映射均保留。"
        "这是多起点局部匹配，不是任意晶体匹配的全局最优证明。",
        "",
        "本次几何合并阈值较严格，受力阈值允许的优化残差可能使同一势阱中的末态仍分为几组。"
        "因此这些是按已声明容差得到的代表组，不把组数解释成不同晶型数或不同能量极小值数。",
        "",
        "C—C最终值逐项记录：严格的1.530–1.540 Å判定与“约1.53 Å”的化学合理性判定分开保存。"
        "优化若略低于1.530 Å不会被悄悄改写成严格区间通过；本次宽松化学检查界限为1.50–1.57 Å。",
        "",
    ]
    if summary.get("directionally_unstable_count"):
        lines += [
            "**0°高对称驻点的额外检查：** 在保持Pnam、晶胞和全部周期键长的情况下，"
            "两链按对称关系分别转动±1°，两个方向都比原结构降低0.06528234 meV/CH₂。"
            "因此该结构虽达到受力阈值，沿此方向仍不稳定，标为`unstable_stationary`并排除出保留池。"
            "原结构、优化日志及两侧单点结果全部保留。此检查不等于对其他结构完成声子或全Hessian验证。",
            "",
        ]
    if summary.get("candidates"):
        lines += ["表中链平面角均相对于+a轴，按等价关系折叠到0–90°。", "",
                  "| 起点角 | 状态 | E (eV/CH₂) | 最大力 (eV/Å) | 终态链角 | 分组 |",
                  "|---:|---|---:|---:|---:|---|"]
        for row in summary["candidates"]:
            def value(key, fmt):
                v = row.get(key)
                return format(v, fmt) if v is not None else "—"
            lines.append(f"| {row['initial_angle_deg']}° | {row['status']} | "
                         f"{value('energy_eV_per_CH2','.8f')} | {value('max_raw_force_eV_A','.6f')} | "
                         f"{value('final_chain_plane_angle_deg','.3f')} | {row.get('cluster_id','—')} |")
    lines += [
        "",
        "预优化得到的是固定胞、指定对称性下的驻点；未做声子/Hessian验证，不能把每个驻点都称为稳定晶型。"
        "新的结构不能直接继承旧P1候选的收敛结论，额外密采样检查见单独记录；缺失时保持待检查。",
        "",
    ]
    if dense:
        lines += [
            "对每个最终代表的同一几何，从2×3×8（576原子）加密到3×4×10（1440原子）做独立单点。"
            "同时检查每CH₂能量变化<1 meV、最大逐原子力矢量变化<0.005 eV/Å，"
            "以及密采样下最大力<0.01 eV/Å。结果只针对这些代表，不自动外推给每个组成员。",
            "",
            "| 代表组 | 候选 | ΔE (meV/CH₂) | 最大ΔF (eV/Å) | 密采样最大力 (eV/Å) | 通过 |",
            "|---|---|---:|---:|---:|---|",
        ]
        for row in dense["results"]:
            lines.append(f"| {row['group_id']} | {row['candidate_id']} | "
                         f"{row['energy_change_meV_per_CH2']:.6f} | "
                         f"{row['max_force_vector_change_eV_A']:.6f} | "
                         f"{row['dense_fmax_eV_A']:.6f} | {'是' if row['passed'] else '否'} |")
        lines += [""]
    lines += [
        "## Tashiro XRD基准的实际用途",
        "",
        "原PDF Table I的59行d(obs)已独立回读核对，按刊载值去重为45个观测峰位。"
        "原文是纤维衍射，不是这里重新取得的完整粉末连续谱。"
        "原文Eq.11–13按Fcalc相关权重拆分重叠强度，因此本次不将这些表列Fobs当作独立完整强度基准，"
        "也没有用精修坐标给候选提供答案。",
        "",
        f"{xrd['candidate_count']}个初态按指定Pnam计算的允许峰位置比较完全相同。"
        f"描述性阈值ΔQ≤{xrd['q_tolerance_A_inv']} Å⁻¹下，"
        f"有{xrd['candidates'][0]['observations_within_tolerance']}/45个观测峰接近允许位置，"
        f"RMS ΔQ={xrd['candidates'][0]['rms_delta_q_A_inv']:.8f} Å⁻¹。"
        "这里采用观测峰到最近允许位置的多对一比较；阈值不是实验标准差，不是精修R因子，也不是验收线。"
        "d(obs)=1.4694 Å是唯一超出这个描述性阈值的项目，保留其残差，没有删除。",
        "",
        "这里计算的是目标Pnam允许的晶格反射位置，不是由原子结构因子得到的强度谱。"
        "偶然更高对称性可能带来的额外消光未纳入该位置表。",
        "",
        "**同一固定晶胞和Pnam下，峰位置无法独立区分setting角，因此本次不产生XRD堆积排名。**"
        "若要以实验强度鉴别堆积，需要可靠的未拆分积分强度或可追溯原始衍射数据。"
        "当前结果可以用于下一阶段DFT校准，但不能宣称已由XRD确认唯一原子结构。",
        "",
        (f"最终{xrd_count}个代表也已重新核对目标Pnam操作和同一套45个峰位；"
         "结果仍只有一套晶胞/空间群允许位置表，没有产生堆积排名。"
         if final_xrd and (xrd_count := final_xrd["candidate_count"]) else
         "最终代表的XRD位置复核尚待完成。"),
        "",
        "## 文件在哪里",
        "",
        "- `POSCAR_list_initial.txt`：16个初态POSCAR的完整路径。",
        "- `POSCAR_list_final.txt`：全部预优化结束后，各保留代表POSCAR的完整路径。",
        "- `candidate_catalog.csv`、`generation_log.jsonl`：288个参数组合的通过/淘汰/重复原因。",
        "- `candidates/cand_XXXX/input.cif`：3个独立位点+8个Pnam操作的CIF；同目录`POSCAR`为完整12原子晶胞。",
        "- `candidates/cand_XXXX/geometry_validation.json`：密度、键长、周期拓扑、最近接触和对称性。",
        "- `candidates/cand_XXXX/relaxation/attempts/`：每个预优化的输入、原生输出、逐步轨迹、力、电子精度复核及终态CIF/POSCAR。",
        "- `optimization_catalog.csv`、`optimization_summary.json`：预优化与合并结果；未完成项目明确保留状态。",
        ("- `representatives_preview/group_XXX/`：当前预览代表；最终代表等待全部计算结束。"
         if summary.get("pending_or_failed") else
         "- `representatives/group_XXX/`：各组代表的CIF/POSCAR及成员对应。"),
        "- `representative_convergence/`：同一代表几何的更密周期采样、原生输出与能量/力变化。",
        "- `plots/initial_candidates_ab.png`：16个初始堆积的ab投影；`plots/optimization_angle_and_energy.png`：起点到终态的角度和相对能量。",
        "- `reference_xrd/`、`xrd_benchmark_initial/`、`xrd_benchmark_final/`：XRD来源与初、末态逐峰位置差。",
        "- `protocol.json`、各`source_snapshot/`、`audits/`：参数、程序版本、原始来源核验与独立检查。",
        "",
        "每套CIF与POSCAR均回读核对晶胞、组成和周期坐标，结果及SHA-256保存在file_validation.json。"
        "CIF使用CH₂作为式量单位，Z=4；完整晶胞为C₄H₈。部分最初生成尝试遇到90°浮点边界错误，"
        "已在量化开始前修复并重做完整几何检查，审计保留在audits/generation_attempt_01.json。",
        "",
        "## 与后续路线的衔接",
        "",
        "本步骤的范围是候选提出、几何排错、低成本局部优化和重复合并。"
        "DFT收敛、PBE-D3(BJ)/optB88-vdW/r²SCAN-D4比较、GFN2误差标定与3σ筛选尚属后续阶段。"
        "10¹–10³是通用候选规模目标；在本次强约束下初始有16个独特候选，"
        "优化后若汇聚为更少组则按实际数量保留，不人为制造重复结构。",
        "",
    ]
    (ROOT/"README_流程与结果.md").write_text("\n".join(lines), encoding="utf-8")
    print(ROOT/"README_流程与结果.md")


if __name__ == "__main__":
    main()
