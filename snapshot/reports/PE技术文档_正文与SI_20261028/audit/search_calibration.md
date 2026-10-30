# PE 候选生成 低成本模型与 GFN2 标定技术审计

审计日期：2026-10-28。本文所有定量内容来自本机已保存协议、结果 JSON、数据 CSV 和历史审计，未新增量子化学计算。本次对16点标定统计进行了独立复算，结果见 `search_calibration_verified.json`。路径根目录均为 `E:/PE_DFT/`。

## 1 本轮项目与状态

| 项目目录 | 目的 | 已完成内容及适用状态 |
|---|---|---|
| calculations/pe_gfn2_first_run_20261020 | 实现核验和单初猜验证 | tblite 0.7.0 首算、力/应力差分与周期平移验证；随后单个P1候选完成周期采样与26步固定胞优化 |
| calculations/pe_pnam_candidates_20261020 | 约束候选生成和局部优化 | 288参数组合、248通过接触过滤、16独特初态均完成GFN2；剔除1个方向不稳驻点，15条终态合并为4代表 |
| calculations/pe_gfn2_dft_calibration_20261021 | 同几何相对能误差标定 | 21套冻结几何、21 GFN2及29 DFT生产单点、29配对；16初态构成经验σ样本，4终态及1不稳对照单列 |
| reports/PE_GFN2_3sigma_20261021 | 标定数据和图表归档 | 67条计算记录索引、29配对、原始快照与复算脚本齐备 |
| calculations/pe_energy_models_benchmark_20261020 | 同16起点DREIDING→DFTB与GFN2公平比较 | 预检和输入归档；2026-10-21已取消，32次计划正式优化的启动数为0。不能用作完成的速度或成功率对照 |

## 2 候选模型和输入先验

主批模型固定a=7.432、b=4.929、c=2.545 Å，α=β=γ=90°，体积93.22927476 Å³。完整周期胞C4H8，12原子、4个CH2、两条沿c轴的无限链。空间群约束为Pnam（No.62，Hall 297，Hall符号 −P 2c 2n），使用8个操作，不混写成另一轴序的Pnma。晶胞与Pnam是Tashiro 1997的实验先验，不能声称本轮独立预测晶胞和空间群。

原子坐标由化学几何独立生成：C—C=1.54 Å，C—H=1.09 Å；c轴周期和C—C共同决定C—C—C角。两碳纯平移周期只兼容完整全反式链；t、g+、g−两扭转9组合中8个含gauche组合因周期闭合解析排除。该排除针对指定短周期，不是对长周期gauche缺陷的能量否定。

setting角是链平面在ab平面的投影相对+a的角度，不是链内扭转角。36个角0–175°、步长5°，乘4个半格横向链轴位置及2个碳镜面相位1/4、3/4，共288组合。第二条链由Pnam关联，角度不独立自由。生成器未读取精修原子坐标或实验setting角；Tashiro d(obs)只作外部峰位对照，未将Fcalc、拆分Fobs用于生成、几何过滤或优化。

直接来源：`calculations/pe_pnam_candidates_20261020/protocol.json`、`input_provenance.json`、`generate_candidates.py`、`audits/conformation_enumeration.json`。原始PDF为 `examples/PE_Tashiro_1997/sources/Tashiro_1997_original_Wiley.pdf`，来源记录SHA-256为4750738ed0c2e29a569d3b402d1a94d46db920f3f5f40f4ba953b0d937ad15f0；晶胞/链轴/波长/Pnam对应印刷1692页、PDF第16页。

## 3 过滤 去重与候选计数

非键接触采用Bondi半径C=1.70、H=1.20 Å，阈值0.85×半径和：C···C 2.89、C···H 2.465、H···H 2.04 Å；排除周期化学图距离1或2条键的近邻，其他跨胞像均检查。另有链间C···C≥3.4 Å条件。

| 环节 | 数量 | 减少原因 |
|---|---:|---|
| 枚举组合 | 288 | 36×4×2有限离散参数空间 |
| 几何通过 | 248 | 40个80–100°组合H···H过近 |
| 初态去重后 | 16 | 232重复；允许晶格点群、原点平移和同元素指派；保留0–75°每5°一个 |
| GFN2预优化完成 | 16 | 全部初态送算，没有粗能量淘汰 |
| 排除方向不稳后 | 15 | 0°高对称驻点沿保持Pnam的±1°两侧均降低0.06528234 meV/CH2 |
| 终态联合合并 | 4组 | 成员数12、1、1、1；能差≤1 meV/CH2、周期RMS≤0.02 Å、最大位移≤0.05 Å同时满足 |

初态晶体学密度均0.999359120 g/cm³，与外部结晶相晶体学参考0.997 g/cm³比为1.002366，位于0.95–1.05范围。固定组成/晶胞使密度完全不随堆积变化，只能用于输入合理性检查。0.997是Martin与Passaglia 1966的晶体学参考，非Tashiro同样品独立实测密度。其页码应221–224；protocol.json末页225是已记录笔误。

4组是规定容差下的代表，不能直接称4种稳定晶型；没有全Hessian或声子证明。几何合并阈值与力停止阈值共同影响组数，同一势阱内残差可造成多个代表。严格C—C 1.530–1.540 Å检查与宽化学安全范围1.50–1.57 Å分开；保留结果约1.529 Å，不能写严格窄区间通过。

直接来源：`candidate_catalog.csv`（本次验证288行）、`optimization_catalog.csv`（16行）、`optimization_summary.json`、`audits/zero_angle_curvature/20261020T083728_400937Z/summary.json`、`audits/candidate_counts_timing_and_historical_models_20261020.md`。

## 4 GFN2 实际实现与数值设置

能量/力引擎为tblite 0.7.0内置GFN2-xTB，ASE 3.29.0负责生成、调用及优化；未运行独立xtb.exe。原胞12原子，Γ重复超胞用于周期采样，不能把复制数写成显式k网格。

| 阶段 | 设置 |
|---|---|
| 候选粗优化 | 2×2×6重复胞，fmax 0.03 eV/Å |
| 候选精优化 | 2×3×8，576原子，fmax<0.01 eV/Å；固定胞，按目标Pnam 8操作投影 |
| 电子精度 | 优化accuracy=0.01；终态冷启动accuracy=0.001核验；电子占据300 K |
| 追加收紧 | 15° cand_0025、40° cand_0065各自从已保存终态继续；2×3×8 fmax=0.003 eV/Å |
| 代表加密 | 同几何3×4×10，1440原子；每CH2能差<1 meV，最大逐原子ΔF<0.005 eV/Å，加密自身fmax<0.01 eV/Å |
| 标定生产 | 冻结几何，2×3×8、accuracy=0.001、max_iterations=1000、4线程、中性闭壳层、电子温度0 K |

追加收紧的触发证据是旧cand_0025在加密采样下fmax=0.0115254200 eV/Å，虽然采样差通过但仍超过力目标；旧结果保留，新几何重新验收，不能拿旧几何密采样结果证明新几何。原胞能量=超胞能量/复制数，原胞力=等价副本力平均；同时保留原始/投影力与副本偏差。300 K仅为电子占据设置，没有300 K原子动力学，也未计算晶格自由能。

| 代表 | 候选 | 初始角/终态角 ° | 预优化能量 eV/CH2 | 原始fmax eV/Å | 加密ΔE meV/CH2 | 加密最大ΔF eV/Å | 加密fmax eV/Å |
|---|---|---|---:|---:|---:|---:|---:|
| group_001 | cand_0025 | 15 / 44.36834 | -86.1563225447 | 0.00165798 | -0.249259 | 0.004287 | 0.003914 |
| group_002 | cand_0113 | 70 / 45.51577 | -86.1562890974 | 0.00688889 | -0.244139 | 0.004351 | 0.008303 |
| group_003 | cand_0041 | 25 / 41.39659 | -86.1561227780 | 0.00499600 | -0.266332 | 0.004099 | 0.008054 |
| group_004 | cand_0081 | 50 / 49.78783 | -86.1557713379 | 0.00634521 | -0.231993 | 0.004500 | 0.008397 |

上表为GFN2预优化代表，不是后续DFT优化终态；完整16行应使用optimization_catalog.csv作为SI。

## 5 标定样本和方法

21冻结几何=16初态+4代表终态+1排除的0°驻点对照。生产DFT：PBE-D3(BJ) 21、optB88-vdW 4、r²SCAN-D4 4；生产GFN2 21，共29能量配对。67条索引包含37 DFT和30 GFN2（含诊断、复用），其中65新计算、2个先前DFT点复用；不应将67全称为本轮新生产点。65=63队列计算+2同批温度pilot。

29个生产DFT统一H档1800 eV、k=4×8×12、gpts=112×72×44；U档2000 eV、5×10×15、128×80×48仅4几何PBE迁移核验。L档1400 eV、3×6×9、96×60×36是前阶段对照。GPAW 26.7.0、gpaw-data 1.0.1 C/H PBE PAW、ASE 3.29.0、dftd3 1.6.0、dftd4 4.2.0。实际PAW来源不是SSSP。DFT能量/密度/本征态阈值为1e−9 eV/electron、1e−9 electron/electron、1e−11 eV²/electron；DFT占据展宽0。每DFT任务4 MPI进程，最多2并发。

同一配对须同坐标同晶胞；各方法绝对零点不同，误差先在各方法内相对R0=final_group_001归一化：

ΔE_i^M = 1000 [E_i^M − E_R0^M] / 4 （meV/CH2）

e_i = ΔE_i^DFT − ΔE_i^GFN2；b=mean(e)；σ=std(e, ddof=1)；预测 ΔEhat=ΔE_GFN2+b。

后续analysis_policy.json明确将主统计解释改为有适用范围的经验σ；冻结计算协议未改。四个连续角度块留出分析保留为补充诊断，但主报告σ是同批16初态样本标准差，不能写作外部预测标准差。没有拟合斜率；R²仅是斜率固定1、均值偏差校正的同批描述量。4终态和1不稳对照不混入样本数。

| 统计量 | 独立复算结果 |
|---|---:|
| n | 16 |
| b meV/CH2 | +4.6383338759 |
| σ meV/CH2 | 1.9206975072 |
| 3σ meV/CH2 | 5.7620925215 |
| 未校正RMSE meV/CH2 | 4.9972645143 |
| 均值校正RMSE meV/CH2 | 1.8597073646 |
| 最大绝对中心残差 meV/CH2 | 2.5264521760 |
| 同批均值校正R² | 0.9831347079 |

适用域为PE、固定Tashiro胞、C4H8、同一全反式链构建、初态setting角0–75°每5°；不能自动外推不同晶胞、长周期缺陷、其他聚合物或4个预优化终态。3σ是经验能量保护带，不代表99.73%正确概率。W低能窗口未指定，自动能量淘汰数0，4代表全保留。16→4来自几何/能量联合合并和方向检查，不是3σ筛掉12个。

直接来源：`calculations/pe_gfn2_dft_calibration_20261021/{protocol.json,analysis_policy.json,samples.json,tables/paired_energy_results.csv,audit/final_delivery_check.json}`；便于论文使用的数据在 `reports/PE_GFN2_3sigma_20261021/data/initial_16_energy_calibration.csv` 和 `report_metrics.json`。

## 6 同几何四终态比较与数值核验

| 代表组 | GFN2 | PBE-D3(BJ) | optB88-vdW | r²SCAN-D4 |
|---|---:|---:|---:|---:|
| 001 | 0 | 0 | 0 | 0 |
| 002 | +0.0334473 | −0.0235036 | +0.0054752 | +0.0663434 |
| 003 | +0.1997666 | +0.4388176 | +0.4581107 | +0.2296739 |
| 004 | +0.5512067 | +0.3085682 | +0.4275698 | +1.2653287 |

单位meV/CH2，相对各方法自己的group_001；几何均为GFN2冻结末态。此表不是各泛函独立优化后最低能结构对照，三泛函细排序不同，不支持唯一获胜晶型。001/002差很小；能量误差小不能代替DFT受力合格。四终态PBE/GFN2逐原子力矢量差RMSE约0.4684–0.4699 eV/Å，最大应力分量差25.7816–26.3572 kbar。

| 核验 | 几何数 | 实测最大变化 | 解释 |
|---|---:|---|---|
| DFT H→U | 4 | 0.00178277 meV/atom绝对能；相对R0能0.000782107 meV/CH2；力1.48521e−5 eV/Å；应力0.00110920 kbar | 四项既定门槛均通过，只是抽查点变化 |
| GFN2 2×3×8→3×4×10 | 6 | 0.462302 meV/CH2；力0.00450357 eV/Å；相对R0变化0.213043 meV/CH2 | 能量/力门槛通过；应力最大变化0.272684 kbar单列，不能宣称全部应力≤0.1 kbar |
| GFN2 0→300 K占据 | 3 | 能量4.26326e−11 meV/CH2；力1.29617e−15 eV/Å；应力4.36629e−12 kbar | 仅电子温度诊断 |

直接来源：`reports/PE_GFN2_3sigma_20261021/data/{endpoint_4_method_energies.csv,forces_stress_29_pair_errors.csv,numerical_gate_checks.csv}`。原始判定在calibration目录audit的三个gate JSON。

## 7 初次接口验证与单候选预研

GFN2最小胞Γ首算能量−339.9276806710 eV/胞、−84.9819201678 eV/CH2，fmax=1.99628809 eV/Å，9轮电子迭代；不能用这个粗采样绝对能与后续加密结果相减解释优化收益。沿链1、2、4、6倍Γ重复，能量依次−84.981920168、−85.996643284、−86.166501659、−86.172374194 eV/CH2；4→6仍差5.872535 meV/CH2。

首算的第11原子y力差分误差4.415e−7 eV/Å、yy应力差分误差3.582e−6 eV/Å³、平移能变化2.274e−13 eV、平移最大力分量变化8.338e−14 eV/Å、总力矢量模5.169e−15 eV/Å，均通过本项目既定门槛。这是选定分量的实现一致性，不是对实验精度证明。

DFTB+25.1/tblite0.3.0虽然通过其官方参考测试，但在该PE初猜的力和应力差分未过本项目门槛；调整差分步长、SCF未解决，原因未定位，后端隔离。后续采用tblite0.7.0通过上述检查。

单P1候选随后用2×3×8进行26步ASE BFGS固定胞优化，maxstep=0.06 Å，fmax目标0.01 eV/Å，无额外空间群约束。同采样初/终态能量−344.557771991 / −344.608121771 eV/胞，下降0.050350 eV/胞=12.587 meV/CH2；fmax从0.605896降到0.009876 eV/Å。更密联合采样ΔE=−0.302357 meV/CH2、maxΔF=0.003617、fmax=0.007675 eV/Å。终态C—C 1.5291067–1.5293844 Å，C—H 1.0920436–1.0928022 Å，最近链间原子距离2.633006 Å。

来源：`calculations/pe_gfn2_first_run_20261020/result_summary.md`，`runs/20261020T061353_138002Z/summary.json`，`convergence_relaxation/runs/20261020T063341_922507Z/{result_report_zh.md,summary.json,protocol.json}`。

## 8 历史搜索与公平比较的边界

| 历史项目 | 后端与数据 | 数量及验收 |
|---|---|---|
| pe_xrd_seed_search_20260922 | 45 d(obs)；OpenMM CUDA DREIDING零电荷→DFTB+/3ob；位置损失与能量Pareto | 64初始+32反馈=96提案；17 invalid，12力场completed、67unconverged；6 DFTB中5完成1失败；5终态固定胞原子力≤0.05 eV/Å，无应力验收；结构对照0/5获严格等价证据，匹配结果indeterminate |
| pe_xrd_rerun_v2_20260922 | CUDA DREIDING C_3/H_刚性链/晶胞探索→DFTB3/3ob-3-1 | 533种子尝试→64有效初始+32反馈=96；90 GPU评价；6 DFTB完整记录；原胞Γ分支0/6联合力/应力合格 |
| pe_xrd_rerun_v2_kcorrected_20260922 | 同6父候选，DFTB显式4×4×8修正，无新GPU搜索 | 6投影力<0.05 eV/Å，但允许应力1090.23–3609.98 bar>50 bar，0/6联合通过；末态6×6×16单点复核均通过采样一致性，最大应力变化1.188417 bar，联合验收结论不变 |
| MASTERS 216网格基准 | 规范记录P1两链网格，LAMMPS+AIREBO 3×5×9 | 216组合→50严格等价代表；本次未重审该原始批次，仅作为历史规范背景 |

旧6名额是每晶格家族3个的预设DFTB预算，不是能量自然筛到6个。旧DFTB不等于本轮GPAW全Kohn–Sham DFT。不同搜索域/模型/收敛标准下，96→6不能与16→4当同一筛选率比较。DREIDING有Lennard-Jones项；旧DFTB3无额外色散，不应称整条DREIDING→DFTB路线无色散。

时间：旧seed 290.520341 s，v2原运行392.170692 s，k修正新增268.445095 s；本轮GFN2 18个attempt墙钟和16393.146634 s，首末跨度3341.462721 s（55.69 min），其中16初始attempt和14912.611444 s、追加2次1480.535190 s。早期cand_0041用4线程，其余2线程。独立密采样7文件按SHA去重为5次独立结果（含旧失败），墙钟和2154.252024 s；±1°两侧单点68.520358 s。逐任务墙钟和、并发日历跨度、CPU时间不可混用，无公平同输入同收敛目标速度比。标定批从最早pilot到结束4.16975087 h，含并发，也不是模型速度基准。

新公平对照原拟16×2=32正式路线，包含DREIDING→DFTB3与GFN2，统一Pnam固定胞、原始力0.01 eV/Å、200总几何步预算、4线程、AB/BA交替。该项目已取消；state.json production_attempts_started=0。因此只能在SI中记录已保存预检/取消状态，不能写完成32次对照，不能自动重启。

直接来源：`calculations/pe_pnam_candidates_20261020/audits/candidate_counts_timing_and_historical_models_20261020.{md,json}`；两v2目录的`run_summary.json`；`reports/pe_xrd_seed_search_20260922/{report_summary.json,archive_audit.json,reference_validation.json}`；`reports/pe_xrd_rerun_v2_20260922/k_sampling_diagnostic/final_states/result.json`；`calculations/pe_energy_models_benchmark_20261020/{README.md,state.json,CANCELLED_BY_USER.json,protocol.json}`。

## 9 正文与SI数据建议

正文适合放：候选数量流程、同几何能量定义、16初态相对能曲线与误差、b/σ/3σ及适用域、4代表三泛函相对能对照。必须同时说明晶胞/Pnam先验、经验σ非外推覆盖保证，以及未执行能窗淘汰。

SI应列：288组合目录及过滤理由；16初态/末态全坐标和POSCAR/CIF；每个优化attempt轨迹/SCF/力；0°±1°反例；4代表加密采样；21冻结结构哈希；29配对的绝对/相对能和力/应力；67全部记录索引（复用和诊断标签）；16完整标定数据；数值门槛与结果；软件和PAW版本/哈希；历史未收敛及取消分支说明；可复算脚本和文件清单。已有数据表如下。

| 可用附件 | 直接路径 |
|---|---|
| 全288候选与16优化结果 | calculations/pe_pnam_candidates_20261020/candidate_catalog.csv；optimization_catalog.csv |
| 结构和轨迹 | calculations/pe_pnam_candidates_20261020/candidates/；representatives/；representative_convergence/ |
| 16点标定 | reports/PE_GFN2_3sigma_20261021/data/initial_16_energy_calibration.csv |
| 4代表同几何方法对比 | reports/PE_GFN2_3sigma_20261021/data/endpoint_4_method_energies.csv |
| 全21 GFN2、29 DFT | reports/PE_GFN2_3sigma_20261021/data/all_21_gfn2_production_energies.csv；all_29_dft_production_energies.csv |
| 配对和力应力 | reports/PE_GFN2_3sigma_20261021/data/all_29_energy_pairs.csv；forces_stress_29_pair_errors.csv |
| 67索引和数值检查 | reports/PE_GFN2_3sigma_20261021/data/all_67_calculation_index.csv；numerical_gate_checks.csv |
| 来源归档 | reports/PE_GFN2_3sigma_20261021/data/source_file_manifest.csv；archive_path_map.csv；raw/calculation_snapshot/；raw/reused_DFT/ |
| 已有图 | reports/PE_GFN2_3sigma_20261021/figures/fig01_initial_angle_energy_residual.png；fig02_endpoint_method_energy.png；fig03_initial_energy_mean_correction.png |
| 原始复算 | reports/PE_GFN2_3sigma_20261021/scripts/analyze.py |

报告不能用后续DFT优化或变胞终态替换这里的GFN2冻结几何后仍沿用本阶段σ。不同阶段终态应使用明确的stage、method、geometry_id和reference_id字段关联。
