# GPAW 环境与 PE 三方法验证及收敛报告

日期：2026-10-20。阶段版本：依据 `protocol.json` 第5版及本次更新时已完成的记录。前一版全文保留为[第4版阶段快照](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/reports/GPAW环境与PE三方法验证及收敛报告_v4快照_20261020.md)。

**最新执行范围：继续当前固定结构复核与晶胞验证工作，高精度参考优化已恢复运行。** 此前把用户“晶胞验收不急”误解为可以停止参考优化，属于执行判断错误，不能写成用户要求中止。参考档于13:22:52 UTC保存后中止的历史保留，原尝试仍为未完成、未通过。候选档自然结束并以退出码0释放资源后，参考档新尝试于13:30:21 UTC从父运行的完整评价0003几何重新求解电子态并启动优化；恢复过程没有中断其他活跃量子计算。候选档的科学验收仍待独立核查，不能据进程结束宣称两档晶胞差 <0.5% 已通过。固定几何阶段可先整理结果，完整晶胞验收仍按原判据执行。

完整参数发布的原有门控保持不变；在所有必需条件满足前，不生成宣称“全部通过”的参数锁。后续如提供仅适用于已验证固定几何范围的参数模板，必须单独写明适用范围和待补验证；本次未编写或发布该模板。正式第3级 GFN2 粗筛与 3σ 误差标定尚未启动，已有的 GFN2 单点、预优化与去重属于准备阶段，不等同于该正式阶段已经开展。

**三方法的代表性力及六分量应力验证已通过，PBE 在固定 Tashiro 参考几何上的截断能、k点和网格扫描也已分别通过。** 当前组合候选为1400 eV、Gamma 中心 k=(3,6,9)、网格96×60×36及 XC stencil=4。这仍是待联合验证的候选设置：三方法×两套固定结构的三精度复核正在进行；晶胞候选档正常结束待核查，参考档已恢复运行，生产参数尚未锁定。

已通过的技术验证包含111个成功单点、27组导数比较，以及另外6个加严电子自洽阈值端点构成的3组复核。它们回答“实现给出的能量、力和应力是否一致”；正式收敛回答“需要多高的数值精度”。三种泛函对真实 PE 的准确程度，以及它们之间的结构/能量分歧，需要在数值参数通过完整验收后再比较。

## 1. 当前完成范围

| 环节 | 当前状态 | 可以支持的结论 |
|---|---|---|
| 独立环境、依赖与 PAW 文件核验 | 已完成 | 三方法使用同一套已归档 C/H PBE PAW；所需 Libxc、libvdwxc、D3/D4 和 MPI 接口可用。 |
| 三方法 PE 单点能量、全部原子力、六分量应力 | 已完成 | 三种实现均返回有限结果且完成电子自洽；DFT 与外加色散贡献分别保存。 |
| 两步长导数比较 | 已完成，27/27组通过 | 每种方法的3个代表性力分量和全部6个应力分量均满足预定门槛。不是全部36个原子力分量逐项差分。 |
| 电子自洽阈值敏感性复核 | 已完成，3/3组通过 | 在指定敏感剪切端点上，自洽阈值收紧10倍后的能量和差分应力足够稳定。 |
| PBE 固定参考几何的截断能、k点和网格扫描 | 三项分别已通过 | 得到1400 eV、k=(3,6,9)、网格96×60×36组合候选；更高精度联合验证仍需完成。 |
| 4/8 MPI 进程一致性 | 指定参考点已通过 | 已检点的能量、力和应力在预定容限内一致；不能外推为所有条件的独立验证。 |
| 三方法及独立 GFN2 几何的收敛迁移检查 | 进行中 | L/M/H三精度、两套冻结结构共18个计划单点，尚不能将 PBE 主扫描直接推广。 |
| 两结构相对能的数值收敛 | 进行中 | 各方法的相对结构能须在三精度下满足≤0.5 meV/atom，目前未形成完整通过结论。 |
| 两档精度变胞优化、晶胞长度差 <0.5% | 候选档正常结束待独立核查；参考档新尝试运行中 | 原参考尝试中止历史保留，新尝试从其完整几何恢复。候选进程退出码0不代替科学验收，两档胞长及冷启动整套验证仍未完成。 |
| 生产参数锁定 | 未锁定 | 只在全部必需检查完成后记录有版本和适用范围的设置。 |
| 正式第3级 GFN2 粗筛与 3σ 标定 | 尚未启动 | 已有 GFN2 单点、预优化和去重为准备工作，本阶段不进入正式粗筛或误差标定。 |

状态依据为[当前方案第5版](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/protocol.json)、[最终导数结果](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/validation_audit/derivative_results_deriv_final.json)、[自洽噪声复核](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/validation_audit/scf_noise_validation.json)及[固定几何收敛明细](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/tables/convergence_summary.json)。后续计算会继续更新各进度表，本报告只对本次列明的已完成证据作结论。

## 2. 环境与统一 PAW 来源

计算运行于 WSL Ubuntu 24.04 的独立环境，采用 GPAW 26.7.0 的 legacy CPU 平面波实现。环境路径为 `/opt/pe_dft/envs/pe-gpaw-20261020`。本轮没有把 GPU 性能或其他软件环境的验收结果当作该实现的验证证据。

| 组件 | 实际版本/实现 | 用途 |
|---|---|---|
| GPAW | 26.7.0，legacy CPU、平面波 | 自洽 DFT、PAW、原子力和晶胞应力 |
| ASE / spglib | 3.29.0 / 2.7.0 | 结构与计算接口；后续对称约束、变胞优化 |
| Libxc | 7.1.2 | 半局域 PBE、optB88 相关组成和 r²SCAN |
| libvdwxc | 0.5.0，已确认应力接口可用 | optB88-vdW 的非局域相关 |
| simple-dftd3 / dftd3-python | 1.6.0 | 周期 D3(BJ)，含 ATM 三体项 |
| dftd4 / dftd4-python | 4.2.0 | 周期 D4，含 ATM 三体项 |
| PAW 数据包 | gpaw-data 1.0.1 | 三方法共用原生 GPAW PBE C/H 文件 |
| Python / NumPy / SciPy | 3.12.15 / 2.5.4 / 1.18.1 | 数值运行环境 |
| OpenMPI | 5.0.11 | 并行运行；显式使用 GPAW 的 cgpaw MPI 后端 |

最初的 GPAW 预编译包缺少需要的绑定，最终使用同版本官方源码重新构建。归档比较了2408个官方 Python/C/C++ 源文件，未发现科学源码修改；131个 conda 下载包的文件哈希与包元数据一致。实际安装为 editable 方式，运行依赖已归档的固定源码目录。初次依赖求解、预编译包和打包失败记录均保留。环境本身的构造检查与后续 PE 数值验证分别归档，详见[环境核验记录](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/environment/environment_verified.json)和[安装说明](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/environment/环境安装与接口验收.md)。

三种方法的实际结果均记录了同一组 PAW 文件路径与 SHA-256：

| 元素 | 价电子 / 冻结芯电子 | 压缩 PAW 文件 SHA-256 |
|---|---:|---|
| C | 4 / 2 | `5a769dd14bb5d462c57fadc4dd8d3b7e70428b90701eca0d20b6f0666495cdfd` |
| H | 1 / 0 | `43da6533cea08dec1ee7861140be503c9945863d6d8c7791e01e55e2161a88b8` |

这是 **GPAW 原生 PBE PAW**，不是原计划的 SSSP UPF 数据，不能宣称获得 SSSP 文件认证或照搬其推荐截断能。C/H 数据具有 meta-GGA 所需的芯层动能密度数组；H 无冻结芯，其相应零值合理。实装 GPAW 的 Libxc 绑定支持并执行关闭 FHC enforcement 的默认设置。这些检查解决接口和数据兼容性，不等于证明 PAW 数据集对所有聚合物和泛函没有误差。[方法与 PAW 审查](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/method_audit/METHOD_AUDIT.md)

## 3. 三种方法如何实现

| 方法 | 本轮实际定义 | 色散设置 |
|---|---|---|
| PBE-D3(BJ) | 自洽 PBE 加 D3 的 BJ 阻尼修正 | PBE 对应阻尼参数；ATM 开启；二体/三体/CN 半径为60/40/60 Bohr。 |
| optB88-vdW | optB88 交换、相应局域相关及 libvdwxc 非局域相关 | 非局域项已包含在 GPAW 的 DFT 能量中，不再叠加 D3 或 D4。 |
| r²SCAN-D4 | 自洽 `MGGA_X_R2SCAN+MGGA_C_R2SCAN` 加 D4 | r²SCAN 对应参数；ATM 开启；二体/三体/CN 半径为60/40/30 Bohr。 |

D3、D4 的二体和三体截断平滑宽度均为0.05 Bohr，三轴周期性均开启。ASE 接口接收的长度单位是 Å，实际运行时经过 Bohr 到 Å 的明确换算。D3 的参数字典显式包含 PBE 方法名和 ATM 状态，未采用仅凭名称推断三体项的做法。完整拟合参数、来源与调用证据见[方法定义](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/method_audit/method_specification.json)。

能量导数统一使用与力一致的电子自由能 `energy_free_eV`。GPAW 的展宽外推能量另存为 `energy_extrapolated_eV`；本轮绝缘 PE 使用零占据展宽，仍保留两个字段以免以后混用。该电子自由能不包含声子、零点能或有限温度晶体自由能。对两种外加色散方法，总自由能、力和应力分别为 DFT 与色散对应量之和。后续优化的组合接口补充了色散模块缺失的 `free_energy` 属性声明，保持同一物理定义。

应力使用 ASE 的拉伸为正约定，顺序为 xx、yy、zz、yz、xz、xy，原始单位 eV/Å³，同时保存 kbar。压缩为正的静水压力为应力迹的负三分之一。C₄H₈ 原胞含12个原子、4个 CH₂；按原子和按重复单元归一化时分别除以12和4。不同泛函的绝对总能量不能直接用来判定哪种方法更准确。

## 4. PE 参考结构与独立候选

主数值扫描采用 Tashiro 1997 实验精修的正交 PE，晶胞为 **a=7.432 Å、b=4.929 Å、c=2.545 Å**，c 为链轴。按原文 Table II 的 C、H(A)、H(B) 三个位点与 Pnam 等效位置展开，得到12原子 C₄H₈ 晶胞。项目已独立核对 Hall 297 的8个操作，与归档结构的最大周期位置差为约2.48×10⁻¹⁵ Å。[参考结构与原文核查](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/reference_audit/README_参考结构与来源.md)

Pnam 和常见标准标签 Pnma 属于第62号空间群的不同轴设置。本轮保留原始 Pnam 轴，显式核对 Hall 297；默认软件输出的 Pnma 字符串本身不能证明发生了换轴。原有 P1 文件头、12原子全晶胞 CIF 保留，并新增正式声明 Pnam #62、Hall 符号 `-P 2c 2n` 的[对称 CIF](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/reference/Tashiro_Pnam.cif)：仅列3个不等价4c位点及8个对称操作，避免把全晶胞原子重复展开。写后由 ASE 重新读取得到12原子，与原 POSCAR 的最大周期位置差为1.77636×10⁻¹⁵ Å，晶胞保持一致；来源哈希、轨道闭合和逐原子配对见[导出审计](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/reference/Tashiro_Pnam.audit.json)。这是文件表达和几何验证，不是重新精修或量子计算。

原文氢位置受 C–H 长度和 H–C–H 角约束，碳位置也有实验不确定度，所以实验精修模型是有来源的比较参考。它不必在 DFT 势能面上受力为零。本轮用于数值校准的精修坐标，也不应与“从重复单元生成候选”的坐标来源混淆。

后续迁移检查使用本项目独立化学构建、经 GFN2 固定胞优化的 group_001（cand_0025）。该候选使用实验晶胞与 Pnam 先验，原子坐标生成没有引用 Table II。它与实验参考坐标不同，不能混进同一条截断能或 k 点差值序列。

| 用途 | 文件 | SHA-256 |
|---|---|---|
| 主扫描实验参考 | [POSCAR](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/reference/POSCAR) | `e54173a343c8a1bac20d048fdeaa99241986ee699e9f4591d65fc0985b76ef73` |
| 独立候选迁移检查 | [POSCAR_GFN2](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/reference/POSCAR_GFN2) | `f3db63cbf993da085a1c5bb4f772eb30d2e4debf4ae80fb474e738db79599bd0` |

## 5. 已完成的能量—力—应力验证

为了检验非零剪切应力，验证采用从参考 PE 施加一般剪切及少量原子位移得到的探针，且关闭点群约化。该结构用于检验导数关系，不用于候选排名或与实验结构拟合。

最终探针使用600.012 eV、k=(1,2,3)、显式粗网格64×40×24及 XC stencil=4。三个代表性力分量为 C 的 atom0-x、H(A) 的 atom4-y、H(B) 的 atom8-z，位移步长0.001和0.002 Å；每个应力分量使用10⁻⁵和2×10⁻⁵两档应变。工程剪切的两个对称非对角元为 ±h/2，应力差分为 `[E(+h)−E(−h)]/(2hV)`，V 为探针中心体积。

每种方法包含1个中心、12个原子位移端点和24个应变端点，合计37点；三方法共111点。每种方法3组力与6组应力，共27组。全部端点核对了同一几何变换、PAW、方法及非扫描参数，并合并全部 MPI rank 的波函数和密度/Hamiltonian G 集合、网格与 k 点记录，未把基组变化当作通过。

下表汇总**总能量对应总力/总应力**的最不利误差，均在两个步长及全部已检分量中取最大值：

| 方法 | 最大力差分误差 (eV/Å) | 力两步长最大差 (eV/Å) | 最大应力差分误差 (kbar) | 应力两步长最大差 (kbar) | 结果 |
|---|---:|---:|---:|---:|---|
| 门槛 | 0.001 | 0.0005 | 0.1 | 0.05 | 四项均须通过 |
| PBE-D3(BJ) | 0.000128498 | 0.0000986308 | 0.00178094 | 0.000774909 | 通过 |
| optB88-vdW | 0.000121077 | 0.0000969985 | 0.00410831 | 0.0000697203 | 通过 |
| r²SCAN-D4 | 0.000129788 | 0.000101631 | 0.000801533 | 0.00000290628 | 通过 |

DFT、外加色散及总量的逐分量结果另列于[最终导数明细表](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/validation_audit/derivative_results_deriv_final.md)和[CSV](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/validation_audit/derivative_results_deriv_final.csv)。optB88 的非局域相关已计入 DFT，“单独色散项为零”只是数据拆分约定。

### 电子自洽阈值敏感性

在三方法的 xy 剪切 ±10⁻⁵ 两个端点，将能量、密度和本征态三项自洽阈值同时收紧10倍。几何、实际 G 集合、PAW、程序与其余参数保持一致，分别检查 DFT 和总量。

| 方法 | 两端最大总自由能变化 (eV/胞) | 总差分应力变化的绝对值 (kbar) | 结果 |
|---|---:|---:|---|
| 门槛 | 1×10⁻⁸ | 0.02 | 两项均须通过 |
| PBE-D3(BJ) | 1.40687×10⁻¹² | 4.15183×10⁻⁷ | 通过 |
| optB88-vdW | 3.97904×10⁻¹³ | 5.00662×10⁻⁷ | 通过 |
| r²SCAN-D4 | 1.53477×10⁻¹² | 7.81521×10⁻⁷ | 通过 |

这是指定端点对自洽阈值收紧的稳定性证据，不是精确自洽误差的真值。完整原子力和六分量应力的变化也已保存，未另加未经约定的验收门槛。[自洽敏感性明细](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/validation_audit/scf_noise_validation.md)

## 6. 误差诊断与数值修正

前期未通过的结果全部保留，修正后的批次单独运行或以显式色散后处理文件归档。没有提高原定误差容限，也没有更换三种泛函。

| 发现 | 诊断与修正 | 证据范围 |
|---|---|---|
| 周期色散截断边界影响微小应变差分 | D3/D4 的两体/三体平滑宽度显式设为0.05 Bohr；D3 的 CN 半径由40增至60 Bohr，BJ 参数和 ATM 保持不变。 | 参考 PE 上只增加 CN 后，色散能变化+0.0121913 meV/atom、最大力矢量变化3.21124×10⁻⁵ eV/Å、最大应力分量变化0.00289920 kbar。旧 DFT 部分可原样保留，再在完全相同几何上重算色散；此过程不是重跑 DFT。 |
| optB88 非零剪切应力存在稳定的网格误差 | 固定其它条件，加密网格，再将 XC 密度梯度 stencil 从2改为4。 | yz 误差依次为0.437885、0.0922938、0.00391963 kbar，支持空间梯度离散误差的解释；不据此宣称 libvdwxc 有程序缺陷。 |
| 600 eV 下微小应变造成密度 G 向量换壳 | 仅依据几何和倒空间集合选择600.012 eV，再逐点核对波函数及密度/Hamiltonian 集合。 | 对该探针、指定网格、k 点和应变范围成立。600.012 eV 是导数检查的技术参数，**不是生产推荐截断能**，不能外推到变胞轨迹。 |

上述变更的完整数据与源码证据见[色散截断敏感性](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/method_audit/dispersion_cutoff_sensitivity.json)、[应力源码诊断](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/method_audit/stress_diagnosis/STRESS_SOURCE_DIAGNOSIS.md)、[探针截断集合检查](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/validation_audit/basis_safe_cutoff.json)及[与原计划的差异说明](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/reports/与原计划的差异及影响.md)。

## 7. 正式收敛：已通过部分与尚待完成部分

### 固定参考几何上的 PBE 主扫描

初始400—800 eV范围未达到能量和应力门槛，扩展至1200 eV后仍未形成稳定末端。原64×40×24网格序列保留为未通过记录。由于 GPAW 平面波密度描述同时随波函数截断改变，且受实空间网格表示上限约束，随后另建96×60×36网格序列，以1000、1200、1400、1600、1800、2000 eV完成扫描；不同网格序列没有拼接成同一条曲线。这里的密度 `pd2` 名义截断为波函数截断的2倍，`pd3` 名义截断为8倍，均受网格上限影响。因此扩展截断与加密网格是相关的数值检查，不能只改变一个名义数值就假定密度表示已经足够。

下表为三条相容序列各自候选设置与已确认稳定末点的差异：

| 扫描 | 候选 → 本序列末点 | 其余主要设置 | ΔE (meV/atom) | 最大力矢量差 (eV/Å) | 最大应力分量差 (kbar) | 结果 |
|---|---|---|---:|---:|---:|---|
| 截断能 | 1400 → 2000 eV | k=(4,8,12)，网格96×60×36 | 0.04189464 | 0.0000474977 | 0.07317446 | 通过 |
| k点 | (3,6,9) → (5,10,15) | 1200 eV，网格96×60×36 | 0.00833273 | 0.000274405 | 0.00236840 | 通过 |
| 实空间网格 | 96×60×36 → 128×80×48 | 1400 eV，k=(3,6,9) | 0.00000524916 | 0.000000176710 | 0.00000318548 | 通过 |

三条序列均至少有三个独立设置，末两段相邻差已通过，而且候选及其后所有点对本序列末点都通过。k点扫描使用显式 Gamma 中心，实际已完成至(5,10,15)；原计划(6,12,18)为可继续扩展的设置，没有把未运行点列为完成。组合候选为 **1400 eV、k=(3,6,9)、gpts=(96,60,36)、stencil=4**。由于这些扫描的固定条件不同，组合候选仍需后述联合加密与方法/几何迁移复核，不能仅由三个单因素结果宣称全部生产参数已锁定。详见[固定几何收敛进度](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/reports/固定几何收敛进度.md)。

在指定 PBE 参考点上，4与8个 MPI 进程的一致性也已通过：自由能差3.55271×10⁻¹³ eV/胞、最大力矢量差3.34459×10⁻¹³ eV/Å、最大应力分量差2.63753×10⁻¹¹ kbar，分别低于10⁻⁶ eV/胞、10⁻⁶ eV/Å及10⁻⁵ kbar门槛。这允许把进程数改变作为经检查的计算安排，但其证据范围仍限于已检方法、几何和设置。[并行一致性验收](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/validation_audit/parallel_invariance.md)

| 正式收敛指标 | 门槛 | 比较方式 |
|---|---|---|
| 能量 | ≤0.5 meV/atom | 相同几何、同一方法和 PAW 的总自由能差 |
| 原子力 | ≤0.005 eV/Å | 先逐原子计算力矢量差，再取最大范数；不是两个最大力数值相减 |
| 应力 | ≤0.1 kbar | 六个分量差的最大绝对值 |
| 两结构相对能的变化 | ≤0.5 meV/atom | 比较同一方法下 `(E_GFN2−E_Tashiro)/12` 在不同精度间的变化 |
| 晶胞长度 | <0.5% | 在接受参数与更高精度下分别优化，以高精度终态长度为分母比较对应轴 |

固定几何分析按方法、完整有序几何、PAW、程序及全部非扫描参数分组。至少需要三个不同设置，末三点的两段相邻差都通过才形成暂定平台；推荐设置及其后所有点还要对末点通过全部对应门槛。最高精度点本身不会自动成为推荐设置。晶胞标准由方案的0.005容限落实为严格小于0.005的验收，临界值不自动通过。

### 正在运行的三方法、两结构、三精度联合复核

Tashiro 与 GFN2 两套冻结结构分别在三种方法下使用以下 L/M/H设置，共计划18个单点。L是待采用的组合候选，M/H联合加密截断能、k点和网格；“L”仅是相对命名，不表示使用低质量数据。

| 精度 | 截断能 (eV) | k点 | 网格 |
|---|---:|---|---|
| L | 1400 | 3×6×9 | 96×60×36 |
| M | 1600 | 4×8×12 | 112×72×44 |
| H | 1800 | 4×8×12 | 112×72×44 |

每一方法、每一几何分别检查总自由能、原子力和应力的变化，同时严格配对两结构，检查相对能 `ΔE结构=(E_GFN2−E_Tashiro)/12`。每种方法至少需要三个相容且递增的精度，末两段及候选相对共同最高精度的变化都须≤0.5 meV/atom。相对结构能为负，只表示该方法和这两套冻结坐标下 GFN2 结构能量更低；不等于优化后晶相稳定性，也不能跨泛函比较绝对能量。当前这些联合与相对能检查仍为进行中，不能用早期的一对结果提前验收。[两结构相对能进度](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/reports/相对结构能收敛.md)

### 晶胞优化与共同高精度冷启动复核

晶胞比较采用 PBE-D3(BJ)，以 Pnam 对称约束同时优化原子和三个轴长。这是对原规范“PBE扫描并观察晶胞”的明确调整：固定几何主扫描保持PBE，变胞阶段选择含色散的既定方法，目的是把色散修正纳入PE链间堆积的晶胞验证。[GPAW官方D3说明](https://gpaw.readthedocs.io/documentation/extensions/d3.html)支持其长程相互作用修正用途。本项目没有据此断言纯PBE必定不结合；此处也不验证另外两方法的平衡晶胞。

候选精度优化首次完整评价及轨迹接口已独立核对；其原任务现已自然结束，退出码0，终态科学验收仍待独立核查。高精度参考优化于2026-10-20 12:59 UTC启动，13:22:52 UTC保存后中止，见[进程记录](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/relaxations/cell_reference/process.json)及[原中止记录](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/relaxations/cell_reference/interruption_record.json)。中止记录确认当时参考档各MPI进程均已退出，状态为 `deferred`、验收为 `not_completed_not_passed`。记录中的后置原因是当时执行方的判断；用户随后明确不接受停止工作，现已恢复参考优化。原记录作为历史保留，不表示用户授权了中止，也不表示科学验收失败或成功终态。

[中止前快照](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/relaxations/cell_reference/before_defer_20261020T132252Z)保存输入、环境与源码快照、4个完整评价（0000—0003）、轨迹和优化器重启资料。没有电子波函数GPW检查点；新尝试 `cell_reference_resume_01` 从完整评价0003对应的轨迹末帧冷启动电子态，并新建优化器。候选档自然退出释放资源后，新尝试已于13:30:21 UTC以4个MPI进程启动，见[恢复运行记录](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/relaxations/cell_reference_resume_01/process.json)。原尝试与新尝试分别保留，不将恢复启动记作优化通过。此前13:23:00 UTC恢复的是GFN2几何的DFT迁移队列协调进程，详见[调度记录](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/scheduling_events.jsonl)；这与参考晶胞优化的13:30:21 UTC恢复是两个不同事件。

晶胞程序保留 Hall297 输入轴，分别记录原始与对称投影后的力和应力，并以原始最大力≤0.002 eV/Å、原始最大应力分量≤0.05 kbar、优化器收敛及目标对称性共同验收。单次优化通过不能代替两档精度的晶胞长度比较。

参考优化已根据[独立网格桥接](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/validation_audit/pbe_reference_grid_bridge.json)选为1800 eV、k=(4,8,12)、gpts=(96,60,36)；[选择记录](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/relaxations/cell_reference/grid_selection_decision.json)及修改前g112输入均保留。这与固定结构迁移H档的g112明确区分。桥接比较同一Tashiro几何的PBE DFT分量，g112/4 MPI对g96/8 MPI得到能量、力、应力差4.25413×10⁻⁶ meV/atom、1.47354×10⁻⁷ eV/Å、3.69086×10⁻⁶ kbar，并引用独立4/8 MPI一致性证据；它不替代最终几何的冷启动复核。

候选优化的[曲线PNG](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/reports/figures/relaxation_cell_candidate.png)、[PDF版](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/reports/figures/relaxation_cell_candidate.pdf)及[图像来源指纹](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/reports/figures/relaxation_cell_candidate.sources.json)已提供。**图中包含范围以其来源指纹为准，图像本身不构成终态独立验收。** 候选任务虽已正常退出，仍不能由曲线宣称两档比较通过。图中轴长变化相对初始实验胞，不能将其幅度替代“两档精度优化终态相差<0.5%”的比较。

两档优化成功后，还须将两个终态分别在**同一更高精度设置下从头初始化电子态**进行单点复核，保留终态全精度轨迹坐标及输入/结果哈希。冷启动单点相对各自优化最后评价的自由能变化须≤0.5 meV/atom、最大力矢量变化≤0.005 eV/Å、最大应力分量变化≤0.1 kbar；同时检查原始力和应力、对称性与来源一致。缺少真实冷启动证据时，晶胞比较即使长度接近也保持待完成状态。两个成功终态将分别追加真正的 Pnam 独立位点 CIF，保留原有全晶胞 P1 文件。

**本阶段尚未锁定生产参数。** 三方法/两结构联合复核、相对结构能收敛、两档变胞比较与共同高精度冷启动复核全部完成后，才可写出带软件、PAW、色散、截断能、k点、网格和适用结构范围的锁定记录。即使数值误差满足门槛，泛函、PAW、色散模型、零点振动及实验温度差异仍需分别考虑。

## 8. 与原计划的差异及影响

| 原计划/设想 | 实际实施 | 对结果解释的影响 |
|---|---|---|
| QE 与 SSSP 校验的 PBE PAW | GPAW 与统一原生 PBE C/H PAW | 保留 PAW 路线，数据集和数值实现改变；必须重新建立收敛，不能使用 SSSP 认证表述。 |
| 直接使用预编译环境 | 官方同版本 GPAW 源码本地构建 | 补齐实际需要的库绑定；科学源码未改，固定源码与二进制均可追溯。 |
| D4 4.3.0 | D4 4.2.0 | 解决与 D3 的底层依赖冲突；按实际版本报告和冻结方法。 |
| QE 截断值或 Ry 示例直接移植 | 以 eV 重新扫描 | 单位可换算，跨 PAW/软件的名义截断值不代表同等精度。 |
| 一次 PBE 收敛后永久沿用 | PBE 主扫描加三方法、网格及独立几何复核 | 避免把一个方法/几何的数值稳定性未经检查地推广。 |
| PBE扫描并观察晶胞 | 固定几何主扫描仍用PBE，两档变胞比较改用PBE-D3(BJ) | 将色散修正纳入PE链间堆积的晶胞验证；该胞长数值验证只覆盖PBE-D3(BJ)，不自动覆盖另外两方法的平衡晶胞。 |
| 400—800 eV初始截断扫描 | 保留失败序列，扩至2000 eV；在96×60×36网格另建相容扫描 | 至1200 eV的应力仍不稳定，且密度倒空间表示受网格上限影响。计算量增加，但验收阈值不变。 |
| 分开检查截断能与k点即可 | 新增网格扫描及三方法×两结构×三精度联合加密 | 单因素平台不足以证明组合候选在每种方法和几何上均适用。 |
| 只检查各结构的绝对总能量变化 | 增加两套冻结结构相对能的三精度收敛 | 直接控制候选间能量差所受数值误差；门槛为0.5 meV/atom。 |
| 两档优化只比较晶胞长度 | 同时要求两终态在共同高精度下冷启动复核能量、力和应力 | 避免以两个尚有数值偏差或电子初态依赖的终态直接锁定参数。 |
| 仅归档 P1 全晶胞 CIF | 额外生成经还原验证的3独立位点 Pnam CIF | 改善晶体学表达；原始文件和原子几何保留，不是另一套结构。 |
| “通过后只剩泛函误差” | 区分数值、PAW、色散、泛函及温度等误差 | 参数收敛不等于实验准确性，也不证明某一方法普遍最好。 |

三种选定方法的身份保持不变。CN 半径、截断平滑、stencil 和探针截断位置的调整均已说明原因和数据影响。详见[完整差异报告](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/reports/与原计划的差异及影响.md)。

## 9. 能量与原始数据入口

以下三个文件是最终验证批次的探针中心结果，便于直接查看已通过的总能量、全部力、六分量应力、DFT/色散分项以及输入和 PAW 来源。它们是导数探针单点，不能当作优化终态或生产排名。

| 方法 | 完整结果 | 逐次电子自洽能量 |
|---|---|---|
| PBE-D3(BJ) | [result.json](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/runs/deriv_final_pbe_d3bj_center_E600.012/result.json) | [scf_iterations.jsonl](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/runs/deriv_final_pbe_d3bj_center_E600.012/scf_iterations.jsonl) |
| optB88-vdW | [result.json](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/runs/deriv_final_optb88_vdw_center_E600.012/result.json) | [scf_iterations.jsonl](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/runs/deriv_final_optb88_vdw_center_E600.012/scf_iterations.jsonl) |
| r²SCAN-D4 | [result.json](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/runs/deriv_final_r2scan_d4_center_E600.012/result.json) | [scf_iterations.jsonl](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/runs/deriv_final_r2scan_d4_center_E600.012/scf_iterations.jsonl) |

读取最终能量时使用 `total.energy_free_eV`；其单位为每晶胞 eV。`total.forces_eV_A` 为全部原子力，`total.stress_kbar` 为六分量应力。`dft` 与 `dispersion` 保存对应分项。自洽迭代文件记录电子迭代过程，不能把其中未收敛的中间能量当最终结果。

正式收敛的能量入口为[完整固定几何 JSON 表](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/tables/convergence_summary.json)及[CSV 表](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/tables/convergence_summary.csv)，三精度两结构配对能量另见[相对结构能表](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/tables/relative_structure_energies.json)。候选晶胞优化的全部电子自洽与结构步记录位于[候选优化目录](E:/PE_DFT/calculations/pe_gpaw_convergence_20261020/relaxations/cell_candidate)；此目录存在不表示优化已成功，须以最终 `result.json` 及独立验收状态判定。

每个运行目录还保留输入结构、参数 JSON、程序与方案快照、GPAW 日志、标准输出/错误以及各 MPI rank 的基组记录。分析器选择成功计算点并记录原始哈希；失败尝试不会被删除，也不会以失败时最后一帧冒充成功结果。正式扫描的进展以固定几何报告及各运行目录为准。

### 本版报告依据的关键文件指纹

| 文件 | SHA-256 |
|---|---|
| `protocol.json`，第5版 | `7724a1fe1fa9ad5979cf78dc5629096c5297cfba9e84204860f80d259a92ff00` |
| `derivative_results_deriv_final.json` | `46f5e8792862737fb4dec779d8b9be6b0168595005f7614ff1fcef702ef64877` |
| `scf_noise_validation.json` | `b97f1e47ef859ecdac86ce245528cc35d67c3cd00126e0fbfb4359fd73567f32` |
| `parallel_invariance.json` | `33721b739857c4e97aac8a7249a4137a24ab1282efed32471a1f4b198814dfad` |
| `convergence_summary.json`，本次读取的阶段快照 | `871e54d627afacc2bc8b4cea70676f01c6143ef53ac03dab38681368dbc6cfed` |
| `reference/Tashiro_Pnam.cif` | `ce855a08a7cc9e2b9327e78635f7320a68d96651b1ee4354639a61d58f371f9c` |
| 第4版阶段报告备份 | `f6068cb1bf08d5a49dd485f0403c1b3a1dafba735edff90efa18f8e8e9725c74` |

其中收敛汇总和相对能表仍会随运行继续更新；表中指纹是本次读取时的版本，不要求未来进度表与旧指纹相同。后续更新晶胞与迁移结果时，应同时更新本报告的阶段状态与依据指纹；本版给出了已通过的 PBE 固定几何候选，尚未给出最终生产参数或完整收敛通过结论。
