# 本轮 PE 量子方法与结果审计

审计日期 2026-10-28。原始计算目录只读；本笔记以实际 JSON/CSV 的最新记录为依据。路径前缀均为 `E:/PE_DFT/`。本笔记未启动新的量子计算。

## 1 可用于技术文档的当前结论

- 2026-10-21 的三方法数值验证已完成并发布参数锁。`pe_gpaw_convergence_20261020/reports/` 内若干旧 Markdown 仍写“进行中/未锁定”，属于阶段快照，不能当成当前状态。当前依据为 `calculations/pe_gpaw_convergence_20261020/completion_state.json` 与 `releases/20261021T043247_514339Z/parameters_lock.json`。
- 四个 GFN2 代表在固定实验晶胞内，以 0.002 eV/Å 原始力门槛精修后形成 2 个容差组；继续到 0.0002 eV/Å 后 4/4 通过且六对全部满足几何与能量条件，合并为 1 组。g4 是记录的代表，但微小能差不可解读为可靠的物理能量排序。
- 2026-10-27 从上述最终固定胞 g4 起步的变胞优化已经完成。独立冷启动，以及单独 k 点、实空间网格、截断能加严复核均通过。它和 10-20/21 从实验参考坐标起步的两档精度变胞验证是不同计算。
- PBE0 单点是用户取消，选定 1 个代表，接受 0 个，无有效最终能量、无 SCF 收敛、无 EXX k 点收敛认证。不可列入已完成方法对比，亦不能把取消解释为科学失败。

## 2 软件与方法

本轮主线运行于 WSL Ubuntu 24.04，GPAW 26.7.0 legacy CPU 平面波实现，ASE 3.29.0，spglib 2.7.0，gpaw-data 1.0.1，Libxc 7.1.2，libvdwxc 0.5.0，simple-dftd3 / dftd3-python 1.6.0，dftd4 / dftd4-python 4.2.0，Python 3.12.15，NumPy 2.5.4，SciPy 1.18.1，OpenMPI 5.0.11。环境 `/opt/pe_dft/envs/pe-gpaw-20261020`；GPAW cgpaw MPI 后端。本轮正式 GPAW 数据不可写作 GPU 结果。

三方法统一用 GPAW 原生 PBE C/H PAW，并非 SSSP 选择的 UPF。C 价电子/冻结芯为 4/2，H 为 1/0。压缩文件 SHA256：C=`5a769dd14bb5d462c57fadc4dd8d3b7e70428b90701eca0d20b6f0666495cdfd`；H=`43da6533cea08dec1ee7861140be503c9945863d6d8c7791e01e55e2161a88b8`。

| 方法 | 实际实现 | 色散参数 |
|---|---|---|
| PBE-D3(BJ) | GPAW 自洽 PBE + D3(BJ) | ATM 开；二体/三体/CN=60/40/60 Bohr；平滑宽度0.05/0.05 Bohr；s6=1,s9=1,alp=14,a1=0.4289,s8=0.7875,a2=4.4407 |
| optB88-vdW | optB88-vdW，libvdwxc 非局域相关 | 非局域项已计入 DFT 能量，不额外叠加 D3/D4 |
| r²SCAN-D4 | MGGA_X_R2SCAN+MGGA_C_R2SCAN + D4 | ATM 开；二体/三体/CN=60/40/30 Bohr；平滑宽度0.05/0.05 Bohr；s6=1,s9=1,alp=16,s8=0.6018749,a1=0.51559235,a2=5.77342911 |

统一：C4H8，12 原子，4 CH2，三轴周期，Gamma 中心 k 网格，16 bands，电子展宽0，XC stencil=4，点群开启，4 MPI ranks。SCF 能量/密度/本征态门槛分别 1e-9 eV/electron、1e-9 electron/electron、1e-11 eV²/electron。自由能字段为电子自由能，加外部色散；不含声子、零点能及有限温度核自由能。ASE 应力拉伸为正，Voigt 次序 xx yy zz yz xz xy，1 eV/Å³约1602.1766 kbar。

来源：`calculations/pe_gpaw_convergence_20261020/releases/20261021T043247_514339Z/parameters_lock.json` 的 physics/environment/fixed_parameters；`method_audit/method_specification.json`；`environment/environment_verified.json`。初始 `calculations/pe_dft_convergence_20261020/protocol_draft.json` 中 QE/Ry 扫描属于未采用草案，不能拼入 GPAW/eV 的正式方法。

## 3 数值验证与参数层级

基本差分收敛门槛：ΔE≤0.5 meV/atom，最大原子力矢量差≤0.005 eV/Å，六应力分量最大差≤0.1 kbar。优化收敛与两档胞长比较另验收：原始 Fmax≤0.002 eV/Å，max|stress|≤0.05 kbar，Pnam 闭合、优化器收敛，以及两档每轴相对胞长差严格<0.5%。

| 用途 | Ecut eV | kmesh | 基础网格 |
|---|---:|---|---|
| 发布候选精度 L | 1400 | 3×6×9 | 96×60×36 |
| 方法/几何联合检查 M | 1600 | 4×8×12 | 112×72×44 |
| 方法/几何联合检查 H、后续固定胞生产精修 | 1800 | 4×8×12 | 112×72×44 |
| 数值校准的高精度参考变胞 | 1800 | 4×8×12 | 96×60×36 |
| 10-27 g4 变胞 | 1800 | 4×8×12 | 112×72×44 |

PBE 在固定 Tashiro 几何的正式截断扫描1000/1200/1400/1600/1800/2000 eV（k4×8×12、g96×60×36），候选1400相对2000差为0.04189464 meV/atom、4.74977e-5 eV/Å、0.07317446 kbar。k 扫描2×4×6、3×6×9、4×8×12、5×10×15（1200 eV、g96），候选k3相对k5差为0.00833273 meV/atom、0.000274405 eV/Å、0.00236840 kbar。网格扫描g96/g112/g128（1400 eV、k3），g96相对g128差为5.24916e-6 meV/atom、1.76710e-7 eV/Å、3.18548e-6 kbar。

上述数值为对已检最高档的差，不是相对精确解的真误差。早期较粗网格400—1200 eV的失败序列仍保留，不与正式序列混接。最终释放门控含8项：环境方法和代码、导数与SCF噪声、PBE截断/k/grid、相对结构能、三方法×两几何×L/M/H、MPI与网格桥接、两档变胞与冷启动、证据文件一致性，均通过。

导数验证共111成功点，每方法37点（1中心+12原子位移端点+24应变端点）；27组=每方法3个代表性力分量+6应力分量，非全部36原子力分量。另6个更严SCF端点/3组复核。以下为总量最坏误差：

| 方法 | 力差分误差 eV/Å | 力两步长差 eV/Å | 应力差分误差 kbar | 应力两步长差 kbar |
|---|---:|---:|---:|---:|
| 门槛 | 0.001 | 0.0005 | 0.1 | 0.05 |
| PBE-D3(BJ) | 0.000128498 | 0.0000986308 | 0.00178094 | 0.000774909 |
| optB88-vdW | 0.000121077 | 0.0000969985 | 0.00410831 | 0.0000697203 |
| r²SCAN-D4 | 0.000129788 | 0.000101631 | 0.000801533 | 0.00000290628 |

注意600.012 eV是避免该导数探针 G 向量进出截断球面的特殊参数，不是生产推荐。D3 CN40→60和XC stencil2→4是有记录的数值修正，不能改写旧结果或声称换了泛函。

三方法×两冻结几何共18个L/M/H点均通过。两冻结几何的相对能定义为(E_GFN2几何−E_Tashiro几何)/12：

| 方法 | L meV/atom | M meV/atom | H meV/atom |
|---|---:|---:|---:|
| PBE-D3(BJ) | -0.7384605469 | -0.7379636240 | -0.7379654694 |
| optB88-vdW | -0.3927948944 | -0.3923833000 | -0.3924992487 |
| r²SCAN-D4 | -1.4438002156 | -1.4441415398 | -1.4443154661 |

这两几何是数值迁移验证样本，不能代表最终四代表精修的跨泛函排序。来源：`tables/convergence_summary.csv`、`tables/relative_structure_energies.csv`、`validation_audit/derivative_results_deriv_final.csv`、`validation_audit/scf_noise_validation.json` 及发布目录完整 evidence。

## 4 从实验参考起点的两档变胞数值验收

| 运行 | 胞长a/b/c Å | 总能 eV/cell | Fmax eV/Å | max应力 kbar | 步数 |
|---|---|---:|---:|---:|---:|
| L cell_candidate | 7.1267168374 / 4.9036884078 / 2.5555648720 | -66.78206623909374 | 0.000195438714 | 0.0287043623 | 33 |
| 参考 cell_reference_resume_02 | 7.1300676744 / 4.9065189830 / 2.5556136076 | -66.78263584715239 | 0.000055841864 | 0.0041785586 | 22 |

相对胞长差 a/b/c=0.0469958661% / 0.0576900903% / 0.0019070018%，最大0.0576901%<0.5%。两个终态在共同高精度下冷启动核验通过。L终态相对自身优化末点的ΔE/ΔF/Δσ=0.0472221254 meV/atom、0.00027207710 eV/Å、0.0745306945 kbar，均在差分阈值内。原先中止及恢复尝试应完整保留，当前采用 cell_reference_resume_02，不能选旧未收敛末帧。

来源：`calculations/pe_gpaw_convergence_20261020/relaxations/{cell_candidate,cell_reference_resume_02}/result.json`、`relaxations/cross_precision_verification.json`、发布 `release_check.json` 的 details.two_relaxations_and_cold_start。

## 5 固定胞四代表精修与加严

固定 a/b/c=7.432/4.929/2.545 Å，Pnam #62 Hall297，FixSymmetry，LBFGS maxstep=0.05 Å/maxsteps=200；H参数；不以应力作为固定胞收敛目标；每终态独立冷启动。

| 代表 来源 | 常规终能 eV/cell | 常规Fmax eV/Å | 加严终能 eV/cell | 加严Fmax eV/Å | 加严setting angle ° |
|---|---:|---:|---:|---:|---:|
| g1 cand_0025 | -66.77560670952218 | 0.00175402638 | -66.77581676505100 | 0.000126121333 | 45.749313010 |
| g2 cand_0113 | -66.77581136546064 | 0.00033043603 | -66.77581676453569 | 0.000160864611 | 45.750785664 |
| g3 cand_0041 | -66.77581666973727 | 0.00073623306 | -66.77581672813125 | 0.000046204010 | 45.734634614 |
| g4 cand_0081 | -66.77581673830392 | 0.00097528203 | -66.77581676510076 | 0.000171616357 | 45.748669392 |

加严唯一科学参数变化是 Fmax 0.002→0.0002 eV/Å，从前一接受终态重建电子态与优化器。常规四代表分成{g4,g3,g2}与{g1}两组；加严后全并一组，完整六对最大RMS=0.000324746926 Å，最大单原子位移=0.000394452672 Å，最大能差=9.24237753e-6 meV/CH2。去重阈值为同一映射同时 RMS≤0.02 Å、max≤0.05 Å、ΔE≤1 meV/CH2；使用周期边界、平移、正交反射、同元素置换，组内所有两两关系通过。未计算声子/Hessian，不证明独立多晶型或全局最低能。

加严后绝对能量与色散分量（eV/cell）：g1 DFT=-66.14952269364146,D3=-0.6262940714095436；g2=-66.14952283116980,-0.6262939333658826；g3=-66.14952208656348,-0.6262946415677678；g4=-66.14952284549459,-0.6262939196061678。

来源：`calculations/pe_dft_representative_refinement_20261022/tables/final_energy_results.csv`、`calculations/pe_dft_tight_force_refinement_20261022/tables/final_energy_results.csv` 与 `pairwise_geometry_energy.csv`。g1/g2最终采用attempt_02，g3/g4为attempt_01。全部评价轨迹、全部SCF步能、每对原子映射分别见两批次tables下 `all_evaluation_energy_trajectory.csv`、`all_scf_iteration_energies.csv`、`pairwise_atom_mappings.csv`。

## 6 最新 g4 变胞结果与独立数值加严

起点为加严固定胞 g4 的 `runs/group_004/attempt_01/final_structure.extxyz`（SHA256=64594fc1266a1ca37fe0c711bf729d28bf5f03949dd81444ed4dd17517ad874d）。开放内部原子和三个独立轴向应变，角度90°不变，外压0；FrechetCellFilter+FixSymmetry+LBFGS，cell mask=[1,1,1,0,0,0]，无等比例缩放或恒体积限制。Fmax≤0.002、max六应力≤0.05门槛。

20个优化步（evaluation_0000—0020）完成，耗时7429.1753 s；最后评价23轮SCF。终态a=7.129505992653368,b=4.908355808680994,c=2.555579019463927 Å，V=89.43032104374859 Å³。E_total=-66.78262060562021，E_DFT=-66.12706181019914，E_D3=-0.6555587954210721 eV/cell。优化末点Fmax=0.000620977796537，max|σ|=0.011392258693 kbar；应力xx/yy/zz=-0.000472430686/+0.011392258693/-0.004767385185 kbar，剪切约1e-11 kbar；Pnam残差8.21e-15 Å。

独立冷电子终态验证passed=true：Fmax=0.000620992240432 eV/Å，max|σ|=0.011392943841 kbar；E_cold−E_relax=-5.25801624e-13 eV/cell；23轮SCF/342.2528 s，未复用密度/波函数。终态结构SHA256=1f440539de6a2fca7910a82110edf83b5f7c8782b7193ac2904eed1934a921ab。

| 加严项目 | 比较基准 | 变化 | ΔE meV/atom | 最大ΔF eV/Å | 最大Δσ kbar | 自身Fmax eV/Å | 自身max应力 kbar |
|---|---|---|---:|---:|---:|---:|---:|
| k_only | 冷启动基线 | k4×8×12→5×10×15 | -0.000220649011 | 1.01447094e-5 | 8.81144566e-5 | 0.000620296133 | 0.01130482938 |
| grid_only | 冷启动基线 | g112×72×44→128×80×48 | 7.78229037e-7 | 5.01790793e-9 | 6.50177226e-7 | 0.000620992015 | 0.01139359402 |
| ecut_only | grid_only | E1800→2000，共同g128 | -0.001482668273 | 4.95043661e-6 | 0.001013139699 | 0.000619885107 | 0.01038045432 |

三点均23轮SCF、正常退出，实际基组集合/哈希/同几何/单因素/差值重算/绝对力应力/差分门槛均独立通过。ecut_only必须与grid_only比，不能直接把同时变化网格和Ecut的基线差称作单因素截断效应。其E分别为-66.78262325340887、-66.78262059628199、-66.78263838830127 eV/cell。

这些是相同终态几何的数值加严单点，未在加严设置重新优化，也未得到无限基组极限。k_only的长wall time包含已核实的Windows休眠，不适宜当作纯计算性能数据。未完成此新变胞g4的三泛函优化或声子、ZPE、温度自由能验证。

来源：`calculations/pe_g4_variable_cell_20261027/protocol.json`、`runs/attempt_01/result.json`、`final_cold_verification.json`、`cold_verification/result.json`、`convergence_checks/protocol.json`、`convergence_checks/independent_review.json`。结构文件位于 `runs/attempt_01/final.cif`、`final_POSCAR`、`final_structure.extxyz`；归档波函数为 `final.gpw`（约803 MB）。

## 7 PBE0 当前状态

计划为固定胞最终去重后每组1代表，PBE0 25%EXX + 方法专属D3(BJ)+ATM，1800 eV/k4×8×12/g112×72×44，16 bands，0电子展宽，SCF门槛同主线；GPAW legacy PW HybridXC，4 MPI，parallel={kpt:1,domain:4,band:1}；只单点，不计算力/应力。

`calculations/pe_pbe0_singlepoint_review_20261023/tables/singlepoint_summary.json` 当前 status=cancelled_by_user、complete=false、selected_count=1、accepted_count=0、rows=[]、SCF_converged=false、valid_final_energy_available=false、EXX_kpoint_convergence_certified=false、termination_is_scientific_failure=false。依据为用户取消记录 `control/cancellation_20261023T112717Z/cancellation_record.json`。技术文档应写“执行后由用户取消，保留日志，无可发表最终数值”；不能填入迭代能量充当终态。即使将来此单点完成，只有1组也不能独立检验候选排序。

## 8 正文及 SI 建议数据分配

正文：四代表收敛到一个容差组；最终固定胞结构及对应实验比较；最新g4变胞结果作为独立的静态晶胞预测结果；关键数值收敛摘要。正文应区分固定实验胞重现原子位置与0压静态变胞预测，不能把固定胞的胞长相符算作预测成功。

SI：全部生成/优化/校准项目主表；完整软件与PAW哈希；三泛函定义；扫描曲线和全部点；导数与SCF噪声；两档变胞及恢复谱系；四代表常规/加严前后完整E、F、应力、角度；六对映射与容差；最新g4坐标/CIF与三种单因素加严；无结果/取消项目日志清单。原始运行日志、波函数、逐SCF能量作为电子附件/数据仓库，正文/SI用索引追溯。

技术边界：已验证的是指定C4H8几何/胞和参数范围。没有声子/Hessian、零点能/热膨胀、所有初态全局搜索保证，也未用最新g4做完整跨泛函几何优化；这些应作为后续研究项目，不能编写成已完成结果。
