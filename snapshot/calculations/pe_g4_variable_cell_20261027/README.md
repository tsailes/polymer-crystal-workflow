# 最终 g4 的 PBE-D3(BJ) 变胞优化

2026-10-27 用户要求立即执行，实际于北京时间 11:14:33 后台启动。

- 起点为 `pe_dft_tight_force_refinement_20261022/runs/group_004/attempt_01/final_structure.extxyz` 的逐字副本。该终态已经独立电子复核通过；原始文件保持不变。
- 不从 Tashiro 原子坐标开始。初始胞是 g4 原有固定实验胞 7.432 / 4.929 / 2.545 Å。
- GPAW 26.7.0 / ASE 3.29.0；PBE-D3(BJ)+ATM；1800 eV、Gamma 中心 4×8×12、基础网格 112×72×44、16 bands、0 eV 电子展宽。
- SCF 能量、密度和本征态残差门槛分别为 1e-9 eV/electron、1e-9、1e-11 eV²/electron。
- 复用已经验收的 FrechetCellFilter + FixSymmetry 实现。开放内部原子和 Pnam 允许的三个独立轴向应变；90°角保持不变，无恒体积或等比例缩放限制；外压目标为零。
- LBFGS maxstep=0.05，单次上限200步；原始最大原子力≤0.002 eV/Å，六应力分量最大绝对值≤0.05 kbar，并要求 Pnam 闭合和优化器收敛。
- 优化达标后由同一协调器自动建立 `cold_verification/`，固定新终态，重新初始化电子并重新检查力和应力。没有复用最终电子密度或波函数。
- 后续分别加密 k 点、提高截断能及多方法单点不由本协调器启动，另由主任务安排。

## 运行记录

- `protocol.json`、`source/provenance.json`：起点和科学设置。
- `control/input_check.json`：仅输入/接口预检，无量子计算。
- `control/windows_launch_receipt.json`：Windows 隐藏启动回执。
- `control/workflow_status.json`：流程当前状态。
- `runs/attempt_01/process.json`：MPI 命令和进程号。
- `runs/attempt_01/gpaw.txt`：电子计算日志。
- `runs/attempt_01/scf_iterations.jsonl`：逐轮 SCF 能量。
- `runs/attempt_01/evaluations.jsonl`：逐几何原始能量、力、应力、晶胞和坐标。
- `runs/attempt_01/trajectory.traj`、`optimizer.log`：轨迹与优化器记录。
- `runs/attempt_01/result.json`：优化终态验收，只有完成后才生成。
- `final_cold_verification.json`：独立终态验收，只有对应步骤完成后才生成。

协调器使用独占文件锁，检测到已有结果时拒绝重复覆盖。若需续算，应建立新 attempt 并保存来源和参数一致性，不放宽验收门槛。

## 首次完整评估与晶胞更新（2026-10-27 11:20 左右）

`runs/attempt_01/evaluation_0000.json` 已完成23轮SCF，耗时347.856秒：

- 起始总能量：−66.77581676509409 eV/cell。
- 原始最大原子力：0.0001717965190 eV/Å。
- 六应力分量：约 +3.85170045、+1.96008656、−12.97435856、0、0、0 kbar。
- 起始 a/b/c：7.432 / 4.929 / 2.545 Å；Pnam闭合残差2.74e−15 Å。

随后优化器已经提交第1个晶胞更新，`evaluation_input_0001.json` 中 a/b/c 为 **7.4082421574 / 4.9209753855 / 2.5725969857 Å**，三个夹角仍90°，正在该新几何上求电子。该记录是正在计算的试探几何，不是已收敛终态。

`optimizer.log` 的 fmax 包含晶胞广义力，不能当作原子最大力；原子力与应力分别以 `evaluation_*.json` 的原始量判定。Windows启动器的132字节stderr仅为WSL localhost代理未镜像提示；MPI计算自身stderr在首次评估时为空。
