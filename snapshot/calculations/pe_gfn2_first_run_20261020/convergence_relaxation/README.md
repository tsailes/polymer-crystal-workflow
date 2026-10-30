# PE 周期采样与固定晶胞优化

已完成运行：`runs/20261020T063341_922507Z`，最终判定 `accepted`。

阅读[中文报告](runs/20261020T063341_922507Z/result_report_zh.md)，其中包括计算解释、参数选择、原子力验收、流程前三个级别的含义和未完成的研究环节。

本阶段从已有独立构建的一个12原子PE周期初猜出发，以GFN2-xTB/tblite 0.7.0完成数值采样检查和固定胞局部优化。选定2×3×8 Γ点重复超胞，优化时原胞12个原子位置自由变化、平移周期和晶胞不变。26步后最大原子力0.009876 eV/Å。终态3×4×10组合加密的最大原子力0.007675 eV/Å，能量变化0.302357 meV/CH₂，最大力向量变化0.003617 eV/Å，均通过预设门槛。

终态两条链的锯齿平面接近平行。它只是一个优化候选，尚未验证为实验正交PE堆积。当前晶胞尺寸采用实验先验；尚未做多初猜搜索、DFT计算或3σ标定。

## 交付文件

- [最终CIF](runs/20261020T063341_922507Z/relaxation_2x3x8/pe_relaxed.cif)
- [最终extxyz](runs/20261020T063341_922507Z/relaxation_2x3x8/pe_relaxed.extxyz)
- [优化曲线](runs/20261020T063341_922507Z/plots/fixed_cell_optimization.png)
- [初末结构比较](runs/20261020T063341_922507Z/plots/initial_relaxed_xy.png)
- `summary.json`：所有数值验收条件；`initial_decision.json`：初态选参依据。
- `relaxation_2x3x8/optimization.traj`：完整27帧轨迹（第0帧加26步）。
- 各单点目录中的 `input.extxyz`、`result.json`、`native_output.log`：原始输入、结果和原生日志。
- `preview_26steps/`：复核尚未全部完成时保存的预览；最终图在 `plots/`。

## 复现

在PowerShell运行本目录的 `run_convergence_relaxation.ps1`。
它使用已经安装的WSL Ubuntu-24.04环境 `/opt/pe_dft/envs/gfn2-pe-20261020`，
依次执行初态采样、受力差分验证、参数判定、固定胞优化和终态复核。
每次生成新时间戳目录；任一验收失败都会停止，不覆盖之前的原始结果。
完整重跑需要多次较大超胞计算。

`reproduce.py`记录全流程依赖顺序；`workflow.py`负责单点和优化；
`assess_results.py`只比较有对应关系的细化计算；`plot_results.py`只绘制已有数据。
这些脚本的语法与本轮实际调用路径已检查；没有为验证启动第二次完整重跑。

方法仍为GFN2-xTB。独立 `xtb` 程序和 `tblite` 库是不同软件入口，本轮所有合格数据均来自指定的tblite版本。上一级目录中隔离的旧DFTB+后端未参与本轮优化。
