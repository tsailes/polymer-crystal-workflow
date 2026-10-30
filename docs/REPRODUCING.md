# 复现说明 / Reproduction guide

本仓库提供原始计算源码和只读证据快照，以及两个可迁移的入口：几何候选再生成、冻结几何 DFT 单点输入准备。轻量检查可以独立执行；完整的 GFN2→DFT 搜索和优化仍需安装量化计算环境并迁移历史调度路径。本次源码整理不代表在另一台机器上完成了全流程量化复算。

以下命令从仓库根目录运行。Python 3.12 是历史运行环境。`snapshot/` 保留原始字节及日期目录；新计算输出写入 `work/` 等新目录，禁止直接在证据快照内重跑。

## 1. 校验快照

只需 Python 标准库：

```bash
python scripts/verify_snapshot.py
```

该检查核对已归档文件和记录的一致性，不启动量化计算。它不能代替新的数值收敛检查。

## 2. 重新生成有限候选网格

创建独立 Python 环境，并安装轻量依赖：

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-analysis.txt
python scripts/reproduce_candidates.py --output work/candidates
```

Windows PowerShell 的环境激活命令为 `.\.venv\Scripts\Activate.ps1`。依赖版本来自原始 `protocol.json` 和 `environment_verified.json`；版本锁定不表示已验证所有操作系统组合。

辅助程序先把原始 `generate_candidates.py` 原样复制到空目录，再用当前 Python 执行，不修改生成算法或科学参数。输出目录已包含文件时拒绝运行；换一个新目录即可重跑。

预期为：288 枚举、248 几何过滤通过、40 拒绝、232 等价重复、16 独立初态。`reproduction_verification.json` 检查这些计数、全部候选的编号/参数/分类和代表对应关系，并记录实际软件版本与源文件哈希。它不要求含软件环境信息的整个新 manifest 与历史 manifest 字节一致。

这个步骤只有化学构造、几何过滤和结构等价性计算；不运行 GFN2、DFT 或实验拟合，也不证明连续构型空间全局穷尽。

## 3. 准备一项冻结几何 DFT 单点复算

```bash
python scripts/prepare_dft_input.py --output work/g4-cold-point
```

此命令使用已验收的**变胞 g4 终态**及其独立冷启动单点输入。它核对结构和原始 runner 的 SHA256，仅把 `structure` 和 `point_runner` 两个路径改成本地文件名；截断能、k 网格、实空间网格、SCF 阈值及其他参数完全保留。此步骤本身不启动 GPAW。

生成文件：

| 文件 | 用途 |
|---|---|
| `input.json` | 迁移路径后的实际单点输入 |
| `structure.extxyz` | 原始变胞终态结构，字节不变 |
| `run_point.py` | 原始 GPAW 单点程序，字节不变 |
| `expected_result.json` | 历史已接受的结果，仅用于对照 |
| `preparation_manifest.json` | 源哈希、迁移字段与准备记录 |

### Linux / WSL 量化环境

历史生产环境：GPAW 26.7.0 legacy CPU plane-wave、gpaw-data 1.0.1 的原生 PBE C/H PAW、ASE 3.29.0、NumPy 2.5.4、SciPy 1.18.1、dftd3 1.6.0、OpenMPI 5.0.11、Libxc 7.1.2。GFN2 阶段另用 tblite 0.7.0；它不属于这个 DFT 单点的依赖。

构建与环境证据见：

- [`snapshot/calculations/pe_gpaw_convergence_20261020/environment/`](../snapshot/calculations/pe_gpaw_convergence_20261020/environment/)
- [`snapshot/calculations/pe_g4_variable_cell_20261027/environment/`](../snapshot/calculations/pe_g4_variable_cell_20261027/environment/)

这些历史安装记录保留 `/opt/pe_dft/...` 路径。它们用于说明已运行过的环境，不是适用于所有机器的直接安装脚本。`pip_freeze_final.txt` 中本地构建地址也不能直接当作通用 pip 安装清单。先准备兼容的 MPI/GPAW/Libxc/dftd3 环境，将该环境的 `python` 和 `mpiexec` 放到 PATH，并将 `GPAW_SETUP_PATH` 指向同版本 PAW 数据目录；确认 C/H 数据哈希与历史结果相同。

在已经准备好上述环境的 Linux shell 中，具体运行命令为：

```bash
export GPAW_MPI_BACKEND=cgpaw
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
mpiexec -n 4 python work/g4-cold-point/run_point.py work/g4-cold-point/input.json
```

这是实际量化计算，会使用计算资源。输入为 H 生产设置：1800 eV、Γ 中心 4×8×12 k 网格、112×72×44 实空间网格、16 bands、零占据展宽、stencil 4、PBE-D3(BJ)+ATM。`GPAW_SETUP_PATH` 和共享库依赖须先配置，不能用不同赝势替代后仍声称是同一复算。

新结果位于 `work/g4-cold-point/result.json`。检查 `status=passed_single_point`、`SCF_converged=true`、软件版本、PAW 数据哈希、结构哈希，再比较能量、逐原子力和应力；跨机器运行不能仅凭进程退出码认定科学复现。历史独立单点的最大原子力为 0.0006209922404323713 eV/Å，最大应力分量为 0.011392943840518 kbar；这两项是参考输出，不是新复算结果。正式验收需遵守技术路线中冻结的误差指标。

## 4. 全流程源码入口及迁移边界

| 阶段 | `snapshot/` 下原始目录 / 入口 | 复现需要注意 |
|---|---|---|
| 几何候选 | `calculations/pe_pnam_candidates_20261020/generate_candidates.py` | 使用上面的安全入口，避免覆盖快照 |
| GFN2 预优化 | 同目录 `relax_candidates.py` | 动态导入相邻 `pe_gfn2_first_run_20261020/convergence_relaxation/workflow.py`；保留相对目录，安装 tblite |
| GPAW 数值验证 | `calculations/pe_gpaw_convergence_20261020/scripts/` | 先完成方法/导数/截断能/网格等验收；调度器有历史绝对路径 |
| 同几何标定 | `calculations/pe_gfn2_dft_calibration_20261021/scripts/` | 冻结几何和配对身份必须一致；全组统计不能用不同几何混配 |
| 固定胞精修 | `calculations/pe_dft_representative_refinement_20261022/scripts/` | 结合 tight-force 目录的输入和终态验收，区分各 attempt |
| 严格固定胞 | `calculations/pe_dft_tight_force_refinement_20261022/` | 以接受的终态与独立验收记录为准 |
| 变胞优化 | `calculations/pe_g4_variable_cell_20261027/scripts/run_relax.py` | 从接受的固定胞 g4 启动；保持 Pnam、仅放开 a/b/c；非本页单点复算 |
| 独立验证 | 同目录 `run_terminal_checks.py`、`cold_verification/`、`convergence_checks/` | 固定终态分别做独立电子初值和更严格数值设置检查 |

历史调度器、队列和部分报告程序包含 `/mnt/e/PE_DFT`、`E:\PE_DFT`、`/opt/pe_dft` 等路径。端到端迁移时必须在新的工作副本中逐一配置路径、环境可执行文件和所需上游文件；移植后对新的输入与代码重新记录哈希。不要静默覆盖原始 provenance，也不要把旧队列当作已经适配到本机的一键工作流。

本仓库未包含大型波函数/重启二进制、软件安装目录、全文文献或所有历史探索分支。已归档的结果与轻量复现检查支持审阅和后续复算准备；新机器完整 GFN2/DFT 重跑、正式论文版本和长期归档标识仍需分别完成。
