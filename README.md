# Polymer crystal reconstruction and refinement workflow

聚合物晶体结构重建、分层优化与独立数值验收的论文配套源码准备仓库。

项目面向聚合物晶体，PE（聚乙烯）是当前已完成验证的典型实例，PA6 已开展单候选 GFN2 探索性先导与问题诊断，尚未形成验收通过的量化主线。具体材料的化学拓扑、对称性、搜索范围和收敛参数需分别确定。

**当前版本：v0.1 准备快照，最新进度同步 2026-10-31。** 本仓库公开保存已完成 PE 主线及 PA6 单候选先导的计算脚本、参数、结构与精选证据，用于后续论文写作和正式代码归档。PE 的 2026-10-30 冻结快照保持不变，PA6 最新证据单独保存于 `snapshot_pa6/`。它尚不是完成独立机器全流程重跑验证的正式论文 release；当前验证范围见[复现说明](docs/REPRODUCING.md)。

English summary: This repository develops a workflow for prior-constrained polymer crystal reconstruction and hierarchical refinement, with polyethylene (PE) as the currently validated case study and PA6 as an exploratory single-candidate pilot. In the PE case, four GFN2 representatives converge to one tolerance group after tight fixed-cell DFT refinement; its representative is subsequently relaxed at zero external pressure under orthorhombic and Pnam constraints. The PA6 pilot has not produced the diagnosed target hydrogen-bond network, and periodic-cell sensitivity remains unresolved. This v0.1 preparation snapshot is not a global crystal-prediction benchmark or a fully independently rerun computational release.

**PA6 当前问题：** 一个 152 原子候选在基础胞收紧优化后原始 Fmax 为 0.004295 eV/Å，目标氢键诊断仍为 0；同几何 b 加倍后 Fmax 为 0.032857 eV/Å，未通过本次 0.005 门槛。初始注册、约束路径和周期处理的影响尚待分离，不能据此认定 GFN2 方法失效。完整证据、限制和待验证问题见 [2026-10-31 PA6 进度](docs/PA6_PROGRESS_20261031.md)。

## 计算流程与当前验证实例

当前 PE 验证实例为含两条沿 c 轴全反式链的 C4H8 周期胞，共 12 原子、4 个 CH2。生成输入采用实验晶胞及 Pnam 对称性先验；候选坐标由化学几何构造，不读取文献精修坐标。

```text
输入与先验冻结
  → 对称性约束候选生成
  → 周期接触检查与等价去重
  → GFN2 预优化与代表确认
  → 同几何 DFT 标定
  → 四代表固定胞 DFT 精修与严格力验收
  → 全体两两比较及终态合并
  → 代表零压正交变胞
  → 独立冷启动、数值加严、实验对照
```

DFT 方法和数值验证是正式标定及精修的前置条件。完整定义见[方法说明](docs/METHODS.md)及原始[PE 技术路线 v1.0](snapshot/reports/PE技术路线_v1_0_20261029/技术路线.md)。

## 当前 PE 实例的已完成结果

| 项目 | 已验收结果 |
| --- | --- |
| 候选数量链 | 288 枚举 → 248 接触通过 → 16 不等价初态 |
| GFN2 | 16 个初态完成优化；排除 1 个方向不稳定驻点；15 个保留终态归为 4 个代表 |
| 同几何标定 | 21 套冻结几何，29 个 DFT 生产配对；主统计只使用 16 个初态 |
| 严格固定胞 DFT | 4/4 代表通过原始 Fmax ≤ 0.0002 eV/Å；全部 6 对通过合并条件，归为 1 个容差组 |
| g4 变胞 | 原始 Fmax ≤ 0.002 eV/Å、六应力分量最大绝对值 ≤ 0.05 kbar；独立冷启动和三项数值加严通过 |

固定胞和变胞是两个不同结果。固定胞为 a/b/c = 7.432/4.929/2.545 Å；变胞为约 7.129506/4.908356/2.555579 Å。固定胞原子精修不等于零压平衡胞，变胞终态也不能沿用固定胞更严的力阈值。

经验 σ = 1.920698、3σ = 5.762093 meV/CH2，仅描述声明的 16 初态样本域。本轮未设置能量窗口，**依据 3σ 自动淘汰的候选数为 0**；四个 GFN2 代表均进入 DFT 精修。

## 仓库入口

| 路径 | 用途 |
| --- | --- |
| [docs/METHODS.md](docs/METHODS.md) | 科学流程、参数、单位及验收规则 |
| [docs/REPRODUCING.md](docs/REPRODUCING.md) | 安装、可执行入口、验证范围与运行边界 |
| [docs/STATUS.md](docs/STATUS.md) | 已完成、取消和未开展工作的明确区分 |
| [docs/PUBLICATION.md](docs/PUBLICATION.md) | 正式论文版本与引用准备 |
| [snapshot/](snapshot/) | 保持项目相对目录结构的原始文件精选快照 |
| [snapshot_pa6/](snapshot_pa6/) | 来自独立 PA6 工作区的单候选先导、终态与诊断证据 |
| [docs/PA6_PROGRESS_20261031.md](docs/PA6_PROGRESS_20261031.md) | 当前 PA6 成果、问题与未执行的下一步对照 |
| [GPT Pro 咨询与本地评议](docs/reviews/20261031-gpt-pro/README.md) | 2026-10-31 已完成的网页咨询、回复原文及建议边界；不构成新增计算 |
| [scripts/verify_snapshot.py](scripts/verify_snapshot.py) | 文件完整性和科学计数核对；不启动量化计算 |
| [scripts/reproduce_candidates.py](scripts/reproduce_candidates.py) | 在独立输出目录复现候选生成 |
| [requirements-analysis.txt](requirements-analysis.txt) | 轻量分析与候选生成依赖 |
| [PA6 计算实验方案](snapshot/reports/PA6计算实验设计_20261029/计算实验方案.md) | 2026-10-29 历史方案；后续单候选进度见 PA6 独立快照 |

快速核对与候选生成：

```bash
python scripts/verify_snapshot.py
python -m pip install -r requirements-analysis.txt
python scripts/reproduce_candidates.py --output work/candidates
```

这些命令不代表已重新完成 GFN2、DFT 或全部实验对照。量化计算所需的外部软件、PAW 数据、环境与历史脚本路径要求请先阅读[复现说明](docs/REPRODUCING.md)。

## 结果适用范围

目前证据支持：在给定 PE 化学拓扑、晶胞、Pnam 先验及有限搜索域内，多起点经分层优化收敛到一个严格固定胞容差组，其代表完成零压正交变胞及数值验收。

未完成周期声子、完整 Hessian、零点能、有限温度核自由能或最终 g4 的跨泛函几何优化。索引 XRD 峰位比较不是连续强度谱拟合；不报告 Rwp、整谱匹配率或全局最低能证明。仓库同时保存 PA6 后续方案及其三分支设计、已完成的单候选 GFN2 先导与诊断；PA6 的新版完整候选库、D1/D14 和 DFT 均未运行，不属于已完成的 PE 计算结果。

仓库公开不等于已授予开源许可；当前未附授权许可证，也未登记论文 DOI、归档 DOI 或正式作者名单。使用和引用状态见[发表准备说明](docs/PUBLICATION.md)。
