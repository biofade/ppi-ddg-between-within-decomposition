# BENCHMARK NOVELTY MATRIX（§25）

> **诚实规则**：除下面**已核实**的两项外，其余非 ours 行仍为 `UNKNOWN`。
> 依 §25，**不为制造新颖性而软化既有工作**；在缺乏核实的情况下填 YES 或 NO 都是错的。

> ⚠️ **本轮 PHASE 1 的核实结果显著削减了 V8 的新颖性主张，见下节「已核实的重叠」。**

---

## 0. 已核实的重叠（**必须先读**）

### R1 — Bushuiev et al. 2024，*Revealing data leakage in protein interaction benchmarks*
[arXiv:2404.10457](https://ar5iv.labs.arxiv.org/html/2404.10457)（ICLR 2024 workshop 系列）

**这篇论文直接覆盖了 V8 的 H4（划分膨胀）核心主张。** 其量化结论：

| 划分策略 | 结构性泄漏率 |
|---|---:|
| PPI-code 划分（如 `3BTD_I_E`） | **86%** |
| PDB-code 划分 | **65%** |
| 按 PDB 提交时间 | 61% |
| **序列相似度划分（MMseqs2, 30% identity）** | **30%** |
| Linclust 聚类 | 37% |
| **SKEMPI v2.0 + PPI-code 划分（Luo et al. 2023）** | **56%** |
| **SKEMPI v2.0 + 领域专家分组（Bushuiev et al. 2024）** | **0%** |

其结论与建议：**元数据与序列相似度划分都不够**，应使用**界面 3D 结构相似度**（iDist / iAlign / TM-align）划分。

**⇒ 因此 V8 必须撤回/放弃以下潜在主张（§4/§71）：**
- ❌ 「首次发现宽松划分会膨胀 PPI ΔΔG 性能」—— **已被 R1 占据**
- ❌ 「首次建议使用非泄漏划分」—— 已被 R1 占据（且 R1 的建议比 V8 的 family-OOD 更强）

**V8 相对 R1 仍然独有的部分**（诚实表述）：
- R1 量化的是**泄漏率**（测试样本有多少有近重复训练样本）；
  V8 量化的是**这种泄漏对三轴能力指标的差异化影响**，尤其是
  `C_B`（跨靶定位信号）从 **+0.43 塌到 ≈0** 的**协方差分解**机制，
  以及「Overall Pearson 膨胀 4.5× 而 Efficiency@10 只膨胀 2.5%」的**不均匀性**。
- R1 未做**模型排名反转**、**目标加权敏感性**、**匹配宇宙**、**能力轴分解**。

### R2 — PPIformer 的划分与污染（§50）

PPIformer（[repo](https://github.com/anton-bushuiev/PPIformer)，MIT，[arXiv:2310.18515](https://arxiv.org/abs/2310.18515)）
**在 SKEMPI v2.0 上微调**，且其作者正是 R1 的作者群 —— 他们**自己**采用了
基于**界面相似度**的非泄漏划分（`Bushuiev et al. 2024`，R1 中报告 0% 泄漏）。

**⇒ 两个后果：**
1. PPIformer **已属于 R1 意义上的「非泄漏基准」阵营**；
   V8 若拿它做对比，**不能**把它当作「被泄漏基准误导」的例子。
2. 但依 §50，**SKEMPI v2.0 微调的 checkpoint 不能自动进入 Tier 1** ——
   必须显式审计其训练家族与 V8 `OUTER_FOLDS_V7_2` 的家族是否重叠。

### R3 — V8 采用的 `seq40` 划分正属于 R1 证明**仍有 30% 泄漏**的那一类

V8 的 `U_MAS` 含 `L3_seq40` 划分（序列 40% 同一性簇）。
R1 证明序列相似度划分（30% identity）**仍有 30% 结构性泄漏**。

**⇒ V8 必须在论文中承认：`seq40` 不是干净的划分**，
而 `L4_family`（家族）也只是 R1 意义上的**近似**；
R1 推荐的**界面相似度划分** V8 **未实现**。这是一个明确的局限。

---

## 1. 新颖性矩阵（未核实项保持 UNKNOWN）

> `papers/*.txt`（12 篇）的机器候选抽取见 `benchmark_landscape.csv`，
> 仅为**候选**，逐条附命中原文，需人工核实。

| method | strict task-specific family-OOD | overall metrics | per-target ranking | explicit inter/within decomposition | calibration axis | selection utility axis | target weighting audit | matched-universe comparison | rank reversal analysis | candidate-set calibration | benchmark-guided remedy case study |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `PreFold-dG` | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN |
| `CATH-ddG` | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN |
| `PPIformer` | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN |
| `HERMES` | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN |
| `Prompt-DDG` | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN |
| `RDE-Network` | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN |
| `BA-DDG` | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN |
| `EBM-DDG` | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN |
| `Light-DDG` | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN |
| `Pythia-PPI` | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN |
| `H3-DDG` | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN |
| `MARS-DDG` | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN |
| `USP-ddG` | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN |
| `ReCLIP` | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN |
| `CR-PPI-ddG (ours)` | YES | YES | YES | YES | YES | YES (但区分度不足, 见 DECISION_ALIGNMENT_REPORT) | YES | PARTIAL (row_id 严格匹配; 外部未审计) | YES | NOT_RUN | NOT_RUN |

## 已可确认的「新颖性已被占据」项（避免 §4 的禁用声明）

以下概念在本项目**自己的历史工作里就已出现**，因此 V8 **不得**声称 first：

- per-target / per-structure Spearman（V6.1/V6.2 已用）
- family-OOD 划分（OUTER_FOLDS_V7_2 已用）
- Top-K prioritization（V6.2 / V7.3 Eff@10 已用）
- baseline/mutation 分解（V7.1 `C_B`/`C_W` 已用）

## V8 可能仍然独有的组合（**待外部核实后才能定稿**）

- 严格 task-specific family-OOD **+** 精确匹配宇宙 **+** 三轴能力分解
  **+** 目标加权敏感性审计 **+** 模型秩反转 **+** decision alignment
  **+** benchmark-guided remedy 案例

> ⚠️ 这只是一个**候选组合**。在外部审计完成之前，
> **不得**在论文中据此声明任何优先级。
