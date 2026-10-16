# V6.1 PHASE 6 — 论文大纲

> **每个 Figure / Table 都标注了数据来源（Phase + CSV 路径）**，符合 V6.1 PHASE 6 要求。
> 本文件是**骨架**，不是完整论文。

---

# Title（候选）

**Within-Target Ranking, Not Global Correlation, Is What Separates Physics-Based and
Sequence-Based PPI ΔΔG Predictors**

（备选，更保守）
**A Same-Protocol Re-Evaluation of PPI ΔΔG Predictors: Global Correlation Masks
Within-Target Ranking Ability**

---

## Abstract（150 words）

Protein–protein interaction (PPI) ΔΔG predictors are increasingly reported with high
global correlation, yet it is unclear whether these numbers translate into the ability
to rank mutations **within a single target** — the actual decision faced in
affinity-maturation campaigns. We re-evaluated 13 published predictor outputs
(5 on the official CATH-ddG held-out split; 8 on the DSSA-PPI RBM benchmark) under one
protocol, using **cluster bootstrap** with the complex as the sampling unit, and
introduced `Efficiency@10` to measure how much of the oracle's top-10 utility a model
recovers. Physics/energy-based methods (Flex-ddG, FoldX) reached macro per-complex
Spearman **0.36–0.42** (CIs excluding zero), significantly exceeding pure
sequence-representation methods (ESM-IF, proteinMPNN; Δ = **+0.26**, Wilcoxon
**p = 0.0038**). We further show that a single model's overall Pearson changes by
**0.86** purely from the train/test split definition.

> ⚠️ **Abstract 的最终数字必须从 PHASE 1/2 的 CSV 重新读取后填入**，不得沿用本骨架中的示例值。

---

## 1. Introduction

- **段落 1 — 重要性**：PPI 突变效应预测（ΔΔG）在抗体亲和力成熟、酶工程、
  变异致病性解释中的核心地位。
- **段落 2 — 现有 SOTA 与它们的高 Pearson**：CATH-ddG 报 held-out CATH
  Overall Pearson **0.6150**（`LITERATURE_EVIDENCE`，btaf228 Table 2）；
  各方法普遍以 global Pearson / Spearman 为主要指标。
- **段落 3 — 但问题未解**：高 global 指标可能主要来自**跨复合物区分**
  （某些靶点整体更负），而非**靶点内排序**；且不同论文使用不同划分，
  数字不可比。我们前期工作（V5.0）实测：同一模型仅换划分，
  Overall Pearson 可差 **+0.50**（V2.2 记录）。
- **段落 4 — 本文贡献**：
  1. 提出 `Efficiency@10`（相对 oracle 的 top-10 效用回收率）；
  2. 在 **同一协议** 下重评 13 个已发表预测输出（cluster bootstrap，单位=complex）；
  3. 用配对检验证明**物理/能量类显著优于纯序列表示类**（Wilcoxon p = 0.0038）；
  4. 用 4 档划分量化 **split inflation**；
  5. 明确区分「数据集不同 → 只可比分布，不可比数值」。

---

## 2. Methods

### 2.1 Datasets

| 数据集 | 来源 | 规模 | 划分 |
|---|---|---|---|
| CATH-ddG held-out | SKEMPI v2.0，CATH 同源超家族 hold-out | 813 突变 / 53 PPI；**过滤后 755 / 21 complex** | 序列同一性 <30%，无 easy mutation |
| DSSA-PPI RBM | deep mutational scanning 风格 | 1,702 行 / 1 分箱 | case study（真值为适应度） |
| PPB-Affinity（我们） | 五源聚合 | 8,207 行；ΔΔG 配对 7,523 | 复合物级 + ≤40% seq-id + 家族级 |

### 2.2 Models evaluated

**CATH 上 5 个**：Flex-ddG、FoldX、DDGPred、ESM-IF、proteinMPNN
**RBM 上 8 个**：DSSA-PPI、PPIFormer、ESM-2-FFN、MuTon、MpbPPI、ProBASS、
SAAMBE-SEQ、FoldX
**我们的基线**：LightGBM 800 树 + 60 物理特征 + pooled ESM-2 PCA128 + 反对称块 + 正反向增广

### 2.3 Metrics

- `Overall Pearson`（全局，SOTA 常报）
- `Macro per-PDB Spearman`（靶点内排序均值，**本文核心**）
- `Hit@10 (ΔΔG<0)` 与 `Base Rate`（绝对命中率 vs 随机期望）
- **`Efficiency@10 = (Random − Model) / (Random − Oracle)`**（本文提出）

### 2.4 Statistical analysis

- **Cluster bootstrap**：抽样单位 = complex，携带该 complex 全部突变行；5000 次
- **配对 Wilcoxon 符号秩**（PHASE 1.3；两类模型作用于**同一批突变体** → 必须配对）
- **Mann-Whitney U** 仅用于**不配对**情形（PHASE 4：跨数据集分布对比）
- **符号审计**：用独立第三方真值（本地 SKEMPI v2.0 原始表）逐步复算 ΔΔG 验证约定
  （Pearson +0.9869、符号一致 94.99%）

---

## 3. Results

### 3.1 类别对比：物理 vs 序列（**Figure 3**）

> **数据来源**：`experiments/v6_1/phase1_ci/category_paired_values.csv`（PHASE 1.3）

- 物理/能量类（Flex-ddG+FoldX）per-complex Spearman 均值 **+0.3871**
- 纯序列类（ESM-IF+proteinMPNN）**+0.1278**
- 配对 Δ = **+0.2593**，95% CI **[+0.1262, +0.3944]**，**Wilcoxon p = 0.003753**，
  效应量 r = **0.6106**（大），方向一致 **14/21**

### 3.2 单模型 CI（**Table 1**）

> **数据来源**：`experiments/v6_1/phase1_ci/bootstrap_ci_models.csv`（PHASE 1.1）

| 模型 | macro Spearman | 95% CI | 排除 0 |
|---|---:|---|---|
| Flex-ddG | +0.4138 | [+0.3294, +0.4993] | ✅ |
| FoldX | +0.3561 | [+0.2648, +0.4407] | ✅ |
| DDGPred | +0.3234 | [+0.2124, +0.4285] | ✅ |
| proteinMPNN | +0.2199 | [+0.0924, +0.3429] | ✅ |
| ESM-IF | +0.0277 | [−0.0805, +0.1332] | ❌ |

### 3.3 Split Inflation（**Figure 4**）

> **数据来源**：`experiments/v6_1/phase2_split_inflation/split_inflation_results.csv`（PHASE 2）

同一模型在 4 档划分下的 Overall Pearson（**数字待 PHASE 2 完成后填入**）：

| 档 | 划分 | Overall Pearson | Macro per-PDB Spearman |
|---|---|---:|---:|
| L1 | mutation-level random | *待填* | *待填* |
| L2 | complex-level | *待填* | *待填* |
| L3 | ≤40% seq-id + complex | *待填* | *待填* |
| L4 | family-level | *待填* | *待填* |
| **L1 − L4** | **inflation** | ***待填*** | *待填* |

> **已知初步读数**（PHASE 2 冒烟测试，L1 单 seed）：
> Overall Pearson **+0.8446** 而 Macro per-PDB Spearman **−0.0257**
> → **随机划分能造出 0.84 的全局相关，却完全不产生靶点内排序能力**。

### 3.4 V5 baseline 同协议对比（**Figure 6**）

> **数据来源**：`experiments/v6_1/phase4_v5_comparison/v5_cath_distribution.csv`（PHASE 4）
> 与 `experiments/v6_1/phase1_ci/intersection_report.md`（PHASE 1.2）

- **交集**：CATH test 与 V5 test 在 PDB 级与复合物级均只有 **8 个**共同 PDB；
  `interface_key`（含链）交集为 **0** → **非同集**，只做分布对比
- V5 baseline median per-complex Spearman vs CATH 物理类（Mann-Whitney p 待 PHASE 4 填入）

---

## 4. Discussion

### 4.1 为什么物理方法在 within-target 上更好

物理/能量函数（FoldX、flex ddG）**逐突变地计算局部环境能量变化**，
其归纳偏置天然对齐「单点替换改变了什么」；
而 pooled 表示把整条链平均，**位点信息被稀释**
（我们前期实测：`||ΔZ||` 与 ΔΔG 全局 Pearson +0.2263，但**组内 Spearman 仅 +0.0215**）。

### 4.2 为什么纯序列表示失效

ESM-IF / proteinMPNN 是**结构条件序列似然**，其输出尺度与 ΔΔG 的物理量纲无直接对齐；
本审计中 ESM-IF 的 macro CI 含 0，proteinMPNN 虽排除 0 但仅 +0.22。

### 4.3 MIPA 的定位

若**不含学习的局部物理能量**（FoldX）就能达到 macro per-PDB Spearman **0.36**，
则任何「学习局部扰动」的表示方法（如我们的 MIPA 假设）
**必须直接与 FoldX 对比**，而非只与 pooled baseline 对比。
这为 V5.1/V7 的 baseline 选择提供了硬性要求。

### 4.4 局限

1. **CATH 过滤后仅 21 个 complex** → CI 宽、macro 统计量抽样分布可能非正态
2. **RBM 只有 1 个分箱**且真值是适应度非 ΔΔG → 仅趋势性支持
3. **CATH 与 V5 交集仅 8 个 PDB** → 跨数据集比较，非同集
4. **未获取 4 个模型的预测文件**（SAMAffinity / Pythia-PPI / MINT / Prompt-DDG）
5. 我们的基线在 `real_dG` 口径代理标签上评估（PHASE 2 的 L1–L4），
   与 CATH 的 `ddG` 口径不同，**两部分的数字不可互相引用**

---

## 5. Conclusion

在统一协议下，**PPI ΔΔG 预测的 within-target 排序能力主要由方法的归纳偏置决定**：
逐突变的物理/能量方法显著优于池化序列表示方法（Wilcoxon p = 0.0038）。
同时，**单一模型的可见性能强烈依赖划分定义**（split inflation）。
因此，报告 ΔΔG 预测性能时应**同时给出 within-target 指标与划分定义**，
否则数字不可比。

---

## Data & Code Availability

| 内容 | 路径 / URL |
|---|---|
| 外部预测数据（CATH） | `https://github.com/ak422/CATH-ddG`（**license 未声明**），`data/SKEMPI2/skempi_v2_with_all_results.csv` |
| 外部预测数据（RBM） | `https://github.com/Agitw/DSSA-PPI`（**MIT**），`validate_on_RBM/RBM_baselines_comparison_ddG.csv` |
| 统一长表 | `experiments/v6_sota_audit/standardized_preds.csv` |
| 指标实现 | `scripts/v6_1/fastmetrics.py`（与 `scripts/v5/v5_common.py` 口径一致，自检差异 0.00e+00） |
| 本阶段全部脚本 | `scripts/v6_1/00–08*.py` |
| 本阶段全部结果 | `experiments/v6_1/` |

---

## 附：图表 ↔ 数据来源对照（PHASE 6 强制要求）

| 图/表 | 内容 | **数据来源（Phase + CSV）** |
|---|---|---|
| **Figure 1** | Overall Pearson vs Macro per-PDB Spearman | PHASE 1.1 `phase1_ci/bootstrap_ci_models.csv` + `phase2_split_inflation/split_inflation_results.csv` |
| **Figure 2** | Efficiency Gap（top-10 效用 vs Random/Oracle） | PHASE 1.1 `bootstrap_ci_models.csv`（utility 列）+ V6.0 `sota_audit_metrics.csv` |
| **Figure 3** | 类别对比 violin（物理 vs 序列 vs V5） | PHASE 1.3 `phase1_ci/category_paired_values.csv` |
| **Figure 4** | Split Inflation（L1→L4） | PHASE 2 `phase2_split_inflation/split_inflation_results.csv` |
| **Figure 5** | Hit@10 vs Base Rate | V6.0 `sota_audit_metrics.csv` + PHASE 2 `split_inflation_results.csv` |
| **Figure 6** | Per-complex Spearman 分布（V5 vs CATH） | PHASE 4 `phase4_v5_comparison/v5_cath_distribution.csv` |
| **Table 1** | 单模型 CI | PHASE 1.1 `bootstrap_ci_models.csv` |
| **Table 2** | Pairwise 显著性 | PHASE 3 `phase3_pairwise/pairwise_results.csv` |
| **Table 3** | 测试集交集 | PHASE 1.2 `phase1_ci/intersection_report.md` |
| **Table 4** | 符号审计 | V6.0 `v6_sota_audit/sign_audit.csv` |
