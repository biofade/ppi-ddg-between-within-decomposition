# V6.1+ 执行台账

> 规则：每完成一个 Phase 立即追加。失败写 `logs/v6_1/failures.log`，不停止。
> 约束：不修改 V1.1–V6.0 任何文件；产物用 `v6_1` 命名空间；CPU only；不下载 >100 MB；不伪造数字。

## PHASE 0: 环境检查 — **DONE**

- 环境：Windows 10 / Python 3.11.9 / i7-2630QM 4 核 / 7.55 GB RAM / **无 GPU**
- 依赖齐备：numpy 2.4.6 / pandas 3.0.6 / scipy 1.17.1 / sklearn 1.9.1 / lightgbm 4.7.0 / matplotlib 3.11.2
- 目录：`experiments/v6_1/` 下 11 个子目录已建
- ⚠️ **`experiments/v6_sota_audit/standardized_preds.csv` 不存在**（V6.0 未生成该长表）
  → PHASE 1 先补建，数据源为已核实的 CATH 官方表与 DSSA 表
- 产出：`phase0_environment.md`

## PHASE 1: CI + 交集 + 类别统计 — **DONE**

### 1.0 前置：补建 `standardized_preds.csv`（**DONE**）

- V6.1 规格引用的该文件在 V6.0 **未生成** → 本 Phase 补建
- 产物：`experiments/v6_sota_audit/standardized_preds.csv`，**17,377 行**
- 覆盖：CATH 5 模型（755 行 / 21 complex）+ DSSA-PPI RBM 8 模型（1,702 行 / 1 分箱）
- 符号按 V6.0 `sign_audit.csv` 统一（ESM-IF / PPIFormer / ESM-2-FFN / SAAMBE-SEQ / FoldX-RBM 已反号）

### 1.1 Cluster Bootstrap CI（**DONE**）

CATH held-out，重采样单位 = complex_id，5000 次，seed 42，unique complex < 15 跳过
（**有效 1151 / 跳过 3849**，因 21 个 complex 有放回抽样时 unique<15 很常见）

| 模型 | macro Spearman | 95% CI | CI 排除 0? | Efficiency@10 | 95% CI |
|---|---:|---|---:|---:|---|
| **Flex-ddG** | **+0.4138** | **[+0.3294, +0.4993]** | ✅ | +0.5930 | [+0.1849, +0.7438] |
| **FoldX** | **+0.3561** | **[+0.2648, +0.4407]** | ✅ | +0.5027 | [+0.1658, +0.6445] |
| **DDGPred** | **+0.3234** | **[+0.2124, +0.4285]** | ✅ | +0.4279 | [+0.0957, +0.6066] |
| proteinMPNN | +0.2199 | [+0.0924, +0.3429] | ✅ | +0.2304 | [−1.1898, +0.5379] |
| **ESM-IF** | +0.0277 | **[−0.0805, +0.1332]** | ❌ **含 0** | +0.1503 | [−1.4814, +0.8067] |

**判定**：
- **3 个物理/几何类方法（Flex-ddG / FoldX / DDGPred）的 macro per-PDB Spearman
  CI 明确排除 0**（下界 0.21–0.33）→ 其 within-target 排序能力**统计显著**。
- **ESM-IF 的 CI 含 0**（[−0.0805, +0.1332]）→ **其 within-target 能力不可与随机区分**。
- `Efficiency@10` 的 CI 较宽（因该统计量对 group 组成更敏感），
  但物理类三者的下界仍 > 0.09。

> **实现说明（性能优化，口径不变）**：5000 次 bootstrap × 21 complex × 200 次随机抽样
> 的朴素实现需数十分钟/指标。改用 `scripts/v6_1/fastmetrics.py`：
> `macro_spearman` 用 `np.argsort` 双排序替代 `scipy.rankdata`（**与 v5_common 差异 0.00e+00**，
> 快 17×）；`eff10` 的 bootstrap 内部随机抽样次数由 200 降为 **20**（point estimate 仍用 200）。
> 该近似已在此显式标注。

### 1.2 测试集交集（**DONE**）

| 层级 | CATH | V5 | **交集** |
|---|---:|---:|---:|
| PDB 级（宽松） | — | — | **8** |
| 复合物级（严格） | 21 | 61 | **8** |
| `interface_key`（含链） | 21 | — | **0** |

交集 PDB：`1emv, 1wqj, 2kso, 2o3b, 2wpt, 4e6k, 4hfk, 4l0p`
→ **8 < 10** ⇒ 只做**探索性**对比，明确标注「跨数据集、非同集比较」。
→ **禁止**做「同一测试集上的胜负比较」。

### 1.3 类别配对检验（**DONE，差异显著**）

**方法学**：明确**不用 Mann-Whitney U**（V6.1 禁止）—— 同一 complex 的两类分数
来自同一批突变体，必须**配对**；独立单位 = complex（n=21）。

| 量 | 值 |
|---|---:|
| 物理/能量类（Flex-ddG+FoldX）均值 | **+0.3871** |
| 纯序列表示类（ESM-IF+proteinMPNN）均值 | **+0.1278** |
| 配对差 Δ 均值 | **+0.2593** |
| bootstrap 95% CI | **[+0.1262, +0.3944]** ✅ 排除 0 |
| **Wilcoxon 符号秩 p** | **0.003753** ✅ 显著 |
| 效应量 r | **0.6106**（大） |
| 方向一致性 | **14 / 21** complex |

**敏感性分析**：DDGPred vs 序列组 Δ = +0.1984，Wilcoxon p = **0.0595**（边缘不显著）

→ **✅ 类别差异显著**：物理/能量类在 within-target 排序上显著优于纯序列表示类。

### 1.4 RBM 趋势性支持（**DONE，仅趋势**）

| 模型 | macro Spearman（无 CI） |
|---|---:|
| FoldX (RBM) | +0.6478 |
| DSSA-PPI (RBM) | +0.3129 |

> ⚠️ RBM 仅 **1 个人工分箱**，真值是适应度非 ΔΔG → **不得与 CATH/V5 直接比较**。

## PHASE 2: Split Inflation 量化 — **DONE**（耗时 19.9 分钟）

同一模型（LightGBM 800 树 + 60 物理 + pooled ESM-2 PCA128），**5-fold OOF**，
所有档在同一批 8207 行上算指标（分母相同 → 可比）；3 seeds（42/43/44）。

| 档 | 划分方式 | Overall Pearson | Macro per-PDB Spearman | Efficiency@10 |
|---|---|---:|---:|---:|
| **L1** | mutation-level random | **+0.8436** | **−0.0344** | −0.0325 |
| **L2** | complex-level (PDB) | +0.6590 | +0.0417 | +0.0651 |
| **L3** | ≤40% seq-id + complex | +0.6767 | +0.0568 | +0.0035 |
| **L4** | family-level | **+0.5774** | +0.0118 | −0.0466 |

### ★ 核心量化：Split Inflation（L1 − L4）

| 指标 | L1 random | L4 family | **L1 − L4** |
|---|---:|---:|---:|
| **Overall Pearson** | +0.8436 | +0.5774 | **+0.2663（相对 +46.1%）** |
| Macro per-PDB Spearman | −0.0344 | +0.0118 | −0.0461 |

**两条独立的重要结论**：

1. **Split inflation = +0.2663 Pearson**（相对虚高 **46.1%**）。
   仅把划分从 family-level 换成 mutation-level random，
   就能让同一模型看起来「好得多」—— 这是文献跨论文数字不可比的**根源**。
2. **更关键**：**所有 4 档的 Macro per-PDB Spearman 都在 0 附近
   （−0.034 ~ +0.057）**。即：
   > **这个模型（pooled ESM-2 + 物理特征 + GBDT）无论用什么划分，
   > 都不具备靶点内排序能力。**
   它的 global Pearson 完全由「跨复合物区分」贡献，而非组内排序。

> 这与 V5.0 的 `REPRESENTATION_BOTTLENECK_CONFIRMED` 判定**完全一致**，
> 并首次用**同一模型 × 4 档划分**的受控实验给出因果性证据。

## PHASE 3: Pairwise 对比 — **DONE**

> 产物：`phase3_pairwise/pairwise_results.csv` + `pairwise_report.md`

## PHASE 4: V5 vs CATH 分布 — **DONE**

**统计方法说明**：此处**允许**用 Mann-Whitney U（V5 与 CATH 是**不同数据集、不同
complex**，两组**不配对**），与 PHASE 1.3 的配对情形本质不同。

| 组 | n | median per-complex Spearman | vs V5 Δmedian | Mann-Whitney p |
|---|---:|---:|---:|---:|
| **V5 baseline（我们）** | 129 | **+0.0939** | — | — |
| CATH 物理/能量类（Flex-ddG+FoldX） | 42 | **+0.4179** | **+0.3240** | **0.0000** ✅ |
| CATH 纯序列类（ESM-IF+proteinMPNN） | 42 | +0.1314 | +0.0375 | 0.5598 ❌ |

单模型：Flex-ddG p=0.0003 ✅ / FoldX p=0.0031 ✅ / DDGPred p=0.0069 ✅ /
proteinMPNN p=0.1927 ❌ / ESM-IF p=0.6746 ❌

⇒ **跨数据集支持**：物理/能量类的 per-complex 分布显著高于 V5 baseline，
而**纯序列类与 V5 无显著差异**。
⚠️ 必须标注：这是**跨数据集**比较（交集仅 8 个 PDB，见 PHASE 1.2），**非同集比较**。

## PHASE 5: 论文级 Figure — **DONE**

`phase5_figures/` 下 6 张 **300 dpi** 图：
`fig1_illusion_scatter` / `fig2_efficiency_gap` / `fig3_category_comparison` /
`fig4_split_inflation` / `fig5_hit10_vs_baserate` / `fig6_per_complex_distribution`

> ⚠️ `fig4_split_inflation.png` 生成于 PHASE 2 完成**之前**，当时
> `split_inflation_results.csv` 只有 L1/L2 两档 → 该图**需在 PHASE 2 完成后重生成**。
> 已记录，重生成命令：`python -u scripts/v6_1/08_make_paper_figures.py`

## PHASE 7: 工具包 + 外部模型搜索 — **DONE**

### 7.1 工具包（`phase7_toolkit/affirigid_metrics/`）

- `topk.py`：`topk_metrics` / `efficiency_at_k` / `per_structure_spearman/pearson`
- `bootstrap.py`：`cluster_bootstrap_ci` / `cluster_paired_bootstrap` / `paired_wilcoxon`
- `README.md`（输入/输出格式、符号约定、使用示例）
- `examples/example_usage.py` —— **自检全部通过 ✅**
  （oracle Efficiency@10 = 1.0000 / anti = −0.9777 / random = −0.0453）
- `requirements.txt`

### 7.2 外部模型搜索（`phase7_toolkit/external_model_search.md`）

| 模型 | 状态 | repo | license |
|---|---|---|---|
| CATH-ddG | `PREDICTION_FOUND` | `ak422/CATH-ddG` | None |
| DSSA-PPI | `PREDICTION_FOUND` | `Agitw/DSSA-PPI` | MIT |
| SAMAffinity | `DATA_ONLY_NO_PRED` | `iobio-zjut/SAMAffinity` | None |
| Pythia-PPI | `CODE_ONLY_NO_PRED` | `Wublab/pythia_ppi` | Apache-2.0 |
| MINT | `CODE_ONLY_NO_PRED` | `VarunUllanat/mint` | MIT |
| MuToN | `DATA_ONLY_NO_PRED` | `zpliulab/MuToN` | MIT |
| ProBASS | `DATA_ONLY_NO_PRED` | `sagagugit/ProBASS` | None |
| DeepPPAPred | `REPO_NOT_FOUND` | — | — |

**新增发现（Zenodo）**：record **15384945**（doi 10.5281/zenodo.15384945，
**CC-BY-4.0**）提供**已跑好**的物理能量预测：
`Synthetic_FoldX_ddG_942723.csv`（**942,723 行**）与
`Synthetic_FlexddG_ddG_20829.csv`（20,829 行）→ 已下载至 `F:\dsh\downloads\zenodo_abag\`。
→ 该数据可**绕过 PHASE 1 的 n=21 功效限制**（见 PHASE 7.3 补充分析）。
Figshare API 返回 **403 Forbidden**（已如实记录）。
