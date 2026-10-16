# V6.3 PHASE 2 + 3 — 单变量与多变量机制分析

- 时间: 2026-10-14 13:50:10
- **分析单位 = target**（N ≈ 203）；禁止 row-level 当独立样本

- 载入: **203** targets × 111 cols

- 可分析数值特征: **107**（非空 ≥50 个 target）

## PHASE 2 — 单变量关联（Spearman + BH-FDR）

### 2.3 必须报告的变量（针对 primary outcome `abs_mean_error`）

| 特征 | 类型 | Spearman ρ | p | BH q | n |
|---|---|---:|---:|---:|---:|
| `log_n_mutations` | DEPLOYABLE | -0.0268 | 0.705 | 0.789 | 203 |
| `n_mutations` | DEPLOYABLE | -0.0268 | 0.705 | 0.789 | 203 |
| `n_unique_mutation_sites` | DEPLOYABLE | -0.1166 | 0.0975 | 0.255 | 203 |
| `site_coverage_all` | DEPLOYABLE | -0.0682 | 0.333 | 0.554 | 203 |
| `site_coverage_interface` | DEPLOYABLE | -0.1396 | 0.0469 | 0.173 | 203 |
| `site_entropy` | DEPLOYABLE | -0.1553 | 0.027 | 0.132 | 203 |
| `mutation_pair_entropy` | DEPLOYABLE | -0.1004 | 0.154 | 0.338 | 203 |
| `alanine_fraction` | DEPLOYABLE | -0.0254 | 0.719 | 0.796 | 203 |
| `single_mut_fraction` | DEPLOYABLE | -0.0491 | 0.487 | 0.66 | 203 |
| `ddg_std` | LABEL_DERIVED | +0.3064 | 8.72e-06 | 0.00015 | 203 |
| `ddg_range` | LABEL_DERIVED | +0.2319 | 0.00087 | 0.0112 | 203 |
| `ddg_IQR` | LABEL_DERIVED | +0.2668 | 0.000119 | 0.00176 | 203 |
| `struct_n_iface_res_lig` | DEPLOYABLE | +0.1448 | 0.0393 | 0.161 | 203 |
| `struct_n_iface_res_rec` | DEPLOYABLE | +0.0833 | 0.238 | 0.453 | 203 |
| `struct_n_contact_pairs` | DEPLOYABLE | +0.1817 | 0.00946 | 0.075 | 203 |
| `struct_dsasa_total` | DEPLOYABLE | +0.0549 | 0.437 | 0.634 | 203 |
| `struct_interface_residue_ratio` | DEPLOYABLE | +0.1654 | 0.0183 | 0.102 | 203 |
| `mutations_per_site` | DEPLOYABLE | +0.1691 | 0.0159 | 0.102 | 203 |
| `mutation_count_gini` | DEPLOYABLE | -0.0459 | 0.515 | 0.681 | 203 |
| `conservative_fraction` | DEPLOYABLE | -0.0384 | 0.586 | 0.724 | 203 |

### Top associations（按 |ρ|）

**DEPLOYABLE**（outcome = `abs_mean_error`）

| 特征 | ρ | p | BH q |
|---|---:|---:|---:|
| `signed_mean_error_M0` | -0.4774 | 5.97e-13 | 1.54e-11 |
| `struct_chainpair_contact_density` | +0.2230 | 0.00138 | 0.0146 |
| `struct_contacts_per_res_total` | +0.1845 | 0.00842 | 0.0723 |
| `struct_n_contact_pairs` | +0.1817 | 0.00946 | 0.075 |
| `pred_mean_ddg` | -0.1789 | 0.0106 | 0.0783 |
| `struct_n_contact_res_pairs` | +0.1757 | 0.0121 | 0.0834 |
| `mutations_per_site` | +0.1691 | 0.0159 | 0.102 |
| `pred_mean_M0` | -0.1669 | 0.0173 | 0.102 |
| `struct_weakest_chainpair_dsasa` | +0.1664 | 0.0188 | 0.102 |
| `struct_interface_residue_ratio` | +0.1654 | 0.0183 | 0.102 |

**LABEL_DERIVED**（outcome = `abs_mean_error`）

| 特征 | ρ | p | BH q |
|---|---:|---:|---:|
| `squared_mean_error` | +1.0000 | 0 | 0 |
| `extreme_fraction` | +0.5940 | 9.55e-21 | 4.92e-19 |
| `signed_mean_error` | -0.4813 | 3.6e-13 | 1.24e-11 |
| `true_mean_ddg` | +0.4130 | 9.11e-10 | 1.88e-08 |
| `ddg_std` | +0.3064 | 8.72e-06 | 0.00015 |
| `ddg_IQR` | +0.2668 | 0.000119 | 0.00176 |
| `ddg_range` | +0.2319 | 0.00087 | 0.0112 |
| `deleterious_fraction` | +0.2226 | 0.00142 | 0.0146 |
| `beneficial_fraction` | -0.1890 | 0.00691 | 0.0647 |
| `within_kendall` | -0.0772 | 0.315 | 0.541 |

### 2.2 分组变量

| 分组变量 | 检验 | 统计量 | p | 组数 | effect size |
|---|---|---:|---:|---:|---:|
| assay 类型 | Kruskal-Wallis | 12.747 | 0.1209 | 9 | ε²=0.0249 |
| `n_mutations ≥30` vs `<30` | Mann-Whitney | 4737.0 | 0.3991 | 63 vs 140 | rank-biserial=0.0741 |

- ≥30 组 |mean error| 中位 **0.6577**（n=63）vs <30 组 **0.8303**（n=140）

## PHASE 3 — 多变量 / 混杂控制

- **A_size_only**（1 特征）：family-CV R² = **-0.0226**，MAE = 0.7573，Spearman(pred,obs) = -0.1805
- **B_library_composition**（6 特征）：family-CV R² = **-0.0213**，MAE = 0.7532，Spearman(pred,obs) = +0.2172
- **C_context**（12 特征）：family-CV R² = **-0.0783**，MAE = 0.7744，Spearman(pred,obs) = +0.1968
- **D_explanatory_only**（9 特征）：family-CV R² = **+0.0655**，MAE = 0.7081，Spearman(pred,obs) = +0.3474

## 3.7 Cross-validation 汇总（GroupKFold by family）

| 模型 | 特征数 | CV R² | CV MAE | CV Spearman |
|---|---:|---:|---:|---:|
| A_size_only | 1 | -0.0226 | 0.7573 | -0.1805 |
| B_library_composition | 6 | -0.0213 | 0.7532 | +0.2172 |
| C_context | 12 | -0.0783 | 0.7744 | +0.1968 |
| D_explanatory_only | 9 | +0.0655 | 0.7081 | +0.3474 |

## 3.8 Partial association（残差化）

**方法**：对 `log_n_mutations` 与 `abs_mean_error` 分别回归掉控制变量，
再对**残差**求 Spearman。同时报告普通 Spearman。

| 控制集 | 普通 ρ | partial ρ | p | n |
|---|---:|---:|---:|---:|
| none | -0.0268 | -0.0268（= 普通）| — | 203 |
| n_unique_sites | -0.0268 | **-0.0132** | 0.852 | 203 |
| diversity | -0.0268 | **+0.1423** | 0.0428 | 203 |
| diversity+alanine | -0.0268 | **+0.1479** | 0.0352 | 203 |
| full_B | -0.0268 | **+0.1520** | 0.0304 | 203 |
| full_C | -0.0268 | **+0.1914** | 0.00622 | 203 |

## Primary mechanism decision

- 普通 Spearman(log_n_mut, |error|) = **-0.0268**
- 控制 diversity 后 |partial ρ| 最大 = **0.1914**（最小 p = 0.00622）
- 衰减比例 = -615.14%

> ### **`INCONCLUSIVE`**
>
> 控制后关联衰减到中间水平，无明确判定。

- 耗时 8s
