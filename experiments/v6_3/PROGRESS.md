# V6.3 执行台账

> 每完成一个 PHASE 立即追加。失败写 `logs/v6_3/failures.log`，不停止。
> 约束：不改 V1.1–V6.2.1；`v6_3` 命名空间；CPU only；不下载；不伪造数字；
> **只说 association，不说 causation**；严格区分 `DEPLOYABLE` / `LABEL_DERIVED`。

## PHASE 0: Provenance / Freeze — **DONE**

- 新建 `scripts/v6_3/`、`experiments/v6_3/`（9 子目录）、`logs/v6_3/`
- 冻结 **12** 个输入文件的 SHA256 → `phase0_provenance/input_manifest.csv`
- 元数据可用性预判：**`source` 全为 SKEMPI v2.0（单一来源）→ `UNAVAILABLE`**
- 产出：`phase0_provenance/provenance.md`

## PHASE 1: Target-level 机制表 — **DONE**

- **203 targets × 111 列**（来自 230 个 matched-set target，剔除 n<2 的 27 个）
- **0 个字段 UNAVAILABLE**；仅 `within_spearman`/`within_kendall` 非空 84.2%
- 修复两处 bug：
  - **BUG 1**：`parse_mutation` 原用三字母码正则，但数据是**单字母三字符码**
    （如 `A_P2A`）→ 导致 10 个 composition 字段全 NaN。改为单字母码优先，保留三字母兼容。
  - **BUG 2**：`STRUCT_COLS` 未过滤字符串列 → 产生 20 个全 NaN 的
    `struct_<字符串列>`。改为仅取 `is_numeric_dtype`。
- 产出：`phase1_target_table/target_mechanism_table.csv`、
  `feature_dictionary.csv`、`missingness_report.md`

## PHASE 2: 单变量关联 — **DONE**

**Primary outcome `abs_mean_error`（N=203）**：

| 特征 | 类型 | ρ | p | BH q |
|---|---|---:|---:|---:|
| `true_mean_ddg` | LABEL_DERIVED | **+0.4130** | 9.1e-10 | ✅ |
| `ddg_std` | LABEL_DERIVED | **+0.3064** | 8.7e-06 | ✅ |
| `ddg_IQR` | LABEL_DERIVED | +0.2668 | 1.2e-04 | ✅ |
| `ddg_range` | LABEL_DERIVED | +0.2319 | 8.7e-04 | ✅ |
| `struct_chainpair_contact_density` | DEPLOYABLE | **+0.2230** | 0.00138 | **0.0146** ✅ |
| `struct_contacts_per_res_total` | DEPLOYABLE | +0.1845 | 0.00842 | 0.0723 |
| `mutations_per_site` | DEPLOYABLE | +0.1691 | 0.0159 | 0.102 |
| **`log_n_mutations`** | DEPLOYABLE | **−0.0268** | **0.705** | — **不显著** |

- 分组检验：`n_mutations ≥30` vs `<30` 的 Mann-Whitney 见 `group_tests.csv`
- **发现并更正**：PHASE 1 字典曾把 `signed_mean_error_M0`/`pred_mean_*` 标为
  `DEPLOYABLE`，实为**模型输出派生量**；已在汇总中剔除说明
  （PHASE 6 白名单从未包含它们）
- 产出：`phase2_univariate/univariate_associations.csv`、`group_tests.csv`、
  `top_associations.md`

## PHASE 3: 多变量 / 混杂控制 — **DONE**

**family-level CV（GroupKFold by family）**：

| 模型 | 特征数 | CV R² | CV MAE | CV Spearman |
|---|---:|---:|---:|---:|
| A_size_only | 1 | **−0.0226** | 0.7573 | −0.1805 |
| B_library_composition | 6 | **−0.0213** | 0.7532 | +0.2172 |
| C_context | 12 | **−0.0783** | 0.7744 | +0.1968 |
| D_explanatory_only | 9 | **+0.0655** | 0.7081 | +0.3474 |

⇒ **所有可部署模型 CV R² ≤ 0**；仅含 label-derived 的 D 为正。

**Partial association（残差化）**：

| 控制集 | 普通 ρ | partial ρ | p |
|---|---:|---:|---:|
| none | −0.0268 | −0.0268 | — |
| n_unique_sites | −0.0268 | −0.0132 | 0.852 |
| diversity | −0.0268 | +0.1423 | 0.043 |
| diversity+alanine | −0.0268 | +0.1479 | 0.035 |
| full_B | −0.0268 | +0.1520 | 0.030 |
| full_C | −0.0268 | +0.1914 | 0.006 |

### 判定：**`INCONCLUSIVE`（但实质：`n_mutations` 不是解释变量）**

理由：单变量关联**本身不显著**（ρ=−0.0268, p=0.705）→ 不可能是 calibration failure
的解释变量。原脚本的自动判据因 `plain_rho` 接近 0 而产生「衰减 −615%」的
无意义比例，故落回 `INCONCLUSIVE`；**人工核定后的实质结论是
`SIZE_IS_PROXY`（且是「无效 proxy」）**。

产出：`phase3_multivariable/` 5 个文件

## PHASE 4: Family / Source 集中度 + Influence — **DONE**

- **`source` = SKEMPI v2.0（唯一值）→ `UNAVAILABLE`**
- Family 集中度：`>=30` 组占 **63 targets / 6112 rows**
- **Leave-one-largest-family-out：现象不消失** → **`NOT_SINGLE_FAMILY_DRIVEN`**
- **Leave-top-targets-out（关键）**：

| 移除 | macro mean-Pearson | **row-w C_B** |
|---|---:|---:|
| (基线) | +0.5568 | **−0.0121** |
| top 1 | +0.5595 | +0.0084 |
| top 5 | +0.5730 | **+0.0990** |
| top 20 | +0.5853 | **+0.1377** |

⇒ **大靶点是「拖累」而非「支撑」row-weighted C_B**

产出：`phase4_family_source/` 4 个文件

## PHASE 5: Weighting 敏感性 — **DONE**（★ 本轮最重要方法学发现）

| weighting | L1 C_B | **L4 C_B** | L4 macro mean-Pearson |
|---|---:|---:|---:|
| **row (n_c/N)** | +0.4323 | **−0.0121** | +0.5568 |
| **equal-target** | +0.7752 | **+0.5475** | +0.5568 |
| sqrt(n) | +0.4498 | +0.1645 | +0.5568 |
| capped 30 | +0.3427 | +0.1081 | +0.5568 |
| log(1+n) | +0.4487 | +0.2183 | +0.5568 |

- **L4 `C_B` 极差 = 0.5596，且符号翻转**
- **macro mean-Pearson 在所有 weighting 下完全相同（+0.5568）**

⇒ **「row-weighted C_B ≈ 0/负」是 weighting 的产物，不是模型能力的性质。**
V6.2/V6.2.1 的「between-target calibration 失效」应改述为
**「在 row-weighting 下 between 贡献被 mutation-rich 靶点稀释至 0」**。

产出：`phase5_influence/weighting_sensitivity.csv`、`weighting_report.md`

## PHASE 6: Candidate-Set Calibration Test — **DONE**

**Task**：`y_target = true_mean_ddg`；**family-level CV**；严格 deployable 白名单

| 模型 | 特征数 | Pearson | Spearman | RMSE | MAE |
|---|---:|---:|---:|---:|---:|
| **T** target-only | 23 | **−0.0438** | −0.0470 | 1.4697 | 0.9723 |
| **TS** target+set | 43 | **+0.2461** | +0.2083 | 1.4375 | 1.0294 |

`≥30` 子集：T **+0.0878** → TS **+0.2614**

**Paired target bootstrap（5000）**：
- **ΔPearson(TS−T) = +0.2877, 95% CI [+0.1237, +0.4504] ✅ 排除 0**
- ΔRMSE = −0.0303, CI [−0.1243, +0.0672] ❌ 含 0
- `≥30` ΔRMSE = −0.0557, CI [−0.2547, +0.1467] ❌ 含 0

### 判定：**`SET_CONTEXT_SUPPORTED`（限于 Pearson/Spearman）**

> ⚠️ 两条必须同时报告的限制：
> 1. **target-only Pearson = −0.044 → 完全无预测力**
> 2. **绝对误差未显著改善**（ΔRMSE CI 含 0；TS 的 MAE 反而略差 1.0294 vs 0.9723）

产出：`phase6_set_calibration/` 4 个文件

## PHASE 7: Mechanism Decision — **DONE**

**MECHANISM = `LIBRARY_COMPOSITION_DRIVEN`**

（唯一满足判据的分支；A 的「`n_mut` 单变量显著」子项**不**满足，
故结论限定为「**mutation-set 组成**」而非「mutation count」）

产出：`phase7_mechanism/mechanism_decision.md`

## PHASE 8: Final Report — **DONE**

见 `V6_3_FINAL_SUMMARY.md`。**Q1/Q2/Q3 三问均已明确回答。**

## 全部 9 个 PHASE 完成（0–8）；无 `failures.log`
