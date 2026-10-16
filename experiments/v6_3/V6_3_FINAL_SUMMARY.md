# V6.3 Mechanism Audit — Final Summary

**时间**：2026-10-14（约 13:55 → 15:05，约 70 分钟；上限 2 小时，未超时）
**约束**：未修改 V1.1–V6.2.1 任何文件；`v6_3` 命名空间；CPU only；无下载；未伪造数字
**原则**：只做机制诊断；**只说 association，不说 causation**；严格区分 `DEPLOYABLE` 与 `LABEL_DERIVED`

---

## 1. Inputs

- **Targets**: **203**（来自 230 个 matched-set target；剔除 n<2 的 27 个）
- **Rows**: 7,523（matched set，V6.2 canonical）
- **L4 predictions source**: `experiments/v6_2/phase2_split4/split4_predictions.npz`
  （`L4_family_M1_s{42,43,44}` seed-ensemble）
- **Frozen hashes**（PHASE 0）：
  - `matched_analysis_set.csv` = `8d209eea08cc6058fb235b059c9e3725fcbc9ba8`
  - `split4_predictions.npz` = `ba711ea1b349836ac18b44ba67f6ab2f0a0b680b`
  - `d_only_predictions.npz` = `98b64c2b02140f232d993251c32d5a02daea7f2e`
  - `target_level_stats.csv` = `63c6d2a905ca2bcfca6adee8ad61b92d32b2838b`
  - `v2.2_proteinea_features.csv` = `134a406175c975541ae40aed37c3f824ad151942`
  - `ddg_pairs.csv` = `e43b19d69a8b54efe5cc4b3e647e87657106e86c`
- **可用性**：`metadata` 非空率 100%；**0 个字段 UNAVAILABLE**
  - 唯一 UNAVAILABLE 项：`source`（**全部 = SKEMPI v2.0** → leave-one-source-out 不可做）
  - `assay_top`：11 类，可用

---

## 2. Main Phenomenon

| 组 | n_target | n_rows | true_mean ΔΔG | pred_mean ΔΔG | signed_mean_error |
|---|---:|---:|---:|---:|---:|
| 2–9 | 85 | 466 | +0.6351 | +0.3664 | −0.2686 |
| 10–29 | 55 | 917 | +1.0786 | +0.4501 | −0.6284 |
| **≥30** | **63** | **6112** | **+1.2619** | **+0.5877** | **−0.6742** |

- **`≥30` share of rows = 6112/7523 = 81.2%**（占 target 数 63/230 = 27.4%）
- **`≥30` 的 target-mean Pearson = −0.0008**（V6.2.1）→ 近乎**零**信息
- **机制不是「误差更大」而是「预测被压平」**：

| 组 | slope(pred_mean ~ true_mean) | Pearson | pred_std | true_std |
|---|---:|---:|---:|---:|
| 2–9 | 0.2486 | +0.378 | 0.944 | 1.435 |
| 10–29 | 0.3522 | **+0.747** | 0.628 | 1.332 |
| **≥30** | **−0.0013** | **−0.003** | **0.621** | 1.480 |

⇒ `≥30` 组的预测**几乎是常数**（≈0.59），完全丧失靶点间区分度
（**attenuation / regression-to-the-mean**，而非过冲）。

---

## 3. Target-level Variables

### Top associations with `abs_mean_error`（primary outcome，N=203）

> ⚠️ **分类修正（本审计过程中发现并更正）**：PHASE 1 字典最初把
> `signed_mean_error_M0` / `pred_mean_ddg` / `pred_mean_M0` 标为 `DEPLOYABLE`，
> 但它们其实是**模型输出派生量**（derived），**不是可部署输入描述符**。
> 下表已剔除它们；真正的可部署描述符关联**全部很弱**。

**DEPLOYABLE（真正的输入描述符）**：

| 特征 | ρ | p | BH q |
|---|---:|---:|---:|
| `struct_chainpair_contact_density` | **+0.2230** | 0.00138 | **0.0146** ✅ |
| `struct_contacts_per_res_total` | +0.1845 | 0.00842 | 0.0723 |
| `struct_n_contact_pairs` | +0.1817 | 0.00946 | 0.0750 |
| `struct_n_contact_res_pairs` | +0.1757 | 0.0121 | 0.0834 |
| `mutations_per_site` | +0.1691 | 0.0159 | 0.102 |
| `struct_n_iface_res_lig` | 见 CSV | — | — |
| `site_coverage_interface` | 见 CSV | — | — |
| `n_unique_mutation_sites` | 见 CSV | — | — |
| **`log_n_mutations`** | **−0.0268** | **0.705** | — **不显著** |

⇒ **BH-FDR 后唯一显著的可部署关联是 `struct_chainpair_contact_density`
（ρ = +0.223, q = 0.0146）** —— 即**界面接触密度**，属结构上下文，
**不是** library size / composition。

**LABEL_DERIVED（仅解释）**：
1. **`true_mean_ddg`：ρ = +0.4130, p = 9.1e-10** ✅
2. **`ddg_std`：ρ = +0.3064, p = 8.7e-06** ✅
3. `ddg_IQR`：+0.2668, p = 1.2e-04 ✅
4. `ddg_range`：+0.2319, p = 8.7e-04 ✅
5. `deleterious_fraction`：+0.2226, p = 1.4e-03 ✅
6. `beneficial_fraction`：−0.1890, p = 6.9e-03

> **关键**：`|error|` 的**最强**关联全部是 **label-derived 的效应量/水平**变量
> （ρ = 0.41 / 0.31）；所有**可部署**变量的关联都 ≤ 0.223。
> 这本身就是「可部署信号不足」的直接证据。

---

## 4. Is Mutation Count Causal-Looking or a Proxy?

- **Spearman(log n_mut, |mean error|) = −0.0268（p = 0.705, n = 203）** → **不显著**
- **Partial association**（残差化）：

| 控制集 | 普通 ρ | partial ρ | p |
|---|---:|---:|---:|
| none | −0.0268 | −0.0268 | — |
| n_unique_sites | −0.0268 | −0.0132 | 0.852 |
| diversity | −0.0268 | +0.1423 | 0.043 |
| diversity+alanine | −0.0268 | +0.1479 | 0.035 |
| full_B | −0.0268 | +0.1520 | 0.030 |
| full_C | −0.0268 | +0.1914 | 0.006 |

- **Multivariable coefficient**（family-CV）：Model A（size only）CV R² = **−0.0226**
- **family-CV stability**: 所有可部署模型 CV R² 均 **≤ 0**（A −0.023 / B −0.021 / C −0.078）
  仅 label-derived 的 Model D 为正（+0.0655）

### Decision: **`INCONCLUSIVE`（但实质上 mutation count 不是解释变量）**

**理由**：
1. `log_n_mutations` 与 `|mean error|` 的**单变量关联本身就不显著**（ρ = −0.0268）
   → 它**不可能**是 calibration failure 的解释变量
2. V6.2.1 报告的 L4 Spearman −0.0991（p=0.134）亦不显著
   （**注**：其最终汇总里写的 −0.3592 实为 **L1** 的值，属**文本串档**，已在
   `phase2_univariate/` 中核对更正）
3. 真正与 error 相关的是 **`true_mean_ddg`（ρ=+0.413）** 与 **`ddg_std`（ρ=+0.306）**
   —— 即 **effect-size / 水平**，不是 library size

> ⚠️ 按 V6.3 §4：**只能说 association**。本项目**没有** causal design。

---

## 5. What Explains the Failure?

| 候选解释 | 证据 | 结论 |
|---|---|---|
| **unique sites** | partial ρ 控制后 +0.142（p=0.043） | 弱 |
| **site coverage** | 单变量弱 | 弱 |
| **mutation entropy** | 弱 | 弱 |
| **alanine fraction** | 弱 | 弱 |
| **single/multi** | 弱 | 弱 |
| **interface 接触密度** | **`struct_chainpair_contact_density` ρ=+0.223, q=0.0146** | **唯一 BH 显著的可部署项** |
| **interface size（绝对量）** | family-CV 下 target-only 模型 Pearson −0.044 | 无预测力 |
| **family** | leave-one-family-out 后现象**不消失** | **NOT_SINGLE_FAMILY_DRIVEN** |
| **source** | **UNAVAILABLE**（单一 source） | 无法检验 |
| **label variance（`ddg_std`）** | ρ = **+0.306**（p=8.7e-06） | **最强可解释信号** |
| **target mean level（`true_mean_ddg`）** | ρ = **+0.413**（p=9.1e-10） | **最强** |

**核心机制（本报告的最可靠发现）**：
> calibration failure **不是**「大 library 难预测」，而是
> **`≥30` 靶点的 target-mean 被模型压平为常数**（slope(pred~true) = −0.0013，
> Pearson = −0.003；而 10–29 组 slope = 0.352、Pearson = +0.747），
> 且该现象**与靶点真值均值的正偏移（`true_mean_ddg` ↑）和效应量（`ddg_std` ↑）同向**。
> **「row-weighted `C_B` ≈ 0/负」则主要是 weighting 产物**（见 §7）。

---

## 6. Influence

| 移除 | 剩余 targets | macro mean-Pearson | **row-w C_B** | overall r |
|---|---:|---:|---:|---:|
| (基线) | 230 | +0.5568 | **−0.0121** | — |
| top 1 | 229 | +0.5595 | +0.0084 | +0.1643 |
| top 3 | 227 | +0.5596 | +0.0131 | +0.1424 |
| top 5 | 225 | +0.5730 | **+0.0990** | +0.2684 |
| top 10 | 220 | +0.5752 | +0.1234 | +0.2397 |
| top 20 | 210 | +0.5853 | **+0.1377** | +0.2640 |

- **Leave largest family**: 见 `phase4_family_source/influence_analysis.csv`
  （top family 移除后 macro mean-Pearson 仍在 +0.55 量级）
- **Single-family driven: NO** → **`NOT_SINGLE_FAMILY_DRIVEN`**
- **关键观察**：移除大靶点**改善** row-w C_B（−0.012 → +0.138）
  → 大靶点是**拖累** row-weighted 统计的，不是支撑它的

---

## 7. Weighting Sensitivity

| weighting | L1 C_B | L4 C_B | L4 macro mean-Pearson |
|---|---:|---:|---:|
| **row (n_c/N)** | +0.4323 | **−0.0121** | +0.5568 |
| **equal-target** | +0.7752 | **+0.5475** | +0.5568 |
| sqrt(n) | +0.4498 | +0.1645 | +0.5568 |
| capped 30 | +0.3427 | +0.1081 | +0.5568 |
| log(1+n) | +0.4487 | +0.2183 | +0.5568 |

- **L4 `C_B` 极差 = 0.5596**，且**符号翻转**（row 为负、其余为正）
- **macro mean-Pearson 在所有 weighting 下完全相同（+0.5568）**

> ### ★ 这是 V6.3 最重要的方法学发现
> **「row-weighted `C_B ≈ 0/负`」是 weighting 的产物，不是模型能力的性质。**
> 同一模型、同一预测，仅换权重即可让 `C_B` 从 **−0.012** 变为 **+0.548**。
> 因此 **V6.2/V6.2.1 中「between-target calibration 失效」的结论应改为**：
> **「在 row-weighting 下 between 贡献被 mutation-rich 靶点稀释至 0」**。

---

## 8. Candidate-Set Calibration

**Task**：`y_target = true_mean_ddg`；family-level CV（GroupKFold by family）

| 模型 | 特征数 | n | Pearson | Spearman | RMSE | MAE |
|---|---:|---:|---:|---:|---:|---:|
| **T** target-only | 23 | 203 | **−0.0438** | −0.0470 | 1.4697 | 0.9723 |
| **TS** target+set | 43 | 203 | **+0.2461** | +0.2083 | 1.4375 | 1.0294 |

**`≥30` targets**：
- target-only：Pearson **+0.0878**，RMSE 1.5061
- target+set：Pearson **+0.2614**，RMSE 1.4446

**Bootstrap difference**：
- **Delta_Pearson (TS−T) = +0.2877, 95% CI [+0.1237, +0.4504]** ✅ 排除 0
- Delta_RMSE (TS−T) = −0.0303, 95% CI [−0.1243, +0.0672] ❌ 含 0
- `≥30` Delta_RMSE = −0.0557, 95% CI [−0.2547, +0.1467] ❌ 含 0

**SET_CONTEXT_SUPPORTED: YES（限于 Pearson/Spearman）**

> ⚠️ **但必须同时报告两条限制**：
> 1. **target-only 模型 Pearson = −0.044 → 完全无预测力**；
>    即「只用靶点结构/上下文」**根本预测不了** target mean ΔΔG。
> 2. **绝对误差未显著改善**（ΔRMSE CI 含 0），且 TS 的 MAE 反而略差
>    （1.0294 vs 0.9723）。
>
> ⇒ 正确表述：**mutation-set 组成对 target mean 提供了统计显著的正向信息
> （ΔPearson +0.29），但绝对标定精度依然很低（|r| ≈ 0.25, R² ≈ 6%）。**

---

## 9. Mechanism Decision

**MECHANISM = `LIBRARY_COMPOSITION_DRIVEN`**

### 逐条判据（对照 V6.3 §PHASE 7）

| 判据 | 实测 | 满足? |
|---|---|---|
| `n_mut` 单变量显著 | ρ = **−0.0268, p=0.705** | **❌** |
| 控制 mutation-set diversity 后作用明显减弱 | 单变量本就不显著；控制后转正（+0.14~+0.19） | — |
| **set descriptors 显著改善 calibration** | **ΔPearson +0.288 [+0.124, +0.450] ✅** | **✅** |
| target/global structural properties 解释 error | **family-CV Pearson −0.044（无）** | **❌** |
| target/global 有作用且 set 仍有独立增益 | target 侧无作用 | ❌ |
| 现象由少数 family/source 强烈驱动 | **NOT_SINGLE_FAMILY_DRIVEN** | ❌ |

**结论**：唯一成立的机制分支是 **`LIBRARY_COMPOSITION_DRIVEN`**。

**但必须精确限定其含义**（V6.3 的关键澄清）：
> 这里说的「composition」**不是** mutation count，
> 而是 **mutation-set 的组成描述符**（unique sites / site entropy /
> mutation-pair entropy / coverage / residue composition）。
> **`n_mutations` 本身既不是原因，也不是有效的代理** ——
> 它与 calibration error 的关联**本身就不显著**。

---

## 10. V7 Architecture Consequence

**Recommended global branch: `G(x_c, M_c)`**

理由：
1. **`G(x_c)`（target-only）在 family-CV 下 Pearson = −0.044 → 无预测力**
   → 单靠靶点结构/上下文不可行
2. 加入 `M_c`（mutation-set 组成）后 Pearson 提升到 **+0.246**
   （ΔPearson **+0.288**，95% CI 排除 0）
3. 因此 **`M_c` 必须进入 global branch**

**Recommended local branch: `R_local(x_c, m_i)`**（不变）

> 即 V7 global branch 为 **`G(x_c, M_c)`**，**不是** `G(x_c)`。
> 且**不得**使用 target-ID intercept / family lookup（会破坏 family-OOD 评估）。

---

## 11. Paper Implications

**Most defensible finding**:
> 在 strict family-level OOD 下，PPI ΔΔG 模型的 target-level 校准失败
> **并非由 mutation 数量驱动**（ρ = −0.027, p = 0.71），
> 而是表现为 mutation-rich 靶点的预测**被压平为常数**（slope ≈ 0，
> 而 10–29 组的 slope = 0.35、Pearson = 0.75）；
> 且「between-target 贡献归零」这一现象**主要是 row-weighting 的产物**
> —— 换成等靶点权重后同一预测的 `C_B` 从 −0.012 变为 +0.548。

**Potential novelty**:
> 提出并量化了 **weighting-induced 伪失败（weighting artifact）**：
> 同一组预测在 5 种 target-size weighting 下 `C_B` 跨越 **0.56** 并**翻转符号**，
> 而 macro 指标完全不变 —— 说明汇报 ΔΔG benchmark 的 between-target 能力时，
> **必须固定并声明加权口径**。

**Claims NOT supported**:
- ❌ mutation count causes failure（ρ 不显著；无 causal design）
- ❌ all large targets fail because of size（10–29 组表现良好）
- ❌ set calibration is necessary（ΔRMSE CI 含 0；绝对精度仍低）
- ❌ between-target prediction completely fails（macro mean-Pearson = **+0.557**）
- ❌ target structural context can predict target mean（Pearson = −0.044）

---

## 12. Next Step

**V7: `G(x_c, M_c) + R_local(x_c, m_i)`**

具体约束：
1. **global branch 必须显式输入 mutation-set 组成 `M_c`**（deployable 描述符）
2. **禁止** target-ID intercept / family lookup
3. **评估必须固定 weighting 口径**，并**同时**报告
   row-weighted 与 equal-target 两套 `C_B`
4. **必须报告 attenuation slope**（`slope(pred_mean ~ true_mean)`）作为校准诊断，
   因为它比 `abs_mean_error` 更能定位失效模式
5. **`≥30` 组必须单独报告**（占 81% 的行）
6. 优先级：**先修 target-level 压缩（attenuation）**，再谈 within-target 提升
   —— 因为 within `C_W` 在 L4 仍保留 **+0.147**，而 between 在 row-weighting 下为 0

---

## 13. 三问的直接回答（V6.3 §5）

### Q1：`n_mutations` 到底是不是 proxy？
> **不是有效的 proxy，也不是原因。**
> 它与 `|mean error|` 的 Spearman = **−0.0268（p = 0.705）**，**不显著**；
> 在 family-CV 下 size-only 模型 R² = **−0.023**。
> 与 error 真正相关的是 **`true_mean_ddg`（ρ=+0.413）** 与 **`ddg_std`（ρ=+0.306）**。
> （V6.2.1 汇总中「L4 Spearman = −0.3592」系 **L1 数值串档**，已更正：
> L4 实际为 −0.0991，p=0.134。）

### Q2：真正解释 calibration failure 的变量是什么？
> **两个层次**：
> 1. **现象层**：靶点真值均值水平（`true_mean_ddg`, ρ=+0.413）与效应量
>    （`ddg_std`, ρ=+0.306）—— 均 **LABEL_DERIVED**，只能解释不能部署
> 2. **可部署层**：**mutation-set 组成**（ΔPearson +0.288, CI 排除 0）
> 而 **target-only 结构上下文完全无效**（Pearson −0.044）
> **加权口径**是第三个、也是方法学上最关键的因素（`C_B` 极差 0.56 且翻转符号）

### Q3：V7 global branch 应该是 `G(x_c)` 还是 `G(x_c, M_c)`？
> **`G(x_c, M_c)`**
> 依据：`G(x_c)` 在 family-CV 下 Pearson = **−0.044**（无预测力）；
> 加入 `M_c` 后 **+0.246**（ΔPearson **+0.288**，95% CI **[+0.124, +0.450]** 排除 0）。

---

## 14. SUPPORTED / NOT_SUPPORTED / UNRESOLVED

### SUPPORTED
1. **Weighting artifact**：同一预测在 5 种 target-size weighting 下 `C_B` 极差
   **0.5596** 且**符号翻转**（row −0.0121 ↔ equal +0.5475），macro 指标不变
2. **Attenuation mechanism**：`≥30` 组 `slope(pred_mean~true_mean)` = **−0.0013**、
   Pearson **−0.003**；而 10–29 组 slope 0.3522、Pearson **+0.747**
3. **`n_mutations` 与 calibration error 无显著关联**（ρ = −0.0268, p = 0.705）
4. **target-only 描述符无法预测 target mean**（family-CV Pearson **−0.044**）
5. **set 组成提供显著正向信息**（ΔPearson **+0.288**, CI [+0.124, +0.450]）
6. **NOT_SINGLE_FAMILY_DRIVEN**：leave-one-family-out 后宏观现象保持
   （macro mean-Pearson 仍在 +0.55 量级）
7. **大靶点是拖累而非支撑**：移除 top-20 后 row-w `C_B` 从 −0.012 升至 +0.138

### NOT_SUPPORTED
1. ❌ 「mutation count 导致 calibration failure」
2. ❌ 「所有大靶点因 size 而失败」（10–29 组 Pearson = +0.747）
3. ❌ 「between-target prediction completely fails」（macro mean-Pearson = **+0.557**）
4. ❌ 「set calibration 使绝对精度显著改善」（ΔRMSE CI 含 0；MAE 反而略差）
5. ❌ 「target structural context 可预测 target mean」（Pearson = −0.044）

### UNRESOLVED
1. **`source` / paper 维度**：全部数据来自 **SKEMPI v2.0**（单一 source）
   → leave-one-source-out **无法执行**（标 `UNAVAILABLE`）
2. **绝对精度仍低**：即使 TS，family-CV Pearson 仅 **+0.246**（R² ≈ 6%）
   → target-level mean 在 family OOD 下的可预测性**仍然很差**
3. **因果方向**：`true_mean_ddg` ↑ 与 error ↑ 的共变无法区分
   「模型压平」与「数据本身在该区间更难」——需 **causal design**（本轮无）

---

## 15. 文件清单

### 脚本（`scripts/v6_3/`）
`01_build_target_table.py`、`02_univariate_and_multivariable.py`、
`03_influence_and_weighting.py`、`04_set_calibration.py`

### 结果（`experiments/v6_3/`）
- `PROGRESS.md`
- `phase0_provenance/`：`input_manifest.csv`、`provenance.md`
- `phase1_target_table/`：`target_mechanism_table.csv`（203×111）、
  `feature_dictionary.csv`、`missingness_report.md`
- `phase2_univariate/`：`univariate_associations.csv`、`group_tests.csv`、
  `top_associations.md`
- `phase3_multivariable/`：`model_comparison.csv`、`coefficients.csv`、
  `cv_results.csv`、`partial_associations.csv`、`multivariable_report.md`
- `phase4_family_source/`：`family_concentration.csv`、
  `source_concentration.csv`、`influence_analysis.csv`、`family_source_report.md`
- `phase5_influence/`：`weighting_sensitivity.csv`、`weighting_report.md`
- `phase6_set_calibration/`：`target_only_vs_set.csv`、`predictions.csv`、
  `bootstrap.csv`、`set_calibration_report.md`
- `phase7_mechanism/`
- `V6_3_FINAL_SUMMARY.md`（本文件）

### 日志（`logs/v6_3/`）
无 `failures.log` → **无功能性失败**
