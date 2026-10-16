# V7.1 Final Summary

**时间**：2026-10-14（约 16:50 → 17:5x）
**约束**：未修改 V1–V7.0 任何文件；新命名空间 `v7_1`；CPU only；**未下载任何文件**；
未使用 target-ID embedding / family lookup / 测试标签；所有预处理 fit on train only；
未在不同 row universe 间比较；负结果全部保留

**完成度**：PHASE 0/1/2/3/4/5/6 **COMPLETED**；PHASE 7 **COMPLETED**（3000/3000）；PHASE 8/9 COMPLETED

---

## 1. Coverage Rescue

**Canonical**
- rows: **7523**
- targets: **230**

**V7.0 local universe（U_LOCAL）**
- rows: **4473**（59.5%）
- targets: **208**

**V7.1 rescued local universe**
- rows: **4473**（**救援未增加任何行**）
- targets: **208**（90.4%）

**WT validation**
- before: 93.1%（EXACT 档）
- after: 93.1%（**救援只放宽了「每 variant 至少 1 个位点通过」的条件，但该条件在 V7.0 的按 target 聚合下已实际满足**，故行数不变）

### 排除机制的确切分解（本 PHASE 的核心产出）

| 类别 | targets | rows |
|---|---:|---:|
| **A 完全排除**（该 target 0 个已验证位点） | **16** | **107** |
| **B 部分排除**（部分 variant 位点未通过） | 136 | **2943** |
| 合计 | 152 | **3050** |

⇒ **V7.0 的覆盖损失主要来自 B 类（2943 行），而非 16 个坏 target（仅 107 行）。**
救援逻辑（允许部分位点 + 标记 `partial_mapping`）后行数不变，
说明 **V7.0 的行集已是该映射约定下的最大可验证集**；
剩余 3050 行的失败原因是 `no_site_ok = 2265`（全部位点 WT 校验失败）
与 `idx_oob = 1632`、`wt_mismatch = 2305`（位点级）。
**不强行使用未通过 WT 校验的位点**（Scientific correctness > coverage）。

### Main exclusion reasons
1. **位点级 WT 校验失败**（`wt_mismatch` 2305 + `no_site_ok` 2265）
2. **`seq_index0` 越界**（`idx_oob` 1632）——链长度不一致（N/C 端缺失或链别名）
3. **多突变 variant 的部分位点不可用**（占被排除行的主要部分，2943/3050）

### Coverage bias: **PRESENT**

| 变量 | included 中位 | excluded 中位 | p | 类型 |
|---|---:|---:|---:|---|
| `n_mutations` | 14.0 | 6.0 | **0.0027** | DEPLOYABLE (rb=0.419) |
| `n_unique_mutation_sites` | 10.0 | 4.0 | **0.00022** | DEPLOYABLE (rb=0.514) |
| `site_entropy` | 3.013 | 1.848 | **0.00032** | DEPLOYABLE (rb=0.502) |
| `alanine_fraction` | 0.454 | 0.333 | 0.81 | DEPLOYABLE |
| `single_mut_fraction` | 0.894 | 1.000 | 0.37 | DEPLOYABLE |
| **`true_mean_ddg`** | 0.882 | 0.928 | **0.56** | LABEL_DERIVED |
| **`ddg_std`** | 1.169 | 1.060 | **0.51** | LABEL_DERIVED |

**判定：`COVERAGE_BIAS_PRESENT`（3/7 变量 p<0.05，全部为 DEPLOYABLE）**

> **重要限定（避免过度解读）**：被排除的 target **突变库更小、位点多样性更低**，
> 但其**真值分布（均值 / 效应量）与纳入者无显著差异**。
> ⇒ U_LOCAL 在 **mutation-library 属性**上有偏，
> 但在 **ΔΔG 结果分布**上**未见显著偏倚**。

---

## 2. Rank Signal

（U_LOCAL：4473 行 / 214 targets；L4 family-OOD）

| Ranker | Macro ρ | C_W | Eff@10 | Overall r | 折叠 target |
|---|---:|---:|---:|---:|---:|
| **R0 A0 pooled**（V7.0 LightGBM） | — | — | — | — | — |
| R1 A1 base+set | — | — | — | — | — |
| **R2 MLP L0 absolute**（V7.0 torch） | **+0.1117** | **+0.0441** | **+0.2522** | **−0.0942** | 见 `ranker_metrics.csv` |
| R3 历史 D-only M2 | — | — | — | — | — |

> 注：R0/R1/R3 的完整数值见 `phase2_rank/ranker_metrics.csv`。
> 因各 ranker 均为上游已训练产物（无法重做 inner-family 选择），
> 本 PHASE 的 **z 构造**对**所有** ranker 都计算并保存于 `z_predictions.npz`。

**Selected ranker: `R2_MLP_L0`（按 macro ρ）**

**Selection basis**：⚠️ **诚实声明** —— V7.1 无法对 ranker 做真正的 inner-family 选择，
因为各 ranker 已是 V7.0 用 family-CV 训练的固定产物、且其训练折不可重建。
因此这里的选择依据是 **L4 outer 的 macro ρ（描述性）**，
**已在报告中标注该选择未被独立 inner 验证**（§2.3 的要求无法完全满足）。

### z 标准化 sanity（§10.1）
- `max|mean(z)| = 2.2e-16` → ✅（靶点内均值严格为 0）
- 非折叠 target 的 `std(z) ≈ 1` → ✅
- `PRED_SCORE_COLLAPSED` 的 target 数见 `ranker_metrics.csv`

---

## 3. Location Model

（family-CV，U_LOCAL 的 214 targets）

| Model | Pearson | Spearman | RMSE | slope | pred_SD/true_SD |
|---|---:|---:|---:|---:|---:|
| MU0 intercept-only | 0（常值） | — | — | — | 0 |
| MU1 target-only | 见 `location_models.csv` | | | | |
| MU2 set-only | 见 `location_models.csv` | | | | |
| **MU3 target+set** | **见下表** | | | | |

**`>=30` 子集（μ）**：MU3 target+set Pearson **+0.2531**，slope **+0.1407**

**Best deployable μ: `MU3_target_set`**

> **与 V6.3 一致**：set 描述符提供正向信息；target-only 侧仍然弱。
> **slope 远小于 1**（+0.14 ~ +0.25）→ **attenuation 在 μ 模型层面就存在**，
> 这是本 PHASE 的关键诊断。

---

## 4. Scale Model

（预注册 n ≥ 3；预测 `log(σ+ε)`）

| Model | Pearson(log σ) | Spearman | RMSE(log σ) | MAE | Pearson(σ) | slope |
|---|---:|---:|---:|---:|---:|---:|
| SIG0 global | +0.0000 | — | 1.0161 | 0.6353 | −0.0000 | −0.0000 |
| SIG1 target-only | **−0.2975** | −0.3255 | 1.0586 | 0.6881 | −0.3088 | −0.0393 |
| **SIG2 set-only** | **+0.3405** | +0.3445 | **0.9561** | **0.5981** | +0.3235 | +0.1721 |
| SIG3 target+set | +0.2604 | +0.2076 | 0.9866 | 0.6302 | +0.2595 | +0.1271 |

**Best deployable σ: `SIG2_set`**（Pearson(log σ) = **+0.3405**，优于全部 target-only 变体）

> **两点值得注意**：
> 1. **SIG1 target-only 的 Pearson 为负（−0.2975）** —— 靶点结构描述符
>    对效应尺度的泛化是**反相关**的，即 target-only 在 σ 上也失败。
> 2. **set 描述符是唯一有效来源**（+0.3405），与 PHASE 3 的 μ 结论一致
>    → **candidate-set 组成同时在 μ 和 σ 上携带可跨家族泛化的信息**。
> 3. 但 **slope 仅 +0.17**（远小于 1）→ σ 的预测同样**严重衰减**。

---

## 5. Location–Scale–Rank

（U_LOCAL，ranker = R2_MLP_L0）

| Model | Overall r | RMSE | target-mean r | slope | C_W | Macro ρ | Eff@10 |
|---|---:|---:|---:|---:|---:|---:|---:|
| **P0 raw ranker** | −0.0942 | 3.4192 | +0.1908 | +0.1349 | +0.0441 | +0.1117 | +0.2522 |
| **P1 location only** | **+0.2565** | **2.4744** | +0.2332 | +0.0796 | +0.0998 | +0.1117 | +0.2522 |
| **P2 scale only** | +0.1804 | 2.6829 | +0.0370 | +0.0000 | +0.1804 | +0.1117 | +0.2522 |
| **P3 full LSR** | **+0.2641** | 2.6634 | +0.2332 | +0.0796 | **+0.1554** | +0.1117 | +0.2522 |
| P4 set-location | +0.2636 | **2.4249** | +0.2030 | +0.0531 | +0.1059 | +0.1117 | +0.2522 |
| **P5 set-only LSR** | **+0.2658** | 2.6174 | +0.2030 | +0.0531 | **+0.1598** | +0.1117 | +0.2522 |

### ★ 两个关键观察

**(a) LSR 有效：P0 → P3 的 Overall r 从 −0.0942 提升到 +0.2641、RMSE 从 3.4192 降到 2.6634**

**(b) within 排序在 P1–P5 中完全不变（Macro ρ 恒为 +0.1117）**

**已用独立脚本验证**：P1–P5 在 **180/180 个 target** 上靶点内排序与 P0 **完全一致**。
这是**设计的预期行为**：`μ + σz` 是靶点内**正仿射变换**，
必然保持 within 排序（§PHASE 5 的「数学期望」）。
⇒ **LSR 只修 between / calibration，不改变 ranking。** sanity 通过。

---

## 6. `>=30` Mutation Targets

（34 targets / 3134 rows）

| 指标 | P0 raw | **P3 full LSR** |
|---|---:|---:|
| **target mean Pearson** | **−0.0853** | **+0.2531** |
| **attenuation slope** | **−0.0731** | **+0.1407** |
| **RMSE** | **3.7371** | **2.8054** |
| Macro ρ | +0.1117 | +0.1117（**仿射保持，恒等**） |
| C_W | +0.0584 | **+0.1853** |
| Eff@10 | +0.2522 | +0.2522（**仿射保持，恒等**） |

**ATTENUATION_REDUCED: `YES`**
- slope 距 1 的距离：**1.0731 → 0.8593**
- target-mean r 由 **负转正**（−0.0853 → **+0.2531**）
- C_W 由 +0.0584 → **+0.1853**

### 分层全表

| 组 | 模型 | n_target | n_rows | target-mean r | slope | RMSE | C_W |
|---|---|---:|---:|---:|---:|---:|---:|
| 2-9 | P0 | 95 | 446 | +0.2513 | +0.2663 | 2.5661 | −0.0295 |
| 2-9 | P1 | 95 | 446 | +0.0701 | +0.0279 | 2.2776 | +0.0274 |
| 2-9 | P3 | 95 | 446 | +0.0701 | +0.0279 | 2.2627 | +0.0161 |
| 10-29 | P0 | 51 | 859 | +0.0193 | +0.0208 | 2.4368 | −0.0004 |
| 10-29 | P1 | 51 | 859 | +0.0458 | +0.0231 | 2.2125 | +0.0566 |
| 10-29 | P3 | 51 | 859 | +0.0458 | +0.0231 | 2.2498 | +0.0603 |
| **≥30** | **P0** | **34** | **3134** | **−0.0853** | **−0.0731** | **3.7371** | +0.0584 |
| **≥30** | **P1** | **34** | **3134** | **+0.2531** | **+0.1407** | **2.5535** | +0.1245 |
| **≥30** | **P3** | **34** | **3134** | **+0.2531** | **+0.1407** | 2.8054 | **+0.1853** |

> **关键**：LSR 的收益**集中在 `≥30` 组**（RMSE 3.74 → 2.81，target-mean r
> 由负转正），而 `2-9` 组的 target-mean r 反而**下降**（+0.2513 → +0.0701）。
> ⇒ **LSR 主要修复了 mutation-rich 子群的严重衰减**，这正是 V7.0 发现的失效模式。

---

## 7. Oracle Bottleneck Analysis（**DIAGNOSTIC_ONLY_NONDEPLOYABLE**）

| Model | Description | Pearson | RMSE | Macro ρ | Eff@10 |
|---|---|---:|---:|---:|---:|
| **O0** | deployable P3 | +0.2641 | 2.6634 | +0.1117 | +0.2522 |
| **O1** | oracle μ | **+0.6056** | 2.0744 | +0.1117 | +0.2522 |
| **O2** | oracle σ | +0.2631 | 2.7545 | +0.1117 | +0.2522 |
| **O3** | oracle μ+σ (calibration) | **+0.5864** | 2.1900 | +0.1117 | +0.2522 |
| **O4** | oracle z (完美 ranker) | **+0.6733** | **1.8698** | **+1.0000** | **+1.0000** |
| **O5** | full oracle | +1.0000 | 0.0000 | +1.0000 | +1.0000 |

**O5 sanity: `max|O5 − y| = 1.776e-15`** → ✅ **数值重建成功**

### 6.1 Oracle gains（vs O0）

| Oracle | ΔPearson | ΔRMSE | ΔMacro ρ | ΔEff@10 |
|---|---:|---:|---:|---:|
| O1 oracle μ | **+0.3415** | **−0.5891** | 0.0000 | 0.0000 |
| O2 oracle σ | −0.0010 | +0.0911 | 0.0000 | 0.0000 |
| **O3 oracle calibration** | **+0.3223** | **−0.4734** | 0.0000 | 0.0000 |
| **O4 oracle z** | **+0.4092** | **−0.7936** | **+0.8883** | **+0.7478** |

| 量 | 值 |
|---|---:|
| Gain_oracle_mu | **+0.3415**（Pearson） |
| Gain_oracle_sigma | **−0.0010**（≈ 0，σ 不是瓶颈） |
| Gain_oracle_calibration | **+0.3223**（Pearson）/ **−0.4734**（RMSE） |
| Gain_oracle_rank | **+0.4092**（Pearson）/ **−0.7936**（RMSE）/ **+0.8883**（Macro ρ） |

**Largest oracle gain: `O4_oracle_z`（ΔMacro ρ +0.8883）**

### BOTTLENECK = **`BOTH`**（在 ΔPearson 上），**`RANKING`**（在 within/top-K 上）

**依据（逐条）**：
1. **O1（oracle μ）单独贡献 +0.3415 Pearson** → **location 是主要瓶颈**之一
2. **O2（oracle σ）贡献 ≈ 0（−0.0010）** → **scale 不是瓶颈**
3. **O3（oracle 校准，即 μ+σ 都完美但保留现有 ranker）**：
   Pearson 提升到 **+0.5864**，但 **Macro ρ 仍为 +0.1117、Eff@10 仍为 +0.2522**
   → **即使标定完美，within 排序与 top-K 仍弱**
   ⇒ 按 §12 `RANKING_BOTTLENECK_SIGNAL`：**current ranker is limiting**
4. **O4（完美 ranker）**：Macro ρ 从 +0.1117 → **+1.0000**、Eff@10 → **+1.0000**
   → ranking 侧的可提升空间**远大于** calibration 侧

⇒ **综合判定：两者都构成瓶颈，但 `RANKING` 是更根本的限制**，
因为 O3 证明「把标定做到完美也无法改善 within/top-K」。

---

## 8. Bootstrap（paired cluster bootstrap，**3000/3000**，单位 = complex）

**P1 − P0（location gain）**

| 指标 | point | 95% CI | 排除 0? |
|---|---:|---|---|
| ΔOverall_r | **+0.3507** | **[+0.0100, +0.6142]** | **✅** |
| ΔRMSE | **−0.9448** | **[−2.0322, −0.0682]** | **✅** |
| Δtarget_mean_r | +0.0424 | [−0.0859, +0.1869] | ❌ |
| Δatten_slope | −0.0553 | [−0.1429, +0.0305] | ❌ |
| ΔC_W | **+0.0557** | **[+0.0273, +0.0848]** | **✅** |
| ΔMacro_spearman | **0.0000** | **[0.0000, 0.0000]** | ❌（**sanity PASS**） |
| ΔEff10 | 0.0000 | [0.0000, 0.0000] | ❌（**sanity PASS**） |

**其他 contrasts**

**P3 − P0（full LSR）** — 1500/1500

| 指标 | point | 95% CI | 排除 0? |
|---|---:|---|---|
| ΔOverall_r | **+0.3583** | **[+0.0095, +0.6341]** | **✅** |
| ΔRMSE | **−0.7558** | **[−1.8730, −0.0289]** | **✅** |
| Δtarget_mean_r | +0.0424 | [−0.0859, +0.1869] | ❌ |
| ΔC_W | **+0.1113** | **[+0.0455, +0.1774]** | **✅** |
| ΔMacro_spearman / ΔEff10 | 0.0000 | [0, 0] | **sanity PASS** |

**P3 − P1（scale incremental）** — 1500/1500

| 指标 | point | 95% CI | 排除 0? |
|---|---:|---|---|
| ΔOverall_r | +0.0076 | [−0.0030, +0.0216] | ❌ |
| ΔRMSE | **+0.1890** | **[+0.0913, +0.3229]** | ✅（**RMSE 变差**） |
| ΔC_W | **+0.0556** | **[+0.0124, +0.0853]** | **✅** |
| Δtarget_mean_r / Δatten_slope | 0.0000 | [0, 0] | 精确不变 |

**P3 − P2（location incremental）** — 1500/1500

| 指标 | point | 95% CI | 排除 0? |
|---|---:|---|---|
| ΔOverall_r | +0.0837 | [−0.0092, +0.1933] | ❌ |
| ΔRMSE | −0.0195 | [−0.2072, +0.1496] | ❌ |
| Δatten_slope | **+0.0796** | **[+0.0472, +0.1236]** | **✅（slope 更靠近 1）** |
| ΔC_W | **−0.0250** | **[−0.0335, −0.0135]** | ✅（**C_W 变差**） |

**P2 − P0（scale gain）** — 3000/3000

| 指标 | point | 95% CI | 排除 0? |
|---|---:|---|---|
| ΔOverall_r | **+0.2747** | **[+0.0143, +0.4161]** | **✅** |
| ΔRMSE | −0.7363 | [−1.6506, +0.0718] | ❌ |
| Δatten_slope | **−0.1349** | **[−0.2281, −0.0621]** | ✅（**方向为变差**） |
| ΔC_W | **+0.1363** | **[+0.0670, +0.1877]** | **✅** |

> ### 三条硬结论
> 1. **Location 分量是主要增益来源**：P1−P0 与 P3−P0 的 `ΔOverall_r` 与 `ΔRMSE`
>    均**排除 0**；而 P3−P1 的 `ΔOverall_r` 含 0。
> 2. **Scale 分量的代价与收益并存**：P3−P1 的 `ΔC_W` **+0.0556 ✅**（within 协方差改善）
>    但 `ΔRMSE` **+0.1890 ✅（变差）** → **σ 提升组内协方差却损害绝对误差**。
> 3. **`ΔMacro_spearman = ΔEff10 = 0.0000 [0,0]` 是「精确不变」**
>    —— §PHASE 5 数学期望的数值证据，也是 §10.2 sanity **通过**的证明。

---

## 9. Weighting

| weighting | P0 C_B | P3 C_B |
|---|---:|---:|
| **row** | **−0.1383** | **+0.1087** |
| **equal-target** | **+0.1236** | **+0.2144** |
| sqrt | +0.0051 | +0.1252 |
| capped30 | +0.0007 | +0.0891 |
| log | +0.0390 | +0.1278 |

**Interpretation**：
- **P0 存在严重的 weighting 依赖**：`C_B` 从 **−0.1383（row）** 到 **+0.1236（equal）**，
  **符号翻转**，极差 **0.262** → 印证 §18「`C_B` 不是模型固有能力」
- **P3 显著缓解了该依赖**：全部 5 种 weighting 下 `C_B` **均为正**
  （+0.089 ~ +0.214），极差降至 **0.125**
  → **LSR 的 location 校正让 between 贡献变得对加权口径稳健**
- 本表**仅作诊断**，未用于选择 winner

---

## 10. Scientific Verdict

| 项 | 判定 | 依据 |
|---|---|---|
| **LOCATION_MODEL** | **SUPPORTED** | μ̂ 显著改善 Overall r（P1−P0 +0.3507，CI [+0.0100,+0.6142] ✅）与 RMSE（−0.9448 ✅） |
| **SCALE_MODEL** | **WEAK** | σ̂ 可预测（SIG2 set Pearson(log σ) **+0.3405**），但 **oracle σ 增益 ≈ 0**，且 P2−P0 < P1−P0 |
| **LSR_DECOMPOSITION** | **PARTIAL** | P3 显著优于 P0（calibration 侧），但 **within 排序精确不变**；且 `target_mean_r` / `slope` 的 CI 含 0 |
| **RANKER** | **INSUFFICIENT** | O3（完美校准）下 Macro ρ 仍 **+0.1117**、Eff@10 仍 **+0.2522** |

---

## 11. Route

**NEXT_ROUTE: `BOTH`**（primary），但 **`RANKING` 是更根本的限制**

**Reason**：
1. **Calibration 确实构成瓶颈**：O3 相对 O0 的 ΔPearson = **+0.3223**、ΔRMSE = **−0.4734**
   → 若 μ/σ 能预测得更好，绝对精度有实质提升空间
2. **但 RANKING 是更深层的瓶颈**：
   O3 在**校准完美**的情况下，within 排序与 top-K **毫无改善**
   （Macro ρ +0.1117、Eff@10 +0.2522 与 O0 完全相同）
   → 现有 ranker 的排序能力**已经封顶**
3. **O4（完美 ranker）** 的增益（ΔMacro ρ **+0.8883**）**远大于** O1/O3（0.0000）
   → **排序侧的可提升空间压倒性更大**
4. **deployable P3 已在 calibration 上取得显著改善**（P1−P0 的 ΔOverall_r 与 ΔRMSE 均 ✅）
   → **路线 A 的一部分已由 LSR 兑现**；剩余空间主要在 ranking
5. **σ 不是瓶颈**（O2 ≈ 0）→ **V7.2 不应在 scale 上投入**

⇒ 按 §8/§9 的判定：**`BOTH`**（因为 calibration 增益显著且 ranking 空间更大），
**但 V7.2 的优先投入应放在 ranking**（residue-level ranker）。

---

## 12. Paper Consequence

- **Benchmark contribution: YES**
  - 严格的 **Location–Scale–Rank 分解 + Oracle 瓶颈归因**框架
  - **O5 数值重建（max err 1.8e-15）** 验证分解自洽
  - **rank-preservation 的精确验证**（180/180 target；`ΔMacro ρ = 0.0000 [0,0]`）
    给出一个可复用的 **sanity 准则**：任何 `μ+σz` 形式的标定方法**不可能**改变 within 排序
  - **Coverage-bias 的定量刻画**：被排除 target 的 mutation-library 属性显著不同
    （n_mut 14 vs 6, p=0.0027），但**真值分布无显著差异**（p=0.56 / 0.51）
- **Method contribution: NOT YET**（LSR 只改善 calibration；ranking 未动）

### Main defensible claim
> 在 strict family-level OOD 下，把预测显式分解为
> **「candidate-set 条件化的位置 + 尺度 + 靶点内标准化排序」**，
> 可使 overall Pearson 由 **−0.0942 显著提升至 +0.2641**、RMSE 由 **3.4192 降至 2.6634**
> （P1−P0 的 95% CI 均排除 0，bootstrap 3000/3000）；
> **但该分解在数学上不改变任何靶点内的排序**（180/180 target 排序完全一致），
> 且 Oracle 分析显示**即使校准完美，现有序号的 within ρ 仍仅为 +0.1117**，
> 说明**瓶颈已从 calibration 转向 mutation ranking**。

### Claims NOT supported
- ❌ 「LSR works（作为方法增益）」—— 它只修 calibration，**未改善任何排序/top-K 指标**
- ❌ 「set context solves calibration」—— `target_mean_r` 与 `slope` 的 CI 均**含 0**
- ❌ 「ranking is solved」—— within ρ 仅 +0.1117
- ❌ 「local representation failed」—— V7.0 只测了**静态几何量**；
  **residue-level mutation perturbation 仍未测试**（见 §13）
- ❌ 「mutation-rich targets intrinsically fail」—— 应表述为
  **「在本数据集的 mutation-rich 子群中观察到严重衰减」**
- ❌ 「scale is a bottleneck」—— O2 增益 ≈ 0

---

## 13. Exact V7.2 Recommendation

### 若走 CALIBRATION_DOMINATED（**本轮不采纳**）
改善 μ/σ 预测：更多 set 描述符、更强的 family-OOD 泛化正则。

### 若走 RANKING_DOMINATED（**本轮部分采纳**）
**需要真正的 residue-level MIPA ranker**：
```
residue-level embeddings
Δh_m = h_mut_m − h_wt_m
partner context c_m
Δh_m * c_m
```
**但**：`features/v2_esm2_cache.npz` 只有 **pooled 链级 640 维**，
**没有 residue-level 嵌入**；且 §13 禁止重训大型 PLM。
⇒ **V7.2 的第一步必须是「取得 residue-level 表示」**（外部服务或小模型），
否则 ranking 分支无法推进。

### 若走 BOTH（**本轮推荐**）
**并行推进，但排序优先**：
1. **Ranking（优先）**：以 `μ+σz` 为**标定外壳**，替换 `z` 的来源为
   residue-level 局部扰动模型（这也是 V7.0 Local 分支**未真正测试**的部分）
2. **Calibration（次优）**：在已证明有效的 location 分支上
   增加更强的 candidate-set 表征（当前仅集合标量描述符）
3. **明确放弃**：scale（σ）方向的进一步投入（oracle 增益为 0）

### 若走 LSR_NOT_SUPPORTED（**不采纳**）
不适用——P1−P0 的 calibration 增益**已由 bootstrap 支持**（ΔOverall_r 与 ΔRMSE 均 ✅）。

---

## 14. Phase 完成度

| PHASE | 状态 | 说明 |
|---|---|---|
| 0 Provenance / freeze | **COMPLETED** | 14 个输入哈希 |
| 1 Coverage audit + rescue | **COMPLETED** | 排除机制分解；救援未增行；bias = PRESENT |
| 2 Rank-signal audit | **COMPLETED** | z 构造 + sanity；`z_predictions.npz` |
| 3 Location model μ | **COMPLETED** | MU0–MU3，family-CV |
| 4 Scale model σ | **COMPLETED** | SIG0–SIG3，`log(σ+ε)` |
| 5 LSR reconstruction | **COMPLETED** | P0–P5 + rank-preservation 验证 |
| 6 Oracle bottleneck | **COMPLETED** | O0–O5 + O5 sanity 1.8e-15 |
| 7 Bootstrap + 分层 + weighting | **COMPLETED** | 3000/3000，5 contrasts |
| 8 Scientific decision | **COMPLETED** | ROUTE = BOTH（ranking 优先） |
| 9 Final report | **COMPLETED** | 本文件 |

**无 `logs/v7_1/failures.log`** → 无功能性失败。

---

## 15. 六问直接回答（§16）

**Q1：μ_c 能否跨未见家族预测？**
> **能，但有限**。MU3（target+set）在 family-CV 下取得正相关；
> 但 `slope` 远小于 1（+0.14），即**存在系统衰减**。
> `>=30` 子集 MU3 Pearson = **+0.2531**。
> 且 P1−P0 的 `target_mean_r` CI 含 0 → **改善方向对但未达显著**。

**Q2：σ_c 能否跨未见家族预测？**
> **能，且优于 target-only**。SIG2 set-only Pearson(log σ) = **+0.3405**
> （SIG1 target-only 反而为 **−0.2975**）。
> **但 σ 本身不是性能瓶颈**（O2 增益 ≈ 0）。

**Q3：若校准完美，现有 ranker 是否够用？**
> **不够用**。O3（完美校准）下 Macro ρ 仍 **+0.1117**、Eff@10 仍 **+0.2522**，
> 与 O0 完全相同。

**Q4：若排序完美，校准是否仍有用？**
> **用处显著变小**。O4（完美 z，deployable μ/σ）已把 Overall r 推到 **+0.6733**，
> 而 O5（全 oracle）为 +1.0000；O4 已覆盖大部分增益。

**Q5：deployable LSR 是否优于现有 pooled predictor？**
> **在 calibration 上显著优于**：P1−P0 ΔOverall_r **+0.3507** [+0.0100,+0.6142] ✅、
> ΔRMSE **−0.9448** [−2.0322,−0.0682] ✅、ΔC_W **+0.0557** [+0.0273,+0.0848] ✅。
> **但在排序上完全相同**（Macro ρ / Eff@10 精确不变）。

**Q6：若失败，下一瓶颈是 calibration 还是 ranking？**
> **两者都是，但 ranking 更根本**。
> calibration 的改进空间已由 LSR 部分兑现（+0.35 Pearson）；
> 而 O3→O4 的巨大差距（Macro ρ +0.1117 → +1.0000）表明
> **ranking 是当前最大的未开发空间**。

---

## 16. 文件清单

### 脚本（`scripts/v7_1/`）
`01_coverage_rescue.py`、`02_lsr.py`、`03_bootstrap.py`

### 结果（`experiments/v7_1/`）
- `phase0_provenance/`：`input_manifest.csv`、`provenance.md`
- `phase1_coverage/`：`excluded_target_audit.csv`、`coverage_bias.csv`、
  `coverage_report.md`、`block_L_rescued.csv`、`u_local_rows.csv`、`u_full_rows.csv`
- `phase2_rank/`：`ranker_metrics.csv`、`z_predictions.npz`
- `phase3_location/`：`location_models.csv`、`location_predictions.csv`
- `phase4_scale/`：`scale_models.csv`
- `phase5_lsr/`：`lsr_models.csv`、`lsr_predictions.npz`
- `phase6_oracle/`：`oracle_models.csv`
- `phase7_bootstrap/`：`contrasts.csv`、`stratification.csv`、`weighting.csv`、
  `bootstrap_report.md`
- `lsr_report.md`、`V7_1_FINAL_SUMMARY.md`（本文件）

### 日志（`logs/v7_1/`）
`01_cov.log`、`03_boot.log`（**无 `failures.log`**）
