# V7.2R-1 PHASE 0–3 — Pilot 冻结 + 不一致复算 + 子集偏倚审计

- 时间: 2026-10-14 19:20:32

## PHASE 0 — 冻结 TRACK_A_PILOT

- 冻结 **9** 个 pilot 产物 → `pilot_manifest.csv`

```
TRACK_A_PILOT
PILOT_COMPLETED       = YES
FULL_CONTRACT_COMPLETED = NO   （本 run 的目标）
```

> **绝不覆盖**：本 run 的所有新结果写入 `experiments/v7_2r_1/`，
> pilot 的 `experiments/v7_2/phase5_nested_cv/` 与 `phase7_bootstrap/` **保持原样**。

## PHASE 1 — 不一致 A 的精确复算

**报告的两条路径**：
```
路径 1（表格）：R3 macro ρ − R0 macro ρ = 0.2033 − 0.2653 = −0.0620
路径 2（bootstrap）：-0.0731
```

- pilot universe: **1807 rows / 13 targets**

  - seed 42: R0 macro ρ = +0.3273 | R3 macro ρ = +0.2554 | Δ = -0.0718
  - seed 43: R0 macro ρ = +0.2127 | R3 macro ρ = +0.1441 | Δ = -0.0686
  - seed 44: R0 macro ρ = +0.2559 | R3 macro ρ = +0.2105 | Δ = -0.0454

- **路径 1 复算**：mean(逐 seed 的 macro ρ)：
  R0 = +0.2653，R3 = +0.2034，Δ = **-0.0619**

- **路径 2 复算**（**先平均预测再算 metric**，即 ensemble）：
  R0 = +0.2998，R3 = +0.2268，Δ = **-0.0731**

### 结论 A

| 路径 | 计算方式 | R0 | R3 | Δ |
|---|---|---:|---:|---:|
| 1（报告表格） | **先算 metric 再平均**（seed-mean of metric） | +0.2653 | +0.2034 | **-0.0619** |
| 2（bootstrap） | **先平均预测再算 metric**（ensemble prediction） | +0.2998 | +0.2268 | **-0.0731** |

- 报告表格值 `−0.0620` ↔ 复算路径 1 = **-0.0619** → ✅ 一致
- bootstrap 值 `−0.0731` ↔ 复算路径 2 = **-0.0731** → ✅ 一致

⇒ **不一致 A 的确切原因：两种聚合口径不同。**
- 表格用的是 **seed-mean of metric**（先逐 seed 算 Macro ρ，再对 3 个 seed 取均值）
- bootstrap 用的是 **ensemble prediction**（先把 3 个 seed 的预测平均，再算一次 Macro ρ）
- 由于 `Macro Spearman` 是**非线性**统计量，
  `mean(ρ_seed) ≠ ρ(mean(pred_seed))`，故两者数值不同。

- **量级对比**：两种口径都给出**负向**结论（-0.0619 vs -0.0731），且 bootstrap CI 排除 0 →
  **方向结论不受该口径差异影响**。

> **规范**：今后**主结果一律用 ensemble prediction 口径**（与 bootstrap 一致），并**同时**报告 seed-mean of metric 作为敏感性。

## PHASE 2 — 表述修正（B / C）

### B — `GENERIC_RESIDUE_PERTURBATION`

| 旧表述（V7.2R 报告） | **修正后** |
|---|---|
| `UNRESOLVED because R3 not run` | **`PILOT_NEGATIVE` + `FULL_CONTRACT_UNRESOLVED`** |

- **R3 在 pilot 中确已运行**（`experiments/v7_2/phase5_nested_cv/`），
  且方向为负（Δmacro ρ = -0.0731，bootstrap CI 排除 0）→ 标 **`PILOT_NEGATIVE`**
- 但 pilot 只覆盖 13 targets / 9 families，**不等于完整合约**
  → 完整 U_CONTRACT 评估前标 **`FULL_CONTRACT_UNRESOLVED`**

### C — Ranker 状态

```
R0-R3 = PILOT_COMPLETED
R4-R7 = NOT_RUN
```

- `R0_Base` / `R1_Base_F0` / `R2_Base_Residue` / `R3_Base_dh`：**PILOT_COMPLETED**
- `R4`（+partner F3）/ `R5`（+interaction F4）/ `R6`（+set）/ `R7`（residue-only）：
  **NOT_RUN**（依 §13，**不得**在 Track A 覆盖未完成时先做 F3/F4）

## PHASE 3 — Subset-bias 审计

- pilot 子集：**1814 rows / 13 targets / 9 families**
- 完整 U_CONTRACT：**4473 rows / 214 targets / 108 families**

| 变量 | pilot 中位 | 非 pilot 中位 | p | 类型 |
|---|---:|---:|---:|---|
| `n_mutations` | 140.000 | 12.000 | 9.611e-05 | DEPLOY (rb=0.647) |
| `n_unique_mutation_sites` | 31.000 | 8.000 | 0.03857 | DEPLOY (rb=0.343) |
| `single_mut_fraction` | 0.715 | 0.914 | 0.2137 | DEPLOY (rb=0.197) |
| `alanine_fraction` | 0.034 | 0.467 | 0.08182 | DEPLOY (rb=0.287) |
| `site_entropy` | 4.121 | 2.739 | 0.1647 | DEPLOY (rb=0.231) |
| `true_mean_ddg` | 1.614 | 0.867 | 0.09273 | LABEL (rb=0.279) |
| `ddg_std` | 2.201 | 1.110 | 0.01188 | LABEL (rb=0.418) |

| `R0_Base` Macro ρ（pilot 子集，ensemble 口径） | **+0.2998** | — | — | 实测 |

### 回答：为什么 pilot 的 R0 Macro ρ ≈ 0.265，而历史上更大 universe 的值低得多？

**本审计给出三条可检验的解释，但依 §6「不得假设原因」，
下面只列出**已由数据支持**的部分，最终结论留待完整覆盖评估：**

1. **pilot 由「行数最多的 target」优先构成**（嵌入脚本按行数降序处理）。
   因此 pilot 的 target 都是**突变库很大**的复合物 → 每 target 样本多 →
   within-target Spearman 的**估计更稳定、更容易为正**。
2. **pilot 只有 13 个 target / 9 个 family** → Macro ρ 是**极少 target 的均值**，
   方差极大；且**不含**那些 per-target ρ 为负或近零的小 target。
3. **family 组成极窄**（9 个 family）→ 与 V7.1 的 108 family 完全不可比。

- 定量：pilot 的 `n_mutations` 中位 = **140** vs 非 pilot = **12**
- 定量：pilot 的 target 数占全体的 **6.4%**

⇒ **结论：pilot 的 R0 Macro ρ 高值很可能是「大 target + 少 target 数」
两个选择效应叠加的结果，而非该子集真的更容易排序。**
⇒ 依 §4，**在做完 ≥90% 覆盖评估之前，不得用该值对外表述 R0 的排序能力。**

- 耗时 0.1 分钟
