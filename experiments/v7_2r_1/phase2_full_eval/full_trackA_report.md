# V7.2R-1 — 完整 Track A 评估（冻结嵌套 family-CV）

- 时间: 2026-10-14 19:45:14
- learner: **LightGBM**（trees=300）；PCA=32（**train-only**）
- 协议: V7.1.1 冻结 `OUTER_FOLDS_V7_2` + inner family CV
- 主结果口径: **ensemble prediction**；敏感性: seed-mean of metric

## PHASE 1 — Track A universe manifest（**冻结**）

| 项 | 值 |
|---|---|
| U_CONTRACT rows | 4473 |
| **Track A covered rows** | **4233** |
| **coverage** | **94.6%** |
| targets | **174** / 214 |
| families | **91** / 108 |
| n>=5 targets（Macro ρ 资格） | **116** |

- **未覆盖 240 行**，涉及 142 个 target。
- **排除原因**：这些复合物在 `mutation_sites` 中**没有任何位点**同时满足 `chain ∈ PDB`、`seq_index0` 在链长内、且 `mapped_WT == expected_WT`（§19：绝不接受 WT mismatch）。
  即覆盖率受**WS 校验约束**限制，**未因提高覆盖率而放宽校验**。

- 特征块维度：base=188, F0=53, F1=640, F2=322
- **universe: rows 4233 / targets 174 / families 91**

## PHASE 2 — 完整 R0–R3 嵌套 family-CV

| R0_Base | s42 | Macro ρ=+0.0504 | Eff@10=+0.3430 | C_W=+0.1027 | r=+0.1451 |
| R1_Base_F0 | s42 | Macro ρ=+0.0222 | Eff@10=+0.4021 | C_W=+0.0844 | r=+0.1259 |
| R2_Base_Residue | s42 | Macro ρ=+0.0322 | Eff@10=+0.3068 | C_W=+0.0924 | r=+0.1534 |
| R3_Base_dh | s42 | Macro ρ=+0.0220 | Eff@10=+0.4180 | C_W=+0.0900 | r=+0.1417 |
| R0_Base | s43 | Macro ρ=+0.0143 | Eff@10=+0.4375 | C_W=+0.0952 | r=+0.1811 |
| R1_Base_F0 | s43 | Macro ρ=-0.0504 | Eff@10=+0.2547 | C_W=+0.0720 | r=+0.1497 |
| R2_Base_Residue | s43 | Macro ρ=-0.0184 | Eff@10=+0.3280 | C_W=+0.0854 | r=+0.1739 |
| R3_Base_dh | s43 | Macro ρ=-0.0239 | Eff@10=+0.4082 | C_W=+0.0817 | r=+0.1550 |
| R0_Base | s44 | Macro ρ=+0.0501 | Eff@10=+0.2296 | C_W=+0.2024 | r=+0.2615 |
| R1_Base_F0 | s44 | Macro ρ=+0.0472 | Eff@10=+0.2168 | C_W=+0.1938 | r=+0.2494 |
| R2_Base_Residue | s44 | Macro ρ=+0.0572 | Eff@10=+0.3697 | C_W=+0.2125 | r=+0.2766 |
| R3_Base_dh | s44 | Macro ρ=+0.0465 | Eff@10=+0.2066 | C_W=+0.2029 | r=+0.2474 |

## 2.1 两种聚合口径（**必须同时报告**）

| Ranker | 口径 | Macro ρ | Eff@10 | C_W | Overall r |
|---|---|---:|---:|---:|---:|
| R0_Base | seed-mean of metric | +0.0383 | +0.3367 | +0.1334 | +0.1959 |
| R0_Base | **ensemble prediction** | +0.0259 | +0.1882 | +0.1484 | +0.2159 |
| R1_Base_F0 | seed-mean of metric | +0.0064 | +0.2912 | +0.1167 | +0.1750 |
| R1_Base_F0 | **ensemble prediction** | -0.0047 | +0.2843 | +0.1314 | +0.1954 |
| R2_Base_Residue | seed-mean of metric | +0.0237 | +0.3348 | +0.1301 | +0.2013 |
| R2_Base_Residue | **ensemble prediction** | +0.0177 | +0.3425 | +0.1479 | +0.2260 |
| R3_Base_dh | seed-mean of metric | +0.0149 | +0.3443 | +0.1249 | +0.1814 |
| R3_Base_dh | **ensemble prediction** | +0.0116 | +0.4176 | +0.1410 | +0.2025 |

## PHASE 3 — 5000/5000 paired target-cluster bootstrap


### R3_Base_dh − R0_Base

| 指标 | point | 95% CI | 有效 | 排除 0? |
|---|---:|---|---:|---|
| Δmacro_spearman | -0.0143 | [-0.0588, +0.0307] | 5000 | ❌ 含 0 |
| Δeff10 | +0.2293 | [-0.2989, +0.2570] | 5000 | ❌ 含 0 |
| ΔC_W | -0.0074 | [-0.0203, +0.0051] | 5000 | ❌ 含 0 |
| Δoverall_r | -0.0134 | [-0.0320, +0.0096] | 5000 | ❌ 含 0 |
| Δrmse | -0.0096 | [-0.0359, +0.0145] | 5000 | ❌ 含 0 |

- **§23 成功判据（冻结）：NOT SUPPORTED ❌**


### R2_Base_Residue − R0_Base

| 指标 | point | 95% CI | 有效 | 排除 0? |
|---|---:|---|---:|---|
| Δmacro_spearman | -0.0083 | [-0.0478, +0.0303] | 5000 | ❌ 含 0 |
| Δeff10 | +0.1542 | [-0.4554, +0.4913] | 5000 | ❌ 含 0 |
| ΔC_W | -0.0005 | [-0.0185, +0.0173] | 5000 | ❌ 含 0 |
| Δoverall_r | +0.0102 | [-0.0519, +0.0634] | 5000 | ❌ 含 0 |
| Δrmse | -0.0737 | [-0.1382, +0.0045] | 5000 | ❌ 含 0 |

- **§23 成功判据（冻结）：NOT SUPPORTED ❌**


### R1_Base_F0 − R0_Base

| 指标 | point | 95% CI | 有效 | 排除 0? |
|---|---:|---|---:|---|
| Δmacro_spearman | -0.0306 | [-0.0703, +0.0088] | 5000 | ❌ 含 0 |
| Δeff10 | +0.0960 | [-0.5262, +0.4136] | 5000 | ❌ 含 0 |
| ΔC_W | -0.0170 | [-0.0277, -0.0072] | 5000 | ✅ |
| Δoverall_r | -0.0205 | [-0.0484, +0.0045] | 5000 | ❌ 含 0 |
| Δrmse | -0.0293 | [-0.0680, +0.0169] | 5000 | ❌ 含 0 |

- **§23 成功判据（冻结）：NOT SUPPORTED ❌**

- 耗时 41.9 分钟
