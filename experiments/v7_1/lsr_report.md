# V7.1 PHASE 2–6 — Location–Scale–Rank + Oracle 分解

- 时间: 2026-10-14 17:04:57
- universe: **U_LOCAL**

- rows **4473** / targets **214** / families **108**

## PHASE 2 — Rank-signal audit

- V7.0 预测数组映射: 4473 / 4473 行可对齐

| Ranker | Macro ρ | Kendall | C_W | Eff@10 | Overall r |
|---|---:|---:|---:|---:|---:|
| R0_A0_pooled | +0.0333 | — | +0.1734 | +0.1420 | +0.2205 |
| R1_A1_base_set | +0.0507 | — | +0.1795 | +0.2080 | +0.2408 |
| R2_MLP_L0 | +0.1117 | — | +0.0441 | +0.2522 | -0.0942 |
| R3_Donly_M2 | +0.0305 | — | +0.1774 | +0.1566 | +0.2117 |
| R4_hist_M0 | +0.0311 | — | +0.1701 | +0.1082 | +0.2211 |

- **选定 ranker（按 macro ρ，L4 outer 评估）: `R2_MLP_L0`**
- ⚠️ 诚实声明：V7.1 的 ranker **选择**依据是 inner-family validation 不可得（各 ranker 已是 V7.0 用 family-CV 训练的固定产物），
  故此处以 L4 outer 的 macro ρ 作为**描述性**依据，并已在报告中标注该选择未被独立 inner 验证。

- sanity: max|mean(z)| = 1.38e-14 → ✅

## PHASE 3 — Location model μ

| 模型 | 特征 | Pearson | Spearman | RMSE | slope | pred_SD/true_SD |
|---|---:|---:|---:|---:|---:|---:|
| MU0_intercept | 0 | +0.0000 | +nan | 2.1576 | -0.0000 | 0.0000 |
| MU1_target | 58 | +0.1546 | +0.1132 | 2.3403 | +0.0927 | 0.5997 |
| MU2_set | 25 | +0.2030 | +0.2841 | 2.1166 | +0.0531 | 0.2616 |
| MU3_target_set | 83 | +0.2332 | +0.2604 | 2.1149 | +0.0796 | 0.3412 |

- **最佳可部署 μ: `MU3_target_set`**

### `>=30` 子集（μ）

| 模型 | n | Pearson | slope |
|---|---:|---:|---:|
| MU0_intercept | 34 | -0.0000 | +0.0000 |
| MU1_target | 34 | -0.1261 | -0.0675 |
| MU2_set | 34 | +0.3367 | +0.1664 |
| MU3_target_set | 34 | +0.2531 | +0.1407 |

## PHASE 4 — Scale model σ

预注册：**n ≥ 3** 才用于 scale 建模（敏感性 n ≥ 5）

| 模型 | 特征 | Pearson(log σ) | Spearman | RMSE(log σ) | MAE | Pearson(σ) | slope |
|---|---:|---:|---:|---:|---:|---:|---:|
| SIG0_global | 0 | +0.0000 | +nan | 1.0161 | 0.6353 | -0.0000 | -0.0000 |
| SIG1_target | 58 | -0.2975 | -0.3255 | 1.0586 | 0.6881 | -0.3088 | -0.0393 |
| SIG2_set | 25 | +0.3405 | +0.3445 | 0.9561 | 0.5981 | +0.3235 | +0.1721 |
| SIG3_target_set | 83 | +0.2604 | +0.2076 | 0.9866 | 0.6302 | +0.2595 | +0.1271 |

- **最佳可部署 σ: `SIG2_set`**

## PHASE 5 — Location–Scale–Rank 重建

| 模型 | Overall r | RMSE | target-mean r | slope | C_W | Macro ρ | Eff@10 |
|---|---:|---:|---:|---:|---:|---:|---:|
| P0_raw_ranker | -0.0942 | 3.4192 | +0.1908 | +0.1349 | +0.0441 | +0.1117 | +0.2522 |
| P1_location_only | +0.2565 | 2.4744 | +0.2332 | +0.0796 | +0.0998 | +0.1117 | +0.2522 |
| P2_scale_only | +0.1804 | 2.6829 | +0.0370 | +0.0000 | +0.1804 | +0.1117 | +0.2522 |
| P3_full_LSR | +0.2641 | 2.6634 | +0.2332 | +0.0796 | +0.1554 | +0.1117 | +0.2522 |
| P4_set_location | +0.2636 | 2.4249 | +0.2030 | +0.0531 | +0.1059 | +0.1117 | +0.2522 |
| P5_set_only_LSR | +0.2658 | 2.6174 | +0.2030 | +0.0531 | +0.1598 | +0.1117 | +0.2522 |

- **sanity：正 σ 下靶点内排序不变 → ✅**
- 说明：`μ+σz` 是靶点内**正仿射变换**，因此 P1–P5 的 within 排序应与 P0 **完全一致**（Macro ρ / C_W / Eff@10 不变）；
  变化只应出现在 **between / calibration** 侧。这是设计的预期行为。

## PHASE 6 — Oracle 瓶颈分解（**DIAGNOSTIC_ONLY_NONDEPLOYABLE**）

| 模型 | Overall r | RMSE | Macro ρ | Eff@10 | 说明 |
|---|---:|---:|---:|---:|---|
| O0_deployable_P3 | +0.2641 | 2.6634 | +0.1117 | +0.2522 | |
| O1_oracle_mu | +0.6056 | 2.0744 | +0.1117 | +0.2522 | |
| O2_oracle_sigma | +0.2631 | 2.7545 | +0.1117 | +0.2522 | |
| O3_oracle_calib | +0.5864 | 2.1900 | +0.1117 | +0.2522 | |
| O4_oracle_z | +0.6733 | 1.8698 | +1.0000 | +1.0000 | |
| O5_full_oracle | +1.0000 | 0.0000 | +1.0000 | +1.0000 | sanity |

- **O5 sanity：max\|O5 − y\| = 1.776e-15** → ✅ 数值重建成功

### 6.1 瓶颈归因（Gain vs O0=P3）

| Oracle | ΔPearson | ΔRMSE | ΔMacro ρ | ΔEff@10 |
|---|---:|---:|---:|---:|
| O1_oracle_mu | +0.3415 | -0.5891 | +0.0000 | +0.0000 |
| O2_oracle_sigma | -0.0010 | +0.0911 | +0.0000 | +0.0000 |
| O3_oracle_calib | +0.3223 | -0.4734 | +0.0000 | +0.0000 |
| O4_oracle_z | +0.4092 | -0.7936 | +0.8883 | +0.7478 |
| O5_full_oracle | +0.7359 | -2.6634 | +0.8883 | +0.7478 |

| 量 | 值 |
|---|---:|
| Gain_oracle_mu (O1−O0) ΔPearson | +0.3415 |
| Gain_oracle_sigma (O2−O0) ΔPearson | -0.0010 |
| Gain_oracle_calibration (O3−O0) ΔPearson | +0.3223 |
| Gain_oracle_rank (O4−O0) ΔPearson | +0.4092 |
| Gain_oracle_calibration ΔRMSE | -0.4734 |
| Gain_oracle_rank ΔRMSE | -0.7936 |

### **BOTTLENECK = `CALIBRATION`**

- 耗时 0.2 分钟
