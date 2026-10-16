# V6.3 PHASE 6 — Candidate-Set Calibration Test

- 时间: 2026-10-14 13:52:52
- 单位 = **target**；**family-level CV**（T 与 TS 用相同 folds）
- **禁止** label-derived 变量进入模型

- 可用 target-only 特征: **23**
- 可用 set 特征: **20**

- 两者都完整的 target: **203** / 203

## 6.5 结果（family-level CV）

| 模型 | 特征数 | n | Pearson | Spearman | RMSE | MAE |
|---|---:|---:|---:|---:|---:|---:|
| T_target_only | 23 | 203 | -0.0438 | -0.0470 | 1.4697 | 0.9723 |
| TS_target+set | 43 | 203 | +0.2461 | +0.2083 | 1.4375 | 1.0294 |

## 6.5b `>=30` mutation targets 单独报告

| 模型 | n | Pearson | Spearman | RMSE | MAE |
|---|---:|---:|---:|---:|---:|
| T_target_only (2-9) | 85 | -0.1265 | -0.1060 | 1.5016 | 1.0247 |
| TS_target+set (2-9) | 85 | +0.1446 | +0.1586 | 1.4981 | 1.1190 |
| T_target_only (10-29) | 55 | -0.1896 | -0.2565 | 1.3746 | 0.8007 |
| TS_target+set (10-29) | 55 | +0.1875 | +0.2228 | 1.3297 | 0.8102 |
| T_target_only (>=30) | 63 | +0.0878 | +0.1537 | 1.5061 | 1.0513 |
| TS_target+set (>=30) | 63 | +0.2614 | +0.1799 | 1.4446 | 1.0999 |

## 6.6 Paired target bootstrap（单位 = target，5000 次）

| 统计量 | 全样本 | `>=30` 子集 |
|---|---:|---:|
| Delta_RMSE (TS-T) | -0.0303 [-0.1243, +0.0672] | — |
| Delta_MAE (TS-T) | +0.0578 [-0.0175, +0.1327] | — |
| Delta_Pearson (TS-T) | +0.2877 [+0.1237, +0.4504] | — |
| Delta_RMSE (TS-T) | — | -0.0557 [-0.2547, +0.1467] |

## 判定

- T  (target-only):  Pearson -0.0438, RMSE 1.4697
- TS (target+set)  : Pearson +0.2461, RMSE 1.4375
- ΔRMSE(TS−T) = **-0.0322**；ΔPearson = +0.2898

> ### **`SET_CONTEXT_NOT_SUPPORTED`**

> ⚠️ 注意：无论 T 还是 TS，若 family-CV Pearson 都很低（≈0），
> 则说明 **target-level mean 在 family OOD 下不可由这些描述符预测**，
> 这才是对 V7 `G(x_c, M_c)` 设计的真正约束。

- 耗时 60s
