# V7.1.1 PHASE 2 — Metric 资格审计

> **核心区分**：`DATA COVERAGE`（数据可用性）vs `METRIC ELIGIBILITY`（某指标的最小 n 门槛）。
> **禁止**把 metric eligibility 写成模型 coverage loss。

数据：LSR predictions，**4473 行 / 214 targets**

| 阈值 | 满足的 target 数 |
|---|---:|
| n ≥ 1 | **214** |
| n ≥ 2 | **180** |
| n ≥ 3 | **164** |
| n ≥ 5 | **130** |
| n ≥ 10 | **85** |

| 指标 | 最小 target n | 合格 target 数 | 行数 | target 列表 hash |
|---|---:|---:|---:|---|
| Overall Pearson | ≥1 | **214** | 4473 | `47d7ea4d06c3` |
| RMSE | ≥1 | **214** | 4473 | `47d7ea4d06c3` |
| C_B / C_W | ≥1 | **214** | 4473 | `47d7ea4d06c3` |
| target-mean Pearson | ≥1 | **214** | 4473 | `47d7ea4d06c3` |
| attenuation slope | ≥3 | **164** | 4407 | `a2e1e2edecb6` |
| Macro per-target Spearman | ≥5 | **130** | 4288 | `b0fee31a976d` |
| Kendall tau-b | ≥2 | **180** | 4439 | `7fc2f74f6910` |
| Pairwise Accuracy | ≥2 | **180** | 4439 | `7fc2f74f6910` |
| Efficiency@10 | ≥5 | **130** | 4288 | `b0fee31a976d` |
| Hit@10 | ≥2 | **180** | 4439 | `7fc2f74f6910` |
| Precision@10 | ≥2 | **180** | 4439 | `7fc2f74f6910` |

## 2.2 复算 V7.1 报告的数字（确认其 target 基数）

- `z_predictions.npz` 内含 z 变量：['z_R0_A0_pooled', 'z_R1_A1_base_set', 'z_R2_MLP_L0', 'z_R3_Donly_M2', 'z_R4_hist_M0']

| V7.1 报告值 | 复算值 | 使用的 target 数 |
|---|---:|---:|
| Macro ρ = 0.1117 | **+0.1117** | **130**（n≥5） |
| Eff@10 = 0.2522 | **+0.2522** | **130**（n≥5） |
| C_W = 0.0441 | **+0.0441** | **214**（全部 target） |

> **结论**：同一个 `P0_raw_ranker` 预测，
> `Macro ρ` 基于 **130** 个 target、
> `C_W` 基于 **214** 个 target。
> 这些差异**完全是 metric eligibility**，**不是**模型 coverage loss。

