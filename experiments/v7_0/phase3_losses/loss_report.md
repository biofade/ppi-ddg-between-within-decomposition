# V7.0 PHASE 3 — Loss Ablation（A1 特征集）

- 时间: 2026-10-14 16:09:10
- 架构: **A1 = Base + Set**（PHASE 2 唯一三指标同向改善者）
- 实现: MLP(128) CPU torch；epochs=60, lr=0.001
- δ_pair = 0.5 kcal/mol；每 target 上限 200 对
- 强制：target 均值**只用本折 train**；无 label-derived 特征

- row universe: **4473** 行 / **214** targets
- 特征: base 188 + S 25

- L4_family unit 数: 108

| L0_absolute | r=-0.0942 | C_W=+0.0441 | macro_ρ=+0.1118 | Eff@10=+0.2522 | | self-centered macro_ρ=+0.1118 |
| L1_centered | r=-0.1004 | C_W=+0.0428 | macro_ρ=+0.1095 | Eff@10=+0.2404 | | self-centered macro_ρ=+0.1095 |
| L2_pairwise | r=-0.0984 | C_W=+0.0425 | macro_ρ=+0.1112 | Eff@10=+0.2323 | | self-centered macro_ρ=+0.1112 |
| L3_target_balanced | r=+0.0065 | C_W=+0.0491 | macro_ρ=+0.0674 | Eff@10=+0.1124 | | self-centered macro_ρ=+0.0674 |

## 汇总

| Loss | Overall r | C_W | Macro ρ | Kendall | Eff@10 | C_B(row) | C_B(eq) |
|---|---:|---:|---:|---:|---:|---:|---:|
| L0_absolute | -0.0942 | +0.0441 | +0.1118 | +0.1079 | +0.2522 | -0.1383 | +0.1236 |
| L1_centered | -0.1004 | +0.0428 | +0.1095 | +0.1059 | +0.2404 | -0.1432 | +0.1216 |
| L2_pairwise | -0.0984 | +0.0425 | +0.1112 | +0.1031 | +0.2323 | -0.1409 | +0.1266 |
| L3_target_balanced | +0.0065 | +0.0491 | +0.0674 | +0.0318 | +0.1124 | -0.0426 | +0.0477 |

- 耗时 13.3 分钟
