# V6.2.1 PHASE 4 — Target-size Stratification

- 时间: 2026-10-14 13:24:46
- 目标：解释 L4 上「macro-target Between Pearson > 0」
  与「row-weighted C_B ≈ 0/负」的并存
- **预注册分层**：A 三分位 / B 固定阈值(2–9,10–29,≥30) / C 连续相关
- **不事后挑阈值**

## 4.3 每 target 统计（plot-ready CSV）

- 写出 **460** 行（2 档 × 460 个 target-档）

## 4.1 分层（A 三分位 / B 固定阈值）

### L1_random

- target 数: **230**；突变数 min/median/max = **1/10/573**

**方案 A_tertile**

| bin | n_target | n_rows | macro mean-Pearson | mean-abs-error | row-w C_B | row-w C_W | Overall r | Macro within-ρ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| low | 82 | 258 | +0.7884 | 1.3715 | +0.5513 | -0.0270 | +0.5243 | -0.0952 |
| medium | 73 | 861 | +0.7541 | 0.7796 | +0.4515 | +0.0925 | +0.5440 | -0.0058 |
| high | 75 | 6404 | +0.8298 | 0.6163 | +0.4150 | +0.2947 | +0.7097 | +0.1882 |

**方案 B_fixed**

| bin | n_target | n_rows | macro mean-Pearson | mean-abs-error | row-w C_B | row-w C_W | Overall r | Macro within-ρ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 2-9 | 112 | 496 | +0.7811 | 1.2466 | +0.4691 | +0.0421 | +0.5112 | -0.0797 |
| 10-29 | 55 | 916 | +0.8514 | 0.6824 | +0.4164 | +0.0909 | +0.5073 | +0.0547 |
| >=30 | 63 | 6111 | +0.8418 | 0.6103 | +0.4167 | +0.2990 | +0.7157 | +0.2016 |

### L4_family

- target 数: **230**；突变数 min/median/max = **1/10/573**

**方案 A_tertile**

| bin | n_target | n_rows | macro mean-Pearson | mean-abs-error | row-w C_B | row-w C_W | Overall r | Macro within-ρ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| low | 82 | 258 | +0.6019 | 1.5166 | +0.4109 | -0.0197 | +0.3912 | -0.0464 |
| medium | 73 | 861 | +0.5431 | 0.9841 | +0.3162 | +0.0728 | +0.3890 | +0.0069 |
| high | 75 | 6404 | +0.0202 | 1.1977 | -0.0672 | +0.1621 | +0.0949 | +0.1013 |

**方案 B_fixed**

| bin | n_target | n_rows | macro mean-Pearson | mean-abs-error | row-w C_B | row-w C_W | Overall r | Macro within-ρ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 2-9 | 112 | 496 | +0.5886 | 1.3973 | +0.3153 | +0.0394 | +0.3546 | -0.0143 |
| 10-29 | 55 | 916 | +0.6885 | 0.9147 | +0.3076 | +0.0713 | +0.3789 | +0.0355 |
| >=30 | 63 | 6111 | -0.0008 | 1.2575 | -0.0720 | +0.1629 | +0.0909 | +0.0923 |

## 4.1C 连续分析：`|target_mean_error|` vs `log(n_mutations)`

| 档 | 模型 | Spearman | p | n_target |
|---|---|---:|---:|---:|
| L1_random | M0 | -0.3625 | 0.0000 | 230 |
| L1_random | M1 | -0.3592 | 0.0000 | 230 |
| L1_random | M2 | -0.3540 | 0.0000 | 230 |
| L4_family | M0 | -0.0873 | 0.1872 | 230 |
| L4_family | M1 | -0.0991 | 0.1340 | 230 |
| L4_family | M2 | -0.1561 | 0.0179 | 230 |

## 判定

- L4 各 bin 的 row-weighted `C_B` 极差: **0.4830**
- L4 各 bin 的 macro mean-Pearson 极差: **0.6893**
- L4 每 target 突变数范围: **1–573**

### 标记：**`TARGET_SIZE_HETEROGENEITY_PRESENT`**

### 对「macro Between > 0 但 row-weighted C_B ≈ 0」的解释

两个统计量**权重不同**，因此可以并存，**不矛盾**：

| 统计量 | 权重 | 回答的问题 |
|---|---|---|
| `macro`-target mean Pearson | **每个 target 等权** | 「一个随机**靶点**的均值能否被预测」 |
| row-weighted `C_B` | **按突变数加权**（`w_c = n_c/N`） | 「一个随机**突变记录**的跨靶点分量」 |

⇒ 当 mutation-rich target 的均值预测**较差/反号**、而多数小 target 尚可时，
**大 target 会主导 row-weighted 版本并把它拉向 0 或负值**，
而等权 macro 版本仍为正。本节的分层表正是用来检验这一机制的。

- 耗时 6s
