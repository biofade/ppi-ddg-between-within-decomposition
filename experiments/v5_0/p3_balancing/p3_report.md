# V5.0 P3 — Complex / Family Balanced Boosting

> 动机（V5 §57）: 单 complex 最多 ~573 pairs，row-wise 训练隐含「573 rows ≈ 573 独立观测」，这与 cluster bootstrap 的结论矛盾。


## 0. Train 内 complex 规模分布（仅 train，未使用 test 信息）

- train complex 数: **118**；pairs: 5265
- 每 complex pairs: 中位 12，均值 44.6，最大 573
- 最大 complex 占比: **10.88%** of train rows
- top-5 complex 的 pairs: [573, 534, 358, 300, 282] → 合计占 **38.9%**

## 1. 加权方案（test n=1129, n_complex=61）

| seed | variant | Pearson | Spearman | macro-PDB ρ | macro-PDB r | RMSE | AUROC |
|---:|---|---:|---:|---:|---:|---:|---:|
| 42 | W0_equal | +0.3951 | +0.1972 | +0.1048 | +0.1307 | 2.1357 | 0.5981 |
| 42 | W1_inv_sqrt | +0.3342 | +0.1051 | +0.0913 | +0.1156 | 2.2120 | 0.5280 |
| 42 | W2_inv_n | +0.1911 | +0.0353 | -0.0309 | -0.0356 | 2.3779 | 0.4913 |
| 43 | W0_equal | +0.3843 | +0.1997 | +0.0905 | +0.1257 | 2.1529 | 0.6046 |
| 43 | W1_inv_sqrt | +0.3270 | +0.0878 | +0.0223 | +0.0342 | 2.2155 | 0.5247 |
| 43 | W2_inv_n | +0.2227 | +0.0351 | +0.0266 | +0.0366 | 2.3322 | 0.4904 |
| 44 | W0_equal | +0.3790 | +0.1640 | +0.0914 | +0.1438 | 2.1566 | 0.5778 |
| 44 | W1_inv_sqrt | +0.3659 | +0.1067 | +0.0221 | +0.0653 | 2.1680 | 0.5259 |
| 44 | W2_inv_n | +0.2381 | +0.0444 | +0.0011 | +0.0058 | 2.3060 | 0.4967 |

## 2. Complex-balanced bagging（先抽 complex，再抽 ≤K mutations）

| seed | K | Pearson | Spearman | macro-PDB ρ | 聚合 |
|---:|---:|---:|---:|---:|---|
> **complex-balanced bagging 已跳过**（`--skip-bagging`）：K=5 试点已观测为显著负结果，且每 member 需完整训练，按 V5 §125「不允许 endless experimentation」不再全量运行。

## 3. Cross-split 稳定性（V5 §65/§66）

> 在 **12** 个随机 complex-level 划分上分别训练与评估（128 树，控制成本），比较 W0 与 W2 的 Pearson **SD**。
> 注意：这是**稳定性诊断**，不是主评估口径；test 不参与任何参数选择。

- **W0_equal**: n_splits=12, mean=+0.1930, **SD=0.2025**, range [-0.1593, +0.5059]
- **W2_inv_n**: n_splits=12, mean=+0.1616, **SD=0.1753**, range [-0.1457, +0.5005]

- SD 变化: 0.2025 → 0.1753（**+13.4%**）
- GATE（V5 §66 要求 SD 下降 ≥20%）: **未满足**

- 耗时 283s；峰值 RAM 0 MB
