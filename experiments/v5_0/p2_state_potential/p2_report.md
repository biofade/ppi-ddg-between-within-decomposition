# V5.0 P2 — Thermodynamic State Potential (`ΔΔG = h(mut) − h(WT)`)

> 假设（V5 §47）: 学 state function 而非 pair 函数，可**自动**获得反对称性与可加性。


## 0. State 行统计

- 唯一变体行: **5108**；train 行 3241，test 行 999，**交集 0**
- 有绝对 ΔG 的变体行: 5108
- 反例检查：若 train/test 行交集 > 0，说明 state 模型会看到 test 变体 → 必须排除
- **state 模型 fit 行（仅 train，剔除任何 test 出现的变体）: 3241**

## 1. 结果（test n=1129, n_complex=61）

| seed | variant | Pearson | Spearman | macro-PDB ρ | RMSE | MAE | AUROC |
|---:|---|---:|---:|---:|---:|---:|---:|
| 42 | P2-0_pair_plain | +0.3674 | +0.1361 | +0.0835 | 2.1939 | 1.6106 | 0.5627 |
| 42 | P2-1_pair_dgaux | +0.4040 | +0.2129 | +0.1133 | 2.1248 | 1.5197 | 0.6122 |
| 42 | P2-2_state_absolute | +0.0797 | +0.1338 | -0.0002 | 2.9898 | 2.1805 | 0.5630 |
| 42 | P2-3_state_centered | -0.0758 | -0.1102 | -0.0297 | 2.5915 | 1.7685 | 0.4416 |
| 43 | P2-0_pair_plain | +0.3434 | +0.1273 | +0.0625 | 2.2254 | 1.6494 | 0.5571 |
| 43 | P2-1_pair_dgaux | +0.3740 | +0.1876 | +0.1073 | 2.1622 | 1.5471 | 0.5934 |
| 43 | P2-2_state_absolute | +0.0343 | +0.0701 | -0.0183 | 2.9859 | 2.1233 | 0.5228 |
| 43 | P2-3_state_centered | -0.1104 | -0.1710 | -0.0483 | 2.5988 | 1.7807 | 0.3984 |
| 44 | P2-0_pair_plain | +0.3497 | +0.1229 | +0.1339 | 2.2142 | 1.6271 | 0.5519 |
| 44 | P2-1_pair_dgaux | +0.3676 | +0.1625 | +0.1392 | 2.1668 | 1.5571 | 0.5792 |
| 44 | P2-2_state_absolute | +0.0642 | +0.1060 | +0.0019 | 2.9535 | 2.1257 | 0.5449 |
| 44 | P2-3_state_centered | -0.0596 | -0.0753 | -0.0115 | 2.5676 | 1.7296 | 0.4542 |

## 2. 配对比较 vs P2-1（cluster bootstrap, 同 seed, 单位 = PDB）

| 对比 | 指标 | Δ | 95% CI | p | n_clusters |
|---|---|---:|---|---:|---:|
| P2-0_pair_plain − P2-1_pair_dgaux | overall_pearson | -0.0284 | (见 paired_comparisons.csv) | 0.1456 | 61 |
| P2-0_pair_plain − P2-1_pair_dgaux | macro_pdb_spearman | -0.0266 | (见 paired_comparisons.csv) | 0.4328 | 61 |
| P2-2_state_absolute − P2-1_pair_dgaux | overall_pearson | -0.3225 | (见 paired_comparisons.csv) | 0.0040 | 61 |
| P2-2_state_absolute − P2-1_pair_dgaux | macro_pdb_spearman | -0.1255 | (见 paired_comparisons.csv) | 0.2412 | 61 |
| P2-3_state_centered − P2-1_pair_dgaux | overall_pearson | -0.4639 | (见 paired_comparisons.csv) | 0.0017 | 61 |
| P2-3_state_centered − P2-1_pair_dgaux | macro_pdb_spearman | -0.1498 | (见 paired_comparisons.csv) | 0.1456 | 61 |

## 3. 一致性指标（V5 §54）

| variant | seed | E_anti (mean\|f(A,B)+f(B,A)\|) |
|---|---:|---:|
| P2-0_pair_plain | 42 | 2.531116 |
| P2-1_pair_dgaux | 42 | 2.298994 |
| P2-2_state_absolute | 42 | 0.000000 |
| P2-3_state_centered | 42 | 0.000000 |
| P2-0_pair_plain | 43 | 2.526158 |
| P2-1_pair_dgaux | 43 | 2.289327 |
| P2-2_state_absolute | 43 | 0.000000 |
| P2-3_state_centered | 43 | 0.000000 |
| P2-0_pair_plain | 44 | 2.505795 |
| P2-1_pair_dgaux | 44 | 2.269591 |
| P2-2_state_absolute | 44 | 0.000000 |
| P2-3_state_centered | 44 | 0.000000 |

> pair 模型的 E_anti 用「反向训练的独立模型」估计（下界，非严格反对称）；
> state 模型的 E_anti = 0 **由构造保证**（V5 §48 的核心优势）。

## 4. GATE D 判定（V5 §124，成功标准 §55）

| variant | 3-seed 平均 ΔPearson vs P2-1 | Δmacro-PDB ρ | 判定 |
|---|---:|---:|---|
| P2-0_pair_plain | -0.0284 | -0.0266 | 未达标准 |
| P2-2_state_absolute | -0.3225 | -0.1255 | 未达标准 |
| P2-3_state_centered | -0.4639 | -0.1498 | 未达标准 |

- 耗时 1478s；峰值 RAM 0 MB
