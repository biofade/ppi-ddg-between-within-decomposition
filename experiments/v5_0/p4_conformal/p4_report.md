# V5.0 P4 — Conformal Top-K Selection（校准风险选候选）

> 目标（V5 §67）: 预算 ~10 clones 时选「最值得做」的 10 个突变体。
> 方法: one-sided conformal 保守下界选择，**不用 LambdaRank/xendcg**（V5 §69）。


## 1. 主模型 train OOF（complex-level GroupKFold）

- fold 0: n_oof=1053, complexes=23, OOF Pearson **-0.4234**
- fold 1: n_oof=1053, complexes=23, OOF Pearson **+0.5212**
- fold 2: n_oof=1053, complexes=24, OOF Pearson **+0.2513**
- fold 3: n_oof=1053, complexes=24, OOF Pearson **+0.0015**
- fold 4: n_oof=1053, complexes=24, OOF Pearson **+0.2514**

- train OOF 残差: mean +1.5928, std 3.1854, 覆盖 5265/5265

## 2. Top-K 结果（test n=1129, n_complex=61）

| seed | selector | α | P@5 | P@10 | EF@1% | EF@5% | Hit10(<0) | Hit10(<-0.5) | Eff@10 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 42 | **plain (无 conformal)** | — | 1.0000 | 1.0000 | 9.9035 | 6.7202 | 1.0000 | 1.0000 | 0.0881 |
| 42 | conformal_global_a10 | | 1.0000 | 1.0000 | 9.9035 | 6.7202 | 1.0000 | 1.0000 | 0.0881 |
| 42 | conformal_global_a20 | | 1.0000 | 1.0000 | 9.9035 | 6.7202 | 1.0000 | 1.0000 | 0.0881 |
| 42 | conformal_global_a25 | | 1.0000 | 1.0000 | 9.9035 | 6.7202 | 1.0000 | 1.0000 | 0.0881 |
| 42 | mondrian_a10 | | 1.0000 | 1.0000 | 9.9035 | 4.9518 | 1.0000 | 1.0000 | 0.0542 |
| 42 | mondrian_a20 | | 1.0000 | 1.0000 | 9.0032 | 4.9518 | 1.0000 | 1.0000 | 0.0537 |
| 42 | mondrian_a25 | | 1.0000 | 0.9000 | 9.0032 | 4.9518 | 0.9000 | 0.9000 | 0.0542 |
| 43 | **plain (无 conformal)** | — | 0.8000 | 0.9000 | 9.0032 | 6.7202 | 0.9000 | 0.9000 | 0.0838 |
| 43 | conformal_global_a10 | | 0.8000 | 0.9000 | 9.0032 | 6.7202 | 0.9000 | 0.9000 | 0.0838 |
| 43 | conformal_global_a20 | | 0.8000 | 0.9000 | 9.0032 | 6.7202 | 0.9000 | 0.9000 | 0.0838 |
| 43 | conformal_global_a25 | | 0.8000 | 0.9000 | 9.0032 | 6.7202 | 0.9000 | 0.9000 | 0.0838 |
| 43 | mondrian_a10 | | 1.0000 | 0.9000 | 9.0032 | 5.3055 | 0.9000 | 0.9000 | 0.0518 |
| 43 | mondrian_a20 | | 1.0000 | 0.9000 | 9.0032 | 5.1286 | 0.9000 | 0.9000 | 0.0688 |
| 43 | mondrian_a25 | | 1.0000 | 0.9000 | 9.0032 | 5.1286 | 0.9000 | 0.9000 | 0.0688 |
| 44 | **plain (无 conformal)** | — | 1.0000 | 1.0000 | 9.9035 | 5.8360 | 1.0000 | 1.0000 | 0.1323 |
| 44 | conformal_global_a10 | | 1.0000 | 1.0000 | 9.9035 | 5.8360 | 1.0000 | 1.0000 | 0.1323 |
| 44 | conformal_global_a20 | | 1.0000 | 1.0000 | 9.9035 | 5.8360 | 1.0000 | 1.0000 | 0.1323 |
| 44 | conformal_global_a25 | | 1.0000 | 1.0000 | 9.9035 | 5.8360 | 1.0000 | 1.0000 | 0.1323 |
| 44 | mondrian_a10 | | 1.0000 | 1.0000 | 9.9035 | 4.7749 | 1.0000 | 1.0000 | 0.0991 |
| 44 | mondrian_a20 | | 1.0000 | 1.0000 | 9.9035 | 4.5981 | 1.0000 | 1.0000 | 0.1013 |
| 44 | mondrian_a25 | | 1.0000 | 0.9000 | 9.0032 | 4.7749 | 0.9000 | 0.9000 | 0.1013 |

## 3. 配对比较（cluster bootstrap, 同 seed, 单位 = PDB）

| 对比 | 指标 | Δ | p | n_clusters |
|---|---|---:|---:|---:|
| conformal_global_a10 − plain | p_at_10 | +0.0000 | — | 61 |
| conformal_global_a10 − plain | ef_1pct | +0.0000 | — | 61 |
| conformal_global_a10 − plain | hit10_ddg_m05 | +0.0000 | — | 61 |
| conformal_global_a20 − plain | p_at_10 | +0.0000 | — | 61 |
| conformal_global_a20 − plain | ef_1pct | +0.0000 | — | 61 |
| conformal_global_a20 − plain | hit10_ddg_m05 | +0.0000 | — | 61 |
| conformal_global_a25 − plain | p_at_10 | +0.0000 | — | 61 |
| conformal_global_a25 − plain | ef_1pct | +0.0000 | — | 61 |
| conformal_global_a25 − plain | hit10_ddg_m05 | +0.0000 | — | 61 |
| mondrian_a10 − plain | p_at_10 | +0.0000 | — | 61 |
| mondrian_a10 − plain | ef_1pct | +0.0000 | — | 61 |
| mondrian_a10 − plain | hit10_ddg_m05 | +0.0000 | — | 61 |
| mondrian_a20 − plain | p_at_10 | +0.0000 | — | 61 |
| mondrian_a20 − plain | ef_1pct | -0.3001 | — | 61 |
| mondrian_a20 − plain | hit10_ddg_m05 | +0.0000 | — | 61 |
| mondrian_a25 − plain | p_at_10 | -0.0667 | — | 61 |
| mondrian_a25 − plain | ef_1pct | -0.6002 | — | 61 |
| mondrian_a25 − plain | hit10_ddg_m05 | -0.0667 | — | 61 |

## 4. Conformal coverage 诊断（V5 §79，在 train OOF 上）

| nominal | empirical coverage (y ≤ μ+q) | q |
|---|---:|---:|
| 90% | **0.8999** | 6.6571 |
| 80% | **0.8002** | 4.0325 |
| 75% | **0.7500** | 3.2169 |

> 若 empirical coverage ≈ nominal，则**校准本身是正确的**（V5 §78）。
> 但注意：常数 q 不改变候选**排序**，因此 global conformal 的 Top-K
> 与 plain 完全相同 —— 这是方法上的必然，不是失败。要真正改变 Top-K，
> q 必须**依赖样本**（Mondrian 或 similarity-aware）。

## 5. GATE E 判定（V5 §124 / §77）

| 判据 | plain | 最好变体 | 是否满足 |
|---|---:|---:|---|
| p_at_10 | 0.9667 | 0.9667 (conformal_global_a10) [+0.0000] | 未满足 |
| ef_1pct | 9.6034 | 9.6034 (conformal_global_a10) [+0.0%] | 未满足 |
| hit10_ddg_m05 | 0.9667 | 0.9667 (conformal_global_a10) [+0.0000] | 未满足 |

- 耗时 249s；峰值 RAM 0 MB
