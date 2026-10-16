# V5.0 GATE A — V4 最强基线复现 (`C+DGaux_fwd+rev`)

- 设备: **CPU only**（4 核 i7-2630QM / 7.55 GB）；LightGBM 固定 800 树，**不用 val 早停**
- seeds: [42, 43, 44]
- 口径来源: `scripts/v4.0/11_dualhead_dg_ddg.py`（只读）


## 1. 自动泄漏断言（V5 §86–§90）

| 检查 | 值 | 期望 | 结论 |
|---|---:|---:|---|
| `pdb_overlap_train_val` | 0 | 0 | ✅ |
| `pdb_overlap_train_test` | 0 | 0 | ✅ |
| `pdb_overlap_val_test` | 0 | 0 | ✅ |
| `family_overlap_train_val` | 0 | 0 | ✅ |
| `family_overlap_train_test` | 0 | 0 | ✅ |
| `family_overlap_val_test` | 0 | 0 | ✅ |
| `split_exclusive` | 1 | 1 | ✅ |
| `split_covers_all` | 1 | 1 | ✅ |
| `pca_fit_rows_from_train_only` | 1 | 1 | ✅ |

-> ✅ **全部泄漏断言通过**
## 2. 数据规模

- 配对: **7523**（forward）；train 5265 / val 1129 / test 1129
- 复合物: train 118 / val 51 / test 61
- 特征: `D`（反对称块）**188** 维；`S`（对称块）**188** 维
- 辅助目标 `ΔG_mut`: mean -11.123, std 2.648

> test 只有 **61 个复合物** -> cluster bootstrap 会在极少 cluster 上重采样，CI 会很宽。这是数据规模的硬限制，如实报告。

## 3. 阶段 1 — 辅助任务 `ΔG_mut` 折外预测（GroupKFold by pdb_id）

- fold 0: train 4212, oof 1053, ΔG_mut Pearson **+0.3216**
- fold 1: train 4212, oof 1053, ΔG_mut Pearson **+0.5823**
- fold 2: train 4212, oof 1053, ΔG_mut Pearson **+0.2887**
- fold 3: train 4212, oof 1053, ΔG_mut Pearson **+0.5347**
- fold 4: train 4212, oof 1053, ΔG_mut Pearson **+0.4567**

- 辅助任务折外 Pearson 均值 **+0.4368**
- test 上辅助预测 Pearson **+0.4684**
- OOF 覆盖: 5265 / 5265 train 行（val/test 的 OOF 不参与训练，仅 test 用全 train 模型）

## 4. 阶段 2 — 主任务 ΔΔG（C+DGaux_fwd+rev）

| seed | Pearson | Spearman | RMSE | MAE | AUROC | macro-PDB ρ | macro-PDB r | P@10 | EF@1% | Hit10 | 耗时 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 42 | +0.3951 | +0.1972 | 2.1357 | 1.5249 | 0.5981 | +0.1048 | +0.1307 | 0.9000 | 9.0032 | 0.9000 | 22s |
| 43 | +0.3843 | +0.1997 | 2.1529 | 1.5364 | 0.6046 | +0.0905 | +0.1257 | 1.0000 | 9.0032 | 1.0000 | 17s |
| 44 | +0.3790 | +0.1640 | 2.1566 | 1.5476 | 0.5778 | +0.0914 | +0.1438 | 1.0000 | 9.9035 | 1.0000 | 18s |

## 5. GATE A 判定（|ΔPearson| > 0.01 即不符）

| seed | V4 记录 | V5 复现 | Δ | 判定 |
|---:|---:|---:|---:|---|
| 42 | +0.3951 | +0.3951 | -0.0000 | ✅ 一致 |
| 43 | +0.3843 | +0.3843 | -0.0000 | ✅ 一致 |
| 44 | +0.3790 | +0.3790 | +0.0000 | ✅ 一致 |

**GATE A: PASS ✅**


## 6. Cluster bootstrap（V5 §12 主推断：抽样单位 = PDB）

| seed | 指标 | point | 95% CI | n_clusters | 说明 |
|---:|---|---:|---|---:|---|
| 42 | overall_pearson | +0.3951 | [+0.0730, +0.6743] | 61 | cluster_bootstrap_effective |
| 42 | macro_pdb_spearman | +0.1048 | [-0.0219, +0.2385] | 61 | cluster_bootstrap_effective |
| 42 | macro_pdb_pearson | +0.1307 | [-0.0063, +0.2712] | 61 | cluster_bootstrap_effective |
| 43 | overall_pearson | +0.3843 | [+0.0357, +0.6690] | 61 | cluster_bootstrap_effective |
| 43 | macro_pdb_spearman | +0.0905 | [-0.0379, +0.2180] | 61 | cluster_bootstrap_effective |
| 43 | macro_pdb_pearson | +0.1257 | [-0.0064, +0.2585] | 61 | cluster_bootstrap_effective |
| 44 | overall_pearson | +0.3790 | [+0.0346, +0.6649] | 61 | cluster_bootstrap_effective |
| 44 | macro_pdb_spearman | +0.0914 | [-0.0335, +0.2171] | 61 | cluster_bootstrap_effective |
| 44 | macro_pdb_pearson | +0.1438 | [+0.0124, +0.2723] | 61 | cluster_bootstrap_effective |

> **n_clusters = 61** —— 这是 V5 统计推断的真实独立样本量（不是 1129）。CI 必然宽，不得用 row-level CI 冒充。

- 耗时 660s；峰值 RAM 0 MB
