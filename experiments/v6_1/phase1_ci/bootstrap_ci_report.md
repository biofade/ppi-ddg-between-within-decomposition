# V6.1 PHASE 1.1 — Cluster Bootstrap CI（CATH held-out）

- 时间: 2026-10-13 23:23:37
- 数据源: `experiments/v6_sota_audit/standardized_preds.csv` (`dataset_name == 'CATH'`)
- 重采样单位: **complex_id**；次数 **5000**；seed **42**
- 若某次迭代 unique complex < **15** → 跳过并计入 `n_skipped`
- 指标实现: 复用 V5 `v5_common`（`per_structure_spearman` / `topk_metrics`）

- CATH 数据: **3761 行 / 21 complex**；模型 **['DDGPred', 'ESM-IF', 'Flex-ddG', 'FoldX', 'proteinMPNN']**

> ⚠️ **样本量警示**：仅 21 个 complex（< 30）。CI 会较宽，且 macro 统计量的抽样分布可能非正态 → **结论需谨慎**（V6.1 全局约束 9）。

## 结果

| 模型 | n 行 | n complex | 指标 | point | 95% CI | SD | 有效次数 | 跳过 | P(>0) |
|---|---:|---:|---|---:|---|---:|---:|---:|---:|
| Flex-ddG | 741 | 21 | macro_pdb_spearman | +0.4138 | [+0.3294, +0.4993] | 0.0435 | 1151 | 3849 | 1.0000 |
| Flex-ddG | 741 | 21 | efficiency_at_10 | +0.5930 | [+0.1849, +0.7438] | 9.9958 | 1151 | 3849 | 0.9800 |
| FoldX | 755 | 21 | macro_pdb_spearman | +0.3561 | [+0.2648, +0.4407] | 0.0462 | 1151 | 3849 | 1.0000 |
| FoldX | 755 | 21 | efficiency_at_10 | +0.5027 | [+0.1658, +0.6445] | 0.3770 | 1151 | 3849 | 0.9896 |
| DDGPred | 755 | 21 | macro_pdb_spearman | +0.3234 | [+0.2124, +0.4285] | 0.0560 | 1151 | 3849 | 1.0000 |
| DDGPred | 755 | 21 | efficiency_at_10 | +0.4279 | [+0.0957, +0.6066] | 0.3778 | 1151 | 3849 | 0.9878 |
| ESM-IF | 755 | 21 | macro_pdb_spearman | +0.0277 | [-0.0805, +0.1332] | 0.0568 | 1151 | 3849 | 0.6768 |
| ESM-IF | 755 | 21 | efficiency_at_10 | +0.1503 | [-1.4814, +0.8067] | 46.9349 | 1151 | 3849 | 0.6412 |
| proteinMPNN | 755 | 21 | macro_pdb_spearman | +0.2199 | [+0.0924, +0.3429] | 0.0655 | 1151 | 3849 | 0.9991 |
| proteinMPNN | 755 | 21 | efficiency_at_10 | +0.2304 | [-1.1898, +0.5379] | 46.9316 | 1151 | 3849 | 0.7637 |

## 判定（CI 是否显著 > 0）

| 模型 | macro Spearman point | CI 下界 > 0? | Efficiency@10 point | CI 下界 > 0? |
|---|---:|---|---:|---|
| Flex-ddG | +0.4138 | ✅ 是 | +0.5930 | ✅ 是 |
| FoldX | +0.3561 | ✅ 是 | +0.5027 | ✅ 是 |
| DDGPred | +0.3234 | ✅ 是 | +0.4279 | ✅ 是 |
| ESM-IF | +0.0277 | ❌ 否（含 0） | +0.1503 | ❌ 否（含 0） |
| proteinMPNN | +0.2199 | ✅ 是 | +0.2304 | ❌ 否（含 0） |

## 参考：V5 baseline 的对应 CI（来自 V5 主推断，抽样单位 PDB，n=61）

| 指标 | point | 95% CI | n_clusters |
|---|---:|---|---:|
| macro per-PDB Spearman (seed42) | +0.1048 | [−0.0219, +0.2385] | 61 |
| macro per-PDB Spearman (seed43) | +0.0905 | [−0.0379, +0.2180] | 61 |
| macro per-PDB Spearman (seed44) | +0.0914 | [−0.0335, +0.2171] | 61 |

> 该三行的原始产物：`experiments/v5_0/p0_benchmark/baseline_repro_report.md`。
> **注意**：V5 与 CATH 数据集不同，此表**仅用于量级对照**，不是同集比较。

## RBM 趋势性支持（V6.1 PHASE 1.4）

| 模型 | n 行 | n 分箱 | macro Spearman（点估计，**无 CI**） |
|---|---:|---:|---:|
| FoldX | 1702 | 1 | +0.6478 |
| DSSA-PPI | 1702 | 1 | +0.3129 |
| SAAMBE-SEQ | 1702 | 1 | +0.2057 |
| PPIFormer | 1702 | 1 | +0.1904 |
| MpbPPI | 1702 | 1 | +0.0846 |
| ESM-2-FFN | 1702 | 1 | +0.0116 |
| ProBASS | 1702 | 1 | -0.0059 |
| MuTon | 1702 | 1 | -0.0202 |

> ⚠️ **仅作趋势性支持**：RBM 只有 **1 个分箱**（人为按 `pos` 连续分块），
> **无法计算 CI**，且其真值 `EXP_fitness` 是**适应度而非 ΔΔG**。
> 该数据集的任何数值**不得**与 CATH / V5 的数字直接比较。

- 耗时 46s
