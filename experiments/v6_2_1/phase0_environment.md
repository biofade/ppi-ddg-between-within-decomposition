# V6.2.1 — PHASE 0 环境与结果冻结

**时间**：2026-10-14

## 新建目录

```
scripts/v6_2_1/
experiments/v6_2_1/
  PROGRESS.md
  phase1_shuffle_audit/
  phase2_tree_audit/
  phase3_d_only/
  phase4_target_size/
  phase5_corrected_stats/
  phase6_route/
  V6_2_1_FINAL_SUMMARY.md   （PHASE 7 产出）
logs/v6_2_1/
```

## 冻结的 SHA256（审计基线，**只读**）

| 文件 | 大小 (B) | SHA256（前 40） |
|---|---:|---|
| `experiments/v6_2/phase2_split4/split4_predictions.npz` | 1,696,726 | `ba711ea1b349836ac18b44ba67f6ab2f0a0b680b` |
| `experiments/v6_2/phase3_decomposition/pearson_decomposition.csv` | 4,048 | `ce9baa9c4fe2cf16023999f22820ea298388d384` |
| `experiments/v6_2/phase3_decomposition/control_results.csv` | 948 | `4edcb024912700418ec35b8a082682643b18bd14` |
| `experiments/v6_2/V6_2_FINAL_SUMMARY.md` | 14,832 | `d036742a48b8083bc5e71605c50dbff341418edd` |
| `scripts/v6_2/02_split4_ddg.py`（M1/M0 训练脚本） | 24,563 | `6a391a9c08f64a372507619f28fc7b0453ae5156` |
| `scripts/v6_2/03_decomposition.py`（分解实现） | 22,697 | `5c74317a1b85954911f026753857b39e7b143c29` |
| `scripts/v5/v5_common.py`（LightGBM 配置源） | 27,748 | `91bee0814f66737b50edadd937d967644bbf3ec7` |
| `scripts/v5/02_baseline_repro.py`（V5 Gate A） | 13,749 | `f3b333c3fe88a828385eae945795150d42df489f` |
| `experiments/v6_2/phase2_split4/matched_analysis_set.csv` | 754,492 | `8d209eea08cc6058fb235b059c9e3725fcbc9ba8` |
| `experiments/v6_2/phase2_split4/split4_results.csv` | 4,253 | `52d812e77c1f0169b48de251c85edaf7cd9826d2` |

## 环境

| 项 | 值 |
|---|---|
| Python | 3.11.9 (`E:\sjtu\download\python.exe`) |
| CPU | Intel i7-2630QM, 4 核 / 8 线程 |
| RAM | 7.55 GB（可用通常 1.3–2.1 GB） |
| GPU | **无**（CPU only） |
| LightGBM | 4.7.0（需 `vcomp140.dll` shim，封装于 `v5_common.load_lgb()`） |

## V6.2.1 要解决的 4 个问题（来自任务书）

| # | 问题 | 处置 PHASE |
|---|---|---|
| 1 | within-target shuffle 与协方差分解**数学上不一致** | PHASE 1 |
| 2 | 报告对 LightGBM 树数出现 **300 vs 800** 冲突 | PHASE 2 |
| 3 | **M2 实际 ≡ M0**，D-only 对照未真正实施 | PHASE 3 |
| 4 | L4 上 macro Between Pearson 为正，但 row-weighted `C_B ≈ 0/负` | PHASE 4 |

## ⚠️ 开始前的一个数学说明（将在 PHASE 1 用代码验证，不预先下结论）

任务书描述「若只打乱每个 target 内部的 prediction↔mutation 对应关系，
则 between component 应保持不变，**shuffle 后的 covariance 不应把原本的 between
component 一并清零**」。

这在数学上是对的：**组内置换保持组均值不变**，故
`C_B_shuffle ≡ C_B_original` 应当成立（在完全相同的 row universe 上）。
因此 V6.2 报告的 `r_shuffle ≈ −0.019` 与 `C_B ≈ +0.430` **无法同时成立**，
除非旧 shuffle 实现并未做「纯组内置换」。

PHASE 1 将**从 `split4_predictions.npz` 重新独立实现**，并用强制断言检验：
`C_B_shuffle == C_B_original` 与 `r_shuffle == C_B_shuffle + C_W_shuffle`。

**注意**：任务书同时要求不得使用结论「global Pearson 全部来自 between-target」
（关键科学纪律 §1）。本 PHASE 0 不预设任何结论。

## 约束确认

- ✅ 不修改 V1.1–V6.2 任何已有文件（本任务只**读** V6.2 产物）
- ✅ 新命名空间 `v6_2_1`
- ✅ CPU only；**不下载任何文件**
- ✅ 优先复用 V6.2 已保存的 predictions（不重训 M0/M1）
- ✅ 判定只依据 source code / saved predictions / row IDs / hashes / formulas / assertions
