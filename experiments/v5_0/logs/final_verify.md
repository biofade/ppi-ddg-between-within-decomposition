# V5.0 交付前最终自检

- 时间: 2026-10-13 21:14:22

## 1. 产物清单

- 必需产物: **37/37** 存在

- 图: **6** 张
  - `fig1_overall_vs_per_structure.png` (61141 B)
  - `fig2_error_vs_similarity.png` (45399 B)
  - `fig3_mipa_ablation_NOT_AVAILABLE.png` (26543 B)
  - `fig4_topk_enrichment.png` (111527 B)
  - `fig5_conformal_calibration.png` (5027 B)
  - `fig6_split_strictness.png` (56557 B)

- 脚本: **21** 个

## 2. 关键数字与 CSV 的一致性核对

| 项 | 报告声称 | CSV 实测 | 一致? |
|---|---:|---:|---|
| Gate A Pearson (seed42) | 0.3951 | 0.3951 | ✅ |
| Gate A Pearson (seed43) | 0.3843 | 0.3843 | ✅ |
| Gate A Pearson (seed44) | 0.379 | 0.3790 | ✅ |
| P4 plain P@10 | 0.9667 | 0.9667 | ✅ |
| P4 plain EF@1% | 9.6034 | 9.6034 | ✅ |
| P4 plain Efficiency@10 | 0.1014 | 0.1014 | ✅ |
| within sel10 mean ΔΔG | 0.3577 | 0.3577 | ✅ |
| within random10 mean ΔΔG | 0.6277 | 0.6277 | ✅ |
| within oracle10 mean ΔΔG | -0.5709 | -0.5709 | ✅ |
| within hit10(ΔΔG<0) | 0.308 | 0.3080 | ✅ |
| P2-2 ΔPearson vs P2-1 | -0.3225 | -0.3225 | ✅ |
| P3 W2 ΔPearson vs W0 | -0.1688 | -0.1688 | ✅ |

**数字一致性: 全部通过 ✅**

## 3. `INVALID_DUE_TO_BUG` 标记文件（证据保留，不得当作有效结果）

- `experiments\v5_0\all_experiments.csv.INVALID_DUE_TO_BUG`
- `experiments\v5_0\final\variant_summary.csv.INVALID_DUE_TO_BUG`
- `experiments\v5_0\p2_state_potential\p2_results.csv.INVALID_DUE_TO_BUG`
- `experiments\v5_0\p3_balancing\p3_results.csv.INVALID_DUE_TO_BUG`
- `experiments\v5_0\p4_conformal\p4_results.csv.INVALID_DUE_TO_BUG`

> 这些文件是 bug 前结果的**证据保留**，已在报告中说明其不可用性。

## 4. 是否修改过任何 V1.1–V4.0 文件（时间戳红线）

- V5 会话起点（最早 V5 产物 mtime）: **2026-10-13 18:43:00**
- ✅ **未修改任何 V1.1–V4.0 文件** —— 全部 V5 产物仅落在 `scripts/v5/`、`experiments/v5_0/`、`data/v5_0/`；V4.0 前序产物（`scripts/v4.0/`、`features/v4.0_*`、`data/v4.0/`、`experiments/v4.0/`）全部早于 V5 会话起点。

## 5. 残留进程

- 由外层 PowerShell 检查（本脚本无法可靠枚举）

