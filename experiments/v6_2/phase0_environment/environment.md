# V6.2 PHASE 0 — 环境与目录

**时间**：2026-10-14

## 目录（已建）

```
experiments/v6_2/
  README.md
  PROGRESS.md
  phase0_environment/
  phase1_provenance/
  phase2_split4/
  phase3_decomposition/
  phase4_ranking_metrics/
  phase5_bootstrap/
  phase6_leave_one_out/
  phase7_leakage/
  phase8_summary/
scripts/v6_2/
logs/v6_2/
```

## 关键输入存在性（**8/8 OK**）

| 文件 | 存在 | 大小 (B) |
|---|---|---:|
| `experiments/v6_1/phase2_split_inflation/split_inflation_results.csv` | ✅ | 969 |
| `data/v4.0/ddg_pairs.csv` | ✅ | 2,238,640 |
| `data/v4.0/ddg_split.json` | ✅ | 82,187 |
| `experiments/v5_0/p0_benchmark/baseline_repro_predictions.npz` | ✅ | 62,577 |
| `features/v2.2_proteinea_features.csv` | ✅ | 11,852,893 |
| `features/v2_esm2_cache.npz` | ✅ | 52,481,097 |
| `experiments/v6_sota_audit/standardized_preds.csv` | ✅ | 1,270,276 |
| `experiments/v5_0/p4_conformal/p4_predictions.npz` | ✅ | 165,225 |

## 运行时

| 项 | 值 |
|---|---|
| OS | Windows 10 (10.0.19045) |
| Python | 3.11.9 (`E:\sjtu\download\python.exe`) |
| CPU | Intel i7-2630QM, 4 核 / 8 线程 |
| RAM | 7.55 GB（可用通常 1.3–2.1 GB） |
| **GPU** | **无**（CPU only，遵守约束 3） |
| 依赖 | numpy 2.4.6 / pandas 3.0.6 / scipy 1.17.1 / sklearn 1.9.1 / lightgbm 4.7.0 / matplotlib 3.11.2 / torch 2.5.1+cpu |
| LightGBM 前置 | 需 `ctypes.CDLL(vcomp140.dll)` shim → 已封装于 `v5_common.load_lgb()` |

## 磁盘纪律

- 全部产物写入 **F:\DSH_Galatea\Affinity_Intelligence**（工作区）
- 本 PHASE **无下载**（V6.2 全程未新增下载）
- 未写入 C 盘任何位置

## 约束确认

- ✅ 不修改 V1.1–V6.1 任何文件（V6.2 只读 V6.1 与 V5 产物）
- ✅ `v6_2` 命名空间（`scripts/v6_2/`、`experiments/v6_2/`、`logs/v6_2/`）
- ✅ CPU only
- ✅ 未下载任何文件（故 <100 MB 约束自然满足）
- ✅ 失败记录到 `logs/v6_2/failures.log`（本次实际**无失败**）
