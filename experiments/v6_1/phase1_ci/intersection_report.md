# V6.1 PHASE 1.2 — 测试集交集分析

- 时间: 2026-10-13 23:16:22

- **CATH held-out test**（`standardized_preds.csv`, N≥5 过滤后）：**21 complex / 21 唯一 PDB**

- **V5 test**（`p4_predictions.npz` test_group）：**61 complex / 61 唯一 PDB**

- V5 test 涉及的 `interface_key`（严格复合物定义）：**83**

## A. PDB 级别交集（宽松）

- CATH 唯一 PDB: **21**
- V5   唯一 PDB: **61**
- **交集: 8**
- Jaccard = 0.1081

交集 PDB（小写，已剥离链）：

| # | PDB |
|---:|---|
| 1 | `1emv` |
| 2 | `1wqj` |
| 3 | `2kso` |
| 4 | `2o3b` |
| 5 | `2wpt` |
| 6 | `4e6k` |
| 7 | `4hfk` |
| 8 | `4l0p` |

## B. 复合物级别交集（严格）

- CATH complex_id 数: **21**
- V5   complex_id 数: **61**
- **交集: 8**

- V5 `interface_key`（含链信息）与 CATH complex_id 的交集: **0**

## C. 结论与对后续 Phase 的影响

> ### ⚠️ 交集偏少（8 个 PDB）→ 只做探索性对比，结论需谨慎

- 耗时 1s
