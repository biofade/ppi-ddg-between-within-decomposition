# V7.3 PHASE 1+2 — U_PHYS_STRICT + 结构环境缓存

- 时间: 2026-10-15 00:25:34

- canonical：4473 行 / 214 targets

  - [10/214] 1BP3: 383 res, 2 chains, 累计 keep=315, 1.2 min
  - [20/214] 1E96: 363 res, 2 chains, 累计 keep=600, 2.2 min
  - [30/214] 1GC1: 920 res, 4 chains, 累计 keep=816, 3.0 min
  - [40/214] 1JTG: 853 res, 4 chains, 累计 keep=1163, 4.0 min
  - [50/214] 1N8Z: 1015 res, 3 chains, 累计 keep=1416, 8.9 min
  - [60/214] 1SBB: 950 res, 4 chains, 累计 keep=1903, 9.9 min
  - [70/214] 1XD3: 606 res, 4 chains, 累计 keep=1960, 10.7 min
  - [80/214] 2B0U: 799 res, 4 chains, 累计 keep=2004, 12.0 min
  - [90/214] 2DVW: 302 res, 2 chains, 累计 keep=2053, 13.1 min
  - [100/214] 2JEL: 520 res, 3 chains, 累计 keep=2230, 14.5 min
  - [110/214] 2QJA: 386 res, 4 chains, 累计 keep=2290, 18.5 min
  - [120/214] 2VN5: 401 res, 4 chains, 累计 keep=2300, 19.2 min
  - [130/214] 3BTD: 279 res, 2 chains, 累计 keep=2492, 20.3 min
  - [140/214] 3IDX: 749 res, 3 chains, 累计 keep=2609, 21.5 min
  - [150/214] 3N85: 1033 res, 3 chains, 累计 keep=2675, 23.1 min
  - [160/214] 3SE8: 779 res, 3 chains, 累计 keep=2785, 25.0 min
  - [170/214] 4BFI: 394 res, 2 chains, 累计 keep=3088, 25.9 min
  - [180/214] 4I77: 534 res, 3 chains, 累计 keep=3164, 28.4 min
  - [190/214] 4NKQ: 691 res, 3 chains, 累计 keep=3210, 29.9 min
  - [200/214] 4X4M: 1377 res, 6 chains, 累计 keep=3310, 32.0 min
  - [210/214] 5K39: 320 res, 2 chains, 累计 keep=3337, 33.2 min

## 1.1 U_PHYS_STRICT（§9）

| 项 | 值 |
|---|---:|
| canonical rows | 4473 |
| **U_PHYS_STRICT rows** | **3351** |
| targets | **166** |
| families | **86** |
| n>=1 targets | **166** |
| n>=2 targets | **146** |
| n>=3 targets | **129** |
| n>=5 targets | **102** |
| n>=10 targets | **65** |

### 排除原因（§10）

| reason | rows |
|---|---:|
| `ambiguous_mapping` | 720 |
| `wt_mismatch` | 257 |
| `missing_chain` | 83 |
| `idx_out_of_bounds` | 62 |
| **合计排除** | **1122** |

- 结构缓存：**214** 个 target → `F:\DSH_Galatea\Affinity_Intelligence\features\v7_3\pdb_environment_cache`
- SASA 成功：**213** / 214（freesasa）

- 耗时 33.4 分钟
