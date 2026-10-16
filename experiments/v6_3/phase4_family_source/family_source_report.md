# V6.3 PHASE 4 + 5 — Family/Source 集中度 + Influence + Weighting

- 时间: 2026-10-14 13:52:12

- L4 M1 ensemble: 7523 行 / 230 targets

## 4.1 Family 集中度

| 集合 | n_target | n_family | rows | top1 family 占比(targets) | top1 family 占比(rows) |
|---|---:|---:|---:|---:|---:|
| >=30 | 63 | 36 | 6112 | 0.302 | 0.006 |
| <30 | 140 | 86 | 1383 | 0.179 | 0.003 |
| ALL | 203 | 109 | 7495 | 0.217 | 0.001 |

- `>=30` 组 top-5 family：
  - `L98`：19 targets / 1324 rows
  - `L65`：3 targets / 175 rows
  - `L89`：2 targets / 417 rows
  - `L74`：2 targets / 149 rows
  - `L59`：2 targets / 102 rows

## 4.2 Source 集中度

- `source_top` 唯一值数: **1** → 值 = `SKEMPI v2.0`
- **判定：`UNAVAILABLE`** —— 全部 target 来自同一 source，无法做 source 分层或 leave-one-source-out
- `assay_top` 唯一值数: **11** → 可用作 assay 分层

| assay | n_target | n_rows | mean abs_err | mean signed_err |
|---|---:|---:|---:|---:|
| `BLI` | 9 | 98 | 0.7589 | -0.0460 |
| `ELISA` | 5 | 290 | 0.9084 | -0.7980 |
| `FL` | 16 | 1362 | 0.8901 | +0.0049 |
| `IASP` | 15 | 424 | 1.6308 | -1.2786 |
| `ITC` | 35 | 452 | 0.7342 | -0.4342 |
| `Other` | 21 | 989 | 1.3868 | -0.7048 |
| `RA` | 3 | 154 | 0.5451 | +0.5451 |
| `SP` | 7 | 568 | 1.2207 | -1.2207 |
| `SPR` | 89 | 3070 | 1.0658 | -0.4192 |

## 4.3 Leave-one-largest-family-out（在 >=30 组内）

| 移除的 family | 剩余 target | 移除 rows | macro mean-Pearson | row-w C_B |
|---|---:|---:|---:|---:|
| `L98` | 211 | 1323 | +0.5753 | -0.0209 |
| `L65` | 227 | 175 | +0.5570 | -0.0146 |
| `L89` | 228 | 417 | +0.5580 | -0.0042 |
| `L74` | 228 | 149 | +0.5580 | -0.0123 |
| `L59` | 228 | 102 | +0.5566 | -0.0122 |
| `L6` | 228 | 659 | +0.5603 | +0.0071 |
| `L14` | 228 | 780 | +0.5556 | -0.0152 |
| `L128` | 228 | 70 | +0.5565 | -0.0145 |

- **基线**：macro mean-Pearson = **+0.5568**，row-w C_B = **-0.0121**
- 移除任一 top family 后 macro mean-Pearson 范围: [+0.5556, +0.5753]
- **单 family 移除是否使现象消失**: 否

## 4.4 Leave-top-targets-out（按 n_mutations 降序）

| 移除 | 剩余 targets | macro mean-Pearson | row-w C_B | overall r |
|---|---:|---:|---:|---:|
| top 1 | 229 | +0.5595 | +0.0084 | +0.1643 |
| top 3 | 227 | +0.5596 | +0.0131 | +0.1424 |
| top 5 | 225 | +0.5730 | +0.0990 | +0.2684 |
| top 10 | 220 | +0.5752 | +0.1234 | +0.2397 |
| top 20 | 210 | +0.5853 | +0.1377 | +0.2640 |

### 4.3 判定

> **`NOT_SINGLE_FAMILY_DRIVEN`**

## PHASE 5 — Weighting 敏感性

> 目的**不是**挑最好看的 weighting，而是回答：
> **performance conclusion 对 target-size weighting 有多敏感？**

| weighting | L4 C_B | L4 C_W | L4 macro mean-Pearson | n_target |
|---|---:|---:|---:|---:|
| row (n_c/N) | -0.0121 | +0.1468 | +0.5568 | 230 |
| equal-target | +0.5475 | +0.0459 | +0.5568 | 230 |
| sqrt(n) | +0.1645 | +0.0941 | +0.5568 | 230 |
| capped 30 | +0.1081 | +0.0843 | +0.5568 | 230 |
| log(1+n) | +0.2183 | +0.0745 | +0.5568 | 230 |

| weighting | L1 C_B | L1 macro mean-Pearson |
|---|---:|---:|
| row (n_c/N) | +0.4323 | +0.8154 |
| equal-target | +0.7752 | +0.8154 |
| sqrt(n) | +0.4498 | +0.8154 |
| capped 30 | +0.3427 | +0.8154 |
| log(1+n) | +0.4487 | +0.8154 |

- **L4 `C_B` 跨 weighting 的极差 = 0.5596**（min -0.0121 / max +0.5475）
- 符号是否随 weighting 变化: **是（结论不稳定）**

- 耗时 11s
