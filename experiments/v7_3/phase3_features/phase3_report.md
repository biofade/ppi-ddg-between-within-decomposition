# V7.3 PHASE 3 — 物理扰动代理特征

- 时间: 2026-10-15 01:06:46

> **术语纪律（§11）**：以下全部是 `fixed-backbone physics-derived perturbation proxies`，**不是**精确自由能。

## 3.1 20-AA 物理属性表（§12）

| 属性 | 来源/标度 |
|---|---|
| `charge` | formal charge at pH 7.4 (His 取 0.1 部分质子化) |
| `hydropathy` | Kyte & Doolittle, J. Mol. Biol. 157:105 (1982) |
| `volume` | side-chain volume (Zamyatnin 1972 / standard tabulation) |
| `polarity` | Grantham, Science 185:862 (1974) polarity scale |
| `donor` | 侧链 H-bond donor 能力（粗粒度计数） |
| `acceptor` | 侧链 H-bond acceptor 能力（粗粒度计数） |
| `aromatic` | flag: F/W/Y/H |
| `cation` | flag: K/R/H |
| `anion` | flag: D/E |
| `aliphatic` | flag: A/V/L/I/M |
| `gly` | flag: G |
| `pro` | flag: P |

→ `features/v7_3/aa_physics_table.csv`（20 行 × 14 列，**在看 outer 结果之前冻结**）

- U_PHYS_STRICT：**3351** 行 / 166 targets

  - [20] 1E96: n_rows=600 累计, 0.5 min
  - [40] 1K8R: n_rows=1180 累计, 0.8 min
  - [60] 1X1X: n_rows=1948 累计, 1.7 min
  - [80] 2J1K: n_rows=2144 累计, 1.8 min
  - [100] 3BT1: n_rows=2492 累计, 2.0 min
  - [120] 3S9D: n_rows=2743 累计, 2.2 min
  - [140] 4KRO: n_rows=3200 累计, 2.5 min
  - [160] 5CYK: n_rows=3329 累计, 2.5 min

- 完成 **166** targets / **3351** rows

## 3.2 物理特征族统计

| family | n_features |
|---|---:|
| `F_AROM` | 6 |
| `F_BURIAL` | 15 |
| `F_ELEC` | 6 |
| `F_ENV_INTRA` | 87 |
| `F_ENV_PARTNER` | 87 |
| `F_EPI_GEOM` | 5 |
| `F_HBOND` | 6 |
| `F_HYD` | 12 |
| `F_IFACE` | 2 |
| `F_MUT` | 48 |
| `F_PACK` | 12 |
| `F_SALT` | 6 |
| `OTHER` | 13 |
| **合计物理特征** | **305** |

- 输出：`features/v7_3/physics_features.parquet` （3351 行 × 613 列）
- 含 `base_*`（历史 Base）与 `bl_*`（STATIC_LOCAL_CONTROL）

- 耗时 2.6 分钟
