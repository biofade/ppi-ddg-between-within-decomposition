# V7.2R-3 PHASE 0–1 — 修复链丢弃缺陷 + 构建含全部伙伴链的无标签包

- 时间: 2026-10-14 22:02:28

- 目标宇宙：**4219** 行 / **174** targets / **91** families（= V7.2R-2 已覆盖集，保证 U_ABC 可比）

## 1.1 链覆盖（修复后）

| 项 | V7.2R-2（错） | V7.2R-3（修复） |
|---|---:|---:|
| (target,chain) 记录 | 233 | **654** |
| 每 target 链数中位 | 1 | **3** |
| 行覆盖 | 4219 | **4219** |

- 无标签检查：**0** 处违规 → ✅ 无

## 1.2 包内容

| file | bytes | sha256(16) | 含标签 |
|---|---:|---|---|
| `canonical_rows.csv` | 62,792 | `021be6411d90c3e4` | **NO** |
| `variants.csv` | 118,146 | `d5d7da3d38516896` | **NO** |
| `mutation_sites.csv` | 218,805 | `162e8cbf56a75d6b` | **NO** |
| `target_chains.csv` | 148,936 | `1db4cef5f39d0954` | **NO** |
| `variant_reconstruction.csv` | 101,187 | `23853573bad97458` | **NO** |
| `PACKAGE_README.md` | 768 | `cfcd9e45623f3c31` | **NO** |

- **ZIP**: `kaggle_input_package.zip`（0.15 MB）
- **SHA256**: `e918d44230dc3fd4e880ab662c5d4b8c568aa2890c75680ad20aa6ca0bb744c0`

- 耗时 17s
