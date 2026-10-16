# V7.4 PHASE 1 — U_RELAX_STRICT + Kaggle 输入包

- 时间: 2026-10-15 09:06:23

- 起点 U_PHYS_STRICT：**3351** 行 / 166 targets

  - 1000/3351 | keep=965 | 0.2 min
  - 2000/3351 | keep=1944 | 0.3 min
  - 3000/3351 | keep=2909 | 0.4 min

## 1.1 U_RELAX_STRICT（§9）

| 项 | 值 |
|---|---:|
| U_PHYS_STRICT rows | 3351 |
| **U_RELAX_STRICT rows** | **3253** |
| targets | **162** |
| families | **85** |
| n>=1 targets | **162** |
| n>=2 targets | **145** |
| n>=3 targets | **126** |
| n>=5 targets | **100** |
| n>=10 targets | **64** |

### 追加排除原因（V7.4 §9）

| reason | rows |
|---|---:|
| `chain_break_near_site_B` | 34 |
| `chain_break_near_site_A` | 24 |
| `insertion_code_at_site` | 19 |
| `chain_break_near_site_U` | 11 |
| `chain_break_near_site_G` | 6 |
| `gap_near_site_C` | 2 |
| `chain_break_near_site_E` | 1 |
| `chain_break_near_site_I` | 1 |
| **额外排除** | **98** |

## 1.2 Kaggle 输入包

| file | bytes | 含标签 |
|---|---:|---|
| `variants.csv` | 86,614 | **NO** |
| `mutation_sites.csv` | 142,095 | **NO** |
| `target_chains.csv` | 49,856 | **NO** |
| `canonical_rows.csv` | 48,521 | **NO** |
| `PACKAGE_README.md` | 925 | **NO** |

- ZIP：`kaggle_input_package.zip`（0.08 MB）
- SHA256：`0e00174cfbe5db01e7d7cebd6301b72ea45e09a8837057876c5fc059dddc75d2`
- **无标签检查：0 处违规 → ✅ 无**
- (target,chain) 记录：**213**；target 数 **162**

- 耗时 0.5 分钟
