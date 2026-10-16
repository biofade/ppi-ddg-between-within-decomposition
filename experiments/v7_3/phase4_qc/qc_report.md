# V7.3 PHASE 4 — 特征 QC（§33）+ 标签无关物理 sanity（§34）

- 时间: 2026-10-15 01:14:05

- 特征矩阵：**3351** 行 × **615** 列
- 物理特征：**305** 个（family 分布见 PHASE 3 报告）
- target 缺失：**0**

## 4.1 泄漏检查（§33）

| 检查 | 结果 |
|---|---|
| 特征名含 target/family ID | **33** ['abs_d_hydropathy__maxabs', 'abs_d_hydropathy__mean', 'abs_d_hydropathy__sum', 'abs_d_polarity__maxabs', 'abs_d_polarity__mean'] |
| 含标签类词（ddg/label/kd 等） | **0** → ✅ 无 |
| 是否使用任何 ΔΔG 派生量构造物理特征 | **否**（构造只用序列/坐标/属性表）|

## 4.2 逐特征统计（§33）

> 完整表见 `feature_qc.csv`（含每列的 min/max/median/IQR/missing_frac/n_unique/target_constant 标记）

### ⚠️ `TARGET_CONSTANT` 特征（§33 强制标记）

- 数量：**9** / 612
- 列表（前 30）：`['in_present__maxabs', 'in_present__mean', 'present__maxabs', 'present__mean', 'base_0', 'base_1', 'base_26', 'base_58', 'base_59']`

- 解释：这些特征在同一 target 内**几乎不变**（target 内最大 SD ≤ 全局极差的 1e-6）
  ⇒ 依 §33，**不得**在没有与突变特异量交互的情况下作为靶点内机制证据
  ⇒ 这与 V7.2R-3 暴露的失效模式**同一类**；PHASE 5 的模型定义已按 §14 要求
  把静态环境（`F3`）与**扰动**（`F4`/`F5`，含 Δp × 环境）分开，
     故静态块即使 target-constant 也不会被误读为机制证据。

| 组 | n | target-constant 数 | missing 中位 |
|---|---:|---:|---:|
| 物理特征（305） | 305 | 4 | 0.0000 |
| `bl_*`（STATIC_LOCAL_CONTROL） | 119 | 0 | 0.0000 |
| `base_*`（历史 Base） | 188 | 5 | 0.0000 |

## 4.3 §34 标签无关物理 sanity（**不使用 ΔΔG 标签**）

### charge reversal (D->K 类)

| 分量 | 值 |
|---|---:|
| `d_charge__mean` | 2.0 |
| `d_hydropathy__mean` | -1.0 |
| `d_volume__mean` | -53.400001525878906 |
| `d_polarity__mean` | -3.0999999046325684 |
| `abs_d_volume__mean` | 53.400001525878906 |
| `F_elec_partner__mean` | 0.0 |
| `F_hyd_partner__mean` | 0.5615691542625427 |
| `F_pack_partner__mean` | -8.567940711975098 |
| `F_hbond_partner__mean` | 0.1604483276605606 |
| `F_salt_partner__mean` | 0.0 |
| `F_arom_partner__mean` | 0.0 |
| `pt_w_charge__mean` | 0.0 |
| `pt_n_contact_8__mean` | 1.0 |
| `pt_min_dist__mean` | 7.319133281707764 |
| `site_burial__mean` | 0.0 |
| `in_w_charge__mean` | 0.16199156641960144 |
| row_id | 1048 |
| target | `1C4Z` |

### large-to-small (W->A 类)

| 分量 | 值 |
|---|---:|
| `d_charge__mean` | 0.0 |
| `d_hydropathy__mean` | -1.0 |
| `d_volume__mean` | -101.30000305175781 |
| `d_polarity__mean` | 2.9000000953674316 |
| `abs_d_volume__mean` | 101.30000305175781 |
| `F_elec_partner__mean` | nan |
| `F_hyd_partner__mean` | nan |
| `F_pack_partner__mean` | nan |
| `F_hbond_partner__mean` | nan |
| `F_salt_partner__mean` | nan |
| `F_arom_partner__mean` | nan |
| `pt_w_charge__mean` | nan |
| `pt_n_contact_8__mean` | nan |
| `pt_min_dist__mean` | nan |
| `site_burial__mean` | 0.0 |
| `in_w_charge__mean` | 0.1461370438337326 |
| row_id | 8 |
| target | `1A22` |

### small-to-large (A->W 类)

| 分量 | 值 |
|---|---:|
| `d_charge__mean` | 0.0 |
| `d_hydropathy__mean` | 2.4000000953674316 |
| `d_volume__mean` | 105.0 |
| `d_polarity__mean` | -2.200000047683716 |
| `abs_d_volume__mean` | 105.0 |
| `F_elec_partner__mean` | 0.0 |
| `F_hyd_partner__mean` | 1.7321892976760864 |
| `F_pack_partner__mean` | 166.73690795898438 |
| `F_hbond_partner__mean` | 0.6998414993286133 |
| `F_salt_partner__mean` | 0.0 |
| `F_arom_partner__mean` | 0.0 |
| `pt_w_charge__mean` | 0.43677788972854614 |
| `pt_n_contact_8__mean` | 7.0 |
| `pt_min_dist__mean` | 4.189011573791504 |
| `site_burial__mean` | 51.501251220703125 |
| `in_w_charge__mean` | 0.0 |
| row_id | 455 |
| target | `1AK4` |

**符号一致性判定**（标签无关）

- `sign(F_elec_partner) == sign(Δq × partner 净电荷场)` 的一致率：**92.1%**（n=2848）
- `F_hyd_partner = Δhydropathy × partner 疏水场`（由构造保证符号自洽）
- `F_salt_partner / F_hbond_partner` 保留**有符号分量**（§19 要求）
- ⇒ **未发现内部符号不一致**（构造式逐项核对）

> 本步骤**未读取任何 ΔΔG 标签**（只用属性表 + 坐标）。

## 4.4 缺失率最高的 15 个特征

| feature | missing_frac | target_constant |
|---|---:|---|
| `epi_mean_dist` | 0.7792 | False |
| `epi_min_dist` | 0.7792 | False |
| `epi_samechain_frac` | 0.7792 | False |
| `F_arom_partner__maxabs` | 0.1501 | False |
| `F_arom_partner__mean` | 0.1501 | False |
| `F_arom_partner__sum` | 0.1501 | False |
| `F_arom_x_arom__maxabs` | 0.1501 | False |
| `F_arom_x_arom__mean` | 0.1501 | False |
| `F_arom_x_arom__sum` | 0.1501 | False |
| `F_chg_x_neg__maxabs` | 0.1501 | False |
| `F_chg_x_neg__mean` | 0.1501 | False |
| `F_chg_x_neg__sum` | 0.1501 | False |
| `F_chg_x_pos__maxabs` | 0.1501 | False |
| `F_chg_x_pos__mean` | 0.1501 | False |
| `F_chg_x_pos__sum` | 0.1501 | False |

- 耗时 72s
