# V7.2R-2 PHASE 4 — Kaggle Track B 特征导入

- 时间: 2026-10-14 21:04:45
- ZIP: `F:\dsh\downloads\v7_2r_2\V7_2R_2_TRACKB_FEATURES.zip`
- 实际 SHA256: `f3f0b8b12cb469668effc2cc15cedf90b6cdee9cf6ac4d8e44ef1966cbd231d9`
- 输入包 SHA256（记录）: `5bf0a25005f675db2f7360139bd8de87d14609011a704302565e6f6402bf537c`

- 解压 184 个条目 → `F:\DSH_Galatea\Affinity_Intelligence\features\v7_2r_2\_unzip`

- 找到 **174** 个 per-target npz
  - 已处理 40/174
  - 已处理 80/174
  - 已处理 120/174
  - 已处理 160/174

## 4.1 导入结果

| 项 | 值 |
|---|---|
| per-target npz | **174** / 174 |
| Track B rows | **4219** |
| unique row_id | **4219** |
| targets | **174** |
| families | **91** |
| dim(r_wt) | **1280** |
| dim(r_mut) | **1280** |
| dim(delta_r) | **1280** |

## 4.2 完整性检查

| 检查 | 结果 |
|---|---|
| 非预期 row_id（不在 canonical） | **0** 个文件 → ✅ 无 |
| 重复 row_id | **0** 个文件 → ✅ 无 |
| 无法映射到 canonical target 的行 | **0** |
| Kaggle `target_id` vs canonical 一致率 | **100.00%**（4219 行） |

> **`target identity` 一律采用 `canonical_target_labels.csv`**（§7），
> Kaggle 侧自报的 `target_id` 仅用于一致性核对。

- 耗时 5s
