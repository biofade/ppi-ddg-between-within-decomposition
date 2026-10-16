# V7.0 PHASE 1b — 修正后的 BLOCK L（局部结构环境，93.1% WT 校验）

- 时间: 2026-10-14 14:05:59
- **正确映射约定**：`链 = ann_chain`，`链内第 seq_index0 个残基`（0-based）
- 实测 WT 一致率 **93.1%**（EXACT 档）；对照 `pdb_resseq`/`ann_pos1`/`numbering_offset` 均 46.9%
- 只使用 **WT 校验通过** 的位点

- canonical rows: **7523**；有映射记录: **6738**
- mapping_status: {'EXACT': 5503, 'AMBIGUOUS': 3563, 'ALIGNMENT_CONFIDENT': 2580}

## 映射统计（逐 pair_row）

| 项 | 数 |
|---|---:|
| n_total | 6738 |
| n_nopdb | 0 |
| n_chain_missing | 0 |
| n_idx_oob | 1632 |
| n_wt_mismatch | 2305 |
| n_ok | 4473 |

### BLOCK L 覆盖
- **rows: 4473 / 7523 = 59.5%**
- **targets: 208 / 230**
- 特征列（除元数据）: **118**

- 写出 `matched_coverage_rows.csv`（PHASE 2 的统一 row universe）

- 耗时 116s
