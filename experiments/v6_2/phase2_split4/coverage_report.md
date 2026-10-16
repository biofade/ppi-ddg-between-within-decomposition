# V6.2 PHASE 2.1 — Canonical Row Universe（matched_analysis_set）

- 时间: 2026-10-14 08:51:17


- `build_pair_features` 产出: **7523** 行（真 ΔΔG forward 配对）
- 唯一复合物 (pdb_id): **230**

## state_id 构造（V6.2 §2.4 要求的 state identity）

- 序列列: ['Ligand Sequences', 'Receptor Sequences']
- 唯一 `state_id`: **6047** / 7523 行
  → 说明每个 `(pdb, WT序列, MUT序列)` 组合基本唯一（符合预期：每对独立）

## family_label（≤40% 序列同一性 + union-find）

- V4 split 提供的 `family_of_complex`: **231** 个映射
- 唯一 family: **112**

## 重复观测 / twin

- `(pdb, mutation)` 重复的行: **4722**（2011 个重复键）
- 本表为 forward-only，故不含 forward/reverse twin；重复键来自同一 (pdb, mutation) 的多条实验记录

## 过滤（四档都可用）

- 剔除 `NaN_ddG`: **0** 行
- 剔除 `missing_mutation`: **0** 行
- 剔除 `missing_pdb_id`: **0** 行
- 剔除 `missing_family_label`: **0** 行
- 剔除 `missing_features`: **0** 行

- **matched_analysis_set: 7523 行**（原始 7523）
- excluded: **0 行**

- 在 matched set 内，同一 `state_id` 出现多行的组数: **1019**（L1 需保证它们同 fold）

## 覆盖统计

| 项 | 值 |
|---|---:|
| 原始 forward 配对 | 7523 |
| **matched_analysis_set** | **7523** |
| excluded | 0 |
| matched 唯一 pdb_id | 230 |
| matched 唯一 family | 112 |
| matched 唯一 state_id | 6047 |
| ddG mean / std | +1.5222 / 2.3956 |
| ddG min / max | -11.4939 / +12.4388 |
| Base Rate (ΔΔG<0) | **0.2363** |

> **四档主比较全部基于 `matched_analysis_set`**；
> 另报 full-coverage sensitivity（V6.2 §2.1）。

- 耗时 26.3s
