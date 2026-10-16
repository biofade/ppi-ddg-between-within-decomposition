# V6.3 PHASE 1 — Target-level 机制表

- 时间: 2026-10-14 13:49:01
- 主键: `pdb_id`（= target）；另含 `family`
- **严格区分 `LABEL_DERIVED` 与 `DEPLOYABLE`**（V6.3 §1.1）

- matched set: **7523** 行 / **230** targets
- ddg_pairs forward 对齐: **7523** 行

- L4 预测可用行: **7523**；M1 ensemble 来自 3 seeds

- 能映射到 v2.2 行的 mutation: **7521 / 7523**（100.0%）
- 可用结构描述符（数值型）: **58** 个

## 结果：**203** 个 target（原始 230）

## 缺失率（非空比例 < 90% 的字段）

| 字段 | 非空 % | 类型 |
|---|---:|---|
| `within_spearman` | 84.2 | LABEL_DERIVED |
| `within_kendall` | 84.2 | LABEL_DERIVED |

- **完全不可用（UNAVAILABLE）字段: 0** 

- 耗时 18s
