# V5.0 — 严格 EXACT 位点子集（GATE B 的最保守口径）

> 筛选条件: `validated_status ∈ {EXACT, ALIGNMENT_CONFIDENT}` **且**
> `pdb_resseq_real > 0` **且** `real_aa == wt_aa` **且** `ann_pos1 == pdb_resseq_real`。
> 目的: 排除 `|offset|` 极大的偶然匹配（GATE B 已实测存在该现象）。

- 输入位点: **7869**

| 条件 | 剩余位点 |
|---|---:|
| + validated_status ∈ {EXACT, ALIGNMENT_CONFIDENT} | 4828 |
| + pdb_resseq_real > 0 | 4828 |
| + real_aa == wt_aa（氨基酸身份一致） | 4828 |
| + ann_pos1 == pdb_resseq_real（offset=0） | 1523 |

- **严格 EXACT 位点: 1523 / 7869 = 19.35%**
- 覆盖的唯一 pair: **1216**
- 覆盖的唯一 PDB: 73
- 覆盖的 side: {'receptor': 928, 'ligand': 595}

> ⚠️ 有 2 条 pair 的 `pair_row` 超出映射后的行数（7523）—— 这些属于 V4 映射阶段被丢弃的 3 对，已排除。

### 收紧到「可做 MIPA」的条件

| 条件 | pair 数 | 涉及 PDB |
|---|---:|---:|
| 该 pair 的**全部** annotation 都严格映射 | **1140** | 92 |
| + 仅**一侧**序列变化（可判定突变链） | **987** | 82 |

> 「两侧都变化」的 pair 无法确定某条 annotation 属于哪条链，因此不入 MIPA 主实验。

### 严格 EXACT 单侧突变 pair 的 split 分布

| split | pair 数 |
|---|---:|
| train | 697 |
| val | 101 |
| test | 189 |

- **test 上可用于 MIPA 的 pair: 189** （占 test 1129 的 16.7%）

## GATE B 最终判定（严格口径）

- 可用于 MIPA 的 pair 覆盖率 = **987 / 7523 = 13.12%**

**GATE B（严格 EXACT）= `SOLVE_MAPPING_FIRST`**

> 结论：按 V5 §124，覆盖率不足 70% ⇒ **不做 full MIPA**；
> 只做**严格 EXACT 子集 pilot**，并且**必须**在相同 intersection rows 上
> 重跑 baseline 做 coverage fairness（V5 §34）。

- 耗时 10s
