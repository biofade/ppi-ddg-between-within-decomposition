# V5.0 — 映射的独立 PDB 验证

> 目的: `04` 解出的 offset 范围 [−425,+645] 对真实 PDB 编号不合理，
> 必须用本地真实结构文件独立验证，不能凭 alignment 自证。

- `mutation_sites.csv` 位点数: **11651**
- 主实验 pair 数: **5574**

- 本地 PDB 文件: **3031**（目录 `F:\DSH_Galatea\Affinity_Intelligence\data\pdb`）

## E1: `seq_index0` 对应的真实 PDB 残基编号 vs 推断的 `pdb_resseq`

- 参与验证的位点（限主实验 pair）: **7869**

| 检查结果 | 条数 |
|---|---:|
| 可在 PDB 中定位（checked） | 7134 |
| **位置+残基都一致（MATCH）** | **568** |
| 不一致（POS_MISMATCH / AA_DIFF） | 6566 |
| seq_index0 越界 | 735 |
| PDB 文件缺失 | 0 |

- **E1 通过率 = 568 / 7134 = 7.96%**

## E3: 同一 pdb 内 offset 是否自洽

- 主实验 pair 覆盖 pdb 数: 226
- 单 pdb 内 offset **唯一**的 pdb: **127** (56.2%)
- 单 pdb 内 offset **冲突**的 pdb: **99**

- |offset| > 100 的位点: **867**（占主实验位点 11.0%）→ 这些极可能是**假匹配**

## 结论与判据（决定主实验是否可用）

**MAPPING_NOT_VALIDATED** — 映射未获 PDB 独立验证 → **不得**基于该映射做 MIPA 结论；应先修 mapping（V5 §124 <70% 分支）。

- 耗时 9s
