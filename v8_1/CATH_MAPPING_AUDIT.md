# CATH-ddG 映射歧义审计（§13.3）

* 归一化键定义：`(PDB4, residue_position, mutant_AA) — chain and WT dropped`
  （**丢弃链与 WT 残基** —— 这正是歧义来源）
* 外部表行数 = 7085；重复归一化键 = **570**

## canonical 行分级

| 类别 | 行数 |
|---|---:|
| `EXACT_UNIQUE` | 1781 |
| `AMBIGUOUS_KEY` | 2048 |
| `UNMATCHED` | 644 |

## 规则

* **主匹配评估只允许 `EXACT_UNIQUE`**；
* `AMBIGUOUS_KEY` **不得**被静默解析到任一行；
* 官方 held-out CATH 测试集自带无歧义行，**单独评估、不做 canonical 重映射**。
