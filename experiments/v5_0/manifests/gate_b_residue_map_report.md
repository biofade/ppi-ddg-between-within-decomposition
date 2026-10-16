# V5.0 GATE B（修订版）— 序列引导的 PDB 残基映射

> **前置负结果**: `04`+`05` 的映射方案未通过独立验证（E1 通过率 7.96%），
> 根因是**从未确定 dataset 序列 ↔ 具体 PDB 链**的对应（错取最长链）。
> 本脚本按 V5 §25「PDB sequence ↔ dataset sequence global alignment」重做。

- 待重映射位点: **7869**（来自 5574 个主实验 pair）
- 本地 PDB: 3031 个

- 涉及 PDB: **226** → 解析中 ...
   ... 50/226
   ... 100/226
   ... 150/226
   ... 200/226
- 解析完成；有 ATOM 蛋白链的 PDB: **226 / 226**

  ... 2000/7869
  ... 4000/7869
  ... 6000/7869

## 1. 修订后逐位点状态

| validated_status | 位点数 | 占比 |
|---|---:|---:|
| `EXACT` | 3798 | 48.27% |
| `AA_MISMATCH` | 1737 | 22.07% |
| `LOW_IDENTITY_CHAIN` | 1102 | 14.00% |
| `ALIGNMENT_CONFIDENT` | 1030 | 13.09% |
| `INDEX_NOT_IN_ALIGNMENT` | 150 | 1.91% |
| `GAP_AT_MUTATION_SITE` | 52 | 0.66% |

## 2. 链 identity 分布（验证「链对应」是否找对）

- identity =1.000: 4962 (63.06%)
- identity 0.95–0.999: 1417 (18.01%)
- identity 0.80–0.95: 1285 (16.33%)
- identity <0.80: 205 (2.61%)

## 3. 独立验证：`pdb_resseq_real` 处真实残基是否 == annotation WT 残基

- 可检查位点: **7719**；残基一致: **5930** （**76.82%**）

## 4. GATE B 判定（修订版）

- 可用位点（EXACT + ALIGNMENT_CONFIDENT）: **4828 / 7869 = 61.35%**
- 覆盖的主实验 pair: **3587**

**GATE B（修订版）= `SOLVE_MAPPING_FIRST`**

- 耗时 14s
