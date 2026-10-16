# V5.0 GATE B — mutation ↔ 残基映射（正式版，含编号 offset 求解）

> 机制: pair 内 WT/mutant 序列 global alignment → 残基自洽性 + **同链唯一 offset**。
> offset 的引入依据: V5 §25 明确列出 `residue numbering offset` 为必须处理的情形；
> `03b` 实测朴素 1-based 位置仅 73.37% 命中，而失配样例呈**恒定偏移**（如 +10），
> 证明这是编号差异而非 construct 不一致。

- forward pair: **7526**

  ... 2000/7526  EXACT=1464 CONF=683
  ... 4000/7526  EXACT=3057 CONF=970
  ... 6000/7526  EXACT=4448 CONF=1808

## 1. 配对级状态

| pair_status | 对数 | 占比 |
|---|---:|---:|
| MAPPED | 5574 | 74.06% |
| FAILED | 1752 | 23.28% |
| PARTIAL | 200 | 2.66% |

## 2. 位点级状态

| mapping_status | 位点数 | 占比 |
|---|---:|---:|
| EXACT | 5506 | 47.26% |
| AMBIGUOUS | 3565 | 30.60% |
| ALIGNMENT_CONFIDENT | 2580 | 22.14% |

### 编号 offset 分布（验证「恒定偏移」假设）

- offset = 0（直接匹配）: **7249** (62.2%)
- offset ≠ 0（编号偏移，已解）: **4402** (37.8%)
- offset 范围: [-425, 645]，非零 offset 的常见取值: 过多

## 3. 失败原因（如实分类，不丢弃）

| reason | 条数 |
|---|---:|
| NO_CONSISTENT_SUBSTITUTION | 1672 |

## 4. GATE B 判定（V5 §124）

| 项 | 值 |
|---|---:|
| 总 pair | 7526 |
| **主实验可用 pair**（EXACT+CONFIDENT 且该 pair 全部 annotation 映射成功） | **5574** |
| **覆盖率** | **74.06%** |
| EXACT-only 子集（敏感性分析用） | 3824 |

**GATE B = `SUBSET_PILOT`**（阈值: ≥90% full / 70–90% subset pilot / <70% 先解决）

- 耗时 33s
