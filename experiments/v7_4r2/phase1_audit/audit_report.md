# V7.4R2 PHASE 1 — 结构准备审计（纯 Python，无需 OpenMM）

- 缺口安全半径 `GAP_SAFE = 3` 残基（位点距最近缺口需 ≥ 3）

- 起点 U_RELAX_STRICT：**3253** 行 / 162 targets

## 1.1 CANONICAL_WT_STRICT（§6 严格结构排除）

| 项 | 值 |
|---|---:|
| 起点 rows | 3253 |
| **CANONICAL_WT_STRICT rows** | **3214** |
| targets | **161** |
| families | **85** |
| n>=1 targets | **161** |
| n>=2 targets | **144** |
| n>=3 targets | **124** |
| n>=5 targets | **100** |
| n>=10 targets | **62** |

### 排除原因（§6）

| reason | rows |
|---|---:|
| `site_near_chain_gap` | 39 |
| **合计排除** | **39** |

## 1.2 逐 target 缺口统计（PDBFixer 会在此插入缺失残基）

- target 数：**162**
- **含 ≥1 链内缺口的 target：110 / 162**
- 缺口总数：**943**
- 缺口分布：median=2 max=65

> V7.4R2 **禁止**填补缺失残基（§5）。缺口附近的位点已按 `GAP_SAFE=3` 安全半径排除。

