# V8 跨宇宙一致性（PHASE 13 收尾）

> **纪律**：三个外部宇宙的 target 定义、行数、模型数都不同
> ⇒ 依 §21 **只比较方法排序，不比较绝对数值**。

---

## 1. 三个外部宇宙

| 宇宙 | 来源 | target 定义 | 行数 | 方法数 |
|---|---|---|---:|---:|
| `U_HERMES_CANONICAL` | HERMES 团队复算官方预测 | `target_id_canonical`（手工家族） | 3822–3842 | 2 |
| `U_EXT_CANONICAL` | CATH-ddG 表 | `target_id_canonical`（单点匹配） | 3756–3819 | 5 |
| `U_EXT_HELDOUT_CATH` | CATH-ddG 官方 held-out | `complex_PPI` | 813 | 5 |

## 2. 方法 × 宇宙 矩阵

**Overall Pearson**

| method | U_EXT_CANONICAL | U_EXT_HELDOUT_CATH | U_HERMES_CANONICAL |
|---|---:|---:|---:|
| `DDGPred` | +0.3357 | +0.6128 | — |
| `Flex-ddG` | +0.2619 | +0.6022 | — |
| `ddG_FoldX` | +0.2727 | +0.4603 | — |
| `proteinMPNN` | +0.1674 | +0.1837 | **+0.1591** |
| `ESM-IF` | −0.1014 | +0.1421 | — |
| `HERMES` | — | — | **+0.1035** |

**Macro per-target Spearman**

| method | U_EXT_CANONICAL | U_EXT_HELDOUT_CATH | U_HERMES_CANONICAL |
|---|---:|---:|---:|
| `ddG_FoldX` | +0.1502 | +0.3584 | — |
| `Flex-ddG` | +0.1320 | +0.4159 | — |
| `DDGPred` | +0.1047 | +0.3262 | — |
| `proteinMPNN` | +0.0919 | +0.2227 | +0.0841 |
| **`ESM-IF`** | **−0.0346** | **−0.0329** | — |
| `HERMES` | — | — | +0.0678 |

---

## 3. ★ 跨宇宙排序一致性（两个 CATH-ddG 宇宙，n=5）

```
Overall r  :  rho = +0.9000   tau = +0.8000
Macro rho  :  rho = +0.9000   tau = +0.8000
```

**⇒ 在同一个外部方法群上，两个**不同 universe**（不同行集合、不同 target 定义、
813 vs 3819 行）给出**几乎相同的相对排序**（ρ=+0.90）。**

这是 V8 的**度量框架在外部数据上的首个稳健性证据**：
方法的相对优劣**不因宇宙构建方式而翻转**（与内部模型的 PHASE 13 发现形成对比 ——
内部模型池的秩相关**会**因 provenance 翻转）。

> ⚠️ **不得过度解读**：n=5 时 ρ=0.90 的 95% 区间仍宽（约 [+0.1, +0.99]）；
> 且这两个宇宙**同源**（都来自 CATH-ddG 表），不是独立复现。

---

## 4. 跨宇宙一致的**稳健**观察

### 4.1 `proteinMPNN` —— 出现在 **3 个宇宙**

| 宇宙 | Overall r | Macro ρ |
|---|---:|---:|
| `U_HERMES_CANONICAL` | +0.1591 | +0.0841 |
| `U_EXT_CANONICAL` | +0.1674 | +0.0919 |
| `U_EXT_HELDOUT_CATH` | +0.1837 | +0.2227 |

**⇒ 三个来源一致给出 Overall r ≈ +0.16 ~ +0.18**（极差仅 0.025）。
这是本项目**最稳健的单一外部事实**。

### 4.2 `ESM-IF` —— **两个宇宙上 Macro ρ 皆为负且几乎相同**

```
U_EXT_CANONICAL    Macro rho = -0.0346
U_EXT_HELDOUT_CATH Macro rho = -0.0329     ← 差 0.0017
```

**⇒ `ESM-IF` 在靶内排序上不携带正信息的结论，跨宇宙高度一致。**
（这比单宇宙结论可信得多。）

---

## 5. 对 §68 的最终影响（未变）

```text
EXTERNAL_MODEL_REPRODUCTION = NOT_TESTABLE
  独立外部方法族 = 6（HERMES, proteinMPNN, ESM-IF, ddG_FoldX, Flex-ddG, DDGPred）
  < 8 门槛
  且跨宇宙可检验的方法仅 5 个（HERMES 只出现在 1 个宇宙）
  ⇒ 无法检验 V8 的核心主张（校准轴 vs 排序轴的能力分化）在外部是否复现
```

**⇒ 依 §69：**

```text
BENCHMARK_VALIDATION    = PARTIALLY_SUPPORTED
BENCHMARK_UTILITY_LEVEL = MORE_DISCRIMINATIVE
```

**本轮**不提升**判定**，但**提升了外部证据的稳健性**（跨宇宙 ρ=+0.90 的排序一致 +
proteinMPNN / ESM-IF 的跨宇宙一致读数）。

---

## 6. 复现命令

```powershell
python -W ignore -u scripts/v8/42_cross_universe_consistency.py
```
产出：`external/cross_universe_external.csv`、
`external/cross_universe_agreement.csv`、`external/cross_universe_meta.json`
