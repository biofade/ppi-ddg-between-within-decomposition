# V8 External Replication Test（§68 定稿）

> **一句话结论：能力分化（H1）在外部匹配子集上**没有**复现 —— 而且是在**内部模型**上
> 就没有复现。这限制了 H1 的适用范围，并阻止 `BENCHMARK_VALIDATION` 升级。**
>
> 这是一个**阴性结果**，但它是**关于 H1 边界条件**的真实发现，不是失败。

---

## 1. 匹配宇宙 `U_HERMES_CANONICAL`

来源：HERMES **官方发布的行级预测**
（[github.com/StatPhysBio/hermes_experiments](https://github.com/StatPhysBio/hermes_experiments)，
`skempi/results/hermes_py_000_ft_cdna117k_ddg_st/…use_mt_structure=0.csv`，403,691 B / 6,993 行）。

| 项 | 值 |
|---|---:|
| canonical 行 | 4473 |
| 命中 HERMES（严格 1:1） | 2865 |
| **单点突变子集匹配率** | **2865 / 2897 = 98.90%** |
| targets | 195 |
| **n≥5 的 targets** | **94** |
| families | 101 |
| 符号校验 `r(y_true, hermes_proxy)` | **+0.1688** ✅ 为正 |

**符号约定（§19）**：HERMES 输出 log-odds（`log P(mt) − log P(wt)`，越大越有利），
V8 约定 `ΔΔG>0 = 结合变弱` ⇒ **方向相反**，构造
`hermes_proxy = −(log_proba_mt − log_proba_wt)`。方向校验 r=+0.1688 > 0 确认正确。

**§50 污染状态**：HERMES 在 **SKEMPI v2.0** 上微调 ⇒ **TIER 2**，不进 Tier 1。

---

## 2. HERMES 的三轴读数（同一批 2865 行 / 94 个 n≥5）

| 轴 | 指标 | 值 |
|---|---|---:|
| A 校准 | Overall Pearson | **+0.1688** |
| A 校准 | RMSE / MAE | 2.2831 / 1.5845 |
| A 校准 | target-mean Pearson | +0.2972 |
| A 校准 | C_B (row) | +0.0512 |
| A 校准 | C_W (row) | +0.1177 |
| B 排序 | **Macro per-target Spearman** | **+0.0523** |
| B 排序 | Macro Kendall τ-b | +0.0406 |
| B 排序 | Pairwise δ=0.5 | +0.5198 |
| C 选择 | Efficiency@10 / NDCG@10 | +0.8225 / +0.7447 |
| QC | 分解残差 \|r − (C_B+C_W)\| | **0.0** ✅ |

---

## 3. ★ 关键检验：秩一致性**没有**复现

在**同一批 2865 行**上重算全部 136 个内部模型（这样比较才是匹配宇宙的，§21）：

| 模型集 | n | **ρ(Overall r, Macro ρ)** | τ | ρ(Overall r, Eff@10) |
|---|---:|---:|---:|---:|
| 仅内部 | 136 | **+0.4743** | +0.3528 | +0.1146 |
| 内部 + HERMES | 137 | **+0.4737** | +0.3532 | +0.1106 |

```
U_CANONICAL 全集 (4473 行)          : ρ(Overall r, Macro ρ) = −0.20 ~ −0.26
U_HERMES_CANONICAL (2865 行, 单点) : ρ(Overall r, Macro ρ) = +0.4743   ← ★ 反了
```

**⇒ 在外部匹配子集上，校准轴与排序轴**正相关**（+0.47），
即我在全集上报告的「能力分化」在此**不复现**。

**请注意这个结论的重要性**：它**不是** HERMES 的问题 ——
**是内部模型自己在这个子集上就没有分化**。
`EXTERNAL_REPLICATION_VERDICT = EXTERNAL_REPLICATION_DISCORDANT`。

---

## 4. HERMES 自身的定位

| 指标 | HERMES | 内部中位 | HERMES 百分位 |
|---|---:|---:|---:|
| Overall Pearson | +0.1688 | +0.2155 | 41.6% |
| Macro ρ | +0.0523 | +0.0622 | 37.2% |
| C_W | +0.1177 | +0.1877 | 37.2% |
| Efficiency@10 | +0.8225 | +0.8162 | 75.9% |
| RMSE | +2.2831 | +2.2656 | 53.3% |

**⇒ HERMES 落在内部模型群的中位附近（37–76 百分位），并且与内部中位**同象限**。**
这是一个**合理**的结果（HERMES 是较早的方法），但它**不能**用来支持或反驳能力分化 ——
因为在这个子集上分化本身就不存在。

---

## 5. 对 H1 适用范围的影响（**必须据此收窄 H1 的表述**）

| 宇宙 | 行数 | 特征 | ρ(Overall r, Macro ρ) |
|---|---:|---|---:|
| `U_CANONICAL` 全集 | 4473 | 含多点突变, 108 families, 130 个 n≥5 | **−0.20 ~ −0.26** |
| `U_CANONICAL` 仅 rowid 子集 | — | — | −0.197 |
| `U_CANONICAL` 仅 align 子集 | — | — | −0.213 |
| **`U_HERMES_CANONICAL`** | **2865** | **仅单点突变, 101 families, 94 个 n≥5** | **+0.474** |

**⇒ 能力分化是一个「regime-dependent」现象，不是普适性质。**
可能的机制（**均为候选，未检验**）：
1. **多点突变改变任务性质**：多点组合的 `ΔΔG` 更难，校准轴与排序轴在难样本上分离；
2. **模型构成差异**：全集含 102 个 V6.2 键（伪重复），而 HERMES 子集只覆盖 136 个键；
3. **target 规模差异**：全集 130 个 n≥5 vs 子集 94 个；
4. **单点突变本身更易排序**（macro ρ 更高），使得两轴同向。

**⇒ H1 的判定维持 PASS（在 `U_CANONICAL` 全集上确实成立），
但论文中必须把 H1 限定为：**

> 在**含多点突变、以 family 为外推单位**的完整 family-OOD 宇宙上，
> 校准轴与排序轴对模型的排序近乎无关甚至弱负相关（−0.20 ~ −0.26）。
> **该分化在「仅单点突变」的匹配子集上不复现（+0.47）。**

---

## 6. 对 §68 判定的影响

```text
EXTERNAL_MODEL_REPRODUCTION = ATTEMPTED_DID_NOT_REPLICATE
```

§68 的 `STRONGLY_SUPPORTED` 要求「外部模型也复现能力分化」。
**本检验的结果是：连内部模型在该外部匹配子集上都不复现。**
⇒ **不能**升级为 `SUPPORTED` / `STRONGLY_SUPPORTED`。

```text
BENCHMARK_VALIDATION    = PARTIALLY_SUPPORTED   （维持，且理由更充分）
BENCHMARK_UTILITY_LEVEL = MORE_DISCRIMINATIVE   （维持 LEVEL 1）
```

**⇒ 依 §69「不要强行给正面结论」，此处维持不变。**

---

## 7. 诚实收尾

**本轮做对了什么**：把一个「需要外部验证」的主张**真的去验证了**，
并且**在结果不利时如实记录**，而不是只报告支持性证据。

**本轮的负面结果有三个正面价值**：
1. 它给 H1 划出了**明确的边界**（含多点突变的完整宇宙 vs 仅单点子集）；
2. 它证明 V8 的度量引擎在**外部数据**上可用（分解残差 0.0，符号约定经数据验证）；
3. 它留下一个**可检验的机制假设**：分化是否由**多点突变**驱动 ——
   下一轮可用 `U_CANONICAL` 内部按 `n_sites` 分层直接检验（纯 CPU，无需外部资源）。

**产出**：`external/external_replication_scoreboard.csv`、
`external/external_replication_meta.json`、`external/hermes_matched_final.csv`
