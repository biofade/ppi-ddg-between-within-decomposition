# ROUND8_FIX_REPORT

对你第八轮指出的 1 BLOCKER + 4 MAJOR + 2 小问题的处理结果。

```
VALIDATION:  PASS = 63    WARNING = 2    FAIL = 0
主稿 15 页 / 补充 12 页，fatal=0，citation undefined=0，reference undefined=0
```

---

## BLOCKER —— orientation 用了 evaluation labels 却声称 label-free ✅ 已按你的方案处理

### 我先去找独立的 score 定义 —— 找不到

我找到了 **CATH-ddG 的上游仓库**（`_audit/x/cath/CATH-ddG-main`），并确认：

- 那里有一份 `skempi_v2.csv`，与我们下载的**逐列完全一致**（`ESM-IF`/`proteinMPNN`/`ddG`
  等列最大差 = 0）；
- 但仓库**不包含**生成这五个打分列的代码。`data/README.md` 只指向一个 Google Drive 数据包，
  没有任何 score-definition 文档；
- 仓库里出现的 ProteinMPNN 代码（`transfer_model.py`、`protein_mpnn_utils.py`）是
  **CATH-ddG 自己用作结构编码器**，不是用来生成表中 `proteinMPNN` 那一列的。

⇒ **我无法从独立来源确定这五列的方向。你的判断成立：orientation 现在实质上由
evaluation correlation 判定。** 原 Methods 的 "it uses no evaluation labels" 是事实错误。

### 处理：改为诚实口径（采纳你的备选方案）

我**没有**假装找到了独立证据，而是把口径改成：

| | 角色 | 依据 |
|---|---|---|
| **RAW（未改动列）** | **primary analysis** | 无需任何方向决定，纯 provenance 驱动 |
| **sign-oriented** | **exploratory sensitivity** | 明确声明方向由 evaluation correlation 判定，因此是 **label-informed** |

Methods 新增/改写的小节（原 `Orientation alignment for cross-method comparison` →
`Score orientation and its limits`）现在明确写：

> We were unable to locate, in the accompanying repository or documentation, a statement of the
> score definition or sign convention used to produce them. We therefore treat the published
> columns as *raw* and make the raw analysis the primary external analysis …
>
> **the orientation multiplier is selected using the evaluation labels**, because no
> label-independent source for it could be found. It is a single sign bit per method, it involves
> no rescaling or recalibration, and it preserves all information content; but it is a
> label-informed decision and we therefore report it as a sensitivity analysis rather than as the
> primary result. The two orientations lead to opposite conclusions about cross-axis agreement,
> and we report both.

### 结论随之调整（两个方向都报）

| | raw（primary） | oriented（sensitivity） |
|---|---|---|
| ρ / τ / p | **+0.90 / +0.80 / 0.083** | **+0.30 / +0.20 / 0.683** |

Abstract / Results / Discussion 全部改为**两向并列**，并把 final claim 收成：

> The external evidence therefore does not support a claim in either direction: with five
> methods, and with the orientation of the published scores undetermined by their own
> documentation, the panel is both underpowered and orientation-dependent, and we draw no
> directional conclusion from it.

这比原来"strongly positively aligned"或"no significant agreement"都更诚实。

---

## MAJOR 1 —— Note 7 未同步 ✅ 已重建为四面板

Note 7 现在有 **(A) raw 可变覆盖 / (A′) oriented 可变覆盖 / (B) raw 共同交集 /
(B′) oriented 共同交集** 四张表，并给出结论：

> restricting all five methods to their common intersection **leaves the raw panel aligned and
> the oriented panel unaligned, exactly as in the full uniquely mapped subsets**. The coverage
> imbalance is therefore not what drives the difference between the two conclusions.

即：**coverage sensitivity 对两个口径都稳健**，这正是你说的"好消息"。

---

## MAJOR 2 —— Methods 的 permutation p ✅ 已改

> … giving $p = 0.083$ for the raw columns and $p = 0.683$ after the exploratory
> sign-orientation.

---

## MAJOR 3 —— "no sign flip" 与 sign-oriented 矛盾 ✅ 已全部清除

主稿 Table 3 表注、Methods 的 semantics 小节、补充材料 Note 8 的 "No transformation" 段、
Table S3 表注 —— **共 4 处** "no sign flip" 全部改为：

> Raw prediction columns are ingested and preserved exactly as published, with no rescaling or
> recalibration … A sign multiplier is applied only for the orientation sensitivity described in
> Methods.

新增**检查 BG**：主稿或补充材料中再现 "no sign flip" 即 FAIL。

---

## MAJOR 4 —— Table S4 实际是 S3 ✅ 已改

主稿与 Note 8 的 `Table~S4` 全部改为 **Table S3**（并修正了 Methods 里引用
"Supplementary Table~S4 alongside each $o_m$" 的错引）。
新增**检查 BH**：`Table S4` 在补充材料中再现即 FAIL。

---

## 小问题 1 —— Figure 6 Panel B 标题 ✅ 已改

`Method-level capability profile (raw values printed)` →
**`(sign-oriented values printed)`**（图内数值确实是 oriented 的）。
新增**检查 BI**。

## 小问题 2 —— common-intersection 能量型 RMSE 被过度清空 ✅ 已恢复

此前我把 `2.3111` 也一并替换，误伤了 (B) 表。现在 (B)/(B′) 表中：

| 方法 | RMSE |
|---|---|
| DDGPred / FoldX / Flex-ddG | **2.1265 / 2.3623 / 2.1620**（真实重算值） |
| ESM-IF / proteinMPNN | `—` |

新增**检查 BJ**：三个能量型方法的 RMSE 必须存在。

---

## 符号复用 ✅ 已改

orientation multiplier 由 `μ_c` 改为 **`o_m`**（method-specific）：

$$s_m^{\mathrm{aligned}} = o_m\, s_m^{\mathrm{raw}},\qquad o_m \in \{+1,-1\}$$

现在全文 `μ_c` **只**用于 decomposition 的 target 均值，`o_m` 只用于方向乘子。
新增**检查 BK**。

---

## 验证器

本轮新增 8 项（BE–BL），并更新 2 项过期期望（AZ、F、BB）：

| 检查 | 内容 |
|---|---|
| BE | orientation 必须声明为 label-informed，不得声称 label-free |
| BF | raw 与 oriented 结果必须**同时**报告 |
| BG | 不得再出现 "no sign flip" |
| BH | 表引用必须为 S3 |
| BI | Figure 6 Panel B 标题为 sign-oriented |
| BJ | common-intersection 能量型 RMSE 保留 |
| BK | orientation 乘子用 `o_m` 而非 `μ_c` |
| BL | Note 7 已同步为 oriented 分析 |

验证器现共 **65 项**（A–BL），`PASS=63 / WARNING=2 / FAIL=0`。

---

## 仍未做

1. 键级分组重训排除 duplicate-key leakage（需重训，Note 6 已声明限制）。
2. 多数 Figure 仍为内联数组（Figure 3 完全数据驱动；Figure 6 数值现由 `o_m` 注释说明来源）。
3. **WARNING 1**：Data/code availability 提交途径需你定。
4. **WARNING 2**：作者/单位/Email/ORCID/Funding —— 仓库无真实 metadata，未编造。
5. 中文版未同步；Figure 1 未改为满页。

---

## 请你判断的一点

按你的要求，我在**找不到独立 score definition** 时选择了"raw 为 primary、oriented 为
label-informed sensitivity、两者都报"这条路。如果你认为更合适的做法是**彻底删掉
orientation 分析**（只保留 raw），或者**把 oriented 作为 primary 但在标题/摘要里就声明
其为 label-informed exploratory**，请告诉我，我可以再改。

我不打算在没有你判断的情况下擅自二选一，因为这两种处理对审稿人的呈现差别很大。
