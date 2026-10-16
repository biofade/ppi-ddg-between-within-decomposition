# V8 Decision-Alignment Report（§38/§39/§42/§43 / H3）

> **术语（§39）**：这是 **decision alignment**，**不是** external validation ——
> 同一批实验标签同时参与了预测指标与效用指标的计算。
>
> **结论先行：H3 FAIL。** 而且失败的方式比"不显著"更重要 —— 它暴露了本 benchmark
> **选择效用轴本身的区分度缺陷**。详见第 4 节。

---

## 1. 模型层面：各指标与选择效用的秩关联（Spearman）

### universe `U_CANONICAL`（n=160 模型）

| predictor | axis | effort@10 | ndcg@10 | hit@10 |
|---|---|---:|---:|---:|
| `overall_pearson` | calibration | −0.007 | +0.095 | −0.258 |
| `rmse` | calibration(neg) | −0.169 | −0.245 | +0.114 |
| `target_mean_pearson` | calibration | +0.071 | −0.043 | −0.229 |
| `macro_spearman` | ranking | +0.188 | +0.379 | −0.181 |
| `macro_kendall` | ranking | +0.159 | +0.356 | −0.184 |
| `pairwise_d05` | ranking | +0.139 | +0.316 | −0.029 |
| `C_W` | ranking | **−0.214** | −0.051 | **−0.389** |

### universe `U_MAS`（n=42 模型）

| predictor | axis | effort@10 | ndcg@10 | hit@10 |
|---|---|---:|---:|---:|
| `overall_pearson` | calibration | +0.892 | +0.867 | +0.259 |
| `rmse` | calibration(neg) | −0.907 | −0.862 | −0.209 |
| `target_mean_pearson` | calibration | +0.869 | +0.875 | +0.340 |
| `macro_spearman` | ranking | +0.729 | +0.822 | +0.452 |
| `macro_kendall` | ranking | +0.601 | +0.728 | +0.606 |
| `pairwise_d05` | ranking | +0.771 | +0.882 | +0.495 |
| `C_W` | ranking | +0.797 | +0.707 | +0.230 |

---

## 2. H3 初步对照（聚合）——看似支持

| universe | n | 校准轴中位 ρ | **排序轴中位 ρ** | 传统全局中位 ρ | 排序轴更优? |
|---|---:|---:|---:|---:|:--:|
| `U_CANONICAL` | 160 | −0.043 | **+0.055** | −0.088 | ✅ |
| `U_MAS` | 42 | +0.340 | **+0.717** | +0.025 | ✅ |

**但这个聚合掩盖了矛盾**：`C_W` 本身属"排序轴"，却对 Eff@10 是**负相关**（−0.214）；
`Hit@10` 对**所有**预测指标都是负相关或近零。

---

## 3. 稳健性检验（分 seed）——H3 崩溃

`h3_robustness.csv` 逐 seed 重算（U_CANONICAL）：

| subset | axis | effort@10 | ndcg@10 | hit@10 |
|---|---|---:|---:|---:|
| all | calibration | −0.007 | −0.043 | −0.229 |
| all | **ranking** | **+0.159** | **+0.356** | −0.181 |
| seed=42 | calibration | +0.005 | +0.063 | −0.114 |
| seed=42 | ranking | +0.042 | +0.211 | −0.431 |
| seed=43 | calibration | **+0.220** | **+0.213** | −0.045 |
| seed=43 | ranking | +0.094 | +0.316 | −0.316 |
| seed=44 | calibration | +0.019 | −0.120 | −0.255 |
| seed=44 | **ranking** | **−0.115** | +0.091 | −0.349 |

**⇒ 在 seed 43 上「校准轴」反而优于「排序轴」；在 seed 44 上排序轴对 Eff@10 变成负值。
跨 seed 符号不稳定 ⇒ 这个"优势"是噪声，不是能力。**

### `Hit@10` 判定为**退化指标**

- `U_MAS`：`hit@10` 均值 0.6577、**sd 仅 0.0149**、只有 9 个不同取值 ⇒ 几乎是常数，无区分度
- `U_CANONICAL`：sd 0.1340 稍有区分度，但对其余指标的负相关说明它的排序与一切相反

⇒ **`Hit@10` 不应作为效用指标**（本报告保留它仅为透明）。

---

## 4. 关键诊断：选择效用轴存在**发现地板**（discovery floor）

R4 的对比本来会显得自相矛盾：

```
按 Overall Pearson 选 -> Eff@10 regret = 0.0141   （近乎最优）
按 Macro Spearman 选  -> Eff@10 regret = 0.6957   （看似很差）
而 模型层面 ρ(Overall r, Eff@10) = −0.007        （看似无关）
```

`12c_selection_diagnostic.py` 给出了解释 —— **`Efficiency@10` 的值域被压缩**：

### U_CANONICAL（160 模型）`Efficiency@10` 分布

```
min=0.6481  p10=0.8036  median=0.8179  p90=0.8306  max=0.8547
range=0.2066   IQR=0.0143   sd=0.0405
```

**中位 0.818，80% 的模型只分布在 p10–p90 = 0.8036–0.8306（宽 0.027）。**

换算成"标准差异"：

| 选择规则 | Eff@10 绝对差 | normalized regret | **标准差异** |
|---|---:|---:|---:|
| `overall_pearson` | 0.0029 | 0.0141 | **0.07 sd** |
| `macro_spearman` | 0.1437 | 0.6957 | **3.54 sd** |

⇒ **normalized regret 被"最小值 0.648"这一个离群模型放大了 50 倍。**
`macro_spearman` 选的模型实际上距离最优仅 0.144 绝对差，落在整体 sd 的 3.5 倍内，
并非"灾难性"选择。

**更重要的含义**：`Macro ρ` 在 160 个模型间跨越 **0.40**（−0.063 → +0.339），
而 `Eff@10` 只跨越 **0.21**，IQR 仅 **0.014**。
**排序能力的巨大差异几乎不改变 Top-K 效用** —— 因为
`Efficiency@10` 是"选中 gain 之和 / 理想前 10 gain 之和"，
当 n≈20–30 且 utility 尺度由名次定义（gain=1..n）时，
**任何还算合理的模型抓到前 10 名中的一大半**，该比值就被压在高位。

### 结论

> **本 benchmark 的 AXIS C（选择效用）在 U_CANONICAL 上区分度不足；
> 因此它无法可靠地支持或反驳 H3。**
> 这不是"Overall Pearson 意外地好"，而是**效用指标本身太钝**。

---

## 5. §42 Model-selection regret（含 §43 oracle 说明）

> oracle-best 为**事后描述性上界**，**不可部署**（§43）。

| universe | selection_rule | selected_model | regret(E@10) | regret(NDCG@10) |
|---|---|---|---:|---:|
| U_CANONICAL | overall_pearson | `L1_random_A1_base_set` | 0.0141 | 0.0142 |
| U_CANONICAL | macro_spearman | `R1_Base_F0_s42` | 0.6957 | 0.5633 |
| U_CANONICAL | C_W | `R0_Base_s42` | 0.6503 | 0.4163 |
| U_CANONICAL | oracle-best | `L1_random_A3_base_local` | — | — |
| U_MAS | overall_pearson | `L1_random_M1_s42` | 0.0893 | 0.1170 |
| U_MAS | macro_spearman | `L1_random_M0_s43` | 0.0679 | 0.0000 |
| U_MAS | oracle-best | `L1_random_M1_s43` | — | — |

**⇒ 依据第 4 节的诊断，这些 regret 数字受"发现地板"污染，
不应被解读为"Overall Pearson 选择了更好的模型"。**

---

## 6. 对 H3 的判定

```text
H3 Selection alignment = FAIL
```

依据：
1. **分 seed 符号不稳定**：seed 43 上校准轴反超，seed 44 上排序轴为负
   ⇒ 聚合上的"排序轴更优"是噪声；
2. **排序轴内部自相矛盾**：`C_W` 对 Eff@10 为 −0.214，与 `macro_spearman` 的 +0.188 方向相反；
3. **`Hit@10` 退化**（U_MAS sd=0.0149）且与所有指标负相关；
4. **选择效用轴存在发现地板**（IQR=0.0143）⇒ 该轴本身无法支撑这一假设的检验。

**⇒ 依 §69「不要强行给正面结论」，H3 判 FAIL，并记录失败原因是"被检验的效用轴区分度不足"，
而非单纯"两类指标无差别"。**
