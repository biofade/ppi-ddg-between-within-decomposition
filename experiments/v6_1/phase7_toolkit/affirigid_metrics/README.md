# affirigid_metrics — 面向 PPI ΔΔG 的「靶点内效用」评估工具包

## 它解决什么问题

PPI ΔΔG 预测论文普遍只报 **全局 Pearson / Spearman**。但真实实验决策是：

> 在**同一个靶点**内，从候选突变体里挑 **top-K** 去做湿实验。

全局指标会把「跨复合物区分能力」和「复合物内排序能力」混在一起。
我们实测（同一模型、同一特征、同一超参，只换划分）：

| 划分 | Overall Pearson | Macro per-PDB Spearman |
|---|---:|---:|
| mutation-level random | **+0.8446** | **−0.0257** |
| complex-level + ≤40% seq-id + family | 见 `experiments/v6_1/` | 见 `experiments/v6_1/` |

**随机划分能造出 0.84 的全局相关，却完全不产生靶点内排序能力。**
本工具包提供能暴露这一点的指标。

---

## 输入格式

一个长表（pandas DataFrame 或 numpy 数组），最少需要三列：

| 列名（示例） | 含义 |
|---|---|
| `complex_id` | **复合物 / 靶点 ID**（分组键；允许字符串） |
| `true_ddG` | 真实 ΔΔG（kcal/mol） |
| `pred_ddG` | 预测 ΔΔG（kcal/mol） |

可选：`model`、`mutation`、`dataset_name`（用于多模型/多数据集批量评估）。

### ⚠️ 符号约定（必读）

**`ΔΔG < 0` = 增强结合（beneficial）**，`ΔΔG > 0` = 削弱结合。
若你的数据使用相反约定（例如把 `−ΔΔG` 存成正向 score），
**请先取负号再传入**，否则所有 Hit@K / Efficiency 都会反号。

验证方式：随机抽一条「已知增强结合」的突变，检查其 `true_ddG < 0`。

---

## 输出格式

`topk_metrics(...)` 返回一个 dict：

| 键 | 含义 |
|---|---|
| `base_rate_best_true` | 真值最负前 10%（「真强结合体」）的占比 = 随机命中期望 |
| `p_at_5` / `p_at_10` / `p_at_20` | Precision@k（预测 top-k 落在真强结合体内的比例） |
| `ef_1pct` / `ef_5pct` / `ef_10pct` | Enrichment Factor（相对随机的富集倍数） |
| `hit5_beneficial` / `hit10_beneficial` | 预测 top-k 中 `ΔΔG < 0` 的比例 |
| `hit10_ddg_m05` / `hit10_ddg_m10` | 预测 top-10 中 `ΔΔG < −0.5` / `< −1.0` 的比例 |
| `sel10_mean_ddg` | 预测 top-10 的真实 ΔΔG 均值（**越负越好**） |
| `oracle10_mean_ddg` | 全池 oracle top-10 的均值（上界参考） |
| **`efficiency_at_10`** | **`(Random − Model) / (Random − Oracle)`，每复合物内计算后 macro 平均** |
| `macro_pdb_spearman` | 每复合物内 Spearman 的 macro 平均（**靶点内排序能力**） |
| `macro_pdb_pearson` | 同上，Pearson 版 |

`cluster_bootstrap_ci(...)` 返回：

| 键 | 含义 |
|---|---|
| `point` | 点估计 |
| `ci_low` / `ci_high` | 95% CI（cluster bootstrap） |
| `sd` | bootstrap 分布标准差 |
| `n_clusters` | 实际使用的复合物数 |
| `n_valid` / `n_skipped` | 有效 / 被跳过的迭代数（`min_groups` 约束导致） |
| `p_gt_zero` | bootstrap 分布中 > 0 的比例（单侧证据强度） |

---

## 安装与依赖

无第三方安装步骤，直接把 `affirigid_metrics/` 目录放到你的 `PYTHONPATH` 即可。

```
numpy>=1.20
scipy>=1.7      # 仅 paired_wilcoxon 需要
pandas>=1.3     # 仅示例脚本需要
```

---

## 使用示例

见 `examples/example_usage.py`（可直接运行，用合成数据自检）。
