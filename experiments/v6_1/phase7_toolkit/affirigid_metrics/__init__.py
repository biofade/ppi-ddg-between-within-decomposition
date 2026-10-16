# -*- coding: utf-8 -*-
"""affirigid_metrics — 面向 PPI ΔΔG 预测的「靶点内效用」评估工具包.

为什么需要这个包
----------------
PPI ΔΔG 预测论文普遍只报 **全局 Pearson/Spearman**。但在真实亲和力成熟实验中，
决策是：**在同一个靶点内，从候选突变体里挑 top-K 去做湿实验**。
全局指标会掩盖组内能力的缺失 —— 我们实测同一模型在随机划分下
`Overall Pearson` 可达 **+0.84**，而 `Macro per-PDB Spearman` 仅 **−0.03**。

本包提供两个互补的指标族：

1. **Macro per-structure Spearman/Pearson**（`topk.per_structure_*`）
   每个复合物内单独算相关性再 macro 平均 → 直接度量组内排序能力。
2. **Efficiency@K**（`topk.efficiency_at_k`）
   `(Random − Model) / (Random − Oracle)`，在每复合物内选 top-K 后 macro 平均。
   0 = 与随机无异，1 = 达到 oracle，**负值 = 比随机更差**。

以及配套的 **cluster bootstrap**（`bootstrap.*`）：把「复合物」当作抽样单位，
避免把同一靶点的多个突变体误当独立样本而低估置信区间。

快速开始
--------
```python
import pandas as pd
from affirigid_metrics import topk, bootstrap

df = pd.read_csv("predictions.csv")
# 必需列：complex_id, true_ddG, pred_ddG
m = topk.topk_metrics(df.true_ddG, df.pred_ddG, df.complex_id)
print(m["macro_pdb_spearman"], m["efficiency_at_10"], m["ef_1pct"])

ci = bootstrap.cluster_bootstrap_ci(
    topk.per_structure_spearman,
    df.true_ddG.values, df.pred_ddG.values, df.complex_id.values,
    n_boot=5000, min_groups=5)
print(ci["point"], ci["ci_low"], ci["ci_high"])
```

**符号约定（重要）**：`ΔΔG < 0` 表示**增强结合（beneficial）**。
若你的数据使用相反约定，请先取负号再传入。
"""

from . import bootstrap, topk  # noqa: F401

__version__ = "0.1.0"
__all__ = ["topk", "bootstrap"]
