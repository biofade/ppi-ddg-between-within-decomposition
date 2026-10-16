# RECOMPUTE_FEASIBILITY

键级分组重训的可行性评估 **及已执行的修复结果**（用户选择方案 B）。

---

## 一、结论摘要

| 问题 | 结论 |
|---|---|
| L1_random / L2_complex / L4_family 键级泄漏 | **0**（划分构造性保证，无需重训） |
| L3_seq40 键级泄漏 | **存在**：1,295 / 1,328 / 1,280 个 key 跨 fold（seed 42/43/44） |
| 代码与冻结 artifact 是否一致 | **不一致**（脚本 12:52 改动，artifact 11:50 生成） |
| 已执行修复 | ✅ 加入并查集合并 obs_key/state_id，重跑 L3，传播到全部下游 |
| 头条数字是否变化 | **未变**（Random 0.6844、Family-OOD 0.1508、Δr −0.5336、归因 81.13%） |
| 唯一变化的已报告数字 | seq40 那一行（原为 complex 的副本） |

---

## 二、键级泄漏的实测证据

直接执行 `scripts/v6_2/02_split4_ddg.py` 的 `load_all()` / `level_units()` / `make_folds()`
（不重训），重建每个层级的 fold 归属并检验「同一 `obs_key` 是否跨 fold」：

**修复前**

| level | units | 跨 fold 的 obs_key (seed 42/43/44) |
|---|---:|---|
| L1_random | 4,592 | 0 / 0 / 0 |
| L2_complex | 230 | 0 / 0 / 0 |
| **L3_seq40** | **1,144** | **1,295 / 1,328 / 1,280** |
| L4_family | 112 | 0 / 0 / 0 |

**根因**：`level_units()` 对 L1 用并查集显式合并 `obs_key` 与 `state_id`；L2 用 `pdb`、
L4 用 `family`，二者都比 `state_id` 粗，因此自动满足；只有 L3 走
`_seq40_cluster_units()`，该函数**按行做 k-mer Jaccard 贪心聚类，从未与 `obs_key` 合并**。

**修复后**

| level | units | 跨 fold 的 obs_key |
|---|---:|---|
| L1_random | 4,592 | 0 |
| L2_complex | 230 | 0 |
| **L3_seq40** | **301** | **0** |
| L4_family | 112 | 0 |

修复方式与 L1 同一手法：对 (cluster, obs_key, state_id) 做并查集合并。
L3 现在是真正独立的层级（与 L2 的行级 fold 一致率约 20–27%）。

---

## 三、代码/artifact 不一致的发现

时间戳血缘：

```
2026-10-14 10:49   split4_predictions.preL3fix.npz   （L3 与 L2 不同）
2026-10-14 11:50   split4_predictions.npz            （L3 与 L2 逐元素完全相同）
2026-10-14 12:52   scripts/v6_2/02_split4_ddg.py     ← 脚本此后被改动
```

冻结 artifact 中 L3 与 L2 的预测 **max|Δ| = 0**，与论文
「The complex-level and seq40 partitions produced identical assignments」一致；
但当前脚本重跑会得到 1,144 units（≠ 230），因此 **L3 那一档曾无法用代码复现**。
本轮的修复同时消除了这个不一致（现在代码与 artifact 都是 301 units 的独立 L3）。

---

## 四、修复后的数值与影响面

重跑命令：
```
python scripts/v6_2/02_split4_ddg.py --levels L3_seq40 --seeds 42 43 44 --trees 800
```
耗时 **33.8 分钟**（CPU）。

| 档 | r | C_B | C_W | macro ρ | RMSE |
|---|---:|---:|---:|---:|---:|
| Random | 0.6844 | +0.4275 | +0.2591 | 0.0707 | 1.7678 |
| Complex | 0.1931 | +0.0233 | +0.1685 | 0.0497 | 2.6412 |
| **seq40（新）** | **0.1961** | **+0.0239** | **+0.1721** | **0.0631** | **2.6911** |
| Family-OOD | 0.1508 | −0.0072 | +0.1580 | 0.0377 | 2.6952 |

**未变的量**：Random 档、Family-OOD 档全部指标；Δr = −0.5336、ΔC_B = −0.4347、
ΔC_W = −0.1011、between-target 归因 **81.13%**、paired 归因、bootstrap CI。

**已同步的下游**：
- `experiments/v6_2/phase2_split4/split4_results.csv`、`split4_predictions.npz`（L3 行/键）
- `v8_1/FROZEN_RESULTS.json` → `H4_controlled_split.seq40`
- `v8_1/CONTROLSPLIT_IDENTITY.json` → `marginal_medians.seq40`
- 主稿 Table 2 的 seq40 行（0.1961 / +0.0239 / +0.1721 / 0.0631 / 2.6911）
- Figure 1 Panel C：三卡 → **四卡**（Random / Complex / seq40 / Family-OOD 各自独立）
- Figure 2：seq40 数据点改为独立值，并**删除** "identical partition" 标注
- 正文中「complex 与 seq40 产生 identical assignments，故不作为独立比较」的表述改为
  「seq40 是独立的序列同一性划分（301 units vs 230），与 complex 仅约 20% 行一致，
  因此分别报告」

---

## 五、防回归

新增验证检查：

| 检查 | 内容 |
|---|---|
| BM | `H4_controlled_split.seq40` 不得再等于 `complex` |
| BN | 最终 PDF 不得再出现 "identical assignments" / "identical partition" / "Complex / seq40" |
| BO | `KEY_LEAKAGE_EXACT.json` 中不得有任何 level 存在跨 fold 的 obs_key |

原备份保留为 `*.prekeyfix.*` / `*.preL3fix.*` 供对照，但不参与任何计算或验证。
