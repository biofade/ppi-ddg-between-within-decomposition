# V7.2R-3 PHASE 5（修复版）— U_ABC + 冻结嵌套 family-CV

- 时间: 2026-10-14 22:25:10
- PPI context npz: **132**

## 5.1 U_ABC = Track A ∩ Track B ∩ PPI context

| universe | rows | targets | families |
|---|---:|---:|---:|
| Track A | 4233 | 174 | 91 |
| Track B (ESM2-650M chain-isolated) | 4219 | 174 | 91 |
| PPI context | 3452 | 132 | 70 |
| **U_ABC** | 3452 | 132 | 70 |

- n>=5：U_ABC **87** vs V7.2R-2 匹配宇宙 **115** → 相对 **75.7%**
- §16 门槛（>=80%）：**COVERAGE_LIMITED ⚠️**

> ⚠️ **`COVERAGE_LIMITED`**：PPI-context 提取只完成 **132/174** target。
> 原因：Kaggle 会话在 132 个 target 后结束（或主动停止）；
> 依 §16 必须明确标注，且结论必须带此限定。
> **可续跑**：重跑 notebook（`RESUME=True`）会跳过已完成 target。

- 维度：base=188 B_dr=1280 Gt_raw=1280 DG_raw=1280
- **U_ABC: rows 3452 / targets 132 / families 70**

> ⚠️ **§14**：`G_tgt` 对同一 target 内所有突变**完全相同**，在加性模型中**不可能**改变靶点内排序；只有 `DRxG`（交互）可能。

## 5.2 冻结嵌套 family-CV

| B2_chain_isolated | s42 | Macro ρ=+0.0357 | Eff@10=+0.1504 | C_W=+0.2186 | Kendall=+0.0150 | r=+0.3364 |
| B3a_static_PPI | s42 | Macro ρ=+0.0413 | Eff@10=+0.1331 | C_W=+0.2074 | Kendall=+0.0529 | r=+0.2287 |
| B3b_dr_x_PPI | s42 | Macro ρ=+0.0385 | Eff@10=+0.1663 | C_W=+0.2045 | Kendall=+0.0684 | r=+0.1964 |
| B4_PPI_perturb | s42 | Macro ρ=+0.0785 | Eff@10=+0.1742 | C_W=+0.2186 | Kendall=+0.0178 | r=+0.2943 |
| B_shuf_PPI | s42 | Macro ρ=+0.0318 | Eff@10=+0.1146 | C_W=+0.1943 | Kendall=+0.0242 | r=+0.1347 |
| B2_chain_isolated | s43 | Macro ρ=+0.0451 | Eff@10=+0.1928 | C_W=+0.2090 | Kendall=+0.0482 | r=+0.3319 |
| B3a_static_PPI | s43 | Macro ρ=+0.0407 | Eff@10=+0.1716 | C_W=+0.2213 | Kendall=+0.1089 | r=+0.2132 |
| B3b_dr_x_PPI | s43 | Macro ρ=+0.0535 | Eff@10=+0.1585 | C_W=+0.2169 | Kendall=+0.1109 | r=+0.2512 |
| B4_PPI_perturb | s43 | Macro ρ=+0.0423 | Eff@10=+0.1638 | C_W=+0.2133 | Kendall=+0.0081 | r=+0.2928 |
| B_shuf_PPI | s43 | Macro ρ=+0.0577 | Eff@10=+0.1851 | C_W=+0.1909 | Kendall=+0.0673 | r=+0.2911 |
| B2_chain_isolated | s44 | Macro ρ=+0.0178 | Eff@10=+0.1818 | C_W=+0.1985 | Kendall=+0.0015 | r=+0.3559 |
| B3a_static_PPI | s44 | Macro ρ=+0.0839 | Eff@10=+0.2022 | C_W=+0.2195 | Kendall=+0.0660 | r=+0.2884 |
| B3b_dr_x_PPI | s44 | Macro ρ=+0.0760 | Eff@10=+0.2236 | C_W=+0.2214 | Kendall=+0.0869 | r=+0.2968 |
| B4_PPI_perturb | s44 | Macro ρ=+0.0746 | Eff@10=+0.2099 | C_W=+0.2080 | Kendall=+0.0635 | r=+0.3249 |
| B_shuf_PPI | s44 | Macro ρ=+0.0151 | Eff@10=+0.1566 | C_W=+0.1611 | Kendall=+0.0660 | r=+0.2649 |

- **泄漏断言**：outer violation = **0**；inner violation = **0**
- **PCA train-only：PASS**（所有块在折内 fit）

### 3-seed 均值（敏感性口径）

| Model | Macro ρ | Eff@10 | C_W | Kendall | Pairwise δ=.5 | Overall r |
|---|---:|---:|---:|---:|---:|---:|
| B2_chain_isolated | +0.0329 | +0.1750 | +0.2087 | +0.0216 | +0.6426 | +0.3414 |
| B3a_static_PPI | +0.0553 | +0.1690 | +0.2160 | +0.0760 | +0.6512 | +0.2434 |
| B3b_dr_x_PPI | +0.0560 | +0.1828 | +0.2143 | +0.0887 | +0.6500 | +0.2481 |
| B4_PPI_perturb | +0.0651 | +0.1827 | +0.2133 | +0.0298 | +0.6580 | +0.3040 |
| B_shuf_PPI | +0.0349 | +0.1521 | +0.1821 | +0.0525 | +0.6403 | +0.2302 |

### ensemble prediction 口径（主口径）

| Model | Macro ρ | Eff@10 | C_W |
|---|---:|---:|---:|
| B2_chain_isolated | +0.0179 | +0.1311 | +0.2187 |
| B3a_static_PPI | +0.0889 | +0.1661 | +0.2318 |
| B3b_dr_x_PPI | +0.0748 | +0.1669 | +0.2325 |
| B4_PPI_perturb | +0.0481 | +0.2056 | +0.2245 |
| B_shuf_PPI | +0.0412 | +0.1847 | +0.2058 |

- 耗时 13.2 分钟
