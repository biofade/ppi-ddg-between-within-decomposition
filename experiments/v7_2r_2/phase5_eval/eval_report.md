# V7.2R-2 PHASE 5 — U_AB 匹配宇宙 + Track B 评估

- 时间: 2026-10-14 21:10:25
- Track B npz: **174**

## 5.1 U_AB（Track A valid ∩ Track B valid）

| universe | rows | targets | families |
|---|---:|---:|---:|
| Track A (full) | 4233 | 174 | 91 |
| Track B (imported) | 4219 | 174 | 91 |
| **U_AB (matched)** | 4219 | 174 | 91 |

- n≥5 target 覆盖：**115/116 = 99.1%**
- §40 门槛（≥80% of Track A reference）：**PASS ✅**

- 维度：base=188 F2_small=321 B_wt=1280 B_dr=1280
- **U_AB universe: rows 4219 / targets 174 / families 91**

### 5.1b NOT_RUN 模型（如实标注，不伪造）

- **B3_strong_dr_partner** → NOT_RUN：无 partner/MINT 表示（partner_repr 缺失，partner_sanity 为空）
- **B4_interaction** → NOT_RUN：无伙伴向量 p，无法构造 Δr⊙p（§8 明确禁止在维度不兼容时强行做 elementwise）

## 5.2 冻结嵌套 family-CV（B0/A2/B1/B2）

| B0_Base | s42 | Macro ρ=+0.0274 | Eff@10=+0.2893 | C_W=+0.1000 | r=+0.1567 |
| A2_smallPLM_dh | s42 | Macro ρ=+0.0162 | Eff@10=+0.4132 | C_W=+0.0857 | r=+0.1275 |
| B1_strong_repr | s42 | Macro ρ=+0.0685 | Eff@10=+0.3349 | C_W=+0.1052 | r=+0.1506 |
| B2_strong_dr | s42 | Macro ρ=+0.0336 | Eff@10=+0.3010 | C_W=+0.0901 | r=+0.1674 |
| B0_Base | s43 | Macro ρ=+0.0064 | Eff@10=+0.3000 | C_W=+0.0942 | r=+0.1905 |
| A2_smallPLM_dh | s43 | Macro ρ=-0.0248 | Eff@10=+0.4270 | C_W=+0.0835 | r=+0.1587 |
| B1_strong_repr | s43 | Macro ρ=+0.0304 | Eff@10=+0.2964 | C_W=+0.0713 | r=+0.0571 |
| B2_strong_dr | s43 | Macro ρ=-0.0086 | Eff@10=+0.3024 | C_W=+0.0701 | r=+0.0967 |
| B0_Base | s44 | Macro ρ=+0.0325 | Eff@10=+0.2028 | C_W=+0.2029 | r=+0.2361 |
| A2_smallPLM_dh | s44 | Macro ρ=+0.0364 | Eff@10=+0.1961 | C_W=+0.2028 | r=+0.2604 |
| B1_strong_repr | s44 | Macro ρ=+0.0528 | Eff@10=+0.3848 | C_W=+0.1973 | r=+0.2350 |
| B2_strong_dr | s44 | Macro ρ=+0.0666 | Eff@10=+0.3907 | C_W=+0.2043 | r=+0.2528 |

### 3-seed 均值

| Ranker | Macro ρ | Eff@10 | C_W | Kendall | Pairwise δ=.5 | Hit@10 | Overall r |
|---|---:|---:|---:|---:|---:|---:|---:|
| A2_smallPLM_dh | +0.0093 | +0.3454 | +0.1240 | +0.0152 | +0.5712 | +0.7351 | +0.1822 |
| B0_Base | +0.0221 | +0.2640 | +0.1324 | +0.0390 | +0.5752 | +0.7333 | +0.1944 |
| B1_strong_repr | +0.0506 | +0.3387 | +0.1246 | +0.0347 | +0.5697 | +0.7329 | +0.1476 |
| B2_strong_dr | +0.0305 | +0.3314 | +0.1215 | +0.0393 | +0.5693 | +0.7357 | +0.1723 |

- 耗时 8.5 分钟
