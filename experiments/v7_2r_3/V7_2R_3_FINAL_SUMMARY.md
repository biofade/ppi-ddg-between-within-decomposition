# V7.2R-3 Final Summary

**生成时间**：2026-10-14
**ZIP**：`F:\dsh\downloads\v7_2r_3\V7_2R_3_PPI_CONTEXT_FEATURES.zip`
**SHA256**：`b7e9c039e9a8fdc29d874193576d9c0bcfcdab317158945c9e86e395393087e8`

---

## 1. Why V7.2R-3 Was Required

V7.2R-2 的 Kaggle notebook 存在一个**真实的表示缺陷**（§1/§29）：

```python
build_variant()   # 只返回 sites_by_chain 中的链
```

`target_chains.csv` 也由 `sit[["target_id","chain"]].drop_duplicates()` 写出，
因此**只含携带突变位点的链**。实测后果：

| 项 | V7.2R-2 实测 | V7.2R-3 修复后 |
|---|---:|---:|
| (target,chain) 记录 | **233** | **654** |
| 每 target 链数中位 | **1** | **3**（min 2, max 24） |
| 至少缺 1 条伙伴链的 target | **150 / 174 = 86.2%** | **0** |
| 丢失伙伴链总数 | **421** | **0** |

⇒ V7.2R-2 测的是 **`STRONG_GENERIC_PLM_PERTURBATION`（链隔离）**，
**不是** `PPI_AWARE_PERTURBATION`。

**修复**：`build_complex()` 使用该 PDB 的**全部蛋白链**，
突变只施加到受影响链，伙伴链保持 WT（§10/§11）。

---

## 2. PPI Representation Semantics

| 项 | 值 |
|---|---|
| Model | **ESM2-650M over FULL complex**（全部链按链 ID 排序拼接，跨链自注意力） |
| level | **`MUTATION_SITE_PPI`** ✅（§6 要求分类） |
| dimensions | **1280**（`r_wt` / `r_mut` / `delta_r` 各 1280） |
| WT | `r_wt`：全复合物，伙伴链 WT |
| Mut | `r_mut`：受影响链突变，伙伴链 WT |
| partner-dependent | **YES**（由 §3 的 sanity 实测确认） |
| 嵌入链数中位 | **3**（V7.2R-2 为 1） |
| labels_used | **False** |
| supervised_heads_used | **False** |

> 依 §6「不从维度推断语义」：本判定来自 **provenance.json**
> （`representation_level: MUTATION_SITE_PPI`、
> `complex_construction: ALL_PROTEIN_CHAINS (partners retained)`）
> **以及** §3 的伙伴扰动实测。

---

## 3. Partner Sanity Test

10 个 label-free 分层 target（按突变数小/中/大）：

| 量 | 值 |
|---|---:|
| `D_true = ‖r(A alone) − r(A + 真伙伴 B)‖` 中位 | **4.0289** |
| `D_wrong = ‖r(A + 真伙伴) − r(A + 错伙伴 B′)‖` 中位 | **2.6557** |

**verdict：`PARTNER_CONDITIONING_CONFIRMED` ✅**

- 伙伴的加入带来**显著**表示变化（4.03 ≫ 数值噪声）
- 且**错误伙伴**带来的变化较小（2.66 < 4.03）→ 表示**对伙伴身份有选择性**
- 依 §7 门槛（`D_true > 1e-3`）→ **允许继续全量提取**

---

## 4. Coverage

| 项 | 值 |
|---|---:|
| rows | **3452** |
| targets | **132** |
| families | **70** |
| n≥5 | **87** |
| 嵌入链数中位 | **3** |

> ⚠️ Kaggle 提取只完成 **132 / 174** target（V7.2R-2 为 174）。
> **未续跑**（不作假设）；依 §16 明示为覆盖受限。

---

## 5. Matched U_ABC

`U_ABC = Track A ∩ Track B(ESM2-650M 链隔离) ∩ PPI context`

| universe | rows | targets | families |
|---|---:|---:|---:|
| Track A | 4233 | 174 | 91 |
| Track B（链隔离） | 4219 | 174 | 91 |
| PPI context | 3452 | 132 | 70 |
| **U_ABC** | **3452** | **132** | **70** |

- n≥5：U_ABC **87** vs V7.2R-2 匹配宇宙 **115** → 相对 **75.7%**
- **§16 门槛（≥80%）：`COVERAGE_LIMITED` ⚠️**

> **所有下列比较都在同一 U_ABC 上重算**（§16 强制）。
> **结论必须带此覆盖限定**；续跑至 ≥80% 后需重算。

---

## 6. Primary Results

**3-seed 均值（敏感性口径）** — U_ABC 3452 rows / 132 targets / 70 families

| Model | Macro ρ | Eff@10 | C_W | Kendall | Pairwise δ=.5 | Overall r |
|---|---:|---:|---:|---:|---:|---:|
| **B2** chain-isolated | +0.0329 | +0.1750 | +0.2087 | +0.0216 | +0.6426 | +0.3414 |
| **B3a** + static PPI | **+0.0553** | +0.1690 | +0.2160 | +0.0760 | +0.6512 | +0.2434 |
| **B3b** + Δr×PPI interaction | **+0.0560** | **+0.1828** | +0.2143 | **+0.0887** | +0.6500 | +0.2481 |
| **B4** + PPI perturbation | **+0.0651** | +0.1827 | +0.2133 | +0.0298 | **+0.6580** | +0.3040 |
| **B_shuf** shuffled PPI | +0.0349 | +0.1521 | +0.1821 | +0.0525 | +0.6403 | +0.2302 |

**ensemble prediction 口径（主口径，与 bootstrap 一致）**

| Model | Macro ρ | Eff@10 | C_W |
|---|---:|---:|---:|
| B2 chain-isolated | +0.0179 | +0.1311 | +0.2187 |
| **B3a** static PPI | **+0.0889** | +0.1661 | **+0.2318** |
| **B3b** interaction | **+0.0748** | +0.1669 | **+0.2325** |
| **B4** PPI perturb | +0.0481 | **+0.2056** | +0.2245 |
| B_shuf shuffled | +0.0412 | +0.1847 | +0.2058 |

> `mean(ρ_seed) ≠ ρ(mean(pred_seed))` —— 两口径**分别报告，未混用**（§14）。
> 两个口径一致显示：**B3a/B3b > B_shuf ≥ B2**。

---

## 7. B3/B4 − B2 Bootstrap

**1500/1500 paired target-cluster bootstrap**（cluster = canonical target）
> §22 要求 5000；本轮因单机 CPU 成本（单次指标评估 0.40 s × 2 × 1500 × 5）降为 1500，**如实标注**。

### `B3b − B2`（Δr×PPI 交互 vs 链隔离）— **核心对比**

| 指标 | point | 95% CI | 排除 0? |
|---|---:|---|---|
| **ΔMacroSpearman** | **+0.0569** | **[+0.0035, +0.1203]** | **✅ 排除 0** |
| ΔEff@10 | +0.0358 | [−0.1228, +0.2077] | ❌ 含 0 |
| ΔC_W | +0.0138 | [−0.0203, +0.0319] | ❌ 含 0 |
| ΔOverall_r | −0.0879 | [−0.1867, +0.0221] | ❌ 含 0 |
| ΔRMSE | +0.0965 | [−0.0756, +0.2727] | ❌ 含 0 |

⇒ **§23 成功判据：`SUPPORTED` ✅**
（ΔMacroSpearman > 0 且 CI 排除 0；ΔEff@10 与 ΔC_W 点估计均为正；无灾难性 RMSE 退化）

### `B3a − B2`（静态 PPI vs 链隔离）

| 指标 | point | 95% CI | 排除 0? |
|---|---:|---|---|
| **ΔMacroSpearman** | **+0.0710** | **[+0.0147, +0.1345]** | **✅ 排除 0** |
| ΔEff@10 | +0.0349 | [−0.0682, +0.2056] | ❌ 含 0 |
| ΔC_W | +0.0131 | [−0.0202, +0.0298] | ❌ 含 0 |

⇒ **§23：`SUPPORTED` ✅**

### `B4 − B2`（PPI 扰动 vs 链隔离）

| 指标 | point | 95% CI | 排除 0? |
|---|---:|---|---|
| ΔMacroSpearman | +0.0302 | [−0.0184, +0.0815] | ❌ 含 0 |
| **ΔEff@10** | **+0.0744** | **[+0.0092, +0.2066]** | **✅ 排除 0** |
| ΔC_W | +0.0057 | [−0.0085, +0.0162] | ❌ 含 0 |

⇒ **§23：`NOT SUPPORTED` ❌**（Macro ρ CI 含 0）。
依 §18/§47：**Top-K 显著改善但整体排序未显著** → 形式上是 **`TOPK_ONLY_SIGNAL`**
（但 §23 冻结判据要求 Macro ρ 先排除 0，故**不能**记为一般排序增益）。

### `B4 − B3a`（§18：扰动 vs 静态）

| 指标 | point | 95% CI | 排除 0? |
|---|---:|---|---|
| ΔMacroSpearman | **−0.0408** | [−0.0986, +0.0143] | ❌ 含 0 |
| ΔEff@10 | +0.0395 | [−0.0596, +0.1334] | ❌ 含 0 |
| ΔRMSE | −0.1158 | [−0.2961, +0.0212] | ❌ 含 0 |

⇒ **不显著** → **突变诱导的 PPI 变化并未超越静态复合物身份**（§24 的 P2 不成立）。

### `Shuf − 真实`（§19 强制对照）

| 指标 | point | 95% CI | 排除 0? |
|---|---:|---|---|
| ΔMacroSpearman | **−0.0477** | [−0.1033, +0.0055] | ❌ **含 0（边界）** |
| ΔEff@10 | +0.0186 | [−0.0766, +0.1066] | ❌ 含 0 |
| ΔC_W | −0.0260 | [−0.0446, +0.0095] | ❌ 含 0 |

> 真实 PPI **优于**打乱 PPI（Ensemble Macro ρ +0.0889 vs +0.0412；
> seed-mean +0.0553 vs +0.0349），方向一致，但 **CI 含 0（上界 +0.0055，几乎触零）**。
> ⇒ 依 §19：**不足以声称生物伙伴机制的特异性**。

---

## 8. True vs Shuffled Partner

| 口径 | B3a 真实 PPI | B_shuf 打乱 PPI | 差 |
|---|---:|---:|---:|
| Ensemble Macro ρ | **+0.0889** | +0.0412 | +0.0477 |
| seed-mean Macro ρ | **+0.0553** | +0.0349 | +0.0204 |
| C_W（ensemble） | **+0.2318** | +0.2058 | +0.0260 |

**方向支持特异性，但统计上未达显著**（ΔMacro ρ CI [−0.1033, +0.0055]）。

---

## 9. Wrong Partner Control

- **表示层**：已做（notebook cell 08）→ `D_true` 4.03 > `D_wrong` 2.66 ✅
  ⇒ 表示**对伙伴身份有选择性**（§20）
- **排名层**：错伙伴表示**未导出**，故排名层的错伙伴对照 **NOT_RUN**。

---

## 10. Scientific Verdict

```
PPI_RESULT = PPI_CONTEXT_NULL          （最终判定，见 §13 的机制判别）
```

> **判定过程如实记录（含一次自我修正）**：
>
> **第一步（PHASE 5/6/7 后）** 曾判 `STATIC_PPI_CONTEXT_GAIN`，依据是：
> - `B3a − B2` ΔMacro ρ = **+0.0710 [+0.0147, +0.1345]** → 显著 ✅
> - `B3b − B2` ΔMacro ρ = **+0.0569 [+0.0035, +0.1203]** → 显著 ✅
> - `B4 − B3a` ΔMacro ρ = −0.0408（含 0）→ 突变诱导的 PPI 变化无增量
>
> **第二步（PHASE 8 机制判别后）修正为 `PPI_CONTEXT_NULL`**：
> 同方差零信息噪声对照显示 `B3a_real − NC_noise` = **+0.0474 [−0.0013, +0.0981]**，
> **CI 含 0** ⇒ 上述「显著增益」与**噪声所产生的增益统计上不可区分**。
> 详见 §13。

### 必须同时声明的三个限定

**(1) 覆盖受限**：`COVERAGE_LIMITED` —— U_ABC 仅 132/174 target，
n≥5 相对覆盖 **75.7% < 80%**（§16 门槛）。

**(2) §14 的目标恒定退化（PHASE 8 已实测判定）**：`G_tgt` 对同一 target 内所有突变
**数值上恒定** —— 实测 target 内逐维 SD **中位 = 0.000e+00**（半数 target 完全无变化），
全局 max 仅 **2.337e-05**（float32 精度量级）。
因此 `B3a` 与 `B2` 的靶点内预测 Spearman 中位仅 **0.774**（而非理论上的 1.0）
—— LightGBM 在**噪声级差异**上找到了分裂点。
**PHASE 8 的噪声对照进一步证明**：把该噪声换成**同方差纯噪声**，
性能相当（+0.0439 vs +0.0553，差异 CI 含 0）。
⇒ **`B3a` 的增益机制已判定为「噪声级伪随机排序」，非 PPI 信息。**

**(3) 负对照显示信号主要来自身份而非扰动**：

| 对照 | 构成 | Macro ρ |
|---|---|---:|
| `NC0_base_only` | base | +0.0092 |
| **`NC1_wt_only`** | base + **`r_wt`（无差分）** | **+0.0593** |
| `NC2_identity` | base + mutation identity | +0.0113 |
| `NC3_dr_norm` | base + `‖Δr‖` | +0.0300 |
| `NC4_dr_shuffled` | base + **同 target 内打乱的 Δr** | +0.0225 |
| `ref_B2_dr` | base + Δr | +0.0329 |

- **NC1（仅 `r_wt`）高于 B2（Δr）** → 与 V7.2R-2 一致：**微弱信号来自结构身份**
- NC2 表明**不是**氨基酸替换身份编码
- NC3（仅幅度）低于完整 Δr → 幅度不是全部
- NC4（同 target 内打乱 Δr）低于真 Δr → Δr 的**逐突变特异性**有少量贡献，但很小


---

## 11. Allowed Claims

- 已定位并修复 V7.2R-2 的链丢弃缺陷（**233 → 654** 条 (target,chain)；
  86.2% 的 target 曾丢伙伴链）
- **首次产出真正的 PPI-conditioned site-level 表示**（`MUTATION_SITE_PPI`），
  并通过 §7 的伙伴 sanity 检验：`D_true` **4.03** > `D_wrong` **2.66**
- 在**同一 U_ABC**、**冻结嵌套 family-OOD** 协议下，
  加入 PPI 上下文使 **Macro ρ 显著提升**
  （`B3a−B2` **+0.0710** [+0.0147,+0.1345]；`B3b−B2` **+0.0569** [+0.0035,+0.1203]）
- **突变诱导的 PPI 变化（`B4−B3a`）未显示增量**（−0.0408，含 0）
- 真实 PPI 在点估计上**优于**打乱 PPI（+0.0889 vs +0.0412），但 **CI 含 0**
- 负对照显示**结构身份（`r_wt`）比突变扰动承担了更多信号**

---

## 12. Forbidden Claims

- ❌ 「PPI conditioning 解决了 family-OOD 排序」—— Macro ρ 绝对值仍 **≤ 0.09**
- ❌ 「partner-specific 生物学机制已被证明」——
  `Shuf−真实` 的 CI **含 0**（§19 明确禁止）
- ❌ 「`B3a` 的增益是真实 PPI 信息」—— §14 的目标恒定退化 + 打乱对照未分离
- ❌ 「protein foundation models do not work」（§51）
- ❌ 「we invented the underlying PPI representation」（§51）
- ❌ 任何不含 `COVERAGE_LIMITED` 限定的推广

---

## 13. Final Representation-Branch Decision

### 13.1 §14 机制判别（决定性收尾实验，`phase8_noise_control/`）

**设计**：构造**零信息但同方差**的对照 `G_noise`——逐 target、逐维以
`G_real` 的 target 内 SD 为尺度生成独立高斯噪声，并施加 **1e-4 噪声地板**
（使其**对 MINT 更有利** → 检验更严格）。

**校准实测（关键）**：

| 量 | 值 |
|---|---:|
| `G_real` target 内逐维 SD **中位** | **0.000e+00**（即半数 target 内**完全无变化**） |
| `G_real` target 内逐维 SD 全局 max | 2.337e-05 |
| `G_noise` target 内逐维 SD 中位 | 9.007e-05 |

| Model | seed-mean Macro ρ |
|---|---:|
| `B2_dr` | +0.0329 |
| **`B3a_real_PPI`** | **+0.0553** |
| `B_shuf_target` | +0.0349 |
| **`NC_noise_PPI`** | **+0.0439** |

**关键对比（paired cluster bootstrap，1500）**：

| contrast | ΔMacro ρ | 95% CI | 排除 0? |
|---|---:|---|---|
| **`B3a_real − NC_noise`** | **+0.0474** | **[−0.0013, +0.0981]** | **❌ 含 0（边界）** |
| `B3a_real − B_shuf_target` | +0.0477 | [−0.0055, +0.1033] | ❌ 含 0 |
| `NC_noise − B2_dr` | +0.0235 | [−0.0027, +0.0502] | ❌ 含 0 |

### 13.2 判定

```
NOISE_CONTROL_SEPARATED     = NO
CONTINUE_REPRESENTATION_ROUTE = NO
```

**理由（三条独立证据，全部指向同一结论）**：

1. **真实 PPI 未能显著优于同方差噪声**（+0.0474，CI 含 0，上界仅 +0.0981）
2. **真实 PPI 也未能显著优于打乱 PPI**（+0.0477，CI 含 0）
3. **纯噪声特征本身就能把 Macro ρ 从 +0.0329 推到 +0.0439**（+0.0235，
   与 B3a 相对 B2 的 +0.0710 同量级的 **1/3**，且方向一致）

⇒ **`B3a` 相对 `B2` 的显著增益（+0.0710）与「同方差零信息噪声」所产生的
增益在统计上不可区分**（§14 的机制问题**未能解决**）。

⇒ 依 §25：**关闭表示路线**，进入 **V7.3 physics / structure-informed ranking**。

> **注**：`B2_dr` 本身相对 `B0_base`（+0.0092）也几乎无增益，
> 且负对照显示信号主要来自 `r_wt`（身份）而非 Δr（扰动）。
> 与 V7.2R-2 的 `STRONG_BACKBONE_NULL` 一致 —— **表示规模与 PPI 条件化都未解决排序瓶颈**。

---

## 14. Next Step

```
STOP REPRESENTATION ROUTE  →  V7.3 physics / structure-informed ranking
```

依 §25，**不要**：试 ESM-3B / ESM-15B / 另一个通用 PLM / 更深 MLP / transformer。
**立即进入 V7.3**，焦点：

```text
physics
structure
local interaction geometry
energy changes
residue contacts
interface disruption
```

**V7.3 应立即利用的既有资产**：

| 资产 | 用途 |
|---|---|
| `canonical_target_labels.csv`（214 target / 108 family） | 权威 target 命名（§7） |
| `OUTER_FOLDS_V7_2` + V7.1.1 嵌套协议 | 冻结的 family-OOD 评估 |
| `data/v5_0/mutation_sites.csv`（93.1% WT 校验通过） | 位点结构映射 |
| `data/pdb/`（230/230 target 覆盖） | 真实坐标 |
| V7.0 BLOCK L（4473 行，118 个静态几何特征） | **已算好的 structure/geometry 特征** |
| MU3 + z 标定外壳 | 标定/排序解耦 |
| `target_mechanism_table.csv`（111 列） | target 级物理/机制描述符 |

**唯一尚未测试的方向**：**能量/物理层面的局部相互作用**
（不是表示规模，不是 PPI 条件化）。

---

## 附：文件清单

### 脚本（`scripts/v7_2r_3/`）
`00_build_notebook.py`、`01_prepare_ppi_input.py`、`02_import_ppi_features.py`、
`03_matched_eval.py`（**已修复维度 bug + 折内 PCA**）、`04_bootstrap.py`（已向量化）、
`05_controls_and_report.py`、**`06_noise_control.py`（决定性判别）**

### Kaggle（`kaggle/v7_2r_3/`）
`V7_2R_3_PPI_Context_Kaggle.ipynb`、`kaggle_input_package.zip`（**654** 条链记录）、
`README_KAGGLE.md`、`kaggle_input_manifest.csv`

### 结果（`experiments/v7_2r_3/`）
`phase1_prepare/prepare_report.md`、`phase4_import/`（`import_report.md`、
`partner_sanity.csv`、`import_meta.json`）、`phase5_eval/`（`ppi_results.csv`、
`ppi_seedmean.csv`、`ppi_ensemble.csv`、`ppi_predictions.npz`、`u_abc_manifest.csv`、
`eval_report.md`）、`phase6_bootstrap/`（`contrasts.csv`、`bootstrap_report.md`）、
`phase7_controls/`（`negative_controls.csv`、`negative_controls.md`）、
**`phase8_noise_control/`（`noise_control_results.csv`、
`noise_control_contrasts.csv`、`noise_control_predictions.npz`、`noise_control.md`）**

**无 `failures.log`。**
