# V7.0 Final Summary

**时间**：2026-10-14（约 13:55 → 16:45，约 2.8 小时）
**约束**：未修改 V1–V6.3 任何文件；新命名空间 `scripts/v7_0`/`experiments/v7_0`/`logs/v7_0`；
CPU only；**未下载任何文件**；未伪造数字；未为得到明确路线而强行解释

**状态**：PHASE 0、1、1b、2、3、5 **COMPLETED**；PHASE 4 **NOT_RUN**（见 §14 阻塞）；PHASE 6 **PARTIAL**

---

## 1. Coverage

| 项 | 值 |
|---|---|
| canonical rows | **7523** |
| targets（canonical） | **230** |
| **统一 row universe（§5 matched coverage）** | **4473 行 / 214 targets / 108 L4-units** |
| 占 canonical 比例 | **59.5%** |
| **Set feature coverage** | 4450 / 4473 行（23 行缺 V6.3 target 表） |
| **Local feature coverage** | **4473 行（100% of universe）** |
| Global coverage | 4450 / 4473 |

**excluded rows（无 BLOCK L）= 3050（40.5%）**，涉及 16 targets。
按 §5，**所有模型在同一 matched-coverage 子集上训练与评估**，未跨 row universe 比较。

### BLOCK L 的构造与一次重要修正（如实记录）

| 步骤 | 结果 |
|---|---|
| 初版用 `(ann_chain, pdb_resseq)` 定位残基 | WT 一致率 **46.8%** → 一度误判「上游映射不可靠」 |
| 系统排查 4 种约定 | `pdb_resseq` / `ann_pos1` / `numbering_offset` 均 **46.9%** |
| **正确约定**：`链 = ann_chain`，**链内第 `seq_index0` 个残基（0-based）** | **WT 一致率 93.1%**（5123/5503） |
| 最终 BLOCK L | **4473 行 / 208 targets，118 特征**（6/8/10 Å 三半径的接触/partner 理化/界面距离/突变身份） |

⇒ **修正的是本项目最初的调用约定，不是上游数据。** 只有 WT 校验通过的位点进入特征。

---

## 2. Architecture Ablation（L4_family，seed ensemble，4473 行 / 214 targets）

| Model | Overall r | C_W | Macro ρ | Kendall | Eff@10 | **Row C_B** | **Equal C_B** | slope |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **A0 Base** | +0.2205 | +0.1734 | +0.0333 | +0.0159 | +0.1420 | +0.0472 | +0.5052 | +0.2148 |
| **A1 +Set** | **+0.2408** | **+0.1795** | **+0.0507** | −0.0024 | **+0.2080** | +0.0613 | +0.5823 | +0.2029 |
| A2 +Global | +0.1416 | +0.1599 | +0.0827 | +0.0469 | +0.2186 | −0.0183 | +0.4320 | +0.1562 |
| A3 +Local | +0.2487 | +0.1720 | +0.0057 | +0.0114 | +0.1484 | +0.0767 | +0.4810 | +0.2012 |
| A4 +Set+Local | +0.2487 | +0.1780 | +0.0508 | +0.0347 | +0.1101 | +0.0707 | +0.5713 | +0.1933 |
| A5 Set+Local | +0.0514 | +0.0603 | +0.0698 | +0.0630 | +0.2143 | −0.0089 | +0.5489 | +0.1667 |
| A6 Global+Set+Local | +0.0348 | +0.0548 | +0.0672 | +0.0548 | +0.1274 | −0.0199 | +0.5361 | +0.1578 |

**L1_random（diagnostic only）**：A0 r=+0.6893 / C_W=+0.2873 / macro ρ=+0.1092；
A1 r=**+0.7729** / C_W=+0.2563 / macro ρ=+0.0807。

> **注意**：A2 与 A5/A6 的 `row C_B` 为**负**而 `equal C_B` 为**正** → §E weighting artifact 再次出现
> （见 §7）。**禁止**只看单一口径下结论。

---

## 3. Architecture Gains（paired cluster bootstrap，**2000/2000**，单位 = complex）

### Set gain（A1 − A0）
- ΔC_W: **+0.0061**, 95% CI [−0.0102, +0.0229] → **含 0 ❌**
- ΔMacro ρ: **+0.0174**, CI [−0.0334, +0.0697] → **含 0 ❌**
- ΔEff@10: +0.0622, CI [−0.6524, +0.4970] → 含 0 ❌
- ΔOverall r: +0.0203, CI [−0.0723, +0.1041] → 含 0 ❌

### Local gain（A3 − A0）
- ΔC_W: −0.0014, CI [−0.0222, +0.0182] → 含 0 ❌
- ΔMacro ρ: **−0.0275**, CI [−0.0668, +0.0087] → **含 0，且点估计为负** ❌
- ΔEff@10: +0.0094, CI [−0.4443, +0.3529] → 含 0 ❌
- ΔOverall r: +0.0281, CI [−0.0095, +0.0638] → 含 0 ❌

### Set+Local（A4 − A0）
- ΔC_W: +0.0046, CI [−0.0170, +0.0285] → 含 0 ❌
- ΔMacro ρ: +0.0175, CI [−0.0331, +0.0727] → 含 0 ❌
- ΔEff@10: +0.0769, CI [−1.1511, +0.7150] → 含 0 ❌
- ΔOverall r: +0.0281, CI [−0.0738, +0.1291] → 含 0 ❌

### Local incremental（A4 − A1）
- ΔC_W: −0.0015, CI [−0.0145, +0.0176] → 含 0 ❌
- ΔMacro ρ: **+0.0001**, CI [−0.0322, +0.0325] → 含 0 ❌
- ΔOverall r: +0.0079, CI [−0.0204, +0.0441] → 含 0 ❌

### §13 Synergy = M(A4) − M(A1) − M(A3) + M(A0)

| 指标 | synergy | 95% CI | CI>0? |
|---|---:|---|---|
| C_W | −0.0001 | [−0.0159, +0.0156] | ❌ |
| Macro ρ | +0.0276 | [−0.0092, +0.0662] | ❌ |
| Eff@10 | +0.0053 | [−0.6191, +0.4351] | ❌ |

**SET_LOCAL_SYNERGY_SUPPORTED: NO** → **`NO_SYNERGY_EVIDENCE`**（0/3 指标 CI>0）

> ### ★ 结论：**没有任何架构在 strict family-OOD 下取得统计显著的提升。**
> 全部 16 个 contrast×指标组合的 95% CI **全部含 0**。

---

## 4. Loss Ablation（A1 特征集，MLP(128) CPU torch，4473 行 / 214 targets）

| Loss | Overall r | C_W | Macro ρ | Kendall | Eff@10 | Row C_B | Equal C_B |
|---|---:|---:|---:|---:|---:|---:|---:|
| **L0 Absolute** | −0.0942 | +0.0441 | **+0.1118** | +0.1079 | **+0.2522** | −0.1383 | +0.1236 |
| L1 +Centered | −0.1004 | +0.0428 | +0.1095 | +0.1059 | +0.2404 | −0.1432 | +0.1216 |
| L2 +Pairwise | −0.0984 | +0.0425 | +0.1112 | +0.1031 | +0.2323 | −0.1409 | +0.1266 |
| **L3 Target-balanced** | **+0.0065** | **+0.0491** | **+0.0674** | +0.0318 | +0.1124 | −0.0426 | +0.0477 |

- **L1（centered）与 L2（pairwise）相对 L0 均无改善**（Macro ρ −0.002 / −0.001）
- **L3（target-balanced）各项更差**（Macro ρ 0.1118 → **0.0674**，Eff@10 0.2522 → 0.1124）
  → **target-balanced 训练在 strict family-OOD 下不成立**（§14 要求的关键 ablation）
- ⚠️ MLP 的 `overall_r` 为**负**（L0 −0.094），与 LightGBM A1（+0.241）不可直接比
  —— 这是 **learner 差异**，不是 loss 差异；loss 消融只在 MLP 内部相互可比

---

## 5. Physics Rank Distillation — **NOT_RUN**

| 项 | 值 |
|---|---|
| 可用 physics teacher | `standardized_preds.csv` 的 FoldX / Flex-ddG |
| **coverage** | **仅 21 / 230 targets（9.1%）**；共 3761 行（CATH） |
| 与本研究 4473 行 universe 的交集 | 远不足 |
| 结论 | **`UNAVAILABLE`** —— 未执行（§15 要求先报 coverage；不重新下载大型数据） |

**PHYSICS_RANK_SUPPORTED: NOT_RUN**（无证据，也无反证）

---

## 6. `>=30` Mutation Targets（§16）

| 组 | 模型 | n_target | n_rows | target-mean Pearson | slope | C_W | Macro ρ | Eff@10 |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 2–9 | A0 | 95 | 446 | +0.4734 | 0.2536 | +0.0046 | −0.0730 | +0.2533 |
| 2–9 | **A1** | 95 | 446 | **+0.5055** | 0.2272 | −0.0060 | −0.0057 | +0.3553 |
| **≥30** | **A0** | **34** | **3134** | +0.0538 | **0.0197** | +0.2151 | +0.1696 | +0.1214 |
| **≥30** | **A1** | 34 | 3134 | +0.0323 | 0.0091 | **+0.2312** | +0.1633 | +0.1322 |
| **≥30** | A3 | 34 | 3134 | **+0.1541** | **0.0608** | +0.2144 | **+0.1869** | **+0.1774** |
| **≥30** | A4 | 34 | 3134 | +0.0692 | 0.0197 | +0.2276 | +0.1827 | +0.1181 |

**观察**：
- `≥30` 组的 target-mean Pearson 仍**很低**（A0 **+0.0538**），slope 仍接近 0（**0.0197**）
  → **attenuation 问题未被任何架构解决**
- 但 `≥30` 组的 **C_W 与 Macro ρ 反而是三组中最高的**（C_W +0.215, Macro ρ +0.170）
  → within 排序在 mutation-rich target 上**并不差**；失效的**只是 target 均值**

---

## 7. Weighting（§E 强制双报）

| 模型 | C_B(row) | C_B(equal) | 差值 | 符号翻转 |
|---|---:|---:|---:|---|
| A0_base | +0.0472 | +0.5052 | 0.4580 | 否 |
| A1_base_set | +0.0613 | +0.5823 | 0.5210 | 否 |
| A3_base_local | +0.0767 | +0.4810 | 0.4043 | 否 |
| A4_base_set_local | +0.0707 | +0.5713 | 0.5005 | 否 |
| **A5_set_local** | **−0.0089** | +0.5489 | 0.5578 | **是** |
| **A6_global_set_local** | **−0.0199** | +0.5361 | 0.5560 | **是** |

**interpretation**：V6.3 的 weighting artifact 在 V7.0 中**再次复现** ——
两个模型的 `C_B` **符号随加权口径翻转**（row 为负、equal 为正），差值达 **0.556–0.558**。
⇒ **汇报 between-target 能力时必须固定并声明加权口径**，否则结论可被加权方式任意翻转。

---

## 8. OOD Retention

| 模型 | C_W(L1) | C_W(L4) | **Retention_W** |
|---|---:|---:|---:|
| A0_base | +0.2873 | +0.1734 | **0.604** |
| A1_base_set | +0.2563 | +0.1795 | **0.700** |
| A3_base_local | +0.2804 | +0.1720 | 0.613 |
| **A4_base_set_local** | +0.2506 | +0.1780 | **0.710** |

| 模型 | Macro ρ(L1) | Macro ρ(L4) | **Retention_ρ** |
|---|---:|---:|---:|
| A0_base | +0.1092 | +0.0333 | 0.305 |
| **A1_base_set** | +0.0807 | +0.0507 | **0.628** |
| A3_base_local | +0.0931 | +0.0057 | 0.061 |
| A4_base_set_local | +0.0646 | +0.0508 | **0.786** |

> denominator 均为正且稳定，故 ratio 可报。
> **A1/A4 的 within 信号 OOD 保留率高于 A0**（Retention_W 0.60→0.70；
> Retention_ρ 0.31→0.63/0.79）—— 这是**方向性**改善，但见 §3：差值 CI 含 0。

---

## 9. Winner

| 项 | 结果 |
|---|---|
| **WINNER_ARCHITECTURE** | **A1（Base + Set）** —— **方向性最优，但未达统计显著** |
| **WINNER_LOSS** | **L0（absolute only）** —— centered / pairwise 无改善；target-balanced 更差 |
| **WINNER_TEACHER** | **NOT_RUN**（teacher coverage 仅 9.1%） |

**按 §20 判定规则（须满足「MacroSpearman 统计显著改善 + 至少一项 C_W/Eff@10 改善」）**：

- A1：Macro ρ **+0.0174**、C_W **+0.0061**、Eff@10 **+0.0622** —— **三项方向全对**
  **但 Macro ρ 的 95% CI 含 0**（[−0.0334, +0.0697]）
- ⇒ **未满足「statistically meaningful improvement」** → **不授予 STRONG METHOD SIGNAL**
- 更不满足 VERY_STRONG（需 Macro ρ ≥ baseline+0.10）或 BREAKTHROUGH（需 ≥0.20）

> **最终标记：`WEAK_DIRECTIONAL_SIGNAL`（不合格为 STRONG）。**
> A1 是唯一在三项 primary 指标上**同向**改善的架构，但**统计上不可区分于 A0**。

---

## 10. Scientific Verdict

| 主张 | 判定 | 依据 |
|---|---|---|
| **Set context** | **NOT_SUPPORTED** | ΔMacro ρ +0.0174 CI [−0.0334,+0.0697] 含 0；ΔC_W +0.0061 含 0 |
| **Mutation-local representation** | **NOT_SUPPORTED** | ΔMacro ρ **−0.0275** CI [−0.0668,+0.0087]（点估计为负）；A4−A1 macro ρ **+0.0001** |
| **Target-balanced training** | **NOT_SUPPORTED** | L3 各项**更差**（Macro ρ 0.1118→0.0674，Eff@10 0.2522→0.1124） |
| **Within/pairwise objective** | **NOT_SUPPORTED** | L1/L2 相对 L0 无改善（ΔMacro ρ −0.002/−0.001） |
| **Physics rank supervision** | **NOT_RUN** | teacher coverage 9.1%，标 `UNAVAILABLE` |
| **Set×Local synergy** | **NOT_SUPPORTED** | 0/3 指标 CI>0 → `NO_SYNERGY_EVIDENCE` |

---

## 11. Main Result

**Family-OOD baseline（A0 = 历史 pooled M0 特征，同 row universe）**
- Overall r: **+0.2205**
- Macro ρ: **+0.0333**
- C_W: **+0.1734**
- Eff@10: **+0.1420**
- Row C_B / Equal C_B: +0.0472 / +0.5052
- atten slope: +0.2148

**Family-OOD V7 best（A1 = Base + Set）**
- Overall r: **+0.2408**
- Macro ρ: **+0.0507**
- C_W: **+0.1795**
- Eff@10: **+0.2080**
- Row C_B / Equal C_B: +0.0613 / +0.5823
- atten slope: +0.2029

**差值均不显著**（全 CI 含 0）；但 **OOD retention 明显改善**
（Retention_W 0.604 → **0.700**；Retention_ρ 0.305 → **0.628**）。

---

## 12. Paper Decision

| 项 | 判定 |
|---|---|
| **BENCHMARK CONTRIBUTION** | **YES** |
| **METHOD CONTRIBUTION** | **NOT YET** |
| **POTENTIAL BIOINFORMATICS-LEVEL** | **NOT YET** |

**Reason**：
- **Benchmark 可写**：① V7.0 给出**严格的 7 架构 × 3 seed × family-OOD 消融**，
  且全部 contrast 坦承不显著；② **weighting artifact 再次复现**（`C_B` 极差 0.556 且翻转），
  这是对 ΔΔG benchmark 汇报规范的有力方法学批评；③ 揭示 `≥30` 组的失效**只在 target 均值**
  （slope 0.02）而 **within 排序并不差**（C_W +0.215、Macro ρ +0.170）；
  ④ 记录了 BLOCK L 映射约定的一次真实修正（46.8% → 93.1%）
- **Method 尚不可写**：Set / Local / loss 变体**均无统计显著增益**，
  synergy 无证据，physics distillation 未跑通。
  按 §20，**不足以支撑方法学论文**。

---

## 13. Next Step

**结论：weak / 无明显增益。按 §19「If weak」分支处理 —— 不要盲目加大模型规模。**

必须诚实定位**哪一分支失败**：

1. **Set context 分支**：方向正确（三指标同向）但效应量太小（ΔMacro ρ +0.017）
   且 CI 宽（±0.05）。**唯一值得继续的方向。**
   下一步：把 set 描述符换成**更细的表征**（按位点分布而非集合标量），
   并**只在 `≥30` 子集上**评估（该组占 70% 的行、且 C_W 最高）。
2. **Local 分支失败**（ΔMacro ρ **为负**）。可能原因：
   - 局部特征是**静态几何量**（半径内接触数、partner 平均理化），
     与 ΔΔG 的关系可能被 pooled base 已含的物理项覆盖
   - **缺少 residue-level embedding**（`esm2_cache` 只有 pooled 640 维）→
     §7 要求的 `Δh_m = h_mut − h_wt` 与 `Δh_m * c_m` **从未实现**
   - 下一步：若要救 Local，需先解决 residue-level 表示（但不许重训大型 PLM）
3. **Loss 分支失败**：pairwise 与 centered 均无效 → **MLP + ranking loss 不是瓶颈**；
   瓶颈更可能在**特征**（与 A1/A3 的结论一致）
4. **Physics 分支未跑通**：需先取得覆盖 230 targets 的 physics 预测（当前 21），
   或改用**已有的 FoldX/物理 60 维特征**做软排序（base 已含 `phys_diff`）
5. **不要再做**：单纯加深网络 / 加树数 / 换 learner —— §2 已证明
   learner 差异（MLP vs LightGBM）量级大于全部架构差异

**V7 主论文假设（§22）本轮未被支持**：
> 「显式利用候选突变集合上下文与突变位点局部相互作用的模型，能否比 pooled ΔΔG 模型
> 更好地完成未见蛋白家族上的突变筛选。」
> → **当前证据：不能显著更好。**

---

## 14. 完成度（COMPLETED / PARTIAL / NOT_RUN / FAILED）

| PHASE | 状态 | 说明 |
|---|---|---|
| 0 冻结 + baseline | **COMPLETED** | 12 个输入哈希；B0 = 历史 M0（复用不重写） |
| 1 特征块 | **COMPLETED** | S(25) / G(58) / L(118) 分块构造 |
| 1b BLOCK L 修正 | **COMPLETED** | 映射约定修正 46.8% → **93.1%**；4473 行 |
| 2 架构消融 | **COMPLETED** | 7 模型 × 2 档 × 3 seed；L4 primary |
| 3 Loss 消融 | **COMPLETED** | L0–L3 四变体（MLP） |
| 4 Physics rank | **NOT_RUN** | teacher coverage **21/230 = 9.1%** → `UNAVAILABLE` |
| 5 Bootstrap | **COMPLETED** | 4 核心 contrast + synergy + weighting + 分层（2000/2000） |
| 6 表征分析 | **PARTIAL** | 已含 §6 分层、§7 weighting、§8 retention；缺 Set tertile / interface 分组细分 |
| 7/8 汇总 | **COMPLETED** | 本文件 |

**无 `logs/v7_0/failures.log`** → 无功能性失败。

---

## 15. 不允许的结论（§21 遵守情况）

本报告**未**使用以下任何表述：
- ❌「Candidate-set composition causes failure」→ 只写「improves predictive signal」且加 CI 限定
- ❌「MIPA is superior」→ A3/A4 **无** bootstrap 支持，故明确写 NOT_SUPPORTED
- ❌「Physics methods are generally superior」→ physics 部分 **NOT_RUN**
- ❌「Overall Pearson improvement = better protein engineering」→ 明确以
  Macro ρ / Eff@10 / C_W 为 primary，且指出 A2 的 r 低但 Macro ρ 高

---

## 16. 文件清单

### 脚本（`scripts/v7_0/`）
`01_features.py`、`02_local_features_v2.py`、`03_ablation.py`、`04_bootstrap.py`、
`05_losses.py`、`06_fast_boot.py`

### 结果（`experiments/v7_0/`）
- `phase0_baseline/`：`input_manifest.csv`、`baseline_definition.md`
- `phase1_features/`：`block_L_local.csv`、`block_S_set.csv`、`block_G_global.csv`、
  `matched_coverage_rows.csv`、`coverage.json`、`feature_report.md`、
  `local_feature_report.md`、`mapping_audit.csv`
- `phase2_ablation/`：`ablation_metrics.csv`、`ablation_predictions.npz`、`ablation_report.md`
- `phase3_losses/`：`loss_metrics.csv`、`loss_predictions.npz`、`loss_report.md`
- `phase4_physics_rank/`：**（空 — NOT_RUN）**
- `phase5_bootstrap/`：`contrasts.csv`、`synergy.csv`、`weighting.csv`、
  `stratification.csv`、`bootstrap_report.md`
- `V7_0_FINAL_SUMMARY.md`（本文件）

### 日志（`logs/v7_0/`）
`03_ablation.log`、`04_boot.log`、`05_loss.log`、`05_loss2.log`、`06_boot.log`
（**无 `failures.log`**）
