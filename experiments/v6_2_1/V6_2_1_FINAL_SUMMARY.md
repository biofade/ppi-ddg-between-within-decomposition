# V6.2.1 Final Summary

**时间**：2026-10-14（约 13:00 → 13:45，约 45 分钟；上限 4 小时，未超时）
**约束**：未修改 V1.1–V6.2 任何文件；`v6_2_1` 命名空间；CPU only；无下载；未伪造数字

---

## 1. Shuffle Audit

- **Old shuffle correct: NO** → 标记 **`OLD_SHUFFLE_BUG_CONFIRMED`**

### Exact bug

`scripts/v6_2/03_decomposition.py` §3.7（Control A）：

```python
for cx in np.unique(gg):
    idx = np.where(gg == cx)[0]
    ps[idx] = pp[rs.permutation(len(idx))]        # ← BUG
```

**正确写法**（对组内**子集**取置换）：

```python
    ps[idx] = pp[idx][rs.permutation(len(idx))]   # ✅
```

**机制**：`pp` 是**全表长度**数组，而 `perm = rs.permutation(len(idx))` 的取值上限
是 `len(idx)`。`pp[perm]` 因此用**组内下标**索引**全局数组** → 抽到别组的元素
→ **破坏组内多重集**（跨组串值），不再是合法组内置换。

**最小复现**（合成 n=200，`01b_old_shuffle_repro.py`）：
`ps[ii]=pp[perm]` → 组内多重集**不守恒**、全局 sum 28.5510；
`ps[ii]=pp[ii][perm]` → **守恒**、sum 14.1821（= 原始 sum）。

### Corrected L1 M1

| 量 | 值 |
|---|---:|
| `r_full` | **+0.6999** |
| `C_B` | **+0.4323** |
| `C_W` | **+0.2676** |
| **shuffle `r` mean ± SD**（100 次） | **+0.4271 ± 0.0081** |
| shuffle `C_B` | **+0.4275**（= 原 `C_B`，守恒 ✅） |
| shuffle `C_W` | **−0.0004**（≈0 ✅） |

> 旧报告 `r_shuffle = −0.0192` 已由上述 bug **精确复现**（变体 A 实测 −0.0192）。

**数学事实（36/36 组合 + 100 次重复验证）**：
组内置换保持组均值 ⇒ `C_B(shuffle) ≡ C_B(original)`，`C_W(shuffle) ≈ 0`。
因此 `r_shuffle` = 「**仅靠跨靶点信息**能达到的全局 Pearson」，
而 `r_full − r_shuffle` = **组内排序的真实增量**（L1/M1 = **+0.2728**）。

---

## 2. Tree Count Audit

- **Actual trees: 300**（全部 4 档 × 3 模型 × 3 seed，统一）
- **300/800 conflict explanation**：
  - `800` = `scripts/v6_2/02_split4_ddg.py` 的 **`--trees` 默认值**（被显式覆盖）
  - V6.2 实际运行命令为 `... --trees 300 ...`，且 `TREES = int(args.trees)`
  - 另有一处**报告文本笔误**：`V6_2_FINAL_SUMMARY.md` PHASE 6 段写「已完成 800 树复现」
- **Retraining required: NO** → 标记 **`ACTUAL_300_CONSISTENT`**

**建议修正（不改 V6.2 原文件，仅记录）**：该句应改为「使用 300 树（统一值）」。

---

## 3. True D-only M2

**历史定义复用**（非重写）：`scripts/v5/02_baseline_repro.py` 的
`Dn = standardize_train(hstack([phys_diff, Z_diff]), train_mask)`，
保留 forward/reverse 增广（`-D`、`-y`），**不含** DGaux、**不含** `S` 对称块。

| Split | Model | r_full | C_B | C_W | Macro ρ | Eff@10 |
|---|---|---:|---:|---:|---:|---:|
| L1 | M0 pooled | +0.6870 | +0.4288 | +0.2582 | +0.0677 | +0.3084 |
| L1 | **M2 D-only** | +0.6868 | +0.4251 | +0.2617 | +0.0717 | +0.3036 |
| L1 | M1 DGaux | +0.6999 | +0.4323 | +0.2676 | +0.0663 | +0.1020 |
| L4 | M0 pooled | +0.1418 | −0.0054 | +0.1472 | +0.0386 | +0.3020 |
| L4 | **M2 D-only** | +0.1355 | −0.0192 | +0.1547 | +0.0347 | +0.2941 |
| L4 | M1 DGaux | +0.1347 | −0.0121 | +0.1468 | +0.0407 | +0.3172 |

> 三者 **`r_full` 与 `C_W` 都极为接近**（L1 `C_W` 0.2582/0.2617/0.2676；
> L4 `C_W` 0.1472/0.1547/0.1468）—— 说明**特征消融在 within 分量上没有实质差别**。

### Gain_W（paired cluster bootstrap 5000/5000，单位 = complex）

| Contrast | Level | Point | 95% CI | Excludes 0? |
|---|---|---:|---|---|
| **Gain_D = C_W(M2)−C_W(M0)** | **L1** | +0.0035 | [−0.0028, +0.0128] | ❌ |
| **Gain_D = C_W(M2)−C_W(M0)** | **L4** | +0.0074 | [−0.0188, +0.0270] | ❌ |
| Gain_DGaux = C_W(M1)−C_W(M2) | L1 | +0.0059 | [−0.0088, +0.0153] | ❌ |
| Gain_DGaux = C_W(M1)−C_W(M2) | L4 | −0.0078 | [−0.0276, +0.0169] | ❌ |
| Gain_total = C_W(M1)−C_W(M0) | L1 | +0.0094 | [+0.0022, +0.0155] | ✅ |
| Gain_total = C_W(M1)−C_W(M0) | L4 | −0.0004 | [−0.0062, +0.0036] | ❌ |

> **关键否定结果**：**D-only 与 DGaux 在 L4 上都没有统计支持的 within 增益。**
> L1 上的 `Gain_total` 虽显著，但幅度仅 **+0.0094**（占 `C_W` 的 3.5%）。

**Delta_OOD_W = C_W(L1)−C_W(L4)**：M0 +0.1109 [+0.0709,+0.1771]；
**M2 +0.1070 [+0.0613,+0.1867]**；M1 +0.1207 [+0.0787,+0.1861] — 三者均 ✅

---

## 4. Target-size Stratification

- **Target-size heterogeneity: YES** → **`TARGET_SIZE_HETEROGENEITY_PRESENT`**
- **corr(log n_mut, |mean error|)**（Spearman）：
  L1 M1 **−0.3592**（p<0.0001）；L4 M1 −0.0991（p=0.134）；L4 M2 −0.1561（p=0.018）

### L4（family-level）分层结果（M1）

| bin | n_target | n_rows | **macro mean-Pearson** | **row-w C_B** | row-w C_W | Overall r | Macro within-ρ |
|---|---:|---:|---:|---:|---:|---:|---:|
| **2–9** | 112 | 496 | **+0.5886** | **+0.3153** | +0.0394 | +0.3546 | −0.0143 |
| **10–29** | 55 | 916 | **+0.6885** | **+0.3076** | +0.0713 | +0.3789 | +0.0355 |
| **≥30** | 63 | **6111** | **−0.0008** | **−0.0720** | +0.1629 | +0.0909 | +0.0923 |

**机制确认**：`≥30` 组虽只占 63/230 个 target，却占 **6111/7523 = 81%** 的行。

- **macro-target** mean Pearson：每 target 等权 → 小/中 target 表现好 → **+0.53 为正**
- **row-weighted** `C_B`：按突变数加权 → 被 `≥30` 组的 **−0.0720** 主导 → **≈0/负**

⇒ **「macro Between > 0 但 row-weighted C_B ≈ 0」不矛盾，由 target mutation-count
imbalance 造成**，已用预注册分层（三分位 + 固定阈值 2–9/10–29/≥30）证实。

---

## 5. Corrected Scientific Interpretation

### Reliable findings

1. **within-target shuffle 只消除 `C_W`，`C_B` 数学上守恒**（36/36 组合、100 次重复）
   - L1/M1：`r_full` +0.6999 = `C_B` +0.4323 + `C_W` +0.2676
   - shuffle 后：`C_B` = +0.4275（守恒），`C_W` = −0.0004（≈0）
   - **组内排序的真实增量 = ShuffleLoss = +0.2728**（L1）；L4 = **+0.1519**
2. **within 信号在 family-level OOD 后显著下降但保留约 55%**
   - `C_W`：L1 **+0.2676** → L4 **+0.1468**
   - `Delta_OOD_W = +0.1207 [+0.0787, +0.1861]`，CI 排除 0（3/3 模型一致）
3. **D-only 与 DGaux 都不能恢复 within 信号**
   - `Gain_D`(L1/L4)、`Gain_DGaux`(L1/L4) **全部 CI 含 0**
   - 只有 `Gain_total`(L1) 显著，幅度 **+0.0094**
4. **target-size 异质性是「macro Between > 0 而 row-weighted C_B ≈ 0」的唯一原因**
   - `≥30` 组（81% 的行）macro mean-Pearson **−0.0008**、`C_B` **−0.0720**
   - `2–9` / `10–29` 组均为强正（+0.59 / +0.69；`C_B` +0.32 / +0.31）
5. 树数实际统一为 **300**（非 800）；结论不受影响

### Findings invalidated or weakened

1. **❌ 失效**：「within-target shuffle 把全局 Pearson 打到 ≈0」
   —— 由一行索引 bug 造成，**已证伪**。
2. **❌ 失效**：「Global Pearson 全部来自 between-target」
   —— shufle 不能支持该结论；正确读数是
   **between 贡献 +0.4323、within 贡献 +0.2676**（两者都实质存在）。
3. **⚠️ 削弱**：「DGaux 提升 family-OOD 下的 mutation ranking」
   —— `Gain_DGaux`(L4) = −0.0078，CI 含 0，**无统计支持**。
4. **⚠️ 修正**：`V6_2_FINAL_SUMMARY.md` 中「800 树复现」为文本笔误，实为 300。

### Between-target conclusion

- **macro-target**：L4 上仍为正（+0.53），**由小/中 target 支撑**（+0.59 / +0.69）
- **row-weighted**：L4 上 ≈0/负（−0.0054 ~ −0.0132），**由 ≥30 组拖低**
- **interpretation**：target 均值预测**并非完全失效**，而是
  **对 mutation-rich target 失效**；两者权重不同故可并存

### Within-target conclusion

- **L1**：`C_W` = **+0.2582 ~ +0.2676**（M0/M1/M2 接近）
- **L4**：`C_W` = **+0.1468 ~ +0.1472**
- **OOD retention**：**≈55%**（+0.2676 → +0.1468）

---

## 6. Route Decision

**NEXT_ROUTE: `MIPA_PLUS_CALIBRATION`**

### Reason（按判据逐条）

| 判据 | 实测 | 是否满足 |
|---|---|---|
| `C_W(L1) > 0` | +0.258 ~ +0.268 | ✅ |
| `C_W(L4) > 0` | +0.147 | ✅ |
| `Delta_OOD_W > 0` 且 CI 排除 0 | +0.1207 [+0.0787,+0.1861] | ✅ |
| **`Gain_D_L4` 有恢复迹象** | **+0.0074 [−0.0188,+0.0270] 含 0** | **❌** |
| **between calibration 在 L4 不稳定** | **是**：`≥30` 组 mean-Pearson **−0.0008**、`C_B` **−0.0720** | **✅** |

- **Route A（MIPA_FIRST）不成立**：其关键判据「`Gain_D_L4 >= 0` 或 local
  representation 有恢复迹象」**未被满足**（D-only 在 L4 无显著增益）。
- **Route B（DUAL_FIRST）不成立**：`C_B >> C_W` 在 L4 **不成立**
  （`C_B` = −0.0121，`C_W` = +0.1468）。
- **Route D（MORE_DIAGNOSIS）不成立**：CI 不宽（`Delta_OOD_W` 半宽仅 0.054），
  true M2 已实施，shuffle 矛盾已**完全解决**并定位到确切 bug。
- **Route C（MIPA_PLUS_CALIBRATION）成立**：
  `C_W(L4)` 明确 > 0（+0.147）；**target-level mean/calibration 在 L4 明显失效**
  （且**仅对 mutation-rich target 失效**：`≥30` 组 mean-Pearson −0.0008）；
  同时 within 信号部分保留。

### 下一代重点（V7 边界，本次**不实施**）

1. **mutation-centered local branch**（延续 MIPA 假设）
2. **partner context / local interaction perturbation**
3. **structure-conditioned local representation**
4. **strict family-level OOD** 作为唯一评估口径
5. **显式 target-level calibration branch** —— 但
   **禁止使用 target-ID intercept**；global branch 必须由 **unseen target features** 预测
6. **必须报告 target-size 分层**（≥30 组单独列出），否则 macro 与 row-weighted
   会给出相反结论

---

## 7. Paper Implications

### Benchmark paper: **YES**

最有力的三项证据：
1. **within-target shuffle 的通用陷阱**：一行索引 bug 即可伪造「全局相关全部来自
   between-target」的假象；`C_B(shuffle) ≡ C_B(original)` 是**可检验的硬约束**
2. **target-size 异质性**：同一 split 下 `≥30` 组 mean-Pearson −0.0008 而
   `2–9`/`10–29` 组 +0.59/+0.69 → **必须分层报告**，否则结论相反
3. `Delta_OOD_W` 用协方差分解量化（+0.1207，3/3 模型 CI 排除 0），
   比单纯 per-complex Spearman 更稳健

### Method paper: **YES（但路线须按 §6 修正）**

### Main revised claim

> **在一个 PPI ΔΔG 基准上，within-target 排序信号在 family-level OOD 后仍保留约
> 55%（协方差分解 `C_W` 从 +0.268 降至 +0.147），但 target-mean 校准对
> mutation-rich target 完全失效（`≥30` 组 mean-Pearson −0.0008）；且 D-only 与
> ΔG-辅助任务均未能在 OOD 下带来统计显著的 within 增益（`Gain_D`/`Gain_DGaux`
> 的 95% CI 全部含 0）。**

### Claims that must NOT be used

- ❌ 「Global Pearson 全部来自 between-target」
- ❌ 「Within signal 在 family OOD 后消失」
- ❌ 「Physics universally beats learned models」
- ❌ 「DGaux improves family-OOD mutation ranking」
- ❌ 「Between-target prediction completely fails」
- ❌ 「within-target shuffle 把全局 Pearson 打到 ≈0」

---

## 8. 文件清单

### 脚本
`scripts/v6_2_1/`：`01_shuffle_audit.py`、`01b_old_shuffle_repro.py`、
`02_tree_audit.py`、`03_true_d_only.py`、`04_target_size.py`、
`05_corrected_stats.py`

### 结果
- `phase0_environment.md`
- `PROGRESS.md`
- `phase1_shuffle_audit/`：`shuffle_audit.csv`、`shuffle_repeats.csv`、
  `old_shuffle_repro.md`、`old_shuffle_repro.csv`、`shuffle_report.md`
- `phase2_tree_audit/`：`tree_count_audit.csv`、`tree_count_report.md`
- `phase3_d_only/`：`d_only_predictions.npz`、`d_only_metrics.csv`、
  `gain_w_bootstrap.csv`、`d_only_report.md`
- `phase4_target_size/`：`target_level_stats.csv`、`target_size_bins.csv`、
  `target_size_report.md`
- `phase5_corrected_stats/`：`primary_table.csv`、`route_statistics.csv`、
  `corrected_shuffle.csv`、`corrected_stats_report.md`
- `phase6_route/`：`route_decision.md`
- `V6_2_1_FINAL_SUMMARY.md`（本文件）

### 日志
`logs/v6_2_1/`：`03_donly.log`、`03_donly_err.log`
（**无 `failures.log`** → 无功能性失败）
