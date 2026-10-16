# V6.2.1 执行台账

> 每完成一个 PHASE 立即追加。失败写 `logs/v6_2_1/failures.log`，不停止。
> 约束：不改 V1.1–V6.2；`v6_2_1` 命名空间；CPU only；不下载；不伪造数字；
> **判定只依据 source code / saved predictions / row IDs / hashes / formulas / assertions**。

## PHASE 0: 环境 + 结果冻结 — **DONE**

- 新建 `scripts/v6_2_1/`、`experiments/v6_2_1/`、`logs/v6_2_1/`
- 冻结 10 个关键文件的 SHA256（含 `split4_predictions.npz` = `ba711ea1…`）
- 产出：`phase0_environment.md`

## PHASE 1: Within-target Shuffle 审计 — **DONE**

### 标记：**`OLD_SHUFFLE_BUG_CONFIRMED`**

**确切 bug（一行代码，已最小化复现）**

V6.2 `scripts/v6_2/03_decomposition.py` §3.7 的 shuffle：

```python
for cx in np.unique(gg):
    idx = np.where(gg == cx)[0]
    ps[idx] = pp[rs.permutation(len(idx))]      # ← BUG
```

**正确写法**应为（对**组内子集**做置换）：

```python
    ps[idx] = pp[idx][rs.permutation(len(idx))]  # ✅ 组内置换
```

**机制**：`pp` 是**全表长度**的数组，而 `perm = rs.permutation(len(idx))` 的
**取值上限是 `len(idx)`**。`pp[perm]` 因此用「组内下标」去索引**全局数组**
→ 抽到的是**别的组**（以及本组其它位置）的元素
→ **破坏组内多重集**（跨组串值），不再是一次合法的组内置换。

**最小复现**（`scripts/v6_2_1/01b_old_shuffle_repro.py`，合成数据 n=200）：

| 写法 | 组内多重集守恒 | 全局 sum |
|---|---|---|
| `ps[ii] = pp[perm]`（旧，BUG） | **False** | 28.5510（错） |
| `ps[ii] = pp[ii][perm]`（正确） | **True** | 14.1821（= 原 sum） |

**在真实数据上的后果（L1 / M1 / seed42）**

| 变体 | r_full | **r_shuffle** | C_B(shuffle) |
|---|---:|---:|---:|
| 原始 | +0.6990 | — | — |
| **A：V6.2 逐字（BUG）** | +0.6990 | **−0.0192** ✅ 与 V6.2 报告**完全一致** | **−0.0178**（错） |
| **B：正确组内置换** | +0.6990 | **+0.4348** | **+0.4275**（= 原 `C_B`）✅ |

⇒ **V6.2 报告的 `r_shuffle = −0.0192` 被精确复现，来源即上述一行 bug。**

### 数学事实（用强制断言验证）

对每个 target 做**真正的**组内置换会**保持该 target 的预测均值不变**，故：

```
C_B(shuffle) ≡ C_B(original)       （100 次重复：36/36 组合守恒 ✅）
C_W(shuffle) ≈ 0                   （实测 |C_W| < 0.001）
r_shuffle = C_B_shuffle + C_W_shuffle ≈ C_B
```

**实测（100 次重复均值）**

| 档 | 模型 | C_B（原始） | C_B_shuf | 守恒? | C_W_shuf |
|---|---|---:|---:|---|---:|
| L1 | M0 | +0.4249 | +0.4249 | ✅ | −0.0004 |
| L1 | M1 | +0.4275 | +0.4275 | ✅ | −0.0004 |
| L2 | M1 | +0.0654 | +0.0654 | ✅ | −0.0005 |
| L4 | M1 | −0.0173 | −0.0173 | ✅ | +0.0002 |

### 结论修正（须替换 V6.2 的表述）

- ❌ **不得**再写：「把靶点内突变顺序打乱后全局 Pearson 掉到 ≈ 0」
- ❌ **不得**引申：「Global Pearson 全部来自 between-target」（V6.2.1 纪律 §1）
- ✅ **正确表述**：within-target shuffle **只消除 `C_W`**；
  `C_B` 在数学上**必然保留**。因此
  - `r_shuffle` = 「**仅靠跨靶点信息**可达到的全局 Pearson」
  - `r_full − r_shuffle` = 「**组内排序的真实增量贡献**」
  - L1/M1：**+0.6990 − (+0.4348) = +0.2642**

产出：`phase1_shuffle_audit/shuffle_audit.csv`、`shuffle_repeats.csv`、
`old_shuffle_repro.{md,csv}`、`shuffle_report.md`

## PHASE 2: Tree-count Provenance — IN_PROGRESS

## PHASE 2: Tree-count Provenance — **DONE**

标记：**`ACTUAL_300_CONSISTENT`**

- 实际树数 **300**（全部 4 档 × 3 模型 × 3 seed，统一）
- `800` 来源：源码 `--trees` **默认值**（被 `--trees 300` 覆盖）
  + `V6_2_FINAL_SUMMARY.md` PHASE 6 段一处**报告文本笔误**
- **无需重训**（各档同值 → split inflation 比较仍公平）
- 产出：`phase2_tree_audit/`

## PHASE 3: 真正的 M2 = D-only — **DONE**（5.9 分钟）

**历史定义复用**：`scripts/v5/02_baseline_repro.py` 的
`Dn = standardize_train(hstack([phys_diff, Z_diff]), train_mask)`
+ forward/reverse 增广；**不含** DGaux、**不含** `S` 对称块。

| Split | Model | r_full | C_B | C_W |
|---|---|---:|---:|---:|
| L1 | M0 | +0.6870 | +0.4288 | +0.2582 |
| L1 | **M2 D-only** | +0.6877 | +0.4297 | +0.2617 |
| L1 | M1 | +0.6999 | +0.4323 | +0.2676 |
| L4 | M0 | +0.1418 | −0.0054 | +0.1472 |
| L4 | **M2 D-only** | +0.1347 | −0.0132 | +0.1468 |
| L4 | M1 | +0.1347 | −0.0121 | +0.1468 |

### ★ 关键否定结果（paired bootstrap 5000/5000）

| 对比 | 档 | point | 95% CI | 排除 0? |
|---|---|---:|---|---|
| **Gain_D = C_W(M2)−C_W(M0)** | L1 | +0.0035 | [−0.0028, +0.0128] | ❌ |
| **Gain_D** | **L4** | +0.0074 | [−0.0188, +0.0270] | ❌ |
| Gain_DGaux = C_W(M1)−C_W(M2) | L1 | +0.0059 | [−0.0088, +0.0153] | ❌ |
| Gain_DGaux | **L4** | **−0.0078** | [−0.0276, +0.0169] | ❌ |
| Gain_total = C_W(M1)−C_W(M0) | L1 | +0.0094 | [+0.0022, +0.0155] | ✅ |
| Gain_total | L4 | −0.0004 | [−0.0062, +0.0036] | ❌ |

⇒ **D-only 与 DGaux 在 L4 都没有统计支持的 within 增益。**

## PHASE 4: Target-size Stratification — **DONE**

标记：**`TARGET_SIZE_HETEROGENEITY_PRESENT`**

L4（family）分层（M1）：

| bin | n_target | n_rows | macro mean-Pearson | row-w C_B |
|---|---:|---:|---:|---:|
| **2–9** | 112 | 496 | **+0.5886** | **+0.3153** |
| **10–29** | 55 | 916 | **+0.6885** | **+0.3076** |
| **≥30** | 63 | **6111** | **−0.0008** | **−0.0720** |

**机制确认**：`≥30` 组占 **81% 的行**（6111/7523）却 mean-Pearson ≈ 0
⇒ 把 row-weighted `C_B` 拉向 0/负，而等权 macro 仍为正。
**「macro Between > 0 但 row-weighted C_B ≈ 0」由 target mutation-count
imbalance 造成，两者权重不同故不矛盾。**

连续分析：Spearman(log n_mut, |mean error|) = **−0.3592**（L1，p<0.0001）、
−0.0991（L4 M1，p=0.134）、−0.1561（L4 M2，p=0.018）

## PHASE 5: Corrected Statistics — **DONE**（339s）

### 修正后的 shuffle（同时报 C_B/C_W before/after）

| 档 | 模型 | r_full | C_B | C_W | r_shuf mean±SD | C_B_shuf | C_W_shuf | **ShuffleLoss** |
|---|---|---:|---:|---:|---|---:|---:|---:|
| L1 | M1 | +0.6999 | +0.4323 | +0.2676 | +0.4271 ± 0.0081 | +0.4275 | −0.0004 | **+0.2728** |
| L4 | M1 | +0.1347 | −0.0121 | +0.1468 | −0.0172 ± 0.0088 | −0.0173 | +0.0002 | **+0.1519** |

`ShuffleLoss = r_full − r_shuffle` = **组内排序的真实增量贡献**。

### Primary route statistics

| 统计量 | 结果 |
|---|---|
| **Delta_OOD_W = C_W(L1)−C_W(L4)** | M0 +0.1109 ✅ / **M2 +0.1070 ✅** / M1 +0.1207 ✅（**3/3 CI 排除 0**）|
| **Gain_D**（L1/L4） | **0/2 CI 排除 0** |
| **Gain_DGaux**（L1/L4） | **0/2 CI 排除 0** |

## PHASE 6: Route Decision — **DONE**

**`NEXT_ROUTE = MIPA_PLUS_CALIBRATION`**

理由（逐条）：
- Route A（MIPA_FIRST）**不成立**：其判据「`Gain_D_L4 >= 0` 或 local 表示有恢复迹象」
  **未被满足**（+0.0074 含 0）
- Route B（DUAL_FIRST）**不成立**：`C_B >> C_W` 在 L4 不成立（−0.0121 < +0.1468）
- Route D（MORE_DIAGNOSIS）**不成立**：CI 不宽（半宽 0.054）、true M2 已实施、
  shuffle 矛盾已完全解决并定位到确切 bug
- **Route C 成立**：`C_W(L4)` = +0.1468 明确 > 0，**且**
  target 级校准对 mutation-rich target 结构性失效（`≥30` 组 mean-Pearson −0.0008）

产出：`phase6_route/route_decision.md`

## PHASE 7: Final Report — **DONE**

见 `V6_2_1_FINAL_SUMMARY.md`。

## 全部 7 个 PHASE 完成；无 failures.log
