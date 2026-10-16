# V6.2.1 PHASE 6 — Route Decision

**判定**：**`NEXT_ROUTE = MIPA_PLUS_CALIBRATION`**

**依据**：PHASE 1–5 的可复现证据（source code / saved predictions / row IDs /
hashes / formulas / assertions），**不含主观判断**。

---

## 1. 判据逐条核对

| Route A（MIPA_FIRST）判据 | 实测 | 满足? |
|---|---|---|
| `C_W(L1) > 0` | **+0.2582 ~ +0.2676**（M0/M2/M1） | ✅ |
| `C_W(L4) > 0` | **+0.1468 ~ +0.1472** | ✅ |
| `Delta_OOD_W > 0` 且 CI 排除 0 | **+0.1207 [+0.0787, +0.1861]**（M1）；M0/M2 亦 ✅ | ✅ |
| `Gain_D_L4 >= 0` 或 local representation 有恢复迹象 | **+0.0074 [−0.0188, +0.0270]** → **含 0** | **❌** |
| between calibration 在 L4 不稳定 | **是**：`≥30` 组 mean-Pearson **−0.0008**、`C_B` **−0.0720** | ✅ |

| Route B（DUAL_FIRST）判据 | 实测 | 满足? |
|---|---|---|
| `C_B >> C_W` 在 L1 **与** L4 都成立 | L1：+0.4323 vs +0.2676（≈1.6×）→ 勉强；**L4：−0.0121 vs +0.1468 → `C_B < C_W`** | **❌** |
| within 很弱 | L4 `C_W` = **+0.147**，非「很弱」 | ❌ |

| Route C（MIPA_PLUS_CALIBRATION）判据 | 实测 | 满足? |
|---|---|---|
| `C_W(L4)` clearly > 0 | **+0.1468**（`Delta_OOD_W` CI 下界 +0.0787 → 明确非 0） | ✅ |
| target-level mean prediction/calibration 在 L4 明显失效 | **是**，且**仅对 mutation-rich target**：`≥30` 组 mean-Pearson **−0.0008** | ✅ |
| 同时 MIPA / D-only 可进一步提升 `C_W` | **否** —— `Gain_D`/`Gain_DGaux` 全部 CI 含 0 | ⚠️ **部分** |

| Route D（MORE_DIAGNOSIS）判据 | 实测 | 满足? |
|---|---|---|
| CI 太宽 | **否**：`Delta_OOD_W` 半宽 **0.054**（相对 point 0.121） | ❌ |
| true M2 没增益 | **是**（+0.0074 含 0） | ✅ |
| corrected shuffle / calibration 出现**新矛盾** | **否**：矛盾已**完全解决**并定位到确切 bug | ❌ |

---

## 2. 为什么是 Route C 而不是 A 或 D

### 不是 A（MIPA_FIRST）
Route A 明确要求「`Gain_D_L4 >= 0` **或** local representation 有恢复迹象」。
实测 **D-only 在 L4 的 within 增益 CI 含 0**（+0.0074 [−0.0188, +0.0270]），
**DGaux 更是负值**（−0.0078 [−0.0276, +0.0169]）。
⇒ 「靠 mutation-centered 表示就能恢复 OOD within 信号」这一前提**未被证据支持**。
（注意：这**不否定** MIPA 方向本身 —— `C_W(L4) = +0.147` 说明
within 信号**确实存在**；只是**当前实现**未能提升它。因此 MIPA 应作为
**其中一个分支**保留，而非唯一路线。）

### 不是 B（DUAL_FIRST）
`C_B >> C_W` 在 L4 **完全不成立**（`C_B = −0.0121 < C_W = +0.1468`）。

### 不是 D（MORE_DIAGNOSIS）
- CI 不宽（半宽 0.054）
- true M2 已实施
- shuffle 矛盾**已完全解决**，且定位到**确切的一行索引 bug**
⇒ 没有「需要更多诊断」的悬置问题。

### 是 C（MIPA_PLUS_CALIBRATION）
两条独立且都**有统计支持**的事实同时成立：
1. **within 信号存在**：`C_W(L4)` = +0.1468，`Delta_OOD_W` CI 排除 0（3/3 模型）
2. **target-level calibration 在 L4 失效，且失效是结构性的**：
   `≥30` 组（占 81% 的行）mean-Pearson **−0.0008**、`C_B` **−0.0720**，
   而 `2–9`/`10–29` 组为 **+0.5886 / +0.6885**（`C_B` +0.3153 / +0.3076）

⇒ 需要**两个分支**：一个修 target 级校准（global/context），
一个继续做 mutation-centered 局部表示（residual）。

---

## 3. V7 边界（本次**不实施**）

**科学问题**：
> 在 strict family-level OOD 下，如何同时 (a) 让 target-level 均值预测对
> mutation-rich target 也成立，(b) 保留并强化 mutation-specific within 信号？

**要求**：
1. **mutation-centered local environment**（延续 MIPA）
2. **partner context / local interaction perturbation**
3. **structure-conditioned local representation**
4. **strict family OOD 作为唯一评估口径**
5. **显式 global calibration branch** ——
   ⚠️ **禁止使用 target-ID intercept**；
   global branch **必须由 unseen target features 预测**
6. **强制 target-size 分层报告**（`≥30` 组必须单独列出），
   否则 macro 与 row-weighted 会给出**相反**结论
7. 评估必须**同时**报告 `C_B`、`C_W`、macro `ρ`、row-weighted 版本
   与 **ShuffleLoss**（= `r_full − r_shuffle`，即组内排序真实增量）

---

## 4. 复现命令

```powershell
python -u scripts/v6_2_1/01_shuffle_audit.py
python -u scripts/v6_2_1/01b_old_shuffle_repro.py
python -u scripts/v6_2_1/02_tree_audit.py
python -u scripts/v6_2_1/03_true_d_only.py
python -u scripts/v6_2_1/04_target_size.py
python -u scripts/v6_2_1/05_corrected_stats.py
```
