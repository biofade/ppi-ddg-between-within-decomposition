# V6.3 PHASE 7 — Mechanism Decision

**MECHANISM = `LIBRARY_COMPOSITION_DRIVEN`**

> 依据 V6.3 §PHASE 7 的判据表，逐条核对（证据全部来自 PHASE 1–6 的可复现输出）。

---

## 1. 判据逐条核对

| Mechanism | 判据 | 实测 | 满足? |
|---|---|---|---|
| **A** `LIBRARY_COMPOSITION_DRIVEN` | `n_mut` 单变量显著 | ρ(log n_mut, \|err\|) = **−0.0268, p=0.705** | ❌ |
| | 控制 mutation-set diversity 后作用明显减弱 | 单变量本就不显著；控制后**转正**（+0.142 ~ +0.191） | — |
| | **set descriptors 显著改善 calibration** | **ΔPearson(TS−T) = +0.2877, 95% CI [+0.1237, +0.4504]** ✅ | **✅** |
| **B** `TARGET_CONTEXT_DRIVEN` | target/global structural properties 解释 error | family-CV **Pearson = −0.0438**（无预测力） | ❌ |
| | set composition 增益弱 | 增益 **+0.2877**（不弱） | ❌ |
| **C** `BOTH_TARGET_AND_SET` | target/global 有作用 **且** set 有独立增益 | target 侧 Pearson = −0.044（**无作用**） | ❌ |
| **D** `FAMILY_OR_SOURCE_DRIVEN` | 由少数 family/source 强烈驱动、LOO 后消失 | **NOT_SINGLE_FAMILY_DRIVEN**；source **UNAVAILABLE** | ❌ |
| **E** `UNRESOLVED` | 各模型 CV 无稳定解释力、无变量可解释 | TS 有显著正向增益；且有明确 attenuation 机制 | ❌ |

---

## 2. 为什么是 A，以及 A 的**精确含义**

### 支持 A 的唯一强证据
**set descriptors 显著改善 calibration**：
- Model T（23 个 target-only 特征）：family-CV Pearson = **−0.0438**
- Model TS（+20 个 set 特征）：family-CV Pearson = **+0.2461**
- **ΔPearson = +0.2877，95% CI [+0.1237, +0.4504]（排除 0）**
- `≥30` 子集：T **+0.0878** → TS **+0.2614**

### 必须精确限定「composition」的含义（V6.3 关键澄清）

> 这里的 **composition 指 mutation-set 的组成描述符**
> （unique sites / site entropy / mutation-pair entropy / coverage /
> residue composition / single-multi 比例），
> **不是** `n_mutations`。

**`n_mutations` 既不是原因，也不是有效代理**：
1. 它与 `|mean error|` 的 Spearman = **−0.0268（p = 0.705）**，**不显著**
2. family-CV 下 size-only 模型 R² = **−0.0226**（无预测力）
3. V6.2.1 报告的 L4 值 −0.0991（p=0.134）亦不显著
   （其最终汇总中的 −0.3592 实为 **L1** 值，属文本串档，已更正）

---

## 3. 真正的机制链（本审计最可靠的发现）

```
观测现象（V6.2.1）：
  L4 family-OOD 下，≥30 mutation targets 的 target-mean Pearson ≈ 0
  且其占 81% 的行 → row-weighted C_B ≈ 0/负

PHASE 5（方法学分解）：
  「C_B ≈ 0/负」主要是 weighting 产物：
     row-weighting  → C_B = −0.0121
     equal-target   → C_B = +0.5475
     （极差 0.5596，符号翻转；而 macro mean-Pearson 恒为 +0.5568）
  ⇒ 删除 top-20 大靶点后 row-w C_B 从 −0.012 **升至** +0.138
     ⇒ 大靶点是**拖累**，不是支撑

PHASE 1/6（机制分解）：
  ≥30 组的失效模式是 **attenuation（压平）**，不是误差变大：
     slope(pred_mean ~ true_mean):
       2–9    = 0.2486  (Pearson +0.378)
       10–29  = 0.3522  (Pearson +0.747)
       ≥30    = −0.0013 (Pearson −0.003)   ← 预测几乎是常数 ≈ 0.59
     pred_std = 0.621 vs true_std = 1.480（严重压缩）

PHASE 2（关联定位）：
  与 |error| 相关的是 **真值水平与效应量**（LABEL_DERIVED）：
     true_mean_ddg  ρ = +0.413  (p=9.1e-10)
     ddg_std        ρ = +0.306  (p=8.7e-06)
  可部署项中唯一 BH 显著者：
     struct_chainpair_contact_density  ρ = +0.223  (q=0.0146)
  而 n_mutations：ρ = −0.0268（不显著）

PHASE 6（可部署性检验）：
  只用 target 结构上下文 → Pearson −0.044（**完全无预测力**）
  加入 mutation-set 组成   → Pearson +0.246（ΔPearson +0.288，CI 排除 0）
```

**结论**：校准失效的**可部署解释**来自 **mutation-set 组成**；
而**最可解释的**变量（真值水平、效应量）是不可部署的。
两者结合 → V7 的 global branch **必须**同时用到 set 组成，
且**不能**指望只靠靶点结构解决。

---

## 4. 对 V7 的直接约束

**Global branch = `G(x_c, M_c)`**（不是 `G(x_c)`）

| 约束 | 依据 |
|---|---|
| 必须输入 mutation-set 组成 `M_c` | `G(x_c)` family-CV Pearson = −0.044；加 `M_c` 后 +0.246（Δ+0.288, CI 排除 0） |
| 禁止 target-ID intercept / family lookup | 会破坏 family-OOD 评估语义 |
| **必须固定并声明 weighting 口径**，同时报 row 与 equal-target 两套 `C_B` | `C_B` 极差 0.5596 且符号翻转 |
| **必须报告 attenuation slope** `slope(pred_mean ~ true_mean)` | 它比 `abs_mean_error` 更能定位失效模式（≥30 组 slope ≈ 0） |
| `≥30` 组必须单独报告 | 占 81% 的行 |
| 优先级：**先修 target-level 压平**，再谈 within 提升 | `C_W` 在 L4 仍保留 +0.147，而 row-w between 为 0 |

**Local branch = `R_local(x_c, m_i)`**（不变）

---

## 5. 复现命令

```powershell
python -u scripts/v6_3/01_build_target_table.py
python -u scripts/v6_3/02_univariate_and_multivariable.py
python -u scripts/v6_3/03_influence_and_weighting.py
python -u scripts/v6_3/04_set_calibration.py
```
