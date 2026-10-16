# V8 External Universe 2：`U_EXT_HELDOUT_CATH`

> **来源**：CATH-ddG 官方划分列 `cath_fold`，`val` = **813 行 / 53 个 PPI**
> （论文：「54 CATH 超家族 / 53 PPI / 813 突变」，序列同一性 < 30%，
> **不含 easy mutation**，最大 TM-score < 0.6）。
>
> **与 V8 的关系**：这是**独立于本项目手工家族划分**的严格外推测试集
> （CATH 超家族口径），且**自带 5 个基线方法的逐突变预测** ⇒ **零推理**。

---

## 1. 宇宙事实

| 项 | 值 |
|---|---:|
| 行数 | **813** |
| PPI（target 单位 = `complex_PPI`） | **53** |
| PDB | 53 |
| **n≥5 的 target** | **21** |
| 单点突变行 | 543 |
| 多点行 | 270（最多 27 点） |
| 全部 813 行均非 `train` | ✅ |

| 基线模型 | 可用行 |
|---|---:|
| `Flex-ddG` | 799 |
| `DDGPred` / `ddG_FoldX` / `ESM-IF` / `proteinMPNN` | 813（全部） |

---

## 2. 三轴读数（102–103 个内部 target 口径不同；此处 21 个 n≥5 target）

| 模型 | n | **Overall r** | RMSE | **Macro ρ** | C_W | Eff@10 |
|---|---:|---:|---:|---:|---:|---:|
| `DDGPred` | 813 | **+0.6128** | 2.1814 | +0.3262 | +0.3070 | +0.8598 |
| `Flex-ddG` | 799 | **+0.6022** | **2.1077** | **+0.4159** | **+0.3484** | **+0.9037** |
| `ddG_FoldX` | 813 | +0.4603 | 2.5507 | +0.3584 | +0.2968 | +0.8963 |
| `proteinMPNN` | 813 | +0.1837 | 3.0788 | +0.2227 | +0.0763 | +0.8697 |
| `ESM-IF` | 813 | +0.1421 | 3.0969 | **−0.0329** | +0.1241 | +0.7637 |

---

## 3. ★ 一个重要且必须谨慎的观察

**在这个官方 held-out CATH 超家族测试集上，外部基线的靶内排序能力
（Macro ρ 最高 **+0.4159**）明显高于 V8 内部模型在 `U_CANONICAL` 上的
（最高 **+0.3388**，中位 +0.05）。**

**但我明确不把这解读为「外部方法更强」**，因为两个宇宙在以下方面**不可比**：

| 维度 | `U_CANONICAL`（内部） | `U_EXT_HELDOUT_CATH`（外部） |
|---|---:|---:|
| 行数 | 4473 | 813 |
| target 定义 | `target_id_canonical`（手工家族） | `complex_PPI`（PPI 单位） |
| 总 target 数 | 214 | 53 |
| **n≥5 的 target** | **130** | **21** |
| 多点突变占比 | 35% | 33% |
| 模型来源 | 本项目 9 个族 | 5 个外部基线 |

**⇒ 依 §21「不得比较不同宇宙的原始数值」，此处只报告事实，不做强弱结论。**
**⇒ 依 §71，不得表述为「physics/统计法优于学习法」或「外部方法更强」。**

一个**中性的**可报告事实：**`ESM-IF` 在两个宇宙上都表现出靶内排序能力接近于零或为负**
（内部 Macro ρ = −0.0346；外部 = −0.0329）。这是**跨宇宙一致**的，
因此比单宇宙结论更可信。

---

## 4. 秩相关：仍然 **NOT_TESTABLE**

```
U_EXT_HELDOUT_CATH (n=5): rho(Overall r, Macro rho) = +0.7000
                          rho(Overall r, Eff@10)    = +0.4000
                          rho(Macro rho, Eff@10)    = +0.9000
U_EXT_CANONICAL    (n=5): rho(Overall r, Macro rho) = +0.7000
                          rho(Macro rho, Eff@10)    = +0.3000
```

**n = 5 时 Spearman 的 95% 区间约 [−0.4, +0.99] ⇒ 无判别力。**

```text
EXTERNAL_MODEL_REPRODUCTION = NOT_TESTABLE
   两个外部宇宙各只有 5 个可评价模型（< 8 门槛）
   ⇒ 仍无法检验「外部模型是否复现能力分化」
   ⇒ §68 STRONGLY_SUPPORTED 条件**不成立**（原因：不可检，非阴性）
```

---

## 5. 尝试扩池的失败记录（诚实）

| 尝试 | 结果 |
|---|---|
| CATH-ddG 的 `HER2.csv`（419 行） | ❌ **无基线预测列**（仅 `ddG` + 划分标签） |
| CATH-ddG 的 `S285.csv`（285 行） | ❌ 同上 |
| `prompt-fold0/1/2` 三列 | ❌ 原样为 `train`/`val` 字串，**不是数值** |
| PPIformer 官方逐突变预测 | ❌ 未发布（需自跑推理；MIT + Zenodo 权重可得） |
| RDE-Network 官方逐突变预测 | ❌ 未发布（Apache-2.0 + 权重可得，需自跑） |
| Prompt-DDG 官方逐突变预测 | ❌ 未发布（MIT + Google Drive 权重） |

**⇒ 达到 ≥8 个外部模型的唯一可行路径是自跑推理**，而三条路线
（PPIformer / RDE-Network / Prompt-DDG）**均声明 CPU 可行**（README 或
`--device` 参数），**无需 GPU**。

**本轮不执行推理**：`pytorch 2.5.1+cpu` 已装，但 CATH-ddG 环境要求
`pytorch 2.3.1 / py3.10 / cudatoolkit 11.3.1`，且本机**可用 RAM 仅约 1.1–5.6 GB**
（CATH-ddG 论文称 A100 上需 ~19.7 GB）⇒ **本地 CPU 推理的可行性未被证明**。

---

## 6. 复现命令

```powershell
python -W ignore -u scripts/v8/40_external_heldout_cath.py
```
产出：`external/external_heldout_cath.csv`、`external/external_heldout_cath_meta.json`
