# V8 External Multi-Model Comparison（PHASE 13 定稿）

> **资源**：CATH-ddG 仓库附带 `data/SKEMPI2/skempi_v2.csv`（1,291,260 B / 7085 行 / 22 列），
> 含 **逐突变预测**：`Flex-ddG, DDGPred, ddG_FoldX, prompt-fold0/1/2, ESM-IF, proteinMPNN`
> ＋真值 `ddG` ＋官方划分 `cath_fold`（5893 train / 813 val）。
> **§9 reuse-first 的最佳情形：零推理、零 GPU。**

---

## 1. 匹配宇宙 `U_EXT_CANONICAL`

**匹配键**：`(PDB4, 位置, 突变氨基酸)`。

> ⚠️ 首轮用完整突变串匹配得 **0 命中** —— 因为 CATH-ddG 的 `Mutation(s)_cleaned`
> 形如 **`LI38G`**（= 链字母 L + `I38G`），而 canonical 的 tok 是 `I38G`。
> **归一化后（丢弃链字母与 wt、只留位置+突变）才匹配上。**
> 这与上一轮的 snapshot gap 同源，也解释了 HERMES × 本表仅 11.4% 命中。

| 项 | 值 |
|---|---:|
| canonical 行 | 4473（严格 1:1） |
| canonical 单点行 | 2897 |
| **至少一个外部模型命中** | **3822 / 4473 = 85.45%** |
| 覆盖 target | **200 / 214** |

| 外部模型 | 命中 | 匹配率 |
|---|---:|---:|
| `Flex-ddG` | 3756 | 0.8397 |
| `DDGPred` | 3785 | 0.8462 |
| `ddG_FoldX` | 3819 | 0.8538 |
| `ESM-IF` | 3819 | 0.8538 |
| `proteinMPNN` | 3819 | 0.8538 |
| `prompt-fold0/1/2` | **0** | **0.0000** ⚠️ |

**⚠️ `prompt-fold*` 三列全为 0 命中**，且其中含 `train`/`val` 字串
（原样为划分标签而非数值）⇒ **该列不是预测值**，不可用。

---

## 2. 三轴读数（每个外部模型，其自身匹配子集内，102–103 个 n≥5 target）

| 模型 | n | Overall r | RMSE | **Macro ρ** | C_W | Eff@10 |
|---|---:|---:|---:|---:|---:|---:|
| `DDGPred` | 3785 | **+0.3357** | 2.5402 | +0.1047 | **+0.1491** | +0.8177 |
| `ddG_FoldX` | 3819 | +0.2727 | 2.6372 | **+0.1502** | +0.1213 | +0.8114 |
| `Flex-ddG` | 3756 | +0.2619 | 2.5649 | +0.1320 | +0.1197 | +0.8150 |
| `proteinMPNN` | 3819 | +0.1674 | 2.9312 | +0.0919 | +0.0729 | +0.8129 |
| `ESM-IF` | 3819 | **−0.1014** | 2.9607 | **−0.0346** | **−0.0592** | +0.7767 |

**可解读的三点：**
1. **物理/统计法（FoldX、Flex-ddG、DDGPred）整体优于结构条件化逆折叠法**
   （ESM-IF、proteinMPNN）——但**这是 6 个方法、1 个数据集上的描述性观察**，
   依 §71 **不得**表述为「physics generally beats learning」。
2. **`ESM-IF` 三项全负**（r、macro ρ、C_W 皆 < 0）⇒ 其预测方向性有问题，值得单独注意。
3. **`Macro ρ` 与 `Overall r` 在本表内同向**（见第 3 节限制）。

---

## 3. ★ 硬限制：**外部模型池太小，无法检验能力分化**

外部模型池 **n = 5**：
```
rho(Overall r, Macro rho) = +0.7000
rho(Overall r, Eff@10)    = +0.7000
rho(Macro rho, Eff@10)    = +0.3000
```

**但 n=5 时 Spearman 秩相关的 95% 区间极宽，无判别力。**
（n=5 的 Spearman 在 ρ=0.7 时的近似 95% CI 约 **[−0.4, +0.99]**。）

```text
EXTERNAL_MODEL_REPRODUCTION = NOT_TESTABLE
   原因：外部模型池 n=5（< 8 的最小可推断门槛）
   ⇒ 无法在外部池内检验「校准轴与排序轴是否分化」
   ⇒ §68 STRONGLY_SUPPORTED 的条件**仍不成立**（但不是因为阴性，而是因为不可检）
```

**同时必须记录**：`U_EXT_CANONICAL` 上外部模型的 `Macro ρ` **全为正**且与 `Overall r`
同向，**这与 V8 在内部模型上观测到的「弱一致/不稳定」并不矛盾，但也不支持它**。

---

## 4. §50 污染审计（逐方法）

| 方法 | 训练数据 | 与 canonical 的重叠风险 |
|---|---|---|
| `Flex-ddG` / `DDGPred` / `ddG_FoldX` | 各自原始论文的数据；**本 CSV 只是他人复算的结果表** | 需查原论文；本表未提供 |
| `ESM-IF` / `proteinMPNN` | **逆折叠模型，非 ΔΔG 监督** ⇒ 用作评分函数 | 无标签重叠 |
| `prompt-fold*` | — | **列不可用（非数值）** |
| **CATH-ddG 自身** | **SKEMPI v2.0**，5893/6706 = **88% 行属其训练集** | **其自身逐突变预测未发布**（补充材料只有 per-protein 表 S5–S9）⇒ 不在本表 |

**⇒ 本表不含任何在 canonical 上训练过的模型**，故 **5 个外部模型可进 TIER 1/2**；
但 `Flex-ddG/DDGPred/FoldX` 的原始训练集需逐个核实（**未做**）。

---

## 5. 对 §68 判定

```text
BENCHMARK_VALIDATION = PARTIALLY_SUPPORTED   （维持）
```
* H3 FAIL（硬条件，未变）
* 外部复现：**从「阴性」修正为「不可检」**（池太小）
* 两者都不满足 SUPPORTED 的条件

**但本轮取得了实质进展**：`U_EXT_CANONICAL` 是 V8 的**第一个多模型外部匹配宇宙**
（5 个可评价外部模型 × 200 target × 85.45% 行覆盖），
使「外部对比」从「没有资源」变为「有资源但模型数不足」。这是可累积的。

---

## 6. 复现命令

```powershell
python -W ignore -u scripts/v8/39_external_multimodel.py
```
产出：`external/external_multimodel_scoreboard.csv`、`external/external_multimodel_meta.json`
