# V8 外部池扩展（PHASE 13 定稿）— **变体可增，方法族不足**

> **本轮结论：外部模型池在「变体」维度上确实扩大了（5 → 12 个实体），
> 但在「独立方法族」维度上仍然不足（7 个族 < 8 门槛）
> ⇒ `EXTERNAL_MODEL_REPRODUCTION` 仍为 `NOT_TESTABLE`，
> 但**限制的性质变了**：从「资源不存在」变为「方法数不足」。**

---

## 1. 新发现的资源

`github.com/StatPhysBio/hermes_experiments/skempi/results/` 目录列举显示 **15 个 HERMES 变体**
＋ **3 个 proteinMPNN 基线**，每目录内均有 `…-use_mt_structure=0.csv`（逐突变预测）。

**本轮下载并评价成功 7 个**（1 个 404）：

| 变体 | 家族 | n | targets≥5 | Overall r | Macro ρ | C_W | Eff@10 |
|---|---|---:|---:|---:|---:|---:|---:|
| `HERMES_py000` | HERMES | 3826 | 103 | +0.1292 | +0.0806 | +0.0967 | +0.8003 |
| `HERMES_py050` | HERMES | 3826 | 103 | +0.1021 | +0.0634 | +0.0900 | +0.8014 |
| `HERMES_bp000` | HERMES | 3842 | 106 | +0.1002 | +0.0721 | +0.0887 | +0.8016 |
| `HERMES_bp050` | HERMES | 3842 | 106 | +0.1048 | +0.0630 | +0.0864 | +0.8033 |
| `proteinMPNN_002` | proteinMPNN | 3842 | 106 | **+0.1843** | +0.0841 | +0.0842 | **+0.8205** |
| `proteinMPNN_020` | proteinMPNN | 3842 | 106 | +0.1591 | +0.0802 | +0.0714 | +0.8229 |
| `proteinMPNN_030` | proteinMPNN | 3842 | 106 | +0.1580 | **+0.0948** | +0.0699 | +0.8192 |

（口径：`proxy = −(log_p_mt − log_p_wt)` 等，与 V8 的 `ΔΔG>0=结合变弱` 对齐。）

---

## 2. ★ 必须同时报告两套读数（防「伪多样性」）

```
按**变体**计 (n=7):  rho(Overall r, Macro rho) = +0.6786   rho(Macro rho, Eff@10) = +0.3214
按**方法族**计 (n=2): HERMES      Overall r=+0.1035  Macro rho=+0.0678  Eff@10=+0.8015
                     proteinMPNN Overall r=+0.1591  Macro rho=+0.0841  Eff@10=+0.8205
```

**⚠️ 这 7 个变体只来自 2 个方法族。**
4 个 HERMES 变体**共享同一架构与 SKEMPI 训练数据**，仅 seed / backbone（py vs bp）
/ 温度不同 —— 它们**不是独立方法**。

**⇒ 若只报「按变体」的 ρ=+0.68，会得到一个由 4 个近重复点驱动的数，
这在方法学上与 PHASE 13 早先发现的「模型池组成伪影」是同一类错误。**
**依 §69 与本文档自己的教训，两套读数并列报告，且以「按方法族」为准。**

---

## 3. 合并后的外部池总账

| 来源 | 实体数 | **独立方法族数** |
|---|---:|---:|
| `U_HERMES_CANONICAL`（本轮） | 7 | **2**（HERMES, proteinMPNN） |
| `U_EXT_CANONICAL`（CATH-ddG 表） | 5 | **5**（Flex-ddG, DDGPred, FoldX, ESM-IF, proteinMPNN） |
| **合计（去重后）** | **12** | **6**（proteinMPNN 重叠） |
| **可检门槛** | ≥8 | **≥8** |
| **判定** | ✅ 达标 | ❌ **未达标** |

> ⚠️ 两组来自**不同宇宙**（不同的行集合与 target 定义），**不得直接合并计算秩相关**。
> 上表只是「资源总量」账，不是「可合并池」。

**⇒ 独立方法族 = 6 < 8 ⇒ 仍 `NOT_TESTABLE`。**

---

## 4. 一个仍然有价值的独立观察

`proteinMPNN` 在**两个不同宇宙**上都被评价过：
* 本轮（`U_HERMES_CANONICAL` 口径，HERMES 团队复算）：Overall r **+0.158 ~ +0.184**
* CATH-ddG 表（`U_EXT_CANONICAL` 口径）：Overall r **+0.1674**
* CATH-ddG held-out（`U_EXT_HELDOUT_CATH`）：Overall r **+0.1837**

**三个独立来源对同一方法给出一致的量级（+0.16 ~ +0.18）**，
且都远低于物理/统计法（FoldX/Flex-ddG/DDGPred 的 +0.26 ~ +0.61）。

**⇒ 这是跨来源一致的、关于 proteinMPNN 作为 ΔΔG 评分器的描述性事实。**
**但依 §71，不得表述为「学到的表示普遍不行」或「physics beats learning」** ——
这里只有 3 个来源、6 个方法族，且 proteinMPNN 本非为 ΔΔG 设计。

---

## 5. 对 §68 与最终判定的影响

```text
EXTERNAL_MODEL_REPRODUCTION = NOT_TESTABLE
   原因（更新）：独立外部方法族 = 6 < 8 门槛
   （已不再是「无资源」；两个宇宙 + 12 个实体已就位，随时可算）

BENCHMARK_VALIDATION    = PARTIALLY_SUPPORTED   （维持；H3 FAIL 为硬条件）
BENCHMARK_UTILITY_LEVEL = MORE_DISCRIMINATIVE   （维持）
```

**⇒ 依 §69，不因本轮扩大了资源而提升判定。**

**边界已明确**：要达到 `SUPPORTED`，需要 **≥2 个额外的独立外部方法**的逐突变预测。
可行来源（均已核实存在但需自跑推理）：
`PPIformer`（CPU 硬阻断）、`RDE-Network`（Apache-2.0，需自跑）、
`Prompt-DDG`（无评分 CLI）。

---

## 6. 复现命令

```powershell
python -W ignore -u scripts/v8/41_expand_external_pool.py
```
产出：`external/external_expanded_pool.csv`、`external/external_expanded_pool_meta.json`
