# V7.4 — 第 4 轮 Kaggle smoke test 交接（两道门）

**时间**：2026-10-15
**状态**：notebook 已就绪；**未跑全量、未生成任何正式 feature**

> 依据用户第 4 轮指令：**不再改物理协议、不跑全量**，
> 只做 smoke test，并以 **Linux/T4 下的 `GB_XML_PROBE` 为最终裁决**。
> 本机 Windows 观察只能当诊断信息。

---

## 1. 本版修正的两个真实缺陷（自查发现）

| # | 缺陷 | 后果 | 修复 |
|---|---|---|---|
| 1 | `residues_within` **被定义两次**（cell 07 的 nm 版 + cell 08 的 Å 版） | 调用哪个取决于执行顺序，属隐患 | 删除 cell 08 的重复定义，现仅 **1 处** |
| 2 | cell 12 用**旧的 3 参签名**调 `minimize(sw, tw, fw.positions)` | 运行必崩；且原实现只在原子数相等时有意义 | 重写为**真正的 §25 几何增量**（配对松弛 + 位点局部几何） |

---

## 2. cell 06 新增的存在性检查（按用户要求）

除了原有的 `rc` / `CustomGBForce present` / energy finite / forces finite /
Context / Minimization，现在**明确输出**：

```text
CustomGBForce count:
GB particles:
System particles:
GB covers all particles:
GB force group:
Force groups:            <- 每个力的力组，后面分离 restraint 与 biological/GB 用
E_pre_min :
E_post_min:
DeltaE_min:
```

### 硬门控（全部满足才 PASS）

| 门 | 条件 |
|---|---|
| 进程 | `rc == 0` |
| GB 力存在 | `CustomGBForce present == True` |
| **GB 力数量** | **`CustomGBForce count >= 1`** |
| **GB 覆盖性** | **`GB particles == System particles`** |
| 能量 | finite |
| 力 | finite |
| Context | 建成 |
| 最小化 | 跑通 |
| **能量不显著上升** | `DeltaE_min <= 0`（若初始能量 >1e6，即明显真空/崩坏）否则 `<= 5% × |E|` |
| provenance | `GB_MODEL = "GBn2_XML"` |

> **为什么必须查「GB 粒子数 == system 粒子数」**：
> 用户指出的第三种情况 —— `CustomGBForce present=True` 但只覆盖部分粒子，
> 或力没有真正进入最终 system。只看 `present` 会漏掉这类错误。

> **为什么必须查 `DeltaE_min` 不上升**：
> 我本地实测过**真空力场**会「跑通」：`Energy = 1.1e14 kJ/mol`、
> `minimization OK`。**「能跑完」绝不等于「协议正确」。**

---

## 3. cell 08b：第二道门 = 数值噪声地板

### 输出的统计量（按用户指定）

```text
WT→WT ΔE    : n / median / MAD / P95 / max
WT→WT RMSD  : n / median / MAD / P95 / max
repeat      : dE median / P95 / max,  RMSD median / P95 / max
WT→Mut ΔE   : n / median / MAD / P95 / max      <- 20 个真实 mutation smoke
WT→Mut RMSD : n / median / MAD / P95 / max
```

### 真正要回答的问题（比值，**不预设绝对阈值**）

```text
R_E = Q95(|ΔE|_WT→WT)  /  median(|ΔE|_WT→Mut)
R_G = Q95(RMSD_WT→WT)  /  median(RMSD_WT→Mut)
```

### 噪声门判定

| 条件 | 判定 |
|---|---|
| `R_E >= 1.0` | **`NUMERICAL_NOISE_TOO_HIGH` → STOP** |
| `0.5 <= R_E < 1.0` | `MARGINAL`（同量级；需人工判断）→ STOP |
| `R_E < 0.5` | OK（identity 噪声明显低于 mutation 信号）→ 可继续 |

未通过则 `raise RuntimeError("NUMERICAL_NOISE_GATE_FAIL")`。
结果同时写入 `noise_gate.json`（含 `GB_MODEL` 与 `PLAT_NAME`）。
**全程不需要 ΔΔG 标签。**

---

## 4. 执行顺序（用户指定）

```text
1. Import 最新 notebook（覆盖）
2. Accelerator = GPU T4 x2
3. Internet ON
4. RUN_SMOKE_ONLY = True        <- 默认已是 True
5. Run All
6. **不要**切全量
```

### 只需先看 cell 06

**情况 A — PASS**：应看到

```text
GB_XML_PROBE
  OpenMM version: 8.6.1.dev-...
  ForceField files:
    amber14/protein.ff14SB.xml
    implicit/gbn2.xml
  Process rc: 0
  Force classes: [... CustomGBForce ...]
  GB-like forces: [...]
  CustomGBForce present: True
  CustomGBForce count: 1
  GB particles: N
  System particles: N
  GB covers all particles: True
  GB force group: 0
  Energy finite: True
  Force finite: True
  Context creation: True
  Minimization: True
  E_pre_min / E_post_min / DeltaE_min: ...
  PASS
```

⇒ **PRIMARY = `ff14SB + GBn2`，无 protocol deviation。**
**不要再考虑 OBC2，也不要 pin 旧版 OpenMM。** 然后才进入 cell 08b。

**情况 B — FAIL**（`rc != 0` 或 `CustomGBForce present = False`）：
即使 energy finite、minimization PASS，也判 **`GBN2_XML_PROTOCOL_FAIL`**，STOP。
本 cell 会主动抛错并打印下一步二选一：
**(a) 钉稳定版 OpenMM 保 GBn2**（**优先**）；(b) OBC2 protocol amendment。

**情况 C — GB 存在但异常**（`GB particles != System particles`、NaN）：
直接 FAIL。

---

## 5. 请回传的内容

```text
GB_XML_PROBE                      （cell 06 整块）
NUM_NOISE_FLOOR_KJMOL             （cell 08b）
E_pre_min / E_post_min
CustomGBForce count
GB particles / System particles
WT→WT identity-control 统计（median/MAD/P95/max）
20 mutation smoke 统计 + R_E / R_G
```

---

## 6. 两道门都过之后

$$ \boxed{\text{冻结协议} \rightarrow \text{开始全量 V7.4}} $$

**在那之前，一条正式 feature 都不生成。**（当前确实为 0 条正式 feature。）

---

## 7. 这三轮暴露的工程原则（记录备查）

$$ \boxed{\text{successful execution} \neq \text{scientifically valid execution}} $$

* 第 1 轮：CUDA platform 缺失 —— 容易发现
* **第 2 轮（最危险）**：`implicitSolvent=None` ⇒ 程序跑了、finite energy 也出来了，
  但科学协议**已经静默变成 vacuum**。若进入 3253 个变体，最后可能得到
  「relaxation 没有帮助」的结论，从而**错误地否定一条真正有效的路线**。
* 第 3 轮：`implicitSolvent` 参数被移除（API 变更，暴露了第 2 轮的静默错误）

⇒ V7.4 现在的质量高于初始版本：**模型存在性、力场存在性、数值稳定性、
identity noise、repeatability** 都已成为硬门控。
