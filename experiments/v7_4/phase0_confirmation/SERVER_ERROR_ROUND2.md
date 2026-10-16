# V7.4 — 连续两次 `Server Error`：证据重估与分步探针

**现象**（用户实测）：cell 00–05 正常（环境报告、协议、**25 行验证集**都打印了），
**cell 06 没有任何输出**，Kaggle 显示 `Error: Server Error ... try again in 30 seconds`。
连续两次一致。

---

## 1. 重估：能确定什么、不能确定什么

### 能确定的

**25 行验证集来自 cell 05。** cell 05 的成功说明：

* `VAR` / `VAL` 构建正常
* **`import pdbfixer` 正常**（cell 05 有 `from pdbfixer import PDBFixer`）

### 关键推理

cell 06（上一版）的**第一条语句**就是 `print("GB_XML_PROBE")`，它排在
`fetch_pdb` **之后**、`subprocess` **之前**。既然这一行**完全没出现**，
两种解释：

| 解释 | 含义 |
|---|---|
| **A** | kernel 死在 `fetch_pdb`（网络/文件系统）或更早 |
| **B** | kernel 根本没死，是 **Kaggle 前端/API 自己的 `Server Error`**（页面级故障），输出没渲染出来 |

**B 不能排除** —— "Please try again in 30 seconds" 是**平台级**措辞，
不是 notebook 异常。cell 00–05 的输出是**崩溃前已渲染**的，与 B 相容。

⇒ **我不能从现有证据断定 cell 06 是罪魁**。上一版把责任归给 cell 06 是**推断，不是结论**。

---

## 2. 因此本轮不再「修 cell 06」，而是**缩小不确定性**

cell 06 被替换为一个**纯诊断探针**：8 个**互相独立、逐级加码**的步骤，
每步独立 `try/except`、**只打印、失败不中断**、`flush=True`。

| step | 内容 | 若不 OK 说明什么 |
|---|---|---|
| 1 | 纯 stdlib 打印 | kernel 是否活着 |
| 2 | `import numpy` | 基础科学栈 |
| 3 | `import openmm` | OpenMM 可用性 + 版本 |
| 4 | `import pdbfixer` | PDBFixer 可导入 |
| 5 | 枚举 OpenMM platform | CUDA 是否可见 |
| 6 | **`fetch_pdb(VAL.iloc[0])`** | **网络 + 文件系统（上次的可疑点）** |
| 7 | `PDBFixer(...)` 只解析 | PDB 解析 |
| 8 | `ForceField(ff14SB, gbn2.xml)` 只加载 | XML 能否装载（**不建 system**） |

**本单元刻意不做重活**：不 `createSystem`、不建 `Context`、不最小化、
**不用 `subprocess`**（已验证 cell 06 内 `subprocess` / `createSystem`
只出现在注释与一行 print 文本里）。

⇒ 即使某步出问题，也只会打印一行 `FAIL`，**不会拖垮 kernel**。

### 输出形态

```text
====================================================================
STEPWISE PROBE START
====================================================================

[STEP 1] python stdlib
   -> OK  python=3.13.15 platform=Linux-...
...
[STEP 6] fetch_pdb (network + filesystem)
   target=1A22
   -> OK  /kaggle/working/v7_4/pdb/1A22.pdb  (2.3s, 218937 bytes)
...
====================================================================
STEPWISE PROBE SUMMARY: {1: 'OK', 2: 'OK', ..., 8: 'OK'}
====================================================================
```

**⇒ 这一次无论怎么失败，SUMMARY 那一行就会告诉我们死在哪一步。**

---

## 3. 请这样做（先排除平台级故障）

1. **先确认 notebook 到底有没有真的在跑**：在 Kaggle 右上角看
   kernel 状态（`Running` / `Idle` / `Restarting`）；如果显示
   `Restarting` 或多次自动重启，说明 **kernel 崩了**，偏向解释 A。
2. 如果只是页面报 `Server Error` 而 kernel 状态仍为 `Idle`/`Running`
   —— 那可能是**平台前端故障**，稍等重试即可，偏向解释 B。
3. 重新 **Import 最新 notebook**（覆盖），保持
   `RUN_SMOKE_ONLY = True`，Run All。
4. 把 **STEPWISE PROBE 那一整块**（含 SUMMARY 行）贴回。
5. 如果连 `STEP 1` 都没出现 ⇒ **与 notebook 内容无关**，
   是平台/kernel 级问题，我们就换策略（例如新建 notebook、改用
   `File > Import` 而非覆盖，或先只跑 cell 06 单格）。

---

## 4. 状态

* **正式 feature 数：0**（仍为 0）
* 未跑全量；未修改物理协议
* 前几版的 GBn2 硬门控设计**保留**（一旦探针定位清楚，就把
  cell 06 换回完整 GBn2 探针）

---

## 5. 历史记录（供对照）

| 版本 | cell 06 形态 | 结果 |
|---|---|---|
| v1 | `implicitSolvent=None`（静默跳过 GBn2） | 跑通但**协议错误**（真空） |
| v2 | `implicitSolvent=app.GBn2` | `ValueError: ... was never used` |
| v3 | 多候选回落（含不存在的 `GBSAOBCForce`） | ImportError |
| v4 | 官方 XML 路线 | `ValueError`（仍带参数） |
| v5 | 官方 XML + 硬门控 | **Server Error（无输出）** |
| v6 | + 子进程流式进度 + 900s 超时 | **Server Error（无输出）** |
| **v7** | **8 步纯诊断探针（无重活/无 subprocess）** | 待测 |
