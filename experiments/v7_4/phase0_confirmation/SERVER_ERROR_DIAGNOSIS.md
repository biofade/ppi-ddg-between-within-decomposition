# V7.4 — `Server Error` 诊断与流式修复（第 5 轮）

**现象**（用户实测，两次一致）：
* cell 00–05 正常：环境报告、协议、**25 行验证集**都打印了
* **cell 06 之后没有 GB_XML_PROBE 输出**
* Kaggle 显示 `Error: Server Error ... Please try again in 30 seconds.`

**结论**：卡在/崩在 **cell 06**，而且我那一版 cell 06 是**静默黑盒** ——
子进程把所有输出缓存在 stdout，只有全部跑完才一次性返回。
于是「出问题」和「卡死」都无法区分，也拿不到任何诊断。

---

## 1. 我这一版 cell 06 的三个真实缺陷

| # | 缺陷 | 后果 |
|---|---|---|
| 1 | 子进程输出**无流式转发** | 崩溃/卡死时父 cell 什么都拿不到 → 完全无法诊断 |
| 2 | **无超时** | 子进程卡住 ⇒ cell 无限等待 ⇒ Kaggle kernel 被杀 ⇒ **Server Error** |
| 3 | 子进程**未加 `-u`** | stdout 块缓冲，进一步降低可见性 |

> 我本机已实测 `PDBFixer.addMissingAtoms()` 在 Windows dev build 上**确定性段错误**。
> Kaggle（Linux）未必相同，但**同一段代码值得设防**。

---

## 2. 修复（`02j_cell06_streaming.py` + 完整性守卫）

1. **子进程 `-u` 无缓冲 + 每步 `print(flush=True)`**，
   同时把每一步 append 到 `gb_probe_progress.log`。
2. **父进程逐行读取并实时转发**（`readline()` 循环），
   因此**即使子进程段错误，失败点也已经打印出来**。
3. **硬超时 `PROBE_TIMEOUT = 900 s`** —— 超时 `proc.kill()`，
   不再拖垮 kernel。
4. 父进程**即使拿不到最终 JSON**，也会：
   * 打印「未收到子进程 JSON 结果（可能崩溃/超时）—— 以上进度即失败点」
   * 兜底把进度文件里未转发的行补打（`[prog]` 前缀）
5. **补加完整性守卫**：先 `findMissingAtoms()` 统计缺原子/残基数，
   **仅当真的缺**才调 `addMissingAtoms()`（本机段错误嫌疑点）；
   结构完整则跳过并记录。
6. 全部硬门控**保留**：`CustomGBForce count>=1`、
   `GB particles == system particles`、energy/forces finite、Context、
   Minimization、`DeltaE_min` 不显著上升。**无 fallback。**

### 新的 cell 06 输出长这样（关键：有实时进度）

```text
GB_XML_PROBE
  ForceField files:
    amber14/protein.ff14SB.xml
    implicit/gbn2.xml
  platform: CUDA
  child timeout: 900s

--- 子进程实时进度 ---
  [   0.00s] start
  [   0.35s] numpy ok
  [   0.61s] openmm ok: 8.6.1.dev-...
  [   0.72s] pdbfixer imported
  [   0.81s] PDBFixer constructed
  [   0.85s] removeHeterogens
  [   0.91s] findMissingResidues
  [   1.20s] findMissingAtoms -> missingAtoms=0 missingResidues=0
  [   1.21s] structure already complete -> skip addMissingAtoms
  [   1.55s] addMissingHydrogens
  [   1.60s] topology atoms=<N>
  [   2.81s] ForceField loaded
  [   2.82s] generators=[...]
  ...
  Process rc: 0
  CustomGBForce count: ...
  GB particles: ... / System particles: ...
  ...
  PASS/FAIL: ...
```

**⇒ 这一次无论成功或失败，我们都会知道确切是哪一步。**

---

## 3. 请这样跑（与上一轮相同，不要开全量）

```text
1. Import 最新 notebook（覆盖）
2. Accelerator = GPU T4 x2
3. Internet ON
4. RUN_SMOKE_ONLY = True      （默认已是）
5. Run All
```

**只需把 cell 06 从 `GB_XML_PROBE` 开始的整块输出贴给我**（含实时进度行）。

### 三种可能的结果与含义

| 观察 | 含义 | 下一步 |
|---|---|---|
| 进度正常走到 `done`，`PASS` | GBn2 官方 XML 路线可用 | **无协议偏离**，进 cell 08b 看噪声门 |
| 进度停在某一步后 `rc != 0` | 该步崩溃 | 由停在哪一步定位 |
| 进度停在某一步后超时（900s） | 该步挂死 | 由停在哪一步定位 |
| 走到 `CustomGBForce count: 0` | GB 力**未注册**（我本机的观察在 Linux 上重现） | 判 `GBN2_XML_PROTOCOL_FAIL` → 走 (a) 钉稳定版 保 GBn2 |

> 注意：即使 `Energy finite: True` / `Minimization: True`，
> 只要 `CustomGBForce count == 0` 或 `GB particles != System particles`，
> 一律 **FAIL**。这正是「程序跑通 ≠ 协议正确」。

---

## 4. 状态

* **正式 feature 数：0**（仍为 0，未生成任何 feature）
* 未跑全量，未修改物理协议
* 上游哈希未变（V7.3 冻结）
