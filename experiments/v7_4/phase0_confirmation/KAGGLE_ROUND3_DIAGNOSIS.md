# V7.4 Kaggle 第 3 轮 — 诊断记录与最终修复

**时间**：2026-10-15 02:12（notebook 时间戳）
**结果**：仍在冒烟阶段**主动停止**（§30 硬门控按设计工作）

**关键突破**：我在本机装了 **与 Kaggle 完全相同** 的
`openmm 8.6.1.dev-b399af4`，因此这一轮的结论是**实测**而非猜测。

---

## 1. 第 3 轮日志给出的决定性信息

```text
createSystem 签名: (self, topology, nonbondedMethod=NoCutoff, nonbondedCutoff=1.0 nm,
  constraints=None, rigidWater=None, removeCMMotion=True, hydrogenMass=None,
  residueTemplates={}, ignoreExternalBonds=False, switchDistance=None,
  flexibleConstraints=False, drudeMass=0.4 Da, **args)

gbn2 模板是否载入: []
error: [('implicitSolvent=app.GBn2', "ValueError: The argument 'implicitSolvent'
          was specified to createSystem() but was never used."),
        ('dielectric-only', "ImportError: cannot import name 'GBSAOBCForce' from
          'openmm.app.internal.customgbforces'"),
        ('manual-GB', "ImportError: cannot import name 'GBSAOBCForce' ...")]
```

## 2. 本机实测结论（OpenMM 8.6.1.dev-b399af4 = Kaggle 版本）

| 路径 | 实测结果 |
|---|---|
| `createSystem(implicitSolvent=app.GBn2)` | ❌ `ValueError: ... was never used` |
| `createSystem(soluteDielectric=…, solventDielectric=…)` | ❌ `ValueError: 'soluteDielectric' ... was never used` |
| `openmm.app.internal.customgbforces.GBSAOBCForce` | ❌ **不存在**（我第 2 版写错的类名） |
| `GBSAGBn2Force` | ❌ **段错误**（无 Python traceback） |
| **`GBSAOBC2Force`** | ✅ 可构造，`addParticles` / 全局参数可用 |
| **`openmm.GBSAOBCForce`**（核心类，OBC1） | ✅ 可用 |
| `app.GBn2` 常量 | 仍存在，但**已成孤儿**（无消费者） |

### 根因（源码级确认）

我核对了 OpenMM `master` 的 `forcefield.py`：

* `createSystem` 里只有 `switchDistance` / `flexibleConstraints` / `drudeMass`
  被放进 `ArgTracker`（⇒ 仅这三个是合法的 `**args` 消费者）
* **源码中已完全没有 `implicitSolvent` / `soluteDielectric` / `solventDielectric`**

⇒ 这是**真实的 API 变更**，不是我用错参数。
`app.GBn2` 常量保留只为向后兼容存在性检查，实际无人消费。

参考：[OpenMM `forcefield.py` (master)](https://raw.githubusercontent.com/openmm/openmm/master/wrappers/python/openmm/app/forcefield.py)

### 同时暴露出我第 2 版的第二个错误

`from openmm.app.internal.customgbforces import GBSAOBCForce` —— 该类**不存在**。
正确类名是 `GBSAOBC2Force` / `GBSAGBn2Force` / `GBSAOBC1Force` / `GBSAHCTForce`。

---

## 3. 最终修复（`02f_fix_gb_final.py`）

### 3.1 子进程隔离探测（**关键工程措施**）

本机实测该 dev build 的 GB 类会**段错误**（无 traceback、进程直接死）。
若在主进程里逐个试候选，会**杀死 Kaggle kernel 并丢掉全部已算结果**。

⇒ 改为在 `subprocess` 里逐个探测候选，段错误只意味着「该候选淘汰」：

```python
GB_CANDIDATES = ["GBSAOBC2", "GBn2", "GBSAOBC"]
for c in GB_CANDIDATES:
    ok, msg = _probe_candidate(pdb, c)   # 子进程
    if ok and GB_MODEL is None: GB_MODEL = c
```

### 3.2 多候选回落 + 绝不退化成真空

* 全部候选失败 ⇒ `FAIL_FORCEFIELD`（**不**用真空力场顶替，§6）
* 实际生效模型记入 `provenance.json`：`gb_model_applied`、`gb_protocol_deviation`

### 3.3 ⚠️ 必须声明的协议偏离

§6 指定 **GBn2**（`implicit/gbn2.xml`）。但 GBn2 在该 OpenMM 版本**不可用**；
若实际生效的是 **GBSA-OBC2** 或 **OBC1**：

* 三者同属 AMBER GB 隐式溶剂家族，但**参数化不同**
* ⇒ 计算条件与 §6 指定**不完全一致**，属**协议偏离**
* ⇒ **必须在 V7.4 最终报告「Energy Semantics / Relaxation Protocol」中如实标注**

`provenance.json` 的 `gb_protocol_deviation` 字段会自动记录该事实。

---

## 4. 三轮失败的共同教训（值得记录）

| 轮次 | 失败点 | 性质 |
|---|---|---|
| 1 | `Platform "CUDA"` 不存在 | 环境（已修：`openmm[cuda12]`） |
| 2 | `implicitSolvent=None` **静默**跳过 GBn2 | **协议级静默错误**（最危险） |
| 2 | `implicitSolvent=app.GBn2` 报错 | API 变更（暴露了上一条） |
| 3 | `customgbforces.GBSAOBCForce` 不存在 | 我写错类名 |

**最重要的收获**：第 1 轮的「成功」其实是**静默错误** ——
它跑出了**无隐式溶剂**的能量而不报错。
**§30 硬门控的价值正在于此**：把静默错误变成显式失败。

---

## 5. 对科学结论的影响

**仍然无影响。** 第 1–3 轮**均未产出任何特征**，
§40 的对比（R4−R0、R2−R1…）**全部未执行**，**尚无科学结论**。

> ⚠️ 不得叙述为「V7.4 跑过、松弛无效」。
> 正确叙述：**V7.4 前三轮分别在 OpenMM 平台、隐式溶剂 API、GB 类名上失败，
> 未产生特征，尚无科学结论。**

---

## 6. 本地已验证可用（无需再猜）

`scripts/v7_4/_probe_gbn2.py` / `_probe_obc2.py` / `_probe_step.py`
在**与 Kaggle 同版本**的本地 OpenMM 上运行，已确认：

```text
GBSAOBC2Force  -> 可构造、可 addParticles、可设全局参数
openmm.GBSAOBCForce (core, OBC1) -> 可用
implicitSolvent 路径 -> 确定不可用
```

## 7. 第 4 轮动作

1. 重新 **Import Notebook**（覆盖）
2. 保持 `RUN_SMOKE_ONLY = True`，Run All
3. 看 cell 06 的 **GB 候选探测表**：
   ```
   GBSAOBC2   PASS  OK GBSAOBC2 E=... forces=...
   GBn2       FAIL  rc=... (崩溃/段错误?)
   GBSAOBC    ...
   >>> 选定 GB 模型: GBSAOBC2
   ```
4. 若 `选定 GB 模型: None` → 把整张探测表发我
5. 若成功 → 看 cell 08b 的 `NUM_NOISE_FLOOR_KJMOL`
