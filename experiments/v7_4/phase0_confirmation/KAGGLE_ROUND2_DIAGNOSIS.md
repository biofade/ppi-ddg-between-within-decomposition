# V7.4 Kaggle 第 2 轮 — 诊断记录与修复

**时间**：2026-10-15 02:06（notebook 时间戳）
**结果**：仍在冒烟阶段**主动停止**（§30 硬门控按设计工作）

---

## 1. 本轮**进展**（相对第 1 轮）

| 项 | 第 1 轮 | **第 2 轮** |
|---|---|---|
| GPU 可见 | ✅ 2× T4 | ✅ 2× T4 |
| OpenMM 安装 | ✅ | ✅ `8.6.1.dev-b399af4` |
| **OpenMM platform** | ❌ `no registered Platform called "CUDA"` | ✅ **`>>> 实际使用 platform: CUDA`** |
| ff14SB / gbn2 文件加载 | ✅ | ✅ 两个都 `OK` |
| 输入包读取 | ✅ | ✅ 3253 行 / 4712 位点 / 213 链 / 162 targets |
| 标签检查 | — | ✅ **无标签列** |
| **参数化** | ✅ PASS（当时 `implicitSolvent=None`） | ❌ **FAIL**（见下） |

`openmm[cuda12]` 安装**成功**，平台列表从
`['Reference','CPU','OpenCL']` 变为
`['Reference','CPU','OpenCL','CPU','CUDA','OpenCL']`，选中 **CUDA**。
⇒ 第 1 轮的根因**已彻底解决**。

---

## 2. 本轮**新**失败（真实日志）

```text
WT 参数化: FAIL FORCEFIELD:ValueError:The argument 'implicitSolvent' was
specified to createSystem() but was never used.
```

### 这是我的修复引入的

- 第 1 轮：`createSystem(..., implicitSolvent=None)` → **能跑，但 GBn2 没生效**（静默错误）
- 第 2 轮：`createSystem(..., implicitSolvent=app.GBn2)` → **OpenMM 8.6.1.dev 不接受**（显式报错）

⇒ 即 **`implicitSolvent` 参数在 OpenMM 8.6 的开发版里已改变语义或已被移除**。
OpenMM 官方文档中 `createSystem` 的隐式溶剂参数在版本间确有变动
（见 [OpenMM API 文档](https://docs.openmm.org/development/api-python/generated/openmm.app.amberprmtopfile.AmberPrmtopFile.html)；
同类报错亦见于 [perses issue #886](https://github.com/choderalab/perses/issues/886)）。

> **重要**：第 2 轮的失败**优于**第 1 轮的"成功"。
> 第 1 轮会静默产出**无隐式溶剂**的能量 —— 那是错误协议。
> 本轮把它变成**显式失败**，正是硬门控的价值。

### 冒烟门控如实打印

```
MUTANT_BUILD_PASS                  = False
FORCEFIELD_PARAMETERIZATION_PASS   = False
WT_MUT_PAIR_PASS                   = False
NO_NAN_PASS                        = True     <- 占位
ENERGY_FINITE_PASS                 = True     <- 占位
MINIMIZATION_CONVERGENCE_PASS      = True     <- 占位
RESTRAINT_ENERGY_SEPARATED         = True     <- 占位
RELAXED_COORDINATES_VALID          = True     <- 占位
CACHE_RELOAD_PASS                  = True     <- 占位
--> RuntimeError: SMOKE FAIL —— 依 §30 停止
```

后续 6 项仍是**占位常量**（它们在最小化真正跑起来之前无法评估）。
**这些占位项必须在第 3 轮改成真实检查** —— 否则会给出虚假的安全感。
已列入本轮的修复清单。

---

## 3. 本轮修复（`02d_fix_solvent.py`）

不再猜测 API，改为**运行时自省 + 多策略回落**，并**绝不静默退化成真空**：

| 策略 | 做法 | 说明 |
|---|---|---|
| **1** | `createSystem(..., implicitSolvent=app.GBn2)` | 旧版 API |
| **2** | `createSystem(..., soluteDielectric=1.0, solventDielectric=78.5)` + 手工挂 GB 力 | 介电常数路线 |
| **3** | `createSystem(...)` 不含溶剂 + 从 `ForceField` 模板 `createForce()` 挂 GB 力 | 最不依赖 API |
| 全失败 | 返回 `FAIL_FORCEFIELD` | **不**退化成真空 |

配套改动：

1. **打印 `createSystem` 的真实签名**（`inspect.signature`）
   ⇒ 下一轮无论成败都能看到确切 API
2. **打印 ForceField 已载入的力模板**（判断 `gbn2.xml` 是否真被解析）
3. `prepare_system` 返回并记录 **`SOLVENT_MODE`**（实际用了哪条策略）
4. **`provenance.json` 记录 `solvent_mode`**（必须如实声明计算条件）
5. **SMOKE GATE 新增 `IMPLICIT_SOLVENT_APPLIED`**（§30 硬门控）
6. 18 code cell / **0 语法错误**

### 占位门控的处理（本轮已暴露的问题）

cell 08 的 6 项占位（`NO_NAN` / `ENERGY_FINITE` / `MINIMIZATION_CONVERGENCE` /
`RESTRAINT_ENERGY_SEPARATED` / `RELAXED_COORDINATES_VALID` / `CACHE_RELOAD`）
**不是**真检查。据 §30「若任一 fail → STOP」的精神，占位常量会误导。
⇒ **第 3 轮必须让 cell 08b 的 identity 对照结果反过来驱动这些门控**
（数值地板为 `None` 或过大 ⇒ 门控 False ⇒ 停止）。
本轮已把 `IMPLICIT_SOLVENT_APPLIED` 加入，其余待第 3 轮一并处理。

---

## 4. 对科学结论的影响

**仍然无影响。** 第 1、2 轮**均未产出任何特征**，
§40 的对比（R4−R0、R2−R1…）**全部未执行**，**尚无科学结论**。

> ⚠️ 不得叙述为「V7.4 跑过、松弛无效」。
> 正确叙述：**V7.4 第 1–2 轮分别在 platform 与隐式溶剂 API 上失败，
> 未产生特征，尚无科学结论。**

---

## 5. 第 3 轮动作

1. 重新 **Import Notebook**（覆盖旧版）
2. 保持 `RUN_SMOKE_ONLY = True`，Run All
3. **把这三行报我**：
   - `createSystem 签名: ...`
   - `溶剂策略: ...`（或 `FAIL` 时的完整 attempts 列表）
   - `>>> 实际使用 platform: ...`（应为 `CUDA`）
4. 若 `FF_OK = True` 且 `SOLVENT_MODE` 非空 → 继续看 cell 08b 的 `NUM_NOISE_FLOOR_KJMOL`
5. 若仍 FAIL → 把 `createSystem 签名` 与 attempts 发我，可**精确**定位 API
