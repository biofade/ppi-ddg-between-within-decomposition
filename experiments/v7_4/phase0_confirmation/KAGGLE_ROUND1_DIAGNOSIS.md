# V7.4 Kaggle 第 1 轮 — 诊断记录与修复

**接收**：`F:\dsh\downloads\v7_4\V7_4_RELAXED_PHYSICS_FEATURES.zip`
**SHA256**：`02f3afe055e809d17fe8078adc2b34264e3fe44fa23392551364c3fe99dbde75`
**大小**：**2,162 B**（异常小）
**内容**：8 个报告文件，**`features/*.npz` = 0 个**

⇒ **这一轮没有产出任何松弛特征**。下面是从返回物**诊断出的确切原因**。

---

## 1. 诊断：什么成功了，什么失败了

### ✅ 成功的部分（有证据）

| 环节 | 证据 | 结论 |
|---|---|---|
| 环境 | `environment_report.json`：2× Tesla T4（14.56 GB 各）、31.35 GB RAM、CUDA (torch) 可用、internet 可用 | **PASS** |
| OpenMM 安装 | 版本 `8.6.1.dev-b399af4`（conda-forge/pypi） | **PASS** |
| PDBFixer 安装 | 已导入 | **PASS** |
| 输入发现 | `smoke_test.csv` 读到了 1A22 | **PASS** |
| **突变体构建** | `smoke_test.csv`：`wt_atoms=6884, mut_atoms=6884`，`shell_wt=5/10/14` | **PASS** |
| **ff14SB + GBn2 参数化** | 同上，6884 原子系统建成 | **PASS** |
| WT/MUT 配对 | 三条记录均 `PASS` | **PASS** |

### ❌ 失败的部分

**`benchmark.csv` 与 `identity_controls.csv` 全部记录**：

```
OpenMMException: There is no registered Platform called "CUDA"
```

⇒ **能量最小化从未执行成功**。这就是 0 个特征文件的直接原因。

**根因**：conda-forge 分发的 `openmm` 包**不含 CUDA platform**（需 CUDA toolkit
与 openmm 一同编译安装）。Kaggle 上 torch 看得到 GPU，**但不代表 OpenMM 也能用**。

**次要**：`RUN_SMOKE_ONLY` 实际为 `True`（cell 10 打印了跳过信息），
所以即使平台可用，本轮也不会产出特征。

---

## 2. 本轮同时暴露的 3 个**代码级** bug（已修）

| # | Bug | 严重性 | 证据 | 修复 |
|---|---|---|---|---|
| 1 | `mm.CustomExternalPotential()` — **该类不存在** | 高 | 被 `if False` 短路掩盖，未抛错 | 改用 `mm.CustomExternalForce` |
| 2 | **限制力从未加入 system** | **高** | `build_restrained_system(...)` 的返回值被丢弃 | 新增 `add_restraints()` 并真正 `system.addForce(...)` |
| 3 | `implicitSolvent=None` | **高** | `createSystem(..., implicitSolvent=None)` ⇒ **GBn2 完全没生效** | 改为 `implicitSolvent=app.GBn2` |
| 4 | `Platform.getPlatformByName("CUDA")` 硬编码 | 高 | 本轮全部异常 | 新增 `get_platform()`：CUDA → OpenCL → CPU → Reference 回落 |

> Bug 1–3 是**协议级**错误：即使平台可用，上一版 notebook 也会跑出
> **无限制、无隐式溶剂**的最小化 —— 那不是 §6/§14 要求的协议。
> 换言之：**幸好在冒烟阶段就失败了。**

### 单位换算（Bug 2 的配套修正）

`CustomExternalForce` 的表达式 `k*((x-x0)^2+(y-y0)^2+(z-z0)^2)` 用 OpenMM
内部单位 **kJ/mol/nm²**，而 §14 的常数以 **kcal/mol/Å²** 给定：

```
1 kcal/mol/Å² = 4.184 / 0.1² = 418.4 kJ/mol/nm²
K_SHELL = 10  kcal/mol/Å² -> 4184  kJ/mol/nm²
K_OUT   = 100 kcal/mol/Å² -> 41840 kJ/mol/nm²
```

**协议数值本身未改**（仍为 10 / 100 kcal/mol/Å²），只把单位换算补上。

---

## 3. 已实施的修复（`02b_fix_notebook.py` + `02c_fix_relaxation.py`）

| 项 | 状态 |
|---|---|
| `get_platform()` 平台回落 + 打印实际 platform | ✅ |
| 尝试 `openmm[cuda12]`，失败自动回落 | ✅ |
| `CustomExternalForce`（正确类名） | ✅ |
| `add_restraints()` 真正挂到 system | ✅ |
| `implicitSolvent=app.GBn2` | ✅ |
| restraint **单独力组 31** ⇒ 能量可分离（§21） | ✅ |
| `benchmark.csv` / `provenance.json` 记录 **实际 platform** | ✅ |
| identity 单元格 f-string 修复 | ✅ |
| 18 个 code cell **语法错误 0** | ✅ |

**保留旧版**：`V7_4_OpenMM_Relaxation_v1_unfixed.ipynb`（便于对照）

---

## 4. 一轮失败对科学结论的影响

**无影响。** 本轮**没有产出任何特征**，因此 §40 的任何对比（R4−R0、R2−R1 等）
**均未执行**。不得把本轮结果解读为「松弛无效」—— 它只是**工程失败**。

> ⚠️ 必须避免的错误叙述：「V7.4 跑过了，松弛没带来增益」。
> 正确叙述：**V7.4 第 1 轮在最小化阶段因 OpenMM 缺少 CUDA platform 而失败，
> 未产生任何特征，尚无科学结论。**

---

## 5. 下一轮的动作

1. 重新上传 `V7_4_OpenMM_Relaxation.ipynb`（已修复版）
2. 保持 `RUN_SMOKE_ONLY = True`，**Run All**
3. 检查 cell 02 打印的 `>>> 实际使用 platform:` —— 期望 `OpenCL`（T4 支持）或 `CUDA`
4. 检查 cell 08b 的 `NUM_NOISE_FLOOR_KJMOL=` —— **必须很小**
5. 若 `identity floor` 很大 → 停止并回报
6. 一切正常 → 改 `RUN_SMOKE_ONLY = False` 跑全量

**若 platform 仍是 `Reference`（最慢）**：把 cell 09 的估算小时数报给我，
我会改成确定性分片以便跨会话续跑。
