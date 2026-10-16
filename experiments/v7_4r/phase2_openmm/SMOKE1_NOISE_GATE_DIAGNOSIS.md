# V7.4R Track A — Colab smoke 第 1 次：噪声门失败的根因诊断与修复

**结果**：`NOISE_GATE = NUMERICAL_NOISE_TOO_HIGH` → **STOP**（§19/§20 按设计触发）

> ✅ **这是门槛的胜利，不是科学结论。**
> 我们在花两小时跑 3253 个变体**之前**就发现协议有缺陷。
> `RELAXATION_NULL` 仍为 **FALSE** —— 松弛从未被有效检验。

---

## 1. 实测数字（用户 Colab 输出）

### 环境（全部通过）
```
platform = CUDA            ✅（openmm[cuda12] 装成功）
OpenMM   = 8.6.1.dev-b399af4
MUTATION_BUILD_PASS = True ✅（20/20 PASS，6884 原子）
自动解压输入包       ✅（发现并解压 colab_input_package.zip）
Drive mount          ✅
```

### GB_XML_PROBE —— **PASS** ✅

```
CustomGBForce count = 1
GB particles = 6884 = System particles   ✅ 覆盖全部粒子
force classes = [..., 'CustomGBForce']   ✅
Energy finite / Forces finite / Context / Minimization 全 True
GB_MODEL = GBn2_XML                       ✅ 无协议偏离
```

**但这一行露出了根因**：

```
[   8.57s] energy=769233143199.849 kJ/mol        <-- 7.69e11 !!
[ 245.69s] minimized E_post=-62874.832 dE=-769233206074.681
```

| 量 | 值 | 归一 | 判断 |
|---|---|---|---|
| `E_pre` | 7.69e11 kJ/mol | **1.12e8 kJ/mol/atom** | **荒谬** |
| `E_post` | −6.29e4 kJ/mol | −9.13 kJ/mol/atom | 合理 |
| 耗时 | **237 s** | | 在「拆弹」，不在找极小 |

### NOISE GATE —— **FAIL** ❌

```
WT->WT dE   : median 5995.5   MAD 3565.0   P95 10703.0   max 10988.6
              （5 次：+10989, -1299, +2563, +9560, +5995 —— 连符号都变）
WT->WT RMSD : median 1.5042   P95 1.5988   max 1.6101
WT->Mut dE  : median 5085.4   P95 22809.1  max 25405.7
WT->Mut RMSD: median 1.4687   P95 1.5387   max 1.5900

R_E = 2.105   (>= 1.0)
R_G = 1.089   (>= 1.0)
=> NUMERICAL_NOISE_TOO_HIGH -> STOP
```

---

## 2. 根因：两个独立缺陷

### 缺陷 1（主因）：起点严重冲突 ⇒ 最小化在「拆弹」

`E_pre = 1.12e8 kJ/mol/atom` **只可能**来自原子严重重叠。
来源是 `PDBFixer` 的建系顺序：

```
removeHeterogens → findMissingAtoms(23) → addMissingAtoms → addMissingHydrogens(7.0)
```

**新增重原子与氢是直接放上去的，完全没有预松弛。**
最小化花 **237 s** 释放 7.7e11 kJ/mol —— 这不是在找局部极小，是在**拆弹**。

⇒ 之后**所有**配对最小化都从**远离任何物理构象**的起点出发：
走的路径依赖数值偶然，**落进不同局部极小**。
这就是 identity（**完全相同的输入**）两次最小化差 1.2e4 kJ/mol、
**符号都变**的原因。

### 缺陷 2：RMSD ≈1.5 Å 与「搬动哪里」无关 ⇒ 不是生物学信号

5 次 identity 覆盖 `shell = 4, 9, 13, 19, 29`，RMSD **全部**落在 1.40–1.61 Å。
**一个 4 残基的壳不可能合理地动 1.5 Å。**

⇒ 这 1.5 Å 是**整体漂移**，不是位点局部响应。
mutation 的 1.4687 Å 只是**同一个漂移** —— 与 mutation 本身无关。
所以 `R_G ≈ 1.09` 完全在预期之内。

---

## 3. 修复：**只改建系步骤，不动任何 §6 协议数值**

保留原样：`ff14SB + implicit/gbn2.xml`、`pH 7.0`、**8 Å 壳**、
骨架 **10** / 壳外 **100** kcal/mol/Å²、`LocalEnergyMinimizer`、`tol 10`。

### A. 两阶段最小化（staged relaxation，标准结构生物学做法）

| 阶段 | 施加 §14 限制？ | 目的 |
|---|---|---|
| **Stage A** | **否** | 清除 `addMissingAtoms`/`addMissingHydrogens` 引入的原子重叠，把体系带到**可最小化**的状态 |
| **Stage B** | **是**（8 Å / 10 / 100） | **真正**的配对目标松弛；`E_pre`/`E_post` **全部取自 Stage B** |

`MIN_MAXITER_STAGE_A = 20000`（原 2000 不够）。

### B. 新增前置硬门：起点能量合理性（新 cell 08c）

```
E_postA_per_atom < MAX_PRE_ENERGY_PER_ATOM = 1000 kJ/mol/atom
```
不通过 ⇒ `STAGE_A_FAIL` ⇒ STOP。**在任何科学特征之前。**

### C. 真正检查收敛

Stage B 后连续两轮能量变化 `< CONVERGENCE_DTOL = 1e-3 kJ/mol` 才判收敛；
不收敛的行标记 `NOT_CONVERGED`（新版已把 `converged` 记入 diag）。

### D. 新增收敛性硬门（cell 11）

```
|dE_identity| median / |E_post| < IDENTITY_REL_TOL = 0.01
```
旧版之所以「通过」，是因为**根本没有这个对照**。
实测旧版该比值约 `5995 / 62875 ≈ 0.095` —— 比阈值高**近 10 倍**。

> ⚠️ `MAX_PRE_ENERGY_PER_ATOM` / `CONVERGENCE_DTOL` / `IDENTITY_REL_TOL`
> 都是**项目工程阈值，不是领域标准**，必须在科学评估前冻结。

---

## 4. 验证

```
cells=20 code=19 语法错误=0
全部 minimize 调用点均解包 5 个返回值 ✅
§6 协议数值未变 ✅（RADIUS_MOVABLE 8.0 / RESTRAINT 10,100 / PH 7.0 / MIN_TOL 10.0）
UPSTREAM_HASH_CHANGE = 0 ✅
```

---

## 5. 下一步

1. 重新上传最新 notebook（覆盖）
2. `RUN_SMOKE_ONLY = True` → Run all
3. **先看新 cell 08c 的 `STAGE_A_GATE_PASS`**：
   * `PASS` ⇒ Stage A 把体系带回了物理合理量级
   * `FAIL` ⇒ 建系/参数化仍有问题，把该 cell 整块发我
4. 再看 cell 11 的 `CONVERGENCE_GATE_PASS` / `R_E` / `R_G` / `NOISE_GATE`

**两道门（Stage A 起点能量、收敛性/噪声）都通过之前，一条正式 feature 都不生成。**

---

## 6. 状态

```
SCIENTIFIC_FEATURES_GENERATED = 0
SCIENTIFIC_EVALUATION        = NOT_RUN
RELAXATION_NULL              = FALSE   <- 松弛从未被有效检验
```
