# V7.4 — 官方 XML GBn2 路线：本地检查结果与最终决策树实现

**时间**：2026-10-15
**触发**：用户指出 OpenMM 8.6 的官方 GBn2 用法是
`ForceField(ff14SB, "implicit/gbn2.xml")` + `createSystem()`（**不传** `implicitSolvent=`），
并要求先测这条路，且**必须在子进程里测、不允许 fallback**。

---

## 1. 诚实的更正：我第 3 轮的探针并没有测到这条路

我的 `_probe_gbn2.py` case B 写的是：

```python
app.ForceField('amber14/protein.ff14SB.xml', 'implicit/gbn2.xml')
ff.createSystem(..., implicitSolvent=app.GBn2)      # <-- 仍然带了参数
```

⇒ 它测的是**「XML + 参数」**，不是用户指出的**「仅 XML」**。
用户是对的：这个区别是关键。已按用户要求补测。

## 2. 本地补测结果（**必须带重要限制**）

### 2.1 本地环境的已知缺陷（先声明）

本机 Windows venv 的 **`PDBFixer.findMissingAtoms()` 会确定性段错误**
（6/6 次，与 GB 无关，numpy 1.26.4 与 2.4.6 都一样）。
⇒ 本地**无法**把体系加氢到可参数化状态
（跳过该步只剩 6000 原子，`createSystem` 会因缺氢报 template 错误，
**与 GB 无关**）。

**因此本地这次补测不构成对 GBn2 XML 路线的判定。** 但其中一项观察
**不依赖 PDBFixer**，下面单独列出。

### 2.2 不依赖 PDBFixer 的关键观察：`gbn2.xml` 没有注册任何力

这是本次最有价值的一条，且完全绕过了 PDBFixer：

```text
files=['amber14/protein.ff14SB.xml']
   generators: ['HarmonicBondGenerator','HarmonicAngleGenerator',
                'PeriodicTorsionGenerator','NonbondedGenerator']
   NonbondedGenerator.createForce sig: (sys, data, nonbondedMethod, nonbondedCutoff, args)
   溶剂相关属性: []
   _templates 中含 GBSA/GB 的键: []
   _patches   中含 GB 的键: []
   _atomClasses 含 GB 的类型数: 0

files=['amber14/protein.ff14SB.xml','implicit/gbn2.xml']
   generators: ['HarmonicBondGenerator','HarmonicAngleGenerator',
                'PeriodicTorsionGenerator','NonbondedGenerator']   <-- 完全相同
   NonbondedGenerator.createForce sig: (同上)
   溶剂相关属性: []
   _templates 中含 GBSA/GB 的键: []
   _patches   中含 GB 的键: []
   _atomClasses 含 GB 的类型数: 0
```

**观察**：加不加 `implicit/gbn2.xml`，`ForceField` 的 **generator / template /
atomClass 完全一致**；`NonbondedGenerator.createForce` 的签名是
`(sys, data, nonbondedMethod, nonbondedCutoff, args)` —— **没有溶剂参数**。

这与官方文档描述一致：**GB 模型确实是通过 XML 装载的**。
`gbn2.xml` 提供的是**原子类型上的半径/scale 参数**，
而**读取这些参数并构造 GB 力的代码路径**必须存在才会生效。

在我这个 `8.6.1.dev-b399af4` 构建上，该路径**看起来不存在**
（这也是 `ArgTracker` 报 `implicitSolvent ... was never used` 的同一个原因）。

> ⚠️ **限制**：这是**本机 Windows 构建**的观察。
> 无法排除 Windows wheel 与 Linux wheel 的差异。
> ⇒ **最终判定必须来自 Kaggle（Linux）的子进程探针**，这正是 notebook 现在做的事。

### 2.3 用于对照的真空基线（说明门控必要性）

无溶剂（`TAG=none`）在**加氢成功**的那一次运行中：

```text
Energy: 114546315411697.3 kJ/mol      <-- 1.1e14，完全非物理
max|F| = 4.933e+16 kJ/mol/nm
Minimization: OK   dE = -1.145e14
```

⇒ **真空力场会「跑通」并给出 finite 能量**。
这再次证明：**「能跑完」绝不能当作「协议正确」**。
故硬门控必须**同时**要求 `CustomGBForce present`，
而不能只看 rc==0 与 energy finite。

---

## 3. 已实现的 notebook 改动（严格按用户指定的决策树）

`cell 06` 现在只做一件事：**测官方 XML GBn2 路线，无 fallback**。

### 3.1 探测的调用（与用户指定一致）

```python
ff = app.ForceField("amber14/protein.ff14SB.xml", "implicit/gbn2.xml")
system = ff.createSystem(topology,
                         nonbondedMethod=app.NoCutoff,
                         constraints=app.HBonds)
```

**不再出现 `implicitSolvent=`**（已验证 notebook 中零残留）。

### 3.2 子进程隔离（保留）

用户明确要求「不要因为官方路线就取消 crash isolation」。
探测在 `subprocess` 中运行；段错误只表现为 `rc != 0`，不会杀死 kernel。

### 3.3 输出用户要求的完整报告

```text
GB_XML_PROBE
  OpenMM version:
  ForceField files:
    amber14/protein.ff14SB.xml
    implicit/gbn2.xml
  Process rc:
  Force classes:
  GB-like forces:
  CustomGBForce present:
  Particle count:
  Energy:
  Force finite:
  Context creation:
  Minimization:
  PASS/FAIL:
```

### 3.4 硬门控（全部满足才 PASS）

| 门 | 条件 |
|---|---|
| 进程 | `rc == 0` |
| GB 力存在 | `CustomGBForce present == True` |
| XML 加载 | `implicit/gbn2.xml` 已载入 |
| 能量 | `energy_finite == True` |
| 力 | `forces_finite == True` |
| Context | 建成成功 |
| 最小化 | 能跑 |
| provenance | `GB_MODEL = "GBn2_XML"` |

### 3.5 失败即 STOP，**不自动 fallback**

```text
STOP: 官方 XML GBn2 路线未通过硬门控。
  依用户指定的决策树，**不自动 fallback 到 OBC2**。
  下一步二选一（需人工决定）：
    (a) 钉一个经 smoke test 确认可工作的稳定 OpenMM 版本以保住 GBn2；
    (b) 正式提出 protocol amendment 改用 OBC2 作为 sensitivity。
```

notebook 中已**移除**所有 OBC2 / 核心类回落代码（已验证零残留）。

---

## 4. 决策树现状

```text
OpenMM 8.6.1 dev
        │
        ▼
ForceField(ff14SB, implicit/gbn2.xml) -> createSystem()      [cell 06 子进程探针]
        │
        ├── PASS  ->  PRIMARY = GBn2_XML，无协议偏离，开始正式 V7.4
        │
        └── FAIL / SEGV
                 -> STOP（不 fallback）
                 -> 人工决定：(a) 钉稳定版保 GBn2  (b) OBC2 protocol amendment
```

**OBC2 的定位**（按用户要求）：`engineering fallback / sensitivity`，
可用于验证工程流水线（mutation builder / WT→WT / repeatability /
minimization / ZIP / resume / energy decomposition），
但**暂不作为 V7.4 全量科学主协议** —— 否则一旦结果变化，
`GBn2 → OBC2` 会成为混入的第二个变量。

---

## 5. 本轮交付

* `kaggle/v7_4/V7_4_OpenMM_Relaxation.ipynb` —— cell 06 已改为官方 XML GBn2 探针
  （18 code cell，**0 语法错误**）
* 备份：`..._v1_unfixed.ipynb`、`..._v2_solventbroken.ipynb`、`..._v3_gbbroken.ipynb`
* 本地诊断脚本：`scripts/v7_4/_probe_gbxml.py`、`_probe_step.py`、`_probe_obc2.py`

## 6. 第 4 轮请回传

只要 `GB_XML_PROBE` 整块输出（cell 06 会完整打印）。
- 若 `PASS/FAIL: PASS` → 继续，看 cell 08b 的 `NUM_NOISE_FLOOR_KJMOL`
- 若 `PASS/FAIL: FAIL` → 该块会给出 rc / forces / CustomGBForce / energy，据此决定 (a) 或 (b)

> 前三轮**均未产出特征**，§40 对比未执行，**尚无科学结论**。
> 不得叙述为「松弛无效」。
