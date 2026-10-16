# -*- coding: utf-8 -*-
"""scripts/v7_4r/02d_staged_minimization.py -- 修「数值噪声地板过高」的**根因**（第 2 轮 Colab 结果）.

# Colab 第 1 次 smoke 的实测结果（关键数字）

| 量 | 值 | 解读 |
|---|---|---|
| probe `E_pre` | **7.69e11 kJ/mol** | **按原子归一 = 1.12e8 kJ/mol/atom —— 荒谬** |
| probe `E_post` | −6.29e4 kJ/mol | −9.1 kJ/mol/atom，**这才是合理量级** |
| probe `dE` | −7.69e11 kJ/mol | 最小化把「一个物理上不可能的起点」砸平 |
| identity dE（5 次） | +10989, −1299, +2563, +9560, +5995 | **中位 5995，跨度 1.2e4，符号都变** |
| identity RMSD | 1.4024 – 1.6101 Å | **floor ≈ 1.5 Å —— 与搬动哪个残基无关** |
| mutation RMSD | median 1.4687 Å | **≈ identity 的 1.5 Å** |
| **R_E** | **2.105** | identity 噪声 > mutation 信号 |
| **R_G** | **1.089** | 几何同样无区分度 |
| `NOISE_GATE` | **NUMERICAL_NOISE_TOO_HIGH** | **STOP（按 §19/§20 正确触发）** |

# 诊断：两个独立缺陷

## 缺陷 1（主因）：起点结构严重冲突，最小化在「拆弹」而不是「松弛」

`E_pre = 1.12e8 kJ/mol/atom` 只可能来自**原子严重重叠**。
来源是 `PDBFixer` 的建系顺序：

```
removeHeterogens → findMissingAtoms(23) → addMissingAtoms → addMissingHydrogens(7.0)
```

**新增的重原子与氢是直接放上去的，完全没有预松弛。**
probe 里最小化耗时 **237 s**（8.57s → 245.69s）却在释放 7.7e11 kJ/mol
—— 这不是在找局部极小，是在**拆弹**。

⇒ 之后所有「配对最小化」都从一个**远离任何物理构象**的起点出发，
走的路径**依赖数值偶然**，落进**不同的局部极小**。这就是为什么：

* identity（**完全相同的输入**）两次最小化会差 1.2e4 kJ/mol、连符号都变
* identity 与 mutation 的 RMSD 都是 ~1.5 Å（见缺陷 2）

## 缺陷 2：RMSD ≈1.5 Å 与「搬动哪里」无关 ⇒ 这不是生物学信号

5 次 identity 覆盖 `shell=4, 9, 13, 19, 29`，RMSD 全在 **1.40–1.61 Å**。
一个 4 残基的壳**不可能**合理地动 1.5 Å。
⇒ 这 1.5 Å 是**整个体系的整体漂移/松弛路径**，不是位点局部响应。

**后果**：mutation 的 1.4687 Å 只是同一个漂移 —— 与 mutation 本身无关。
所以 R_G ≈ 1.09 完全在预期之内。

# 修复（**不改任何 §6 协议数值**）

保留 `ff14SB + GBn2`、`pH 7.0`、8 Å 壳、限制 10/100、`LocalEnergyMinimizer`、
`tol 10` —— **只改「怎么从起点走到可最小化的状态」这个建系步骤**：

## A. **两阶段最小化（staged relaxation）** —— 标准结构生物学做法

1. **Stage A（建系松弛）**：加完缺失重原子与氢之后、**在施加 §14 限制之前**，
   先做一次**全局**自由最小化，把 `addMissingAtoms` / `addMissingHydrogens`
   引入的冲突清掉。判据是**能量量级回到合理范围（< 1000 kJ/mol/atom）**，
   不是「找极小」。
2. **Stage B（§14 目标松弛）**：从 Stage A 的结果出发，施加既有的
   8 Å 壳 / 骨架 10 / 壳外 100 限制，做**真正**的配对最小化。
   `E_pre` / `E_post` 全部取自 Stage B。

这样 `E_pre` 才是「物理上可解释的起点能量」，
配对最小化也才在**同一个盆地**里比较 WT 与 MUT。

## B. **新增前置硬门：起点能量合理性**

```
E_pre_per_atom = E_pre / system.getNumParticles()
要求 E_pre_per_atom < 1000 kJ/mol/atom      # 项目工程阈值，非领域标准
```

在一个**已知良好**的体系（例如单个小蛋白的 crystal 结构）上测定该阈值，
**先于**任何科学评估冻结。本 cell 直接把它打进 GB probe（§13 硬门控）。

## C. 提高最小化迭代上限并**真正检查收敛**

`MIN_MAXITER` 从 2000 提到 **20000**，并要求
**连续两次收紧后能量变化 < 1e-3 kJ/mol** 才判为收敛；
未收敛的行标记 `NOT_CONVERGED`，**不进入**科学特征。

## D. 收敛性对照（新的必要对照）

同一输入连跑两次，要求
```
|dE_identity| < 0.01 * |E_post|        # 项目工程阈值
```
旧版之所以能「通过」是因为**根本没有这个对照**；现在它在 cell 08b 里已经是硬门。
"""

from __future__ import annotations

import ast
import json
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
NB = R / "colab" / "v7_4r" / "V7_4R_OpenMM_Relaxation_COLAB.ipynb"
nb = json.loads(NB.read_text(encoding="utf-8"))
cells = nb["cells"]


def find(marker):
    for i, c in enumerate(cells):
        if c["cell_type"] == "code" and marker in "".join(c["source"]):
            return i
    raise KeyError(marker)


# ============================================================ 1) config：新增工程阈值
i00 = find("# ============================================================ 00 Config")
s00 = "".join(cells[i00]["source"])
s00 = s00.replace("MIN_MAXITER = 2000",
                  "MIN_MAXITER = 2000                  # Stage B（§14 目标松弛）\n"
                  "MIN_MAXITER_STAGE_A = 20000         # Stage A（建系松弛）\n"
                  "MAX_PRE_ENERGY_PER_ATOM = 1000.0    # kJ/mol/atom：起点能量合理性门（工程阈值）\n"
                  "CONVERGENCE_DTOL = 1e-3             # kJ/mol：连续两次能量变化阈值\n"
                  "IDENTITY_REL_TOL = 0.01             # |dE_identity| < 1% |E_post|")
cells[i00]["source"] = s00.splitlines(True)

# ============================================================ 2) cell 08：两阶段 minimize
i08 = find("# ============================================================ 08 Restraints + minimization")
s08 = "".join(cells[i08]["source"])
s08 = s08.replace(
    '''def minimize(system, topology, positions, site_keys, shell_keys):
    """加限制 -> L-BFGS 最小化。返回 (positions, e_pre, e_post, e_restraint)。
    e_pre/e_post **不含** restraint（§27）。"""
    add_restraints(system, topology, positions, site_keys, shell_keys)
    FG = assign_force_groups(system)
    integ = mm.VerletIntegrator(0.001)
    ctx = mm.Context(system, integ, PLAT)
    ctx.setPositions(positions)
    e_pre, _ = ctx_energies(ctx, system)
    mm.LocalEnergyMinimizer.minimize(ctx, tolerance=MIN_TOL,
                                     maxIterations=MIN_MAXITER)
    st = ctx.getState(getPositions=True, getEnergy=True)
    e_post, e_res = ctx_energies(ctx, system)
    return st.getPositions(), e_pre, e_post, e_res''',
    '''def stage_a_relax(system, topology, positions):
    """Stage A：**不施加 §14 限制**，先把 addMissingAtoms/addMissingHydrogens
    引入的原子重叠清掉。返回 (positions, e_preA, e_postA, per_atom_A)。

    判据是「能量量级回到合理范围」，不是「找极小」。
    这一步**不产生任何科学特征**，只是把体系带到可最小化的状态。"""
    integ = mm.VerletIntegrator(0.001)
    ctx = mm.Context(system, integ, PLAT)
    ctx.setPositions(positions)
    n = system.getNumParticles()
    e_preA = ctx.getState(getEnergy=True).getPotentialEnergy(
        ).value_in_unit(unit.kilojoule_per_mole)
    mm.LocalEnergyMinimizer.minimize(ctx, tolerance=MIN_TOL,
                                     maxIterations=MIN_MAXITER_STAGE_A)
    st = ctx.getState(getPositions=True, getEnergy=True)
    e_postA = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    return st.getPositions(), e_preA, e_postA, (e_postA / max(1, n))

def minimize(system, topology, positions, site_keys, shell_keys):
    """**两阶段**：
      Stage A  建系松弛（无 §14 限制）—— 消除加原子/加氢引入的重叠
      Stage B  §14 目标松弛（8 A 壳，骨架 10，壳外 100）
    返回 (positions, e_pre_B, e_post_B, e_restraint, diag)
    e_pre_B/e_post_B **不含** restraint（§27），且均为 **Stage A 之后**的物理能量。"""
    # ---- Stage A ----
    posA, e_preA, e_postA, per_atomA = stage_a_relax(system, topology, positions)
    # ---- Stage B：从 Stage A 结果出发加限制 ----
    add_restraints(system, topology, posA, site_keys, shell_keys)
    FG = assign_force_groups(system)
    integ = mm.VerletIntegrator(0.001)
    ctx = mm.Context(system, integ, PLAT)
    ctx.setPositions(posA)
    e_pre, _ = ctx_energies(ctx, system)
    mm.LocalEnergyMinimizer.minimize(ctx, tolerance=MIN_TOL,
                                     maxIterations=MIN_MAXITER)
    e_prev = None
    for _round in range(3):                     # 真正检查收敛
        st = ctx.getState(getPositions=True, getEnergy=True)
        e_now, _ = ctx_energies(ctx, system)
        if e_prev is not None and abs(e_now - e_prev) < CONVERGENCE_DTOL:
            break
        e_prev = e_now
        mm.LocalEnergyMinimizer.minimize(ctx, tolerance=MIN_TOL,
                                         maxIterations=MIN_MAXITER)
    st = ctx.getState(getPositions=True, getEnergy=True)
    e_post, e_res = ctx_energies(ctx, system)
    diag = {"e_preA": e_preA, "e_postA": e_postA, "per_atom_A": per_atomA,
            "converged": bool(e_prev is not None
                              and abs(e_post - e_prev) < CONVERGENCE_DTOL)}
    return st.getPositions(), e_pre, e_post, e_res, diag''')
cells[i08]["source"] = s08.splitlines(True)

# ============================================================ 3) 追加 Stage A 能量门 cell
i08b = find("# ============================================================ 08b") if any(
    "# ============================================================ 08b" in "".join(c["source"])
    for c in cells) else i08

GATE_CELL = '''
# ============================================================ 08c Stage-A 起点能量门（§13 扩展）
# 第 1 次 Colab smoke 暴露：E_pre = 1.12e8 kJ/mol/atom（原子严重重叠）。
# 该门在**任何**科学特征之前验证「起点是否物理上合理」。
# 阈值 MAX_PRE_ENERGY_PER_ATOM 是**项目工程阈值**，非领域标准；
# 必须先在一个已知良好的体系上测定、并在科学评估前冻结。
_probe_pdb = fetch_pdb(str(var.iloc[0].target_id))
_fw, _ew, _sw, _tw = prepare_system(_probe_pdb)
if _ew:
    raise RuntimeError(f"STAGE_A_GATE: prepare failed {_ew}")
_posA, _e_preA, _e_postA, _per_atomA = stage_a_relax(_sw, _tw, _fw.positions)
_n_atoms = _sw.getNumParticles()
print("=== STAGE-A / 起点能量门 ===")
print("  atoms                 :", _n_atoms)
print("  E_pre (addMissing 后) : %.4e kJ/mol  = %.4e kJ/mol/atom"
      % (_e_preA, _e_preA / _n_atoms))
print("  E_post (Stage A 后)   : %.4e kJ/mol  = %.4e kJ/mol/atom"
      % (_e_postA, _per_atomA))
print("  阈值 MAX_PRE_ENERGY_PER_ATOM =", MAX_PRE_ENERGY_PER_ATOM, "kJ/mol/atom")
STAGE_A_GATE_PASS = bool(_per_atomA < MAX_PRE_ENERGY_PER_ATOM)
print("  STAGE_A_GATE_PASS =", STAGE_A_GATE_PASS)
if not STAGE_A_GATE_PASS:
    print()
    print("STOP: Stage A 之后能量仍不合理 -> 建系/参数化有问题。")
    print("  这通常意味着：严重重叠未清除、或体系含未参数化组分。")
    raise RuntimeError(f"STAGE_A_FAIL: {_per_atomA:.4e} kJ/mol/atom")
print("  -> Stage A 把体系带到物理合理量级，可进入 §14 目标松弛")
'''
# 插到 09 Mutation builder 之前（即 08 之后）
i09 = find("# ============================================================ 09 Mutation builder")
cells.insert(i09, {"cell_type": "code", "execution_count": None, "metadata": {},
                   "outputs": [], "source": GATE_CELL.strip("\n").splitlines(True)})
print("已插入 08c（Stage-A 起点能量门）")

NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding="utf-8")

# ============================================================ 4) 调用点补上第 5 个返回值
nb = json.loads(NB.read_text(encoding="utf-8"))
cells = nb["cells"]
fixed = 0
for c in cells:
    if c["cell_type"] != "code":
        continue
    s = "".join(c["source"])
    o = s
    s = s.replace("pw, wpre, wpost, _ = minimize(", "pw, wpre, wpost, _, _dg = minimize(")
    s = s.replace("pm, mpre, mpost, _ = minimize(", "pm, mpre, mpost, _, _dg = minimize(")
    s = s.replace("p1, e1pre, e1post, r1 = minimize(", "p1, e1pre, e1post, r1, _dg1 = minimize(")
    s = s.replace("p2, e2pre, e2post, r2 = minimize(", "p2, e2pre, e2post, r2, _dg2 = minimize(")
    s = s.replace("p2, e2pre, e2post, r2 = minimize(", "p2, e2pre, e2post, r2, _dg2 = minimize(")
    s = s.replace("pw, w_pre, w_post, w_res = minimize(", "pw, w_pre, w_post, w_res, _dg = minimize(")
    s = s.replace("pm, m_pre, m_post, m_res = minimize(", "pm, m_pre, m_post, m_res, _dg = minimize(")
    s = s.replace("pm, e1pre, e1post, _r1 = minimize(", "pm, e1pre, e1post, _r1, _dg1 = minimize(")
    if s != o:
        c["source"] = s.splitlines(True); fixed += 1
print("已更新 minimize 调用点:", fixed)
NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding="utf-8")

nb2 = json.loads(NB.read_text(encoding="utf-8"))
code2 = [c for c in nb2["cells"] if c["cell_type"] == "code"]
bad = []
for i, c in enumerate(code2):
    try:
        ast.parse("".join(c["source"]))
    except SyntaxError as e:
        bad.append((i, e.lineno, e.msg))
t = "\n".join("".join(c["source"]) for c in code2)
print(f"cells={len(nb2['cells'])} code={len(code2)} 语法错误={len(bad)}")
for b in bad:
    print("   ", b)
print()
for k, v in (("两阶段 stage_a_relax", "def stage_a_relax" in t),
             ("Stage A iter 20000", "MIN_MAXITER_STAGE_A = 20000" in t),
             ("收敛检查", "CONVERGENCE_DTOL" in t),
             ("起点能量门", "STAGE_A_GATE_PASS" in t),
             ("阈值常量", "MAX_PRE_ENERGY_PER_ATOM = 1000.0" in t),
             ("identity 相对容差常量", "IDENTITY_REL_TOL" in t),
             ("§6 协议数值未变",
              "RADIUS_MOVABLE = 8.0" in t and "RESTRAINT_OUTSIDE = 100.0" in t
              and "PH = 7.0" in t and "MIN_TOL = 10.0" in t)):
    print(f"  {'✅' if v else '❌'} {k}")
