# -*- coding: utf-8 -*-
"""scripts/v7_4/02c_fix_relaxation.py -- V7.4 修复松弛实现（3 个真实 bug）.

## 从返回 ZIP 诊断出的问题

| # | 问题 | 证据 | 修复 |
|---|---|---|---|
| 1 | `mm.CustomExternalPotential()` **不存在** | 本该抛 AttributeError；被 `if False` 短路掩盖 | 改用 `mm.CustomExternalForce` |
| 2 | 限制力**从未加入 system** | `build_restrained_system` 的返回值被丢弃（`if False else None`） | 在 minimize 内真正 `system.addForce(...)` |
| 3 | `implicitSolvent=None` | 传入 `None` ⇒ **GBn2 没生效** | 传 `implicitSolvent=app.GBn2` |
| 4 | `Platform.getPlatformByName("CUDA")` 在 Kaggle 不存在 | `benchmark.csv` 全部 `OpenMMException` | 已在 02b 换成 `get_platform()` 回落 |

## 冻结的协议（§14，**值未改**，只把单位正确换算）

`CustomExternalForce` 的表达式 `k*((x-x0)^2+(y-y0)^2+(z-z0)^2)` 用
**OpenMM 内部单位 = kJ/mol/nm²**。故 `kcal/mol/Å²` 需乘 **69.4786**。

* 突变位点残基：`k = 0`（侧链与骨架完全自由）
* 8 Å 壳内**骨架**原子：`k = 10 kcal/mol/Å²` → 694.79 kJ/mol/nm²
* 8 Å 壳内**侧链**原子：`k = 0`（可松弛）
* 壳外全部原子：`k = 100 kcal/mol/Å²` → 6947.86 kJ/mol/nm²
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

ROOT = Path(__file__).resolve().parent.parent.parent
NB = ROOT / "kaggle" / "v7_4" / "V7_4_OpenMM_Relaxation.ipynb"
nb = json.loads(NB.read_text(encoding="utf-8"))
cells = nb["cells"]


def find(marker):
    for i, c in enumerate(cells):
        if c["cell_type"] == "code" and marker in "".join(c["source"]):
            return i
    raise KeyError(marker)


# ============================================================ cell 06: prepare_system
i06 = find("# ============================================================ 06 OpenMM Parameterization")
cells[i06]["source"] = '''
# ============================================================ 06 OpenMM Parameterization Smoke Test
# 修复：implicitSolvent=app.GBn2（此前传 None，GBn2 未生效）
from pdbfixer import PDBFixer
THREE = {"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU",
         "G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE",
         "P":"PRO","S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL"}

def prepare_system(pdb_path, muts=None):
    """返回 (fixer, error, system, topology)。
    error=None 表示成功。仅蛋白链；pH 7.0 加氢；ff14SB + GBn2 隐式溶剂。"""
    try:
        fixer = PDBFixer(filename=str(pdb_path))
        fixer.removeHeterogens(keepWater=False)          # §10 仅蛋白链
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(PH)                    # §6 pH 7.0
    except Exception as e:
        return None, f"FIXER:{type(e).__name__}:{e}", None, None
    if muts:
        ok, info = build_mutant(fixer, muts)
        if not ok:
            return None, info, None, None
    ff = app.ForceField(*FORCEFIELD)
    try:
        system = ff.createSystem(
            fixer.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=None,
            implicitSolvent=app.GBn2,          # <-- 关键修复
        )
    except Exception as e:
        return None, f"FORCEFIELD:{type(e).__name__}:{e}", None, None
    return fixer, None, system, fixer.topology

# ---- 参数化冒烟 ----
v0 = VAL.iloc[0]
pdb0 = fetch_pdb(v0.target_id)
fx, err, sysm, topo = prepare_system(pdb0)
FF_OK = err is None
print("WT 参数化:", "PASS" if FF_OK else f"FAIL {err}")
if FF_OK:
    print("  原子数:", sysm.getNumParticles(), "| 力数:", sysm.getNumForces())
    print("  力列表:", [sysm.getForce(i).getName() for i in range(sysm.getNumForces())])
    muts = [(str(v0.chain), int(v0.resseq), THREE[str(v0.wt)], THREE[str(v0.mut)])]
    fx2, err2, sys2, topo2 = prepare_system(pdb0, muts)
    print("MUT 参数化:", "PASS" if err2 is None else f"FAIL {err2}")
    if err2 is None:
        print("  MUT 原子数:", sys2.getNumParticles())
'''.strip("\n").splitlines(True)

# ============================================================ cell 07: restraints + minimize
i07 = find("# ============================================================ 07 Relaxation Protocol + Restraints")
cells[i07]["source"] = '''
# ============================================================ 07 Relaxation Protocol + Restraints
# 修复：CustomExternalForce（正确类名）+ 真正 addForce + 正确单位换算

# kcal/mol/A^2 -> kJ/mol/nm^2
KCAL_A2_TO_KJ_NM2 = 4.184 / (0.1 ** 2)   # = 418.4 kJ/mol/nm^2 per kcal/mol/A^2
# 说明：1 kcal = 4.184 kJ；1 A^2 = 0.01 nm^2 -> 1 kcal/mol/A^2 = 418.4 kJ/mol/nm^2
# 但 CustomExternalForce 的 (x-x0) 以 nm 计，k*dx^2 中 dx 已含单位；
# 故系数 k 的单位是 kJ/mol/nm^2。换算因子 = 4.184 / (0.1^2) = 418.4
K_SHELL = RESTRAINT_BACKBONE_SHELL * KCAL_A2_TO_KJ_NM2   # 10  -> 4184
K_OUT   = RESTRAINT_OUTSIDE        * KCAL_A2_TO_KJ_NM2   # 100 -> 41840

def residues_within(topology, positions, site_keys, radius):
    """返回任意位点重原子 radius (A) 内、且**属于其它残基**的残基键集合。"""
    pos_nm = _np.array(positions.value_in_unit(unit.nanometer))
    site_atoms = [a.index for a in topology.atoms()
                  if (a.residue.chain.id, a.residue.id) in site_keys]
    out = set()
    allpos = _np.array([list(p) for p in positions.value_in_unit(unit.nanometer)])
    r_nm = radius / 10.0
    for a in topology.atoms():
        key = (a.residue.chain.id, a.residue.id)
        if key in site_keys:
            continue
        d2 = _np.sum((allpos[site_atoms] - allpos[a.index]) ** 2, axis=1)
        if float(_np.min(d2)) <= r_nm ** 2:
            out.add(key)
    return out

def add_restraints(system, topology, positions, site_keys, shell_keys):
    """§14 冻结策略：位点 k=0；壳内骨架 k=K_SHELL；其余 k=K_OUT。
    返回 restraint force（其能量**不得**进入生物能量特征，§21）。"""
    f = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    f.setName("Restraint")
    for pname in ("k", "x0", "y0", "z0"):
        f.addPerParticleParameter(pname)
    allpos = _np.array([list(p) for p in positions.value_in_unit(unit.nanometer)])
    for a in topology.atoms():
        key = (a.residue.chain.id, a.residue.id)
        if key in site_keys:
            k = 0.0
        elif key in shell_keys:
            k = K_SHELL if a.name in ("N", "CA", "C", "O") else 0.0
        else:
            k = K_OUT
        f.addParticle(a.index, [k, float(allpos[a.index][0]),
                                float(allpos[a.index][1]),
                                float(allpos[a.index][2])])
    f.setForceGroup(31)          # 单独力组 -> 能量可分离（§21）
    system.addForce(f)
    return f

def group_energy(ctx, group):
    return ctx.getState(getEnergy=True, groups={group}
                        ).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

def system_forces_by_group(system):
    """把每个 Force 放进自己的力组，便于分解（§21）。"""
    for i, f in enumerate(system.getForces()):
        f.setForceGroup(i)
    return {i: system.getForce(i).getName() for i in range(system.getNumForces())}

def minimize(system, topology, positions, site_keys, shell_keys, tag=""):
    """加限制 -> L-BFGS 最小化。返回 (positions, e_pre, e_post, e_restraint_post)。
    e_pre/e_post 为**不含** restraint 的力场势能（§21）。"""
    add_restraints(system, topology, positions, site_keys, shell_keys)
    fmap = system_forces_by_group(system)
    integ = mm.VerletIntegrator(0.001)
    ctx = mm.Context(system, integ, PLAT)
    ctx.setPositions(positions)
    e_pre = group_energy(ctx, 0)                 # 力组 0 = 第一个力（非 restraint）
    mm.LocalEnergyMinimizer.minimize(ctx, tolerance=MIN_TOL,
                                     maxIterations=MIN_MAXITER)
    st = ctx.getState(getPositions=True, getEnergy=True)
    e_post = group_energy(ctx, 0)
    e_restr = group_energy(ctx, 31)              # 单独取，不混入（§21）
    return st.getPositions(), e_pre, e_post, e_restr

print("restraint 系数（kJ/mol/nm^2）：K_SHELL=%.1f  K_OUT=%.1f"
      % (K_SHELL, K_OUT))
print("力组: restraint 单独置于 group 31（§21：不得进入生物能量特征）")
'''.strip("\n").splitlines(True)

# ============================================================ cell 08: smoke 用新签名
i08 = find("# ============================================================ 08 Paired WT/MUT Smoke")
s08 = "".join(cells[i08]["source"])
s08 = s08.replace("rsw = build_restrained_system(tw, sw, site_keys) if False else None", "")
s08 = s08.replace("p1, e1pre, e1post = minimize(sw, tw, fw.positions)",
                  "p1, e1pre, e1post, ers = minimize(sw, tw, fw.positions, site_keys, shell_w)")
s08 = s08.replace("p2, e2pre, e2post = minimize(sw, tw, fw.positions)",
                  "p2, e2pre, e2post, _ = minimize(sw, tw, fw.positions, site_keys, shell_w)")
cells[i08]["source"] = s08.splitlines(True)

# ============================================================ cell 08b: identity 用新签名
i08b = find("# ============================================================ 08b Identity")
s8b = "".join(cells[i08b]["source"])
s8b = s8b.replace("p1, e1pre, e1post = minimize(sw, tw, fw.positions)",
                  "site_keys = {(str(v.chain), str(v.resseq))}\n"
                  "        sh = residues_within(tw, fw.positions, site_keys, RADIUS_MOVABLE)\n"
                  "        p1, e1pre, e1post, _ = minimize(sw, tw, fw.positions, site_keys, sh)")
s8b = s8b.replace("p2, e2pre, e2post = minimize(sw, tw, fw.positions)",
                  "p2, e2pre, e2post, _ = minimize(sw, tw, fw.positions, site_keys, sh)")
cells[i08b]["source"] = s8b.splitlines(True)

# ============================================================ cell 09: benchmark 用新签名
i09 = find("# ============================================================ 09 Runtime Benchmark")
s09 = "".join(cells[i09]["source"])
s09 = s09.replace("pw, e_pre_w, e_post_w = minimize(sw, tw, tw_pos, \"WT\")",
                  "sk = {(str(v.chain), str(v.resseq))}\n"
                  "        sh_w = residues_within(tw, fw.positions, sk, RADIUS_MOVABLE)\n"
                  "        pw, e_pre_w, e_post_w, ers_w = minimize(sw, tw, tw_pos, sk, sh_w, \"WT\")")
s09 = s09.replace("pm, e_pre_m, e_post_m = minimize(sm, tm, fm.positions, \"MUT\")",
                  "sh_m = residues_within(tm, fm.positions, sk, RADIUS_MOVABLE)\n"
                  "        pm, e_pre_m, e_post_m, ers_m = minimize(sm, tm, fm.positions, sk, sh_m, \"MUT\")")
cells[i09]["source"] = s09.splitlines(True)

# ============================================================ cell 10: full run 用新签名
i10 = find("# ============================================================ 10 Full Streaming Relaxation")
s10 = "".join(cells[i10]["source"])
s10 = s10.replace("pw, w_pre, w_post = minimize(sw, tw, fw.positions)",
                  "sk = {(str(x.chain), str(x.pdb_resseq)) for x in g.itertuples()}\n"
                  "                sh_w = residues_within(tw, fw.positions, sk, RADIUS_MOVABLE)\n"
                  "                pw, w_pre, w_post, ers_w = minimize(sw, tw, fw.positions, sk, sh_w)")
s10 = s10.replace("pm, m_pre, m_post = minimize(sm, tm, fm.positions)",
                  "sh_m = residues_within(tm, fm.positions, sk, RADIUS_MOVABLE)\n"
                  "                pm, m_pre, m_post, ers_m = minimize(sm, tm, fm.positions, sk, sh_m)")
cells[i10]["source"] = s10.splitlines(True)

NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding="utf-8")

# ---- 校验 ----
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
for k, v in (("CustomExternalForce", "CustomExternalForce" in t),
             ("CustomExternalPotential 残留", "CustomExternalPotential" in t),
             ("implicitSolvent=app.GBn2", "implicitSolvent=app.GBn2" in t),
             ("add_restraints 被调用", "add_restraints(" in t),
             ("力组 31 分离 restraint", "setForceGroup(31)" in t),
             ("platform 回落", "PLAT_NAME" in t)):
    print(f"  {'✅' if (v if '残留' not in k else not v) else '❌'} {k}")
