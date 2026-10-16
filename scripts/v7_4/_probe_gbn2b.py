# -*- coding: utf-8 -*-
"""scripts/v7_4/_probe_gbn2b.py -- 找出 OpenMM 8.6.1.dev 下**可用**的 GB 路径."""

from __future__ import annotations

import inspect
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

import openmm
import openmm.app as app
from pdbfixer import PDBFixer

print("OpenMM:", openmm.version.version)
f = PDBFixer(filename=r"data/pdb/1A22.pdb")
f.removeHeterogens(keepWater=False)
f.findMissingResidues(); f.findMissingAtoms(); f.addMissingAtoms()
f.addMissingHydrogens(7.0)
topo = f.topology
print("原子数:", topo.getNumAtoms())

print()
print("=== NonbondedGenerator 是否支持 implicitSolvent ===")
ff = app.ForceField("amber14/protein.ff14SB.xml")
gen = None
for g in ff.getGenerators():
    if type(g).__name__ == "NonbondedGenerator":
        gen = g
        break
print(" NonbondedGenerator 找到:", gen is not None)
if gen is not None:
    try:
        print(" createForce 签名:", str(inspect.signature(gen.createForce))[:400])
    except Exception as e:
        print(" 签名不可读:", e)
    attrs = [a for a in dir(gen) if "olvent" in a or "ielectric" in a or "GB" in a]
    print(" 溶剂相关属性:", attrs)

print()
print("=== 直接给 NonbondedGenerator 设 implicitSolvent ===")
for name in ("implicitSolvent",):
    if gen is not None and hasattr(gen, name):
        try:
            setattr(gen, name, app.GBn2)
            s = ff.createSystem(topo, nonbondedMethod=app.NoCutoff, constraints=None)
            print("  设 gen.%s=GBn2 后 createSystem ->" % name,
                  [s.getForce(i).getName() for i in range(s.getNumForces())])
        except Exception as e:
            print("  设 gen.%s 失败: %s: %s" % (name, type(e).__name__, str(e)[:140]))
    else:
        print("  gen 无属性", name)

print()
print("=== 手工加 GBSAGBn2Force（正确类）===")
from openmm.app.internal.customgbforces import GBSAGBn2Force
try:
    s = app.ForceField("amber14/protein.ff14SB.xml").createSystem(
        topo, nonbondedMethod=app.NoCutoff, constraints=None)
    gb = GBSAGBn2Force()
    gb.setName("GBn2")
    # CustomGBForce 的介电常数是 global parameter
    names = [gb.getGlobalParameterName(i) for i in range(gb.getNumGlobalParameters())]
    print("  GBn2 全局参数:", names)
    for i, n in enumerate(names):
        if "olvent" in n.lower():
            gb.setGlobalParameterDefaultValue(i, 78.5)
            print("   设 %s = 78.5" % n)
        if "olute" in n.lower():
            gb.setGlobalParameterDefaultValue(i, 1.0)
            print("   设 %s = 1.0" % n)
    s.addForce(gb)
    print("  加力成功 ✅ 力:", [s.getForce(i).getName() for i in range(s.getNumForces())])
except Exception as e:
    import traceback
    print("  失败 ❌", type(e).__name__, str(e)[:200])
    traceback.print_exc()

print()
print("=== GBSAOBCForce (core) 是否可用 ===")
try:
    g2 = openmm.GBSAOBCForce()
    g2.setSolventDielectric(78.5)
    g2.setSoluteDielectric(1.0)
    print("  core GBSAOBCForce ✅ 可用（OBC1 模型）")
except Exception as e:
    print("  core GBSAOBCForce ❌", e)

print()
print("=== GBSAOBC2Force（8.6 内部类）===")
try:
    from openmm.app.internal.customgbforces import GBSAOBC2Force
    g3 = GBSAOBC2Force()
    n3 = [g3.getGlobalParameterName(i) for i in range(g3.getNumGlobalParameters())]
    print("  GBSAOBC2Force ✅ 全局参数:", n3)
except Exception as e:
    print("  GBSAOBC2Force ❌", e)
