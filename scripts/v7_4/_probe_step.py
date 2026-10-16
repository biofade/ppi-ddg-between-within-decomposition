# -*- coding: utf-8 -*-
"""scripts/v7_4/_probe_step.py -- 单步探测（避免无 Python traceback 的段错误）。

用法: python _probe_step.py <step>
步骤: sig | gen | setattr | manual | obc2
**不**用 PDBFixer（避免 C 层崩溃），改用手写 PDBFile 以便聚焦 API。
"""

from __future__ import annotations

import inspect
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

import openmm
import openmm.app as app

STEP = sys.argv[1] if len(sys.argv) > 1 else "sig"
print("STEP =", STEP, "| OpenMM", openmm.version.version, flush=True)

try:
    from pdbfixer import PDBFixer
    f = PDBFixer(filename=r"data/pdb/1A22.pdb")
    f.removeHeterogens(keepWater=False)
    f.findMissingResidues(); f.findMissingAtoms(); f.addMissingAtoms()
    f.addMissingHydrogens(7.0)
    topo = f.topology
    print("原子数:", topo.getNumAtoms(), flush=True)
except Exception as e:  # noqa: BLE001
    print("PDBFixer 失败:", e, flush=True)
    sys.exit(2)

if STEP == "sig":
    for g in app.ForceField("amber14/protein.ff14SB.xml").getGenerators():
        print(" gen:", type(g).__name__, flush=True)
        try:
            print("   createForce:", str(inspect.signature(g.createForce))[:600], flush=True)
        except Exception as e:  # noqa: BLE001
            print("   (no sig)", e, flush=True)
        print("   溶剂属性:", [a for a in dir(g)
                              if "olvent" in a or "ielectric" in a or a == "implicitSolvent"],
              flush=True)

elif STEP == "setattr":
    ff = app.ForceField("amber14/protein.ff14SB.xml")
    gen = next((g for g in ff.getGenerators()
                if type(g).__name__ == "NonbondedGenerator"), None)
    print("gen:", type(gen).__name__, flush=True)
    ok = hasattr(gen, "implicitSolvent")
    print("has implicitSolvent:", ok, flush=True)
    if ok:
        setattr(gen, "implicitSolvent", app.GBn2)
        s = ff.createSystem(topo, nonbondedMethod=app.NoCutoff, constraints=None)
        print("forces:", [s.getForce(i).getName() for i in range(s.getNumForces())],
              flush=True)

elif STEP == "manual":
    from openmm.app.internal.customgbforces import GBSAGBn2Force
    s = app.ForceField("amber14/protein.ff14SB.xml").createSystem(
        topo, nonbondedMethod=app.NoCutoff, constraints=None)
    print("base forces:", [s.getForce(i).getName() for i in range(s.getNumForces())],
          flush=True)
    gb = GBSAGBn2Force()
    gb.setName("GBn2")
    names = [gb.getGlobalParameterName(i) for i in range(gb.getNumGlobalParameters())]
    print("GBn2 全局参数:", names, flush=True)
    for i, n in enumerate(names):
        if "olvent" in n.lower():
            gb.setGlobalParameterDefaultValue(i, 78.5)
        if "olute" in n.lower():
            gb.setGlobalParameterDefaultValue(i, 1.0)
    s.addForce(gb)
    print("加力成功:", [s.getForce(i).getName() for i in range(s.getNumForces())],
          flush=True)

elif STEP == "obc2":
    from openmm.app.internal.customgbforces import GBSAOBC2Force
    gb = GBSAOBC2Force()
    names = [gb.getGlobalParameterName(i) for i in range(gb.getNumGlobalParameters())]
    print("GBSAOBC2Force 全局参数:", names, flush=True)

elif STEP == "core":
    gb = openmm.GBSAOBCForce()
    gb.setSolventDielectric(78.5)
    gb.setSoluteDielectric(1.0)
    print("core GBSAOBCForce OK", flush=True)

print("STEP DONE", flush=True)
