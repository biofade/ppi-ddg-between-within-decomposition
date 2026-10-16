# -*- coding: utf-8 -*-
"""scripts/v7_4r/_local_stageA_test.py -- 本地实测 Stage A 行为（绕开本机 PDBFixer 段错误）.

本机 `PDBFixer.findMissingAtoms()` 段错误 ⇒ **不用 PDBFixer**。
改用 OpenMM 的 `app.PDBFile` + `Modeller.addHydrogens`（对**完整**结构可用）。
目的：在提交 Colab 之前，先验证
  1) 起点能量是否也出现 e11 量级的荒谬值
  2) Stage A（无限制最小化）能否把能量带回合理量级
  3) identity 二次最小化的**可复现性**（这才是噪声门的核心）
"""

from __future__ import annotations

import sys
import time

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

import numpy as np
import openmm as mm
import openmm.app as app
from openmm import unit

PLAT = mm.Platform.getPlatformByName("CPU")
PH = 7.0
PDB = sys.argv[1] if len(sys.argv) > 1 else r"data/pdb/1A22.pdb"


def p(*a):
    print(*a, flush=True)


p("OpenMM:", mm.version.version, "| platform CPU")
pdbf = app.PDBFile(PDB)
p("PDBFile OK: atoms=%d chains=%d" % (pdbf.topology.getNumAtoms(),
                                      pdbf.topology.getNumChains()))

modeller = app.Modeller(pdbf.topology, pdbf.positions)
modeller.addHydrogens(app.ForceField("amber14/protein.ff14SB.xml"), pH=PH)
p("addHydrogens OK: atoms=%d" % modeller.topology.getNumAtoms())

ff = app.ForceField("amber14/protein.ff14SB.xml", "implicit/gbn2.xml")
system = ff.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff,
                         constraints=app.HBonds)
names = [system.getForce(i).getName() for i in range(system.getNumForces())]
p("forces:", names)
gbcount = sum(1 for i in range(system.getNumForces())
              if isinstance(system.getForce(i), mm.CustomGBForce))
p("CustomGBForce count =", gbcount)
n = system.getNumParticles()


def energy(pos, tol=None, iters=0):
    integ = mm.VerletIntegrator(0.001)
    ctx = mm.Context(system, integ, PLAT)
    ctx.setPositions(pos)
    if iters:
        mm.LocalEnergyMinimizer.minimize(ctx, tolerance=tol, maxIterations=iters)
    st = ctx.getState(getPositions=True, getEnergy=True)
    e = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    return st.getPositions(), e


p("")
p("=== 起点能量（addHydrogens 之后，未最小化）===")
t0 = time.time()
pos0, e0 = energy(modeller.positions)
p("  E_start = %.6e kJ/mol  = %.6e kJ/mol/atom   (%.1fs)"
  % (e0, e0 / n, time.time() - t0))

p("")
p("=== Stage A：无限制最小化（模拟新版 stage_a_relax）===")
for iters in (100, 500, 2000, 20000):
    t0 = time.time()
    posA, eA = energy(modeller.positions, tol=10.0, iters=iters)
    p("  iters=%-6d E=%.6e kJ/mol = %.6e kJ/mol/atom  (%.1fs)"
      % (iters, eA, eA / n, time.time() - t0))

p("")
p("=== 可复现性：同一输入连跑两次（噪声门核心）===")
posA1, eA1 = energy(modeller.positions, tol=10.0, iters=20000)
posA2, eA2 = energy(modeller.positions, tol=10.0, iters=20000)
a1 = np.array(posA1.value_in_unit(unit.angstrom))
a2 = np.array(posA2.value_in_unit(unit.angstrom))
rmsd = float(np.sqrt(np.mean(np.sum((a1 - a2) ** 2, axis=1))))
p("  run1 E = %.6e" % eA1)
p("  run2 E = %.6e" % eA2)
p("  |dE|   = %.6e kJ/mol   (rel = %.3e)" % (abs(eA2 - eA1),
                                              abs(eA2 - eA1) / max(1e-12, abs(eA2))))
p("  RMSD   = %.6f A" % rmsd)
p("")
p("=== 判读 ===")
p("  E_start/atom 是否 e8 量级？  %s" % ("是（与 Colab 一致）" if e0 / n > 1e6 else "否"))
p("  最小化后是否回到 ~10 kJ/mol/atom？ %s"
  % ("是" if abs(eA2) / n < 100 else "否"))
p("  两次是否可复现（rel<1e-2）？ %s"
  % ("是" if abs(eA2 - eA1) / max(1e-12, abs(eA2)) < 1e-2 else "否"))
