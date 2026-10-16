# -*- coding: utf-8 -*-
"""scripts/v7_4/_probe_obc2.py -- 精确隔离 GBSAOBC2Force 的段错误点."""

from __future__ import annotations

import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

import numpy as np

def p(*a):
    print(*a, flush=True)

import openmm as mm
import openmm.app as app
from openmm import unit
from pdbfixer import PDBFixer
from openmm.app.internal.customgbforces import GBSAOBC2Force

p("import 完成")
f = PDBFixer(filename=r"data/pdb/1A22.pdb")
f.removeHeterogens(keepWater=False)
f.findMissingResidues(); f.findMissingAtoms(); f.addMissingAtoms()
f.addMissingHydrogens(7.0)
topo = f.topology
n = topo.getNumAtoms()
p("topology 原子数:", n)

ff = app.ForceField("amber14/protein.ff14SB.xml")
system = ff.createSystem(topo, nonbondedMethod=app.NoCutoff, constraints=None)
p("base system 力:", [system.getForce(i).getName()
                      for i in range(system.getNumForces())])

nbf = next(system.getForce(i) for i in range(system.getNumForces())
           if system.getForce(i).getName() == "NonbondedForce")

p("--- A) GBSAOBC2Force() 构造 + createForce ---")
gb = GBSAOBC2Force()
p("  构造 OK")
try:
    f2 = gb.createForce()
    p("  createForce OK ->", type(f2).__name__)
except Exception as e:  # noqa: BLE001
    p("  createForce 失败:", type(e).__name__, e)

p("--- B) addParticles(n) ---")
try:
    gb.addParticles(n)
    p("  addParticles OK, numParticles =", gb.getNumParticles())
except Exception as e:  # noqa: BLE001
    p("  addParticles 失败:", type(e).__name__, e)

p("--- C) setParticleParameters(0, charges, radii, scales) ---")
charges = np.zeros(n); radii = np.full(n, 0.15); scales = np.full(n, 0.8)
for i in range(n):
    q, sig, eps = nbf.getParticleParameters(i)
    charges[i] = q.value_in_unit(unit.elementary_charge)
try:
    gb.setParticleParameters(0, charges, radii, scales)
    p("  setParticleParameters OK")
except Exception as e:  # noqa: BLE001
    p("  setParticleParameters 失败:", type(e).__name__, e)

p("--- D) 全局参数 ---")
try:
    names = [gb.getGlobalParameterName(i) for i in range(gb.getNumGlobalParameters())]
    p("  全局参数:", names)
    for i, nm in enumerate(names):
        if "olvent" in nm.lower():
            gb.setGlobalParameterDefaultValue(i, 78.5)
        if "olute" in nm.lower():
            gb.setGlobalParameterDefaultValue(i, 1.0)
    p("  设默认值 OK")
except Exception as e:  # noqa: BLE001
    p("  全局参数失败:", type(e).__name__, e)

p("--- E) system.addForce(gb) ---")
try:
    system.addForce(gb)
    p("  addForce OK -> ", [system.getForce(i).getName()
                            for i in range(system.getNumForces())])
except Exception as e:  # noqa: BLE001
    p("  addForce 失败:", type(e).__name__, e)

p("--- F) Context 建 + 能量 ---")
try:
    integ = mm.VerletIntegrator(0.001)
    ctx = mm.Context(system, integ, mm.Platform.getPlatformByName("CPU"))
    ctx.setPositions(f.positions)
    e = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole)
    p("  能量 = %.1f kJ/mol" % e)
except Exception as e:  # noqa: BLE001
    p("  Context/能量失败:", type(e).__name__, e)

p("ALL DONE")
