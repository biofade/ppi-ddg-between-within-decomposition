# -*- coding: utf-8 -*-
"""scripts/v7_4/_probe_gbxml.py -- V7.4 官方 XML GBn2 路线探测（**无 fallback**）.

只测一条路径，完全按用户指定的写法：

    ff = app.ForceField("amber14/protein.ff14SB.xml", "implicit/gbn2.xml")
    system = ff.createSystem(topology, nonbondedMethod=app.NoCutoff,
                             constraints=app.HBonds)

**硬门控**（全部满足才 PASS）：
  rc == 0
  | CustomGBForce 存在于 system
  | implicit/gbn2.xml 成功加载
  | 能量 finite
  | forces finite
  | Context 创建成功
  | 最小化能跑

**不允许 fallback**：失败即 FAIL，绝不自动跳 OBC2。

必须作为**子进程**运行（该 dev build 的 GB 代码有段错误风险）。
"""

from __future__ import annotations

import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

import numpy as np

TAG = sys.argv[1] if len(sys.argv) > 1 else "gbn2"
PDB = sys.argv[2] if len(sys.argv) > 2 else r"data/pdb/1A22.pdb"
# TAG 决定单独加载 GB XML 还是仅 ff14SB（对照）
GB_XML = {"gbn2": "implicit/gbn2.xml",
          "obc2": "implicit/obc2.xml",
          "obc1": "implicit/obc1.xml",
          "none": None}.get(TAG, "implicit/gbn2.xml")


def p(*a):
    print(*a, flush=True)


p("GB_XML_PROBE")
p("TAG: %s  (gb xml = %s)" % (TAG, GB_XML))
p("")

import openmm as mm
import openmm.app as app
from openmm import unit

p("OpenMM version: %s" % mm.version.version)
files = ["amber14/protein.ff14SB.xml"] + ([GB_XML] if GB_XML else [])
p("ForceField files:")
for f in files:
    p("  %s" % f)

# --- 建拓扑（加氢）---
try:
    from pdbfixer import PDBFixer
    fx = PDBFixer(filename=PDB)
    fx.removeHeterogens(keepWater=False)
    # NOTE: findMissingAtoms()/addMissingAtoms() segfault on this Windows dev build
    # (independent of GB model). 1A22 already has all heavy atoms, so we only add H.
    fx.addMissingHydrogens(7.0)
    topo = fx.topology
    pos = fx.positions
    p("PDBFixer: OK  atoms=%d" % topo.getNumAtoms())
except Exception as e:  # noqa: BLE001
    p("PDBFixer 失败: %s: %s" % (type(e).__name__, e))
    sys.exit(1)

# --- 官方 XML 路线（不传 implicitSolvent）---
p("")
p("--- ForceField(%s) -> createSystem(NoCutoff, HBonds) ---" % ", ".join(files))
system = None
try:
    ff = app.ForceField(*files)
    # gbn2.xml 是否真的加载：看 generator / 模板 / 全局参数痕迹
    try:
        gen_names = [type(g).__name__ for g in ff.getGenerators()]
    except Exception as e:  # noqa: BLE001
        gen_names = ["<err %s>" % e]
    p("ForceField generators: %s" % gen_names)
    system = ff.createSystem(topo, nonbondedMethod=app.NoCutoff,
                             constraints=app.HBonds)
    p("createSystem: OK")
except Exception as e:  # noqa: BLE001
    p("createSystem 失败: %s: %s" % (type(e).__name__, str(e)[:300]))
    p("")
    p("PASS/FAIL: FAIL")
    sys.exit(1)

# --- 力清单 ---
names = [system.getForce(i).getName() for i in range(system.getNumForces())]
p("Force classes: %s" % names)
try:
    gb_present = any(isinstance(system.getForce(i), mm.CustomGBForce)
                     for i in range(system.getNumForces()))
except Exception:  # noqa: BLE001
    gb_present = False
gb_named = [n for n in names if "GB" in n.upper() or "BORN" in n.upper()]
p("CustomGBForce present: %s" % gb_present)
p("GB-like forces: %s" % gb_named)
p("Particle count: %d" % system.getNumParticles())
try:
    p("Global parameters: %s" % [system.getGlobalParameterName(i)
                                 for i in range(system.getNumGlobalParameters())])
except Exception:  # noqa: BLE001
    pass

# --- 全局参数默认值（看溶剂介电常数是否 78.5）---
p("")
p("--- Context + 能量 ---")
ctx = None
e = np.nan
try:
    integ = mm.VerletIntegrator(0.001)
    ctx = mm.Context(system, integ, mm.Platform.getPlatformByName("CPU"))
    ctx.setPositions(pos)
    st = ctx.getState(getEnergy=True, getForces=True)
    e = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    f = st.getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer)
    fmax = float(np.max(np.abs(np.asarray(f))))
    p("Context creation: OK")
    p("Energy: %.3f kJ/mol" % e)
    p("Force finite: %s (max|F| = %.3e kJ/mol/nm)"
      % (bool(np.all(np.isfinite(np.asarray(f)))), fmax))
except Exception as ex:  # noqa: BLE001
    p("Context/能量失败: %s: %s" % (type(ex).__name__, str(ex)[:200]))
    p("")
    p("PASS/FAIL: FAIL")
    sys.exit(1)

# --- 最小化 ---
p("")
p("--- Minimization ---")
try:
    e_pre = e
    mm.LocalEnergyMinimizer.minimize(ctx, tolerance=10.0, maxIterations=2000)
    st2 = ctx.getState(getEnergy=True)
    e_post = st2.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    p("Minimization: OK  E_pre=%.1f  E_post=%.1f  dE=%.1f" % (e_pre, e_post, e_post - e_pre))
except Exception as ex:  # noqa: BLE001
    p("Minimization 失败: %s: %s" % (type(ex).__name__, str(ex)[:200]))
    p("")
    p("PASS/FAIL: FAIL")
    sys.exit(1)

# --- 硬门控 ---
ok = (system is not None and gb_present and np.isfinite(e)
      and np.all(np.isfinite(np.asarray(f))))
p("")
p("=== HARD GATE ===")
p("  rc                     = 0")
p("  CustomGBForce present  = %s" % gb_present)
p("  gbn2.xml loaded        = %s" % bool(GB_XML))
p("  energy finite          = %s" % bool(np.isfinite(e)))
p("  forces finite          = %s" % bool(np.all(np.isfinite(np.asarray(f)))))
p("  context created        = %s" % (ctx is not None))
p("  minimisation ran       = True")
p("")
p("PASS/FAIL: %s" % ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
