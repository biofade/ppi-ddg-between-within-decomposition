# -*- coding: utf-8 -*-
"""scripts/v7_4/_probe_gbn2.py -- 本地复现 Kaggle 的 GBn2 参数化路径（一次性诊断工具）.

本地已装 **与 Kaggle 完全相同** 的 OpenMM 8.6.1.dev-b399af4。
目的：确定 `implicitSolvent=app.GBn2` 到底能不能用，以及真正失败在哪一步。
"""

from __future__ import annotations

import sys
import traceback

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

import openmm
import openmm.app as app
from openmm.app.internal.customgbforces import GBSAGBn2Force

print("OpenMM:", openmm.version.version)

try:
    from pdbfixer import PDBFixer
    HAVE_FIXER = True
    print("PDBFixer: 可用")
except ImportError as e:
    HAVE_FIXER = False
    print("PDBFixer: 不可用", e)

if not HAVE_FIXER:
    sys.exit(3)

PDB = r"data/pdb/1A22.pdb"
f = PDBFixer(filename=PDB)
f.removeHeterogens(keepWater=False)
f.findMissingResidues()
f.findMissingAtoms()
f.addMissingAtoms()
f.addMissingHydrogens(7.0)
print("加氢后原子数:", f.topology.getNumAtoms())
print()

CASES = [
    ("A) 仅 ff14SB + implicitSolvent=app.GBn2",
     ["amber14/protein.ff14SB.xml"], dict(implicitSolvent=app.GBn2)),
    ("B) ff14SB+gbn2.xml + implicitSolvent=app.GBn2",
     ["amber14/protein.ff14SB.xml", "implicit/gbn2.xml"],
     dict(implicitSolvent=app.GBn2)),
    ("C) 仅 ff14SB，不传溶剂",
     ["amber14/protein.ff14SB.xml"], {}),
    ("D) 仅 ff14SB + solute/solventDielectric",
     ["amber14/protein.ff14SB.xml"],
     dict(soluteDielectric=1.0, solventDielectric=78.5)),
]
for label, ffargs, kw in CASES:
    try:
        ff = app.ForceField(*ffargs)
        s = ff.createSystem(f.topology, nonbondedMethod=app.NoCutoff,
                            constraints=None, **kw)
        names = [s.getForce(i).getName() for i in range(s.getNumForces())]
        print(f"{label:52s} PASS  forces={names}")
    except Exception as e:  # noqa: BLE001
        print(f"{label:52s} FAIL  {type(e).__name__}: {str(e)[:180]}")

print()
print("=== 正确的 GBn2 力类是什么 ===")
gb = GBSAGBn2Force()
print("  GBSAGBn2Force() 构造成功:", type(gb).__name__)
meths = [m for m in dir(gb) if "ielectric" in m or "et" == m[-3:] or "ett" in m]
print("  可用设置方法候选:", [m for m in dir(gb)
                             if not m.startswith("_")][:40])
try:
    gb.setSolventDielectric(78.5)
    print("  setSolventDielectric: OK")
except Exception as e:  # noqa: BLE001
    print("  setSolventDielectric 失败:", e)
try:
    print("  getSolventDielectric:", gb.getSolventDielectric())
except Exception as e:  # noqa: BLE001
    print("  getSolventDielectric 失败:", e)
