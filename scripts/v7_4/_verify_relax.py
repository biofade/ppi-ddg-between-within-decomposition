# -*- coding: utf-8 -*-
"""scripts/v7_4/_verify_relax.py -- 端到端验证 V7.4 松弛栈（本地 OpenMM 8.6.1.dev）.

验证内容（与 notebook 的 cell 06/07/08b 相同的代码路径）：
  1. PDBFixer 建 WT 体系（ff14SB）
  2. **GBSAOBC2Force** 隐式溶剂（替代已不可用的 GBn2）
  3. CustomExternalForce 位置限制（§14 冻结策略）
  4. L-BFGS 最小化 + **identity(WT→WT) 数值噪声地板**（§31）
  5. 突变体构建 + 配对最小化
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
from pdbfixer import PDBFixer
from openmm.app.internal.customgbforces import GBSAOBC2Force

PH = 7.0
RADIUS_MOVABLE = 8.0
K_SHELL_KCAL = 10.0
K_OUT_KCAL = 100.0
KCAL_A2_TO_KJ_NM2 = 4.184 / (0.1 ** 2)
K_SHELL = K_SHELL_KCAL * KCAL_A2_TO_KJ_NM2
K_OUT = K_OUT_KCAL * KCAL_A2_TO_KJ_NM2
MIN_TOL = 10.0
MIN_MAXITER = 2000

THREE = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
         "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
         "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
         "Y": "TYR", "V": "VAL"}

FIXER_CACHE: dict = {}


def build_topology(pdb_path, muts=None):
    """建加氢拓扑（含可选突变）。注意 PDBFixer 有全局状态，故每次新建。"""
    f = PDBFixer(filename=pdb_path)
    f.removeHeterogens(keepWater=False)
    f.findMissingResidues()
    f.findMissingAtoms()
    f.addMissingAtoms()
    f.addMissingHydrogens(PH)
    if muts:
        chains = {c.id: c for c in f.topology.chains()}
        for (ch, rs, wt3, mut3) in muts:
            tgt = None
            for res in chains[ch].residues():
                if res.id == str(rs):
                    tgt = res
                    break
            if tgt is None:
                return None, f"residue {ch}:{rs} not found"
            if tgt.name != wt3:
                return None, f"WT mismatch {ch}:{rs} got {tgt.name} want {wt3}"
            tgt.name = mut3
    return f, None


def make_system(topology):
    """ff14SB + GBSAOBC2Force（GBn2 在 8.6.1.dev 不可用）。"""
    ff = app.ForceField("amber14/protein.ff14SB.xml")
    system = ff.createSystem(topology, nonbondedMethod=app.NoCutoff,
                             constraints=None)
    nbf = next(system.getForce(i) for i in range(system.getNumForces())
               if system.getForce(i).getName() == "NonbondedForce")
    gb = GBSAOBC2Force()
    gb.setName("GBSAOBC2")
    gb.setNonbondedMethod(GBSAOBC2Force.NoCutoff)
    n = topology.getNumAtoms()
    gb.addParticles(n)
    charges = np.zeros(n)
    radii = np.zeros(n)
    scales = np.zeros(n)
    for i in range(n):
        q, sig, eps = nbf.getParticleParameters(i)
        charges[i] = q.value_in_unit(unit.elementary_charge)
        radii[i] = 0.15
        scales[i] = 0.8
    gb.setParticleParameters(0, charges, radii, scales)
    names = [gb.getGlobalParameterName(i) for i in range(gb.getNumGlobalParameters())]
    for i, nm in enumerate(names):
        if "olvent" in nm.lower():
            gb.setGlobalParameterDefaultValue(i, 78.5)
        if "olute" in nm.lower():
            gb.setGlobalParameterDefaultValue(i, 1.0)
    system.addForce(gb)
    return system


def residues_within(topology, positions, site_keys, radius):
    pos = np.array([list(p) for p in positions.value_in_unit(unit.nanometer)])
    site_atoms = [a.index for a in topology.atoms()
                  if (a.residue.chain.id, a.residue.id) in site_keys]
    if not site_atoms:
        return set()
    r2 = (radius / 10.0) ** 2
    out = set()
    for a in topology.atoms():
        key = (a.residue.chain.id, a.residue.id)
        if key in site_keys:
            continue
        if float(np.min(np.sum((pos[site_atoms] - pos[a.index]) ** 2, axis=1))) <= r2:
            out.add(key)
    return out


def add_restraints(system, topology, positions, site_keys, shell_keys):
    f = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    f.setName("Restraint")
    for pname in ("k", "x0", "y0", "z0"):
        f.addPerParticleParameter(pname)
    pos = np.array([list(p) for p in positions.value_in_unit(unit.nanometer)])
    for a in topology.atoms():
        key = (a.residue.chain.id, a.residue.id)
        if key in site_keys:
            k = 0.0
        elif key in shell_keys:
            k = K_SHELL if a.name in ("N", "CA", "C", "O") else 0.0
        else:
            k = K_OUT
        f.addParticle(a.index, [k, float(pos[a.index][0]),
                                float(pos[a.index][1]), float(pos[a.index][2])])
    f.setForceGroup(31)
    system.addForce(f)
    return f


def energies(ctx, system):
    """分离 restraint（group 31）与生物力场能量（group 0..30）。"""
    sys_grp = {0}
    for i, f in enumerate(system.getForces()):
        if f.getForceGroup() != 31:
            sys_grp.add(i)
    e_bio = ctx.getState(getEnergy=True, groups=sys_grp
                         ).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    e_res = ctx.getState(getEnergy=True, groups={31}
                         ).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    return e_bio, e_res


PLAT = mm.Platform.getPlatformByName("CPU")


def minimize(system, topology, positions, site_keys, shell_keys):
    add_restraints(system, topology, positions, site_keys, shell_keys)
    for i, f in enumerate(system.getForces()):
        if f.getForceGroup() != 31:
            f.setForceGroup(0)
    integ = mm.VerletIntegrator(0.001)
    ctx = mm.Context(system, integ, PLAT)
    ctx.setPositions(positions)
    e_pre, _ = energies(ctx, system)
    mm.LocalEnergyMinimizer.minimize(ctx, tolerance=MIN_TOL,
                                     maxIterations=MIN_MAXITER)
    st = ctx.getState(getPositions=True, getEnergy=True)
    e_post, e_restr = energies(ctx, system)
    return st.getPositions(), e_pre, e_post, e_restr


PDB = r"data/pdb/1A22.pdb"
print("=== 1) WT 构建 + 参数化 ===")
t0 = time.time()
fwt, err = build_topology(PDB)
print("  build:", "OK" if err is None else err, "| 原子:", fwt.topology.getNumAtoms(),
      "| %.1fs" % (time.time() - t0))
sysw = make_system(fwt.topology)
print("  力:", [sysw.getForce(i).getName() for i in range(sysw.getNumForces())])
print("  PASS ✅" if any(sysw.getForce(i).getName() == "GBSAOBC2"
                        for i in range(sysw.getNumForces())) else "  FAIL ❌")

print()
print("=== 2) 单位点 WT 最小化 + identity 数值地板（§31）===")
site = {("A", "10")}
shell = residues_within(fwt.topology, fwt.positions, site, RADIUS_MOVABLE)
print("  壳内残基数:", len(shell))
t1 = time.time()
pw1, ep1, eq1, er1 = minimize(sysw, fwt.topology, fwt.positions, site, shell)
print("  WT relax #1: E_pre=%.1f  E_post=%.1f  dE=%.1f  restraint=%.2f  (%.1fs)"
      % (ep1, eq1, eq1 - ep1, er1, time.time() - t1))
# 第二次：完全相同的输入 -> 重复性（§32）
fwt2, _ = build_topology(PDB)
sysw2 = make_system(fwt2.topology)
shell2 = residues_within(fwt2.topology, fwt2.positions, site, RADIUS_MOVABLE)
t2 = time.time()
pw2, ep2, eq2, er2 = minimize(sysw2, fwt2.topology, fwt2.positions, site, shell2)
dE = eq2 - eq1
a1 = np.array(pw1.value_in_unit(unit.angstrom))
a2 = np.array(pw2.value_in_unit(unit.angstrom))
rmsd = float(np.sqrt(np.mean(np.sum((a1 - a2) ** 2, axis=1))))
print("  WT relax #2: E_post=%.1f  (%.1fs)" % (eq2, time.time() - t2))
print("  >>> identity/重复性: |dE| = %.4f kJ/mol   RMSD = %.4f Å" % (abs(dE), rmsd))

print()
print("=== 3) 突变体构建 + 配对最小化 ===")
mut = [("A", 10, "F", "A")]
fm, err = build_topology(PDB, mut)
print("  MUT build:", "OK" if err is None else err,
      "| 原子:", fm.topology.getNumAtoms() if fm else None)
if fm is not None:
    sysm = make_system(fm.topology)
    shellm = residues_within(fm.topology, fm.positions, site, RADIUS_MOVABLE)
    t3 = time.time()
    pm, epm, eqm, erm = minimize(sysm, fm.topology, fm.positions, site, shellm)
    print("  MUT relax: E_pre=%.1f E_post=%.1f dE=%.1f restraint=%.2f (%.1fs)"
          % (epm, eqm, eqm - epm, erm, time.time() - t3))
    print("  >>> F10A: dE_total(mut-wt) = %.2f kJ/mol ; dE_relax_cost = %.2f"
          % (eqm - eq1, (eqm - epm) - (eq1 - ep1)))

print()
print("=== 结论 ===")
print("  GBn2 在 OpenMM 8.6.1.dev 不可用（implicitSolvent 参数已移除，")
print("  且 GBSAGBn2Force 会段错误）-> 改用 **GBSAOBC2Force**（OBC2 隐式溶剂）")
