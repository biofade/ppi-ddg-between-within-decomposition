# -*- coding: utf-8 -*-
"""scripts/v7_4/02g_gbxml_probe.py -- cell 06 改为**官方 XML GBn2 路线**探测（无 fallback）.

按用户指定的决策树实现：

```
OpenMM 8.6.1 dev
   -> ForceField(ff14SB, implicit/gbn2.xml) -> createSystem(...)
        PASS -> PRIMARY = GBn2, 无协议偏离
        FAIL/SEGV -> STOP（不自动跳 OBC2）
```

cell 06 打印用户要求的完整 `GB_XML_PROBE` 报告：
    OpenMM version / ForceField files / Process rc / Force classes /
    CustomGBForce present / Particle count / Energy / Force finite /
    Context creation / Minimization / PASS-FAIL

**硬门控**：rc==0、存在 `CustomGBForce`、gbn2.xml 成功加载、能量 finite、
forces finite、Context 建成、最小化能跑、provenance 记 `GB_MODEL=GBn2_XML`。

**绝不允许 fallback**（OBC2/核心类都不再自动启用）。
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


PROBE_SRC = r'''
# ============================================================ 06 Official-XML GBn2 probe (NO FALLBACK)
# 用户指定的唯一主路径：
#     ForceField("amber14/protein.ff14SB.xml", "implicit/gbn2.xml")
#     -> createSystem(topology, nonbondedMethod=NoCutoff, constraints=HBonds)
# 不再使用 implicitSolvent= 参数（OpenMM 8.6 的 createSystem 已不消费它）。
# **硬门控 + 无 fallback**：失败即 STOP，绝不静默改用 OBC2（那会引入第二个变量）。
import subprocess as _sp, tempfile, os as _os, json as _json

GB_PROBE_CHILD = r"""
import sys, json
import numpy as np
import openmm as mm
import openmm.app as app
from openmm import unit
out = {}
def emit():
    print("@@JSON@@" + json.dumps(out), flush=True)
try:
    from pdbfixer import PDBFixer
    fx = PDBFixer(filename=sys.argv[1])
    fx.removeHeterogens(keepWater=False)
    fx.findMissingResidues()
    fx.findMissingAtoms()
    fx.addMissingAtoms()
    fx.addMissingHydrogens(7.0)
    topo, pos = fx.topology, fx.positions
    out["atoms"] = topo.getNumAtoms()
except Exception as e:
    out["error"] = "PDBFIXER:%s:%s" % (type(e).__name__, str(e)[:200]); emit(); sys.exit(1)
FILES = ["amber14/protein.ff14SB.xml", "implicit/gbn2.xml"]
out["files"] = FILES
out["openmm"] = mm.version.version
try:
    ff = app.ForceField(*FILES)
    out["generators"] = [type(g).__name__ for g in ff.getGenerators()]
    system = ff.createSystem(topo, nonbondedMethod=app.NoCutoff,
                             constraints=app.HBonds)
except Exception as e:
    out["error"] = "CREATESYSTEM:%s:%s" % (type(e).__name__, str(e)[:250]); emit(); sys.exit(1)
names = [system.getForce(i).getName() for i in range(system.getNumForces())]
out["forces"] = names
out["gb_forces"] = [n for n in names if "GB" in n.upper() or "BORN" in n.upper()]
try:
    out["customgb"] = bool(any(isinstance(system.getForce(i), mm.CustomGBForce)
                               for i in range(system.getNumForces())))
except Exception:
    out["customgb"] = False
out["particles"] = system.getNumParticles()
try:
    ctx = mm.Context(system, mm.VerletIntegrator(0.001),
                     mm.Platform.getPlatformByName("CPU"))
    ctx.setPositions(pos)
    st = ctx.getState(getEnergy=True, getForces=True)
    out["energy"] = float(st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    F = np.asarray(st.getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer))
    out["forces_finite"] = bool(np.all(np.isfinite(F)))
    out["energy_finite"] = bool(np.isfinite(out["energy"]))
    out["context"] = True
    mm.LocalEnergyMinimizer.minimize(ctx, tolerance=10.0, maxIterations=2000)
    out["minimized"] = True
    out["energy_post"] = float(ctx.getState(getEnergy=True).getPotentialEnergy()
                              .value_in_unit(unit.kilojoule_per_mole))
except Exception as e:
    out["error"] = "CONTEXT:%s:%s" % (type(e).__name__, str(e)[:200]); emit(); sys.exit(1)
emit()
"""
fd, _p = tempfile.mkstemp(suffix=".py"); _os.close(fd)
with open(_p, "w", encoding="utf-8") as _fh:
    _fh.write(GB_PROBE_CHILD)

print("GB_XML_PROBE")
print("  OpenMM version: (见子进程输出)")
print("  ForceField files:")
print("    amber14/protein.ff14SB.xml")
print("    implicit/gbn2.xml")
_pdb0 = fetch_pdb(VAL.iloc[0].target_id)
_r = _sp.run([sys.executable, _p, str(_pdb0)], capture_output=True, text=True, timeout=1800)
_rc = _r.returncode
INFO = {}
for _ln in (_r.stdout or "").splitlines():
    if _ln.startswith("@@JSON@@"):
        try:
            INFO = _json.loads(_ln[len("@@JSON@@"):])
        except Exception:
            pass
print("  Process rc: %d%s" % (_rc, "  (崩溃/段错误?)"
                              if _rc not in (0, 1) else ""))
if INFO:
    print("  OpenMM version: %s" % INFO.get("openmm"))
    print("  Force classes: %s" % INFO.get("forces"))
    print("  GB-like forces: %s" % INFO.get("gb_forces"))
    print("  CustomGBForce present: %s" % INFO.get("customgb"))
    print("  Particle count: %s" % INFO.get("particles"))
    print("  Energy: %s kJ/mol" % INFO.get("energy"))
    print("  Force finite: %s" % INFO.get("forces_finite"))
    print("  Context creation: %s" % INFO.get("context"))
    print("  Minimization: %s" % INFO.get("minimized"))
    if INFO.get("error"):
        print("  ERROR: %s" % INFO["error"])
elif _r.stderr:
    print("  stderr tail: %s" % _r.stderr.strip().splitlines()[-1][:200])

GB_XML_OK = bool(
    _rc == 0
    and INFO.get("customgb") is True
    and INFO.get("energy_finite") is True
    and INFO.get("forces_finite") is True
    and INFO.get("context") is True
    and INFO.get("minimized") is True
)
GB_MODEL = "GBn2_XML" if GB_XML_OK else None
print("  PASS/FAIL: %s" % ("PASS" if GB_XML_OK else "FAIL"))
try:
    _os.remove(_p)
except Exception:
    pass

if not GB_XML_OK:
    print()
    print("STOP: 官方 XML GBn2 路线未通过硬门控。")
    print("  依用户指定的决策树，**不自动 fallback 到 OBC2**。")
    print("  下一步二选一（需人工决定）：")
    print("    (a) 钉一个经 smoke test 确认可工作的稳定 OpenMM 版本以保住 GBn2；")
    print("    (b) 正式提出 protocol amendment 改用 OBC2 作为 sensitivity。")
    print("  已完成的诊断信息如上（rc / forces / CustomGBForce / energy / ctx / min）。")
    raise RuntimeError("GB_XML_PROBE_FAIL —— 依决策树 STOP（未做模型替换）")

print("  协议状态: 无偏离（PRIMARY = GBn2，官方 XML 装载）")

def prepare_system(pdb_path, muts=None):
    """建 WT/MUT 体系（ff14SB + implicit/gbn2.xml，官方 XML 路线）。"""
    try:
        fixer = PDBFixer(filename=str(pdb_path))
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(PH)
    except Exception as e:
        return None, f"FIXER:{type(e).__name__}:{e}", None, None
    if muts:
        ok, info = build_mutant(fixer, muts)
        if not ok:
            return None, info, None, None
    try:
        ff = app.ForceField(*FORCEFIELD)
        system = ff.createSystem(fixer.topology, nonbondedMethod=app.NoCutoff,
                                 constraints=app.HBonds)
    except Exception as e:
        return None, f"FORCEFIELD:{type(e).__name__}:{e}", None, None
    global SOLVENT_MODE
    SOLVENT_MODE = "implicit-solvent:GBn2_XML"
    return fixer, None, system, fixer.topology

SOLVENT_MODE = None

# ---- 参数化冒烟（WT / MUT 各一次）----
v0 = VAL.iloc[0]
pdb0 = fetch_pdb(v0.target_id)
fx, err, sysm, topo = prepare_system(pdb0)
FF_OK = err is None
print()
print("WT 参数化:", "PASS" if FF_OK else "FAIL")
if not FF_OK:
    print("  error:", err)
else:
    print("  溶剂:", SOLVENT_MODE)
    print("  原子数:", sysm.getNumParticles(), "| 力数:", sysm.getNumForces())
    print("  力列表:", [sysm.getForce(i).getName() for i in range(sysm.getNumForces())])
    muts = [(str(v0.chain), int(v0.resseq), THREE[str(v0.wt)], THREE[str(v0.mut)])]
    fx2, err2, sys2, topo2 = prepare_system(pdb0, muts)
    print("MUT 参数化:", "PASS" if err2 is None else f"FAIL {err2}")
    if err2 is None:
        print("  MUT 原子数:", sys2.getNumParticles())
'''

i06 = find("# ============================================================ 06 OpenMM Parameterization")
cells[i06]["source"] = PROBE_SRC.strip("\n").splitlines(True)

# provenance：明确 GB_MODEL=GBn2_XML
i_pv = find("# ============================================================ 15 Export")
s15 = "".join(cells[i_pv]["source"])
s15 = s15.replace('"gb_model_applied": GB_MODEL,',
                  '"gb_model_applied": GB_MODEL,\n'
                  '        "gb_model_source": "official ForceField XML (implicit/gbn2.xml)",')
cells[i_pv]["source"] = s15.splitlines(True)

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
for k, v in (("官方 XML 路线（无 implicitSolvent）",
              'ForceField(*FILES)' in t and "implicitSolvent=app.GBn2" not in t),
             ("子进程隔离", "GB_PROBE_CHILD" in t),
             ("硬门控 CustomGBForce", "customgb" in t),
             ("硬门控 energy/forces finite", "forces_finite" in t),
             ("失败即 STOP 不 fallback", "GB_XML_PROBE_FAIL" in t),
             ("无 OBC2 自动回落",
              "GBSAOBC2Force" not in t and "mm.GBSAOBCForce" not in t),
             ("GB_MODEL=GBn2_XML", 'GB_MODEL = "GBn2_XML"' in t),
             ("provenance 记 XML 来源", "official ForceField XML" in t)):
    print(f"  {'✅' if v else '❌'} {k}")
