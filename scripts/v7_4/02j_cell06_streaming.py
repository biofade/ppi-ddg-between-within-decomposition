# -*- coding: utf-8 -*-
"""scripts/v7_4/02j_cell06_streaming.py -- cell 06 改为**流式可观测**（修 Server Error）.

## 问题（用户实测）

Kaggle 跑到 cell 05 之后、cell 06 就 **Server Error**（两次一致）。
cell 06 之前是**静默黑盒**：子进程把所有输出缓存在 stdout，
只有全部跑完才一次性返回。于是：
* 出任何问题 -> 父 cell 什么都拿不到
* 子进程卡住   -> cell 一直等 -> Kaggle kernel 超时/被杀 -> `Server Error`

## 修复

1. **子进程逐步写进度文件**（`gb_probe_progress.log`），父进程**实时读取并转发**。
   这样即使子进程段错误/卡死，我们也已看到「死在哪一步」。
2. **加硬超时**（默认 900 s）—— 卡死不再拖垮 kernel。
3. **子进程不再 stdout 缓冲**：每步 `print(flush=True)` 并同时 append 到进度文件。
4. 父进程即使拿不到最终 JSON，也把**已收到的进度**打印出来。
5. 保留全部硬门控（`CustomGBForce count>=1`、`GB particles==system particles`、
   finite、Context、Minimization、`DeltaE_min`），**无 fallback**。
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


SRC = r'''
# ============================================================ 06 Official-XML GBn2 probe (NO FALLBACK)
# 官方路线：ForceField(ff14SB, "implicit/gbn2.xml") -> createSystem()，**不传** implicitSolvent。
# 子进程隔离 + **流式进度**（防静默黑盒导致 Server Error）+ 硬超时。
# **无 fallback**：失败即 STOP，绝不静默改用 OBC2。
import subprocess as _sp, tempfile, os as _os, json as _json, time as _time

PROBE_TIMEOUT = 900          # s：卡死不再拖垮 kernel
_tmpdir = tempfile.mkdtemp()
PROG = _os.path.join(_tmpdir, "gb_probe_progress.log")
CHILD = _os.path.join(_tmpdir, "gb_probe_child.py")

CHILD_SRC = r"""
import sys, json, os, time
PROG = sys.argv[2]
def step(msg):
    line = "[%7.2fs] %s" % (time.time() - T0, msg)
    print(line, flush=True)
    try:
        with open(PROG, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")
    except Exception:
        pass
T0 = time.time()
out = {}
def emit():
    print("@@JSON@@" + json.dumps(out), flush=True)
step("start")
try:
    import numpy as np
    step("numpy ok")
    import openmm as mm
    import openmm.app as app
    from openmm import unit
    step("openmm ok: %s" % mm.version.version)
    out["openmm"] = mm.version.version
except Exception as e:
    out["error"] = "IMPORT:%s:%s" % (type(e).__name__, str(e)[:200]); emit(); sys.exit(1)
try:
    from pdbfixer import PDBFixer
    step("pdbfixer imported")
    fx = PDBFixer(filename=sys.argv[1])
    step("PDBFixer constructed")
    fx.removeHeterogens(keepWater=False); step("removeHeterogens")
    fx.findMissingResidues();             step("findMissingResidues")
    miss = sum(len(v) for v in (fx.missingAtoms or {}).values())
    step("missingAtoms=%d missingResidues=%d"
         % (miss, sum(len(v) for v in (fx.missingResidues or {}).values())))
    fx.findMissingAtoms();                 step("findMissingAtoms")
    fx.addMissingAtoms();                  step("addMissingAtoms")
    fx.addMissingHydrogens(7.0);           step("addMissingHydrogens")
    topo, pos = fx.topology, fx.positions
    out["atoms"] = topo.getNumAtoms()
    step("topology atoms=%d" % out["atoms"])
except Exception as e:
    out["error"] = "PDBFIXER:%s:%s" % (type(e).__name__, str(e)[:250]); emit(); sys.exit(1)
FILES = ["amber14/protein.ff14SB.xml", "implicit/gbn2.xml"]
out["files"] = FILES
try:
    ff = app.ForceField(*FILES); step("ForceField loaded")
    out["generators"] = [type(g).__name__ for g in ff.getGenerators()]
    step("generators=%s" % out["generators"])
    system = ff.createSystem(topo, nonbondedMethod=app.NoCutoff,
                             constraints=app.HBonds)
    step("createSystem OK")
except Exception as e:
    out["error"] = "CREATESYSTEM:%s:%s" % (type(e).__name__, str(e)[:250]); emit(); sys.exit(1)
names = [system.getForce(i).getName() for i in range(system.getNumForces())]
out["forces"] = names
out["force_groups"] = {names[i]: system.getForce(i).getForceGroup()
                       for i in range(system.getNumForces())}
out["gb_forces"] = [n for n in names if "GB" in n.upper() or "BORN" in n.upper()]
step("forces=%s" % names)
gb_idx = []
try:
    for i in range(system.getNumForces()):
        if isinstance(system.getForce(i), mm.CustomGBForce):
            gb_idx.append(i)
except Exception:
    gb_idx = []
out["customgb_count"] = len(gb_idx)
out["customgb"] = bool(len(gb_idx) >= 1)
out["gb_particles"] = (system.getForce(gb_idx[0]).getNumParticles() if gb_idx else None)
out["gb_force_group"] = (system.getForce(gb_idx[0]).getForceGroup() if gb_idx else None)
out["particles"] = system.getNumParticles()
out["gb_covers_all"] = bool(gb_idx and out["gb_particles"] == out["particles"])
step("CustomGBForce count=%s gb_particles=%s system_particles=%s"
     % (out["customgb_count"], out["gb_particles"], out["particles"]))
try:
    ctx = mm.Context(system, mm.VerletIntegrator(0.001),
                     mm.Platform.getPlatformByName(sys.argv[3]))
    step("Context created on %s" % sys.argv[3])
    ctx.setPositions(pos); step("positions set")
    st = ctx.getState(getEnergy=True, getForces=True)
    out["energy"] = float(st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    F = np.asarray(st.getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer))
    out["forces_finite"] = bool(np.all(np.isfinite(F)))
    out["energy_finite"] = bool(np.isfinite(out["energy"]))
    out["context"] = True
    step("energy=%.3f kJ/mol finite=%s forces_finite=%s"
         % (out["energy"], out["energy_finite"], out["forces_finite"]))
    out["E_pre_min"] = out["energy"]
    mm.LocalEnergyMinimizer.minimize(ctx, tolerance=10.0, maxIterations=2000)
    out["minimized"] = True
    out["E_post_min"] = float(ctx.getState(getEnergy=True).getPotentialEnergy()
                              .value_in_unit(unit.kilojoule_per_mole))
    out["DeltaE_min"] = out["E_post_min"] - out["E_pre_min"]
    step("minimized E_post=%.3f dE=%.3f" % (out["E_post_min"], out["DeltaE_min"]))
except Exception as e:
    out["error"] = "CONTEXT:%s:%s" % (type(e).__name__, str(e)[:250]); emit(); sys.exit(1)
emit()
step("done")
"""

with open(CHILD, "w", encoding="utf-8") as _fh:
    _fh.write(CHILD_SRC)

print("GB_XML_PROBE")
print("  ForceField files:")
print("    amber14/protein.ff14SB.xml")
print("    implicit/gbn2.xml")
print("  platform: %s" % PLAT_NAME)
print("  child timeout: %ds" % PROBE_TIMEOUT)
print()
print("--- 子进程实时进度 ---")

_pdb0 = fetch_pdb(VAL.iloc[0].target_id)
proc = _sp.Popen([sys.executable, "-u", CHILD, str(_pdb0), PROG, PLAT_NAME],
                 stdout=_sp.PIPE, stderr=_sp.STDOUT, text=True, bufsize=1)
_t0 = _time.time()
_lines = []
_timed_out = False
while True:
    if proc.poll() is not None:
        break
    if _time.time() - _t0 > PROBE_TIMEOUT:
        _timed_out = True
        proc.kill()
        break
    # 转发子进程 stdout（-u 已无缓冲）
    ln = proc.stdout.readline()
    if ln:
        _lines.append(ln.rstrip("\n"))
        if not ln.startswith("@@JSON@@"):
            print("  " + ln.rstrip("\n"))
    else:
        _time.sleep(0.2)
try:
    _rest = proc.stdout.read() or ""
    for ln in _rest.splitlines():
        _lines.append(ln)
        if not ln.startswith("@@JSON@@"):
            print("  " + ln)
except Exception:
    pass
_rc = proc.returncode
# 兜底：把进度文件里父进程还没看到的内容补打
try:
    with open(PROG, "r", encoding="utf-8") as fh:
        _prog = fh.read().strip().splitlines()
    _seen = set(_lines)
    for ln in _prog:
        if ln and ln not in _seen:
            print("  [prog] " + ln)
except Exception:
    pass

INFO = {}
for _ln in _lines:
    if _ln.startswith("@@JSON@@"):
        try:
            INFO = _json.loads(_ln[len("@@JSON@@"):])
        except Exception:
            pass

print()
print("  Process rc: %d%s" % (_rc, "  (超时被 kill)" if _timed_out else ""))
if INFO:
    print("  OpenMM version: %s" % INFO.get("openmm"))
    print("  Force classes: %s" % INFO.get("forces"))
    print("  GB-like forces: %s" % INFO.get("gb_forces"))
    print("  CustomGBForce present: %s" % INFO.get("customgb"))
    print("  CustomGBForce count: %s" % INFO.get("customgb_count"))
    print("  GB particles: %s" % INFO.get("gb_particles"))
    print("  System particles: %s" % INFO.get("particles"))
    print("  GB covers all particles: %s" % INFO.get("gb_covers_all"))
    print("  GB force group: %s" % INFO.get("gb_force_group"))
    print("  Force groups: %s" % INFO.get("force_groups"))
    print("  Energy finite: %s" % INFO.get("energy_finite"))
    print("  Force finite: %s" % INFO.get("forces_finite"))
    print("  Context creation: %s" % INFO.get("context"))
    print("  Minimization: %s" % INFO.get("minimized"))
    print("  E_pre_min : %s kJ/mol" % INFO.get("E_pre_min"))
    print("  E_post_min: %s kJ/mol" % INFO.get("E_post_min"))
    print("  DeltaE_min: %s kJ/mol" % INFO.get("DeltaE_min"))
    if INFO.get("error"):
        print("  ERROR: %s" % INFO["error"])
else:
    print("  未收到子进程 JSON 结果（可能崩溃/超时）—— 以上进度即失败点。")

_dE = INFO.get("DeltaE_min")
GB_XML_OK = bool(
    _rc == 0
    and not _timed_out
    and INFO.get("customgb") is True
    and INFO.get("customgb_count", 0) >= 1
    and INFO.get("gb_covers_all") is True
    and INFO.get("energy_finite") is True
    and INFO.get("forces_finite") is True
    and INFO.get("context") is True
    and INFO.get("minimized") is True
    and (_dE is not None)
    and ((_dE <= 0.0) if (_dE is not None and (INFO.get("energy") or 0.0) > 1e6)
         else (_dE <= 0.05 * abs(INFO.get("energy") or 1.0)))
)
GB_MODEL = "GBn2_XML" if GB_XML_OK else None
print("  PASS/FAIL: %s" % ("PASS" if GB_XML_OK else "FAIL"))

if not GB_XML_OK:
    print()
    print("STOP: 官方 XML GBn2 路线未通过硬门控（GBN2_XML_PROTOCOL_FAIL）。")
    print("  依决策树，**不自动 fallback 到 OBC2**。")
    print("  下一步二选一（需人工决定）：")
    print("    (a) 钉一个经 smoke test 确认可工作的稳定 OpenMM 版本以保住 GBn2；")
    print("    (b) 正式提出 protocol amendment 改用 OBC2 作为 sensitivity。")
    print("  上面的实时进度已指出失败发生在哪一步。")
    raise RuntimeError("GBN2_XML_PROTOCOL_FAIL —— 依决策树 STOP（未做模型替换）")

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
print()
v0 = VAL.iloc[0]
pdb0 = fetch_pdb(v0.target_id)
fx, err, sysm, topo = prepare_system(pdb0)
FF_OK = err is None
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

i06 = find("# ============================================================ 06 Official-XML GBn2 probe")
cells[i06]["source"] = SRC.strip("\n").splitlines(True)
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
for k, v in (("流式进度文件", "gb_probe_progress.log" in t),
             ("硬超时", "PROBE_TIMEOUT" in t),
             ("逐行转发子进程输出", "readline()" in t),
             ("子进程 -u 无缓冲", '"-u", CHILD' in t),
             ("超时即 kill", "proc.kill()" in t),
             ("拿不到 JSON 也打印进度", "未收到子进程 JSON" in t),
             ("硬门控 count>=1", "customgb_count" in t),
             ("硬门控 particles 相等", "gb_covers_all" in t),
             ("DeltaE_min 门控", "DeltaE_min" in t),
             ("无 OBC2 fallback", "GBSAOBC2Force" not in t),
             ("失败即 STOP", "GBN2_XML_PROTOCOL_FAIL" in t)):
    print(f"  {'✅' if v else '❌'} {k}")
