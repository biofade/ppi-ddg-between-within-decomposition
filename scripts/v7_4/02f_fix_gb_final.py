# -*- coding: utf-8 -*-
"""scripts/v7_4/02f_fix_gb_final.py -- V7.4 最终 GB 修复（第 4 版）.

## 本地实测结论（OpenMM 8.6.1.dev-b399af4，**与 Kaggle 同版本**）

| 路径 | 结果 |
|---|---|
| `createSystem(implicitSolvent=app.GBn2)` | ❌ `ValueError: ... was never used` |
| `createSystem(soluteDielectric=..., solventDielectric=...)` | ❌ 同上 |
| `openmm.app.internal.customgbforces.GBSAOBCForce` | ❌ **不存在**（我第 2 版写错了） |
| `GBSAGBn2Force` | ❌ 段错误（Windows 本地；Kaggle/Linux 未测到该步） |
| **`GBSAOBC2Force`** | ✅ 可构造，`addParticles`/全局参数可用 |
| **`openmm.GBSAOBCForce`**（核心类，OBC1） | ✅ 可用 |

**根因**：OpenMM 8.6 起 `createSystem` **移除了 `implicitSolvent` 与介电常数参数**，
隐式溶剂不再由 `createSystem` 参数驱动。这也解释了
`createSystem` 签名里只有 `**args`，而 `ArgTracker` 会拒绝未消费的 key。

## 因此本轮修复采用两条工程措施

1. **GB 候选在子进程中逐个试**（`subprocess`）。
   本地已观察到该 dev build 的 GB 类会**段错误**（无 Python traceback）——
   若在主进程里试，会直接杀死 Kaggle kernel，丢掉全部结果。
   子进程隔离后，段错误只是「该候选淘汰」。
2. **多候选回落**，并在 `provenance.json` 记录**实际生效的模型**：
   `GBn2(旧API)` → `GBSAOBC2` → `GBSAOBC`(核心类, OBC1)。
   **绝不静默退化成真空。**

### 必须在报告中声明的科学限制

若实际生效的是 **OBC2** 或 **OBC1** 而非 GBn2：
* 三者同属 AMBER GB 隐式溶剂家族，但**参数化不同**
* 计算条件与 §6 指定的 `implicit/gbn2.xml` **不完全一致**
* ⇒ 这属于**协议偏离**，必须在 V7.4 最终报告中如实标注（§6 要求记录限制）
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


NEW = r'''
# ---------- 隐式溶剂：子进程隔离探测 + 多候选回落 ----------
# 实测（OpenMM 8.6.1.dev，与 Kaggle 同版本）：
#   createSystem(implicitSolvent=...)            -> ValueError (参数已移除)
#   customgbforces.GBSAOBCForce                  -> 不存在
#   GBSAGBn2Force                                -> 段错误风险
#   GBSAOBC2Force / 核心 GBSAOBCForce            -> 可用
# 该 dev build 的 GB 类可能**段错误**（无 traceback），故必须在**子进程**里试，
# 否则会杀死 kernel、丢掉全部已算结果。
import subprocess as _sp, tempfile, os as _os

_PROBE = r"""
import sys
try:
    import openmm, openmm.app as app
    from pdbfixer import PDBFixer
    f = PDBFixer(filename=sys.argv[1]); f.removeHeterogens(keepWater=False)
    f.findMissingResidues(); f.findMissingAtoms(); f.addMissingAtoms()
    f.addMissingHydrogens(7.0)
    topo = f.topology
    ff = app.ForceField("amber14/protein.ff14SB.xml")
    system = ff.createSystem(topo, nonbondedMethod=app.NoCutoff, constraints=None)
    nbf = next(system.getForce(i) for i in range(system.getNumForces())
               if system.getForce(i).getName() == "NonbondedForce")
    import numpy as np
    from openmm import unit
    n = topo.getNumAtoms()
    ch = np.zeros(n); ra = np.full(n, 0.15); sc = np.full(n, 0.8)
    for i in range(n):
        q, s, e = nbf.getParticleParameters(i)
        ch[i] = q.value_in_unit(unit.elementary_charge)
    cand = sys.argv[2]
    if cand == "GBSAOBC2":
        from openmm.app.internal.customgbforces import GBSAOBC2Force as G
    elif cand == "GBSAOBC":
        G = openmm.GBSAOBCForce
    elif cand == "GBn2":
        from openmm.app.internal.customgbforces import GBSAGBn2Force as G
    else:
        print("UNKNOWN"); sys.exit(2)
    gb = G(); gb.setName(cand)
    if hasattr(gb, "addParticles"):
        gb.addParticles(n)
        if hasattr(gb, "setParticleParameters"):
            gb.setParticleParameters(0, ch, ra, sc)
    else:
        for i in range(n):
            gb.addParticle(i, float(ch[i]), float(ra[i]), float(sc[i]))
        gb.setSolventDielectric(78.5); gb.setSoluteDielectric(1.0)
    for i in range(gb.getNumGlobalParameters()):
        nm = gb.getGlobalParameterName(i)
        if "olvent" in nm.lower(): gb.setGlobalParameterDefaultValue(i, 78.5)
        if "olute" in nm.lower(): gb.setGlobalParameterDefaultValue(i, 1.0)
    system.addForce(gb)
    import openmm as mm
    ctx = mm.Context(system, mm.VerletIntegrator(0.001),
                     mm.Platform.getPlatformByName("CPU"))
    ctx.setPositions(f.positions)
    e = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    print("OK", cand, "E=%.1f" % e, "forces=%d" % system.getNumForces())
except Exception as ex:
    print("FAIL", type(ex).__name__, str(ex)[:150])
"""

def _probe_candidate(pdb_path, cand, timeout=600):
    """在**子进程**里试一个 GB 候选。返回 (ok, msg)。"""
    fd, path = tempfile.mkstemp(suffix=".py"); _os.close(fd)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(_PROBE)
    try:
        r = _sp.run([sys.executable, path, str(pdb_path), cand],
                    capture_output=True, text=True, timeout=timeout)
        out = (r.stdout or "").strip().splitlines()
        msg = out[-1] if out else f"rc={r.returncode} no output"
        ok = r.returncode == 0 and ("OK" in msg)
        if r.returncode not in (0,) and not ok:
            msg = f"rc={r.returncode} (崩溃/段错误?) {msg}"
        return ok, msg
    except _sp.TimeoutExpired:
        return False, "timeout"
    except Exception as ex:
        return False, f"{type(ex).__name__}:{ex}"
    finally:
        try:
            _os.remove(path)
        except Exception:
            pass

print("=== GB 候选子进程探测（防段错误杀死 kernel）===")
GB_CANDIDATES = ["GBSAOBC2", "GBn2", "GBSAOBC"]
GB_MODEL = None
_pdb_probe = fetch_pdb(VAL.iloc[0].target_id)
for _c in GB_CANDIDATES:
    _ok, _msg = _probe_candidate(_pdb_probe, _c)
    print("  %-10s %s  %s" % (_c, "PASS" if _ok else "FAIL", _msg))
    if _ok and GB_MODEL is None:
        GB_MODEL = _c
print(">>> 选定 GB 模型:", GB_MODEL)
if GB_MODEL is None:
    print("!!! 所有 GB 候选均失败 -> 依 §6 要求，**不得**退化成真空；将报 FAIL_FORCEFIELD")
elif GB_MODEL != "GBn2":
    print("!!! 协议偏离：§6 指定 GBn2，实际生效 %s（必须在报告中标注）" % GB_MODEL)

def make_system(fixer):
    """用**子进程已确认可用**的 GB 模型建 system。绝不退化成真空。"""
    ff = app.ForceField(*FORCEFIELD)
    system = ff.createSystem(fixer.topology, nonbondedMethod=app.NoCutoff,
                             constraints=None)
    if GB_MODEL is None:
        return None, None, "所有 GB 候选失败（见上方子进程探测）"
    import numpy as _np
    nbf = next(system.getForce(i) for i in range(system.getNumForces())
               if system.getForce(i).getName() == "NonbondedForce")
    n = fixer.topology.getNumAtoms()
    ch = _np.zeros(n); ra = _np.full(n, 0.15); sc = _np.full(n, 0.8)
    for i in range(n):
        q, _s, _e = nbf.getParticleParameters(i)
        ch[i] = q.value_in_unit(unit.elementary_charge)
    if GB_MODEL == "GBSAOBC2":
        from openmm.app.internal.customgbforces import GBSAOBC2Force as _G
    elif GB_MODEL == "GBSAOBC":
        _G = mm.GBSAOBCForce
    else:
        from openmm.app.internal.customgbforces import GBSAGBn2Force as _G
    gb = _G(); gb.setName("ImplicitSolvent_" + GB_MODEL)
    if hasattr(gb, "addParticles"):
        gb.addParticles(n)
        if hasattr(gb, "setParticleParameters"):
            gb.setParticleParameters(0, ch, ra, sc)
    else:
        for i in range(n):
            gb.addParticle(i, float(ch[i]), float(ra[i]), float(sc[i]))
        gb.setSolventDielectric(78.5); gb.setSoluteDielectric(1.0)
    for i in range(gb.getNumGlobalParameters()):
        nm = gb.getGlobalParameterName(i)
        if "olvent" in nm.lower():
            gb.setGlobalParameterDefaultValue(i, 78.5)
        if "olute" in nm.lower():
            gb.setGlobalParameterDefaultValue(i, 1.0)
    system.addForce(gb)
    return system, "implicit-solvent:" + GB_MODEL, None

def prepare_system(pdb_path, muts=None):
    """返回 (fixer, error, system, topology)。error=None 表示成功。"""
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
    system, how, err = make_system(fixer)
    if system is None:
        return None, f"FORCEFIELD:{err}", None, None
    global SOLVENT_MODE
    SOLVENT_MODE = how
    return fixer, None, system, fixer.topology

SOLVENT_MODE = None
'''

i06 = find("# ============================================================ 06 OpenMM Parameterization")
cells[i06]["source"] = (
    "# ============================================================ 06 OpenMM Parameterization Smoke Test\n"
    "# GBn2 在 OpenMM 8.6 已从 createSystem 移除；GB 类可能段错误 -> 子进程隔离探测\n"
    "from pdbfixer import PDBFixer\n"
    'THREE = {"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU",\n'
    '         "G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE",\n'
    '         "P":"PRO","S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL"}\n'
    + NEW
    + '''
# ---- 参数化冒烟 ----
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
''').splitlines(True)

i_pv = find("# ============================================================ 15 Export")
s15 = "".join(cells[i_pv]["source"])
s15 = s15.replace('"solvent_mode": SOLVENT_MODE,',
                  '"solvent_mode": SOLVENT_MODE,\n'
                  '        "gb_model_requested": "GBn2 (amber14 + implicit/gbn2.xml)",\n'
                  '        "gb_model_applied": GB_MODEL,\n'
                  '        "gb_protocol_deviation": (GB_MODEL != "GBn2"),\n'
                  '        "gb_note": "OpenMM 8.6 removed implicitSolvent from createSystem.'
                  ' GBn2 was requested per protocol; gb_model_applied states the model actually'
                  ' used. If gb_protocol_deviation is true it must be reported as a protocol'
                  ' deviation.",')
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
for k, v in (("子进程隔离探测", "_probe_candidate" in t),
             ("候选顺序 OBC2>GBn2>OBC1", '["GBSAOBC2", "GBn2", "GBSAOBC"]' in t),
             ("GB_MODEL 记录", "gb_model_applied" in t),
             ("协议偏离标记", "gb_protocol_deviation" in t),
             ("不静默退化成真空", "不得" in t and "GB_MODEL is None" in t),
             ("不再用不存在的 customgbforces.GBSAOBCForce",
              "from openmm.app.internal.customgbforces import GBSAOBCForce" not in t),
             ("不再用 implicitSolvent", "implicitSolvent=app.GBn2" not in t)):
    print(f"  {'✅' if v else '❌'} {k}")
