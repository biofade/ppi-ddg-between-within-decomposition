# -*- coding: utf-8 -*-
"""scripts/v7_4/02k_stepwise_probe.py -- 用**渐进式温和探针**替换 cell 06.

## 为什么换方案

用户连续两次 `Server Error`，且**连 cell 06 的第一行 `GB_XML_PROBE` 都没打印**。
cell 06 的第一条语句是 `print("GB_XML_PROBE")`，它在 `fetch_pdb` 之后、
`subprocess` 之前 —— 也就是说：

* 要么死在 **cell 06 之前**（`fetch_pdb`？）
* 要么 `print` 的输出**根本没被刷新到 notebook**

结论：**不能再依赖「一次跑完整条链」**。必须先确定 kernel 到底死在哪一步。

## 方案：8 个**互相独立、逐级加码**的温和步骤

每一步：
* 独立 `try/except`，**失败不抛异常**（不中断后续步骤）
* **只打印，不做重活**
* `flush=True` 保证输出立刻可见
* 不依赖前一步成功

| step | 内容 | 验证 |
|---|---|---|
| 1 | 纯 stdlib 打印 | kernel 是否活着 |
| 2 | `import numpy` | 基础科学栈 |
| 3 | `import openmm` | OpenMM 可用性 + 版本 + platforms |
| 4 | `import pdbfixer` | PDBFixer 可导入 |
| 5 | 枚举 platform | CUDA 是否可见 |
| 6 | `fetch_pdb(VAL.iloc[0])` | **网络 + 文件系统**（上次的可疑点） |
| 7 | `PDBFixer(...)` 构造 | 解析 PDB |
| 8 | `ForceField(ff14SB, gbn2.xml)` | XML 能否装载（**不建 system**） |

**step 6/7/8 都不做重活**：不 `createSystem`、不 `Context`、不最小化、不 subprocess。
这样即使某步有问题，也只是打印一行 FAIL，**不会拖垮 kernel**。
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
# ============================================================ 06 STEPWISE DIAGNOSTIC PROBE
# 目的：确定 kernel 究竟死在哪一步。每一步独立 try/except、只打印、失败不中断。
# **不做重活**：不 createSystem、不 Context、不最小化、不 subprocess。
import sys as _sys, time as _time, traceback as _tb

def _step(n, label):
    print("\n[STEP %d] %s" % (n, label), flush=True)
    _sys.stdout.flush()

def _ok(msg=""):
    print("   -> OK  %s" % msg, flush=True); _sys.stdout.flush()

def _fail(e):
    print("   -> FAIL %s: %s" % (type(e).__name__, str(e)[:300]), flush=True)
    _sys.stdout.flush()

STEP_RESULTS = {}
print("=" * 68, flush=True)
print("STEPWISE PROBE START", flush=True)
print("=" * 68, flush=True)

# ---- 1 stdlib ----
_step(1, "python stdlib")
try:
    import platform as _pl
    _ok("python=%s platform=%s" % (_sys.version.split()[0], _pl.platform()))
    STEP_RESULTS[1] = "OK"
except Exception as e:
    _fail(e); STEP_RESULTS[1] = "FAIL"

# ---- 2 numpy ----
_step(2, "import numpy")
try:
    import numpy as _np
    _ok("numpy=%s" % _np.__version__)
    STEP_RESULTS[2] = "OK"
except Exception as e:
    _fail(e); STEP_RESULTS[2] = "FAIL"

# ---- 3 openmm ----
_step(3, "import openmm")
try:
    import openmm as _mm
    _ok("openmm=%s" % _mm.version.version)
    STEP_RESULTS[3] = "OK"
except Exception as e:
    _fail(e); STEP_RESULTS[3] = "FAIL"

# ---- 4 pdbfixer ----
_step(4, "import pdbfixer")
try:
    import pdbfixer as _pf
    _ok("pdbfixer imported")
    STEP_RESULTS[4] = "OK"
except Exception as e:
    _fail(e); STEP_RESULTS[4] = "FAIL"

# ---- 5 platforms ----
_step(5, "enumerate OpenMM platforms")
try:
    _names = []
    for _i in range(_mm.Platform.getNumPlatforms()):
        _names.append(_mm.Platform.getPlatform(_i).getName())
    _ok("platforms=%s" % _names)
    print("   CUDA available: %s" % ("CUDA" in _names), flush=True)
    STEP_RESULTS[5] = "OK"
except Exception as e:
    _fail(e); STEP_RESULTS[5] = "FAIL"

# ---- 6 fetch PDB（网络 + FS）----
_step(6, "fetch_pdb (network + filesystem)")
_pdb_path = None
try:
    _v0 = VAL.iloc[0]
    print("   target=%s" % _v0.target_id, flush=True)
    _t0 = _time.time()
    _pdb_path = fetch_pdb(_v0.target_id)
    _ok("%s  (%.1fs, %d bytes)" % (_pdb_path, _time.time() - _t0,
                                   _pdb_path.stat().st_size))
    STEP_RESULTS[6] = "OK"
except Exception as e:
    _fail(e); STEP_RESULTS[6] = "FAIL"

# ---- 7 PDBFixer 构造（只解析，不补原子）----
_step(7, "PDBFixer construct (parse only, no addMissingAtoms)")
_fx = None
try:
    if _pdb_path is None:
        print("   skipped (step 6 failed)", flush=True)
        STEP_RESULTS[7] = "SKIP"
    else:
        _t0 = _time.time()
        _fx = _pf.PDBFixer(filename=str(_pdb_path))
        print("   constructed (%.1fs)  atoms=%d chains=%d"
              % (_time.time() - _t0, _fx.topology.getNumAtoms(),
                 _fx.topology.getNumChains()), flush=True)
        _ok()
        STEP_RESULTS[7] = "OK"
except Exception as e:
    _fail(e); STEP_RESULTS[7] = "FAIL"

# ---- 8 只加载 ForceField（不 createSystem）----
_step(8, "ForceField(ff14SB, implicit/gbn2.xml) load only")
try:
    import openmm.app as _app
    _t0 = _time.time()
    _ff = _app.ForceField("amber14/protein.ff14SB.xml", "implicit/gbn2.xml")
    print("   loaded (%.1fs)" % (_time.time() - _t0), flush=True)
    try:
        _gens = [type(g).__name__ for g in _ff.getGenerators()]
        print("   generators=%s" % _gens, flush=True)
    except Exception as _e2:
        print("   generators unreadable: %s" % _e2, flush=True)
    _ok()
    STEP_RESULTS[8] = "OK"
except Exception as e:
    _fail(e); STEP_RESULTS[8] = "FAIL"

print()
print("=" * 68, flush=True)
print("STEPWISE PROBE SUMMARY: %s" % STEP_RESULTS, flush=True)
print("=" * 68, flush=True)
_sys.stdout.flush()

# 关键判据：kernel 能活到哪一步
if STEP_RESULTS.get(3) != "OK":
    raise RuntimeError("OpenMM 不可用 -> 无法继续（见 STEP 3）")
if STEP_RESULTS.get(8) != "OK":
    print()
    print("NOTE: STEP 8 失败 = 官方 XML 无法装载 -> 与 GB 力是否存在无关，先把装载问题解决。",
          flush=True)
print()
print(">>> 探测结束。请把上面整块 SUMMARY 贴回。", flush=True)
print(">>> 本单元**未**建 system / 未建 Context / 未最小化 / 未用 subprocess。", flush=True)
'''

i06 = find("# ============================================================ 06")
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
for k, v in (("8 个分步探针", "STEPWISE PROBE START" in t),
             ("每步独立 try/except", t.count("except Exception as e:") > 5),
             ("无 subprocess", "subprocess" not in t.split("06 STEPWISE")[1][:4000]),
             ("不 createSystem", "createSystem" not in t.split("06 STEPWISE")[1][:4000]),
             ("不建 Context", "mm.Context" not in t.split("06 STEPWISE")[1][:4000]),
             ("SUMMARY 打印", "STEPWISE PROBE SUMMARY" in t),
             ("第1步即 flush 打印", 'print("=" * 68, flush=True)' in t)):
    print(f"  {'✅' if v else '❌'} {k}")
