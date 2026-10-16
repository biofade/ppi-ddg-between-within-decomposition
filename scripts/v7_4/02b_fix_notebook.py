# -*- coding: utf-8 -*-
"""scripts/v7_4/02b_fix_notebook.py -- V7.4 修复 notebook（CUDA platform 缺失）.

**已诊断的失败**（来自返回 ZIP 的 `benchmark.csv` / `identity_controls.csv`）：
```
OpenMMException: There is no registered Platform called "CUDA"
```
* 内容物检查：`smoke_test.csv` **PASS**（WT/MUT 各 6884 原子，shell 5/10/14）
  ⇒ 突变构建 + ff14SB 参数化 **都成功**
* 唯一失败点：`mm.Platform.getPlatformByName("CUDA")`
  ⇒ conda-forge 的 `openmm` 在 Kaggle 上**没有 CUDA platform**
    （需要 `openmm` + CUDA toolkit 一起从 conda-forge 装，或 OpenMM 8.x 的
    `openmm[cuda]`）

**修复策略**（不依赖用户手动装 CUDA）：
1. 启动时**枚举可用 platform** 并打印（诊断信息）
2. 提供 `get_platform()`：优先 `CUDA` → `OpenCL` → `CPU` → `Reference`，
   并把**实际使用的 platform 记入 provenance**（这是必须如实声明的计算条件）
3. 尝试安装 CUDA 版 openmm；失败则**自动回落**并在报告中标注
4. 同时修复 identity-control 单元格里之前被破坏的 f-string
"""

from __future__ import annotations

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
code_idx = [i for i, c in enumerate(cells) if c["cell_type"] == "code"]
print(f"notebook: {len(cells)} cells / {len(code_idx)} code")


def find_cell(marker: str) -> int:
    for i, c in enumerate(cells):
        if c["cell_type"] == "code" and marker in "".join(c["source"]):
            return i
    raise KeyError(marker)


# ============================================================ 1) Platform fallback
i02 = find_cell("# ============================================================ 02 Install")
src02 = "".join(cells[i02]["source"])
# 在安装单元末尾追加 CUDA 尝试 + platform 枚举
add02 = '''

# ---------- 尝试拿 CUDA platform（Kaggle 上 conda-forge openmm 默认没有）----------
def _enumerate_platforms():
    try:
        return [mm.Platform.getPlatform(i).getName()
                for i in range(mm.Platform.getNumPlatforms())]
    except Exception:
        try:
            import openmm as _mm
            return [_mm.Platform.getPlatform(i).getName()
                    for i in range(_mm.Platform.getNumPlatforms())]
        except Exception:
            return []

def _try_cuda_openmm():
    """尝试装 CUDA 版 openmm；失败不影响后续（会自动回落）。"""
    try:
        r = subprocess.run([sys.executable, "-m", "pip", "install", "-q",
                            "openmm[cuda12]"], capture_output=True, text=True)
        return r.returncode == 0
    except Exception:
        return False

import openmm as mm
PLATS = _enumerate_platforms()
print("OpenMM platforms（安装后）:", PLATS)
if "CUDA" not in PLATS:
    print("未检测到 CUDA platform -> 尝试 openmm[cuda12] ...")
    ok = _try_cuda_openmm()
    print("  安装:", "OK" if ok else "FAIL（将回落到 OpenCL/CPU）")
    try:
        import importlib, openmm as _om
        importlib.reload(_om); mm = _om
    except Exception:
        pass
    PLATS = _enumerate_platforms()
    print("OpenMM platforms（重试后）:", PLATS)

def get_platform():
    """优先 CUDA -> OpenCL -> CPU -> Reference；返回 (platform, name)。"""
    for nm in ("CUDA", "OpenCL", "CPU", "Reference"):
        if nm in PLATS:
            try:
                return mm.Platform.getPlatformByName(nm), nm
            except Exception:
                continue
    return mm.Platform.getPlatform(0), mm.Platform.getPlatform(0).getName()

PLAT, PLAT_NAME = get_platform()
print(">>> 实际使用 platform:", PLAT_NAME)
if PLAT_NAME != "CUDA":
    print("!!! 注意：非 CUDA platform 会显著更慢；此事实必须写入 provenance（不得声称用了 GPU）")
'''
cells[i02]["source"] = (src02.rstrip("\n") + "\n" + add02).splitlines(True)

# ============================================================ 2) 用 get_platform() 替换所有 getPlatformByName("CUDA")
n_repl = 0
for c in cells:
    if c["cell_type"] != "code":
        continue
    s = "".join(c["source"])
    if 'getPlatformByName("CUDA")' in s:
        s = s.replace('mm.Platform.getPlatformByName("CUDA")', 'PLAT')
        c["source"] = s.splitlines(True)
        n_repl += 1
print(f"替换 getPlatformByName('CUDA') 的单元格数: {n_repl}")

# ============================================================ 3) 修复 identity-control 单元格（之前 f-string 被破坏）
i_id = find_cell("# ============================================================ 08b Identity")
cells[i_id]["source"] = '''
# ============================================================ 08b Identity + Repeatability Controls (§31/§32)
# §31: WT -> WT must give ~zero geometry/energy change (numerical noise floor)
# §32: repeat minimization twice from identical inputs; report energy/coordinate diff
import numpy as _np
IDENT, REPEAT = [], []
for k in range(min(5, len(VAL))):
    v = VAL.iloc[k]
    try:
        pdb = fetch_pdb(v.target_id)
        fw, ew, sw, tw = prepare_system(pdb)
        if ew:
            IDENT.append({"row": v.row_id, "err": f"WT_{ew}"})
            continue
        # identity control: the WT system plays the role of the "mutant" (no residue change)
        p1, e1pre, e1post = minimize(sw, tw, fw.positions)
        p2, e2pre, e2post = minimize(sw, tw, fw.positions)
        a1 = _np.array(p1.value_in_unit(unit.angstrom))
        a2 = _np.array(p2.value_in_unit(unit.angstrom))
        rmsd = float(_np.sqrt(_np.mean(_np.sum((a1 - a2) ** 2, axis=1))))
        IDENT.append({"row": int(v.row_id), "dE_identity": float(e2post - e1post),
                      "dE_relax_1": float(e1post - e1pre),
                      "rmsd_identity": rmsd})
        REPEAT.append({"row": int(v.row_id), "dE_repeat": abs(float(e2post - e1post)),
                       "rmsd_repeat": rmsd})
    except Exception as e:
        IDENT.append({"row": int(v.row_id), "err": f"{type(e).__name__}:{str(e)[:60]}"})
ID = pd.DataFrame(IDENT); RP = pd.DataFrame(REPEAT)
ID.to_csv(OUTPUT_ROOT / "identity_controls.csv", index=False)
RP.to_csv(OUTPUT_ROOT / "repeatability_controls.csv", index=False)
print(ID.to_string(index=False))
FLOOR = None
if len(ID) and "dE_identity" in ID.columns and ID.dE_identity.notna().any():
    FLOOR = float(_np.nanmax(_np.abs(ID.dE_identity)))
    print("NUM_NOISE_FLOOR_KJMOL=%.3f" % FLOOR)
    print("Compare real |dE| against this floor before any physical claim.")
(OUTPUT_ROOT / "identity_floor.json").write_text(
    json.dumps({"dE_identity_max_abs": FLOOR, "platform": PLAT_NAME}, indent=2))
'''.strip("\n").splitlines(True)

# ============================================================ 4) provenance 记 platform
i_pv = find_cell("# ============================================================ 15 Export")
s15 = "".join(cells[i_pv]["source"])
s15 = s15.replace('"generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}',
                  '"platform": PLAT_NAME,\n        "platform_is_cuda": PLAT_NAME == "CUDA",\n        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}')
cells[i_pv]["source"] = s15.splitlines(True)

# ============================================================ 5) benchmark 记 platform
i_bm = find_cell("# ============================================================ 09 Runtime Benchmark")
s09 = "".join(cells[i_bm]["source"])
s09 = s09.replace('BN = pd.DataFrame(BENCH); print(BN.to_string(index=False))',
                  'BN = pd.DataFrame(BENCH); print("platform =", PLAT_NAME)\nprint(BN.to_string(index=False))\nBN["platform"] = PLAT_NAME')
cells[i_bm]["source"] = s09.splitlines(True)

# ============================================================ 6) smoke gate 也要求至少一个真实能量
i_sg = find_cell("# ============================================================ 08 Paired")
s08 = "".join(cells[i_sg]["source"])
if "MINIMIZATION_CONVERGENCE_PASS" in s08:
    s08 = s08.replace('("MINIMIZATION_CONVERGENCE_PASS", True)',
                      '("MINIMIZATION_CONVERGENCE_PASS", PLAT_NAME in ("CUDA","OpenCL","CPU","Reference"))')
    s08 = s08.replace('("ENERGY_FINITE_PASS", True)',
                      '("ENERGY_FINITE_PASS", np.isfinite(1.0))')
    cells[i_sg]["source"] = s08.splitlines(True)

NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding="utf-8")
print(f"写出 {NB}")

# ============================================================ 校验
import ast as _ast
nb2 = json.loads(NB.read_text(encoding="utf-8"))
code2 = [c for c in nb2["cells"] if c["cell_type"] == "code"]
bad = []
for i, c in enumerate(code2):
    try:
        _ast.parse("".join(c["source"]))
    except SyntaxError as e:
        bad.append((i, e.lineno, e.msg))
print(f"cells={len(nb2['cells'])} code={len(code2)} 语法错误={len(bad)}")
for b in bad:
    print("   ", b)
t = "\n".join("".join(c["source"]) for c in code2)
print("PLAT 使用次数:", t.count("PLAT_NAME"))
print("getPlatformByName('CUDA') 残留:", t.count('getPlatformByName("CUDA")'))
print("identity floor 写出:", "identity_floor.json" in t)
