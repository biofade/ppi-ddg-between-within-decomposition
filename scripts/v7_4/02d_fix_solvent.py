# -*- coding: utf-8 -*-
"""scripts/v7_4/02d_fix_solvent.py -- V7.4 修复隐式溶剂（OpenMM 8.6 API 差异）.

## 第 2 轮 Kaggle 失败（真实日志）

```text
WT 参数化: FAIL FORCEFIELD:ValueError:The argument 'implicitSolvent' was
specified to createSystem() but was never used.
```

`platform = CUDA` **已经修好**（`openmm[cuda12]` 安装成功）。
现在唯一阻塞点是 `ForceField.createSystem(..., implicitSolvent=app.GBn2)`
在 **OpenMM 8.6.1.dev** 上不被接受。

## 为什么不猜

`implicitSolvent` 在不同 OpenMM 版本中的接受方式不同：
* 旧版本：`createSystem(..., implicitSolvent=app.GBn2)` 直接可用
* 某些版本：只接受 `implicitSolvent=<GBSAOBCForce 模板>`，
  或要求溶剂参数**已由 ForceField 载入**并改用
  `soluteDielectric` / `solventDielectric`

⇒ 本脚本改成**运行时自省 + 多策略回落**，并把每次尝试的真实异常
  写进 notebook 输出，这样下一轮无论成功或失败都能看到确切原因。

策略顺序：
1. `implicitSolvent=app.GBn2`（旧版 API）
2. 不传 `implicitSolvent`，改传 `soluteDielectric` / `solventDielectric`
3. 构造 system **不含**溶剂，再手工 `addForce` 一个由 gbn2.xml
   模板 `createForce()` 生成的 GB 力
4. 全部失败 → 明确报 `FAIL_FORCEFIELD`（**不**静默退化成真空）
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


NEW_SOLVENT_BLOCK = '''
# ---------- 隐式溶剂：运行时自省 + 多策略回落（OpenMM 版本差异）----------
import inspect as _inspect

def _solvent_signature():
    """打印 createSystem 的真实签名，供下一轮诊断。"""
    try:
        sig = _inspect.signature(app.ForceField.createSystem)
        return str(sig)
    except Exception as e:
        return f"<signature unavailable: {e}>"

SIG = _solvent_signature()
print("createSystem 签名:", SIG)

def _ff_templates(ff):
    """列出 ForceField 里已载入的力模板（判断 gbn2.xml 是否真被载入）。"""
    names = []
    for attr in ("_forces", "_ffforces"):
        obj = getattr(ff, attr, None)
        if isinstance(obj, dict):
            names += list(obj.keys())
    return sorted(set(names))

def _make_gb_force(ff):
    """从已载入模板构造一个 GB 力（策略 3 用）。"""
    from openmm.app.internal.customgbforces import GBSAOBCForce as _G
    for key in ("GBSAOBCForce", "GBn2", "GBn2Force"):
        tpl = None
        for attr in ("_forces", "_ffforces"):
            d = getattr(ff, attr, None)
            if isinstance(d, dict) and key in d:
                tpl = d[key]; break
        if tpl is not None:
            try:
                return tpl.createForce()
            except Exception:
                continue
    return None

def make_system(fixer):
    """按策略 1->2->3 建 system；返回 (system, how, error)。
    **绝不**静默退化成真空（§6 要求 GBn2）。"""
    ff = app.ForceField(*FORCEFIELD)
    base = dict(nonbondedMethod=app.NoCutoff, constraints=None)
    attempts = []

    # 策略 1：旧版 API
    try:
        s = ff.createSystem(fixer.topology, implicitSolvent=app.GBn2, **base)
        return s, "implicitSolvent=app.GBn2", None
    except Exception as e:
        attempts.append(("implicitSolvent=app.GBn2", f"{type(e).__name__}:{e}"))

    # 策略 2：介电常数参数
    try:
        s = ff.createSystem(fixer.topology, soluteDielectric=1.0,
                            solventDielectric=78.5, **base)
        gb = _make_gb_force(ff)
        if gb is not None:
            s.addForce(gb)
            return s, "solute/solventDielectric + manual GB force", None
        attempts.append(("dielectric-only", "no GB template found"))
    except Exception as e:
        attempts.append(("dielectric-only", f"{type(e).__name__}:{e}"))

    # 策略 3：不含溶剂建 system，再挂 GB 力
    try:
        s = ff.createSystem(fixer.topology, **base)
        gb = _make_gb_force(ff)
        if gb is not None:
            s.addForce(gb)
            return s, "manual GB (ff template)", None
        attempts.append(("manual-GB", "no GB template found"))
    except Exception as e:
        attempts.append(("manual-GB", f"{type(e).__name__}:{e}"))

    return None, None, attempts

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
print("gbn2 模板是否载入:", _ff_templates(app.ForceField(*FORCEFIELD))[:20])
'''

# ------- 替换 cell 06 的 prepare_system -------
i06 = find("# ============================================================ 06 OpenMM Parameterization")
cells[i06]["source"] = (
    "# ============================================================ 06 OpenMM Parameterization Smoke Test\n"
    "# 修复：GBn2 的 OpenMM 版本差异 -> 自省 + 多策略回落（不静默退化成真空）\n"
    "from pdbfixer import PDBFixer\n"
    'THREE = {"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU",\n'
    '         "G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE",\n'
    '         "P":"PRO","S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL"}\n'
    + NEW_SOLVENT_BLOCK
    + '''
# ---- 参数化冒烟 ----
v0 = VAL.iloc[0]
pdb0 = fetch_pdb(v0.target_id)
fx, err, sysm, topo = prepare_system(pdb0)
FF_OK = err is None
print("WT 参数化:", "PASS" if FF_OK else "FAIL")
if not FF_OK:
    print("  error:", err)
if FF_OK:
    print("  溶剂策略:", SOLVENT_MODE)
    print("  原子数:", sysm.getNumParticles(), "| 力数:", sysm.getNumForces())
    print("  力列表:", [sysm.getForce(i).getName() for i in range(sysm.getNumForces())])
    muts = [(str(v0.chain), int(v0.resseq), THREE[str(v0.wt)], THREE[str(v0.mut)])]
    fx2, err2, sys2, topo2 = prepare_system(pdb0, muts)
    print("MUT 参数化:", "PASS" if err2 is None else f"FAIL {err2}")
    if err2 is None:
        print("  MUT 原子数:", sys2.getNumParticles())
''').splitlines(True)

# ------- provenance 记溶剂策略 -------
i_pv = find("# ============================================================ 15 Export")
s15 = "".join(cells[i_pv]["source"])
s15 = s15.replace('"platform": PLAT_NAME,',
                  '"platform": PLAT_NAME,\n        "solvent_mode": SOLVENT_MODE,')
cells[i_pv]["source"] = s15.splitlines(True)

# ------- smoke gate：加溶剂策略门 -------
i_sg = find("# ============================================================ 08 Paired")
s08 = "".join(cells[i_sg]["source"])
s08 = s08.replace('("FORCEFIELD_PARAMETERIZATION_PASS", FF_OK),',
                  '("FORCEFIELD_PARAMETERIZATION_PASS", FF_OK),\n'
                  '             ("IMPLICIT_SOLVENT_APPLIED", bool(SOLVENT_MODE)),')
cells[i_sg]["source"] = s08.splitlines(True)

NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding="utf-8")

# ---- 校验 ----
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
for k, v in (("自省 createSystem 签名", "_solvent_signature" in t),
             ("策略1 implicitSolvent=app.GBn2", "implicitSolvent=app.GBn2" in t),
             ("策略2 介电常数", "soluteDielectric" in t),
             ("策略3 手工 GB force", "_make_gb_force" in t),
             ("不静默退化成真空", "绝不" in t and "SOLVENT_MODE" in t),
             ("gate 加 IMPLICIT_SOLVENT_APPLIED", "IMPLICIT_SOLVENT_APPLIED" in t),
             ("provenance 记 solvent_mode", '"solvent_mode": SOLVENT_MODE' in t)):
    print(f"  {'✅' if v else '❌'} {k}")
