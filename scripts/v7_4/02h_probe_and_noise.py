# -*- coding: utf-8 -*-
"""scripts/v7_4/02h_probe_and_noise.py -- 补足 cell 06 的存在性检查 + cell 08b 的噪声统计.

按用户第 4 轮要求（**不改物理协议、不跑全量**）：

## cell 06 追加（防止「协议错了但程序跑通」）
* `CustomGBForce count`（>= 1）
* `GB force getNumParticles()` == `System getNumParticles()`
* `GB force getForceGroup()`
* 每个力的 `getForceGroup()`（后面要分离 restraint 与 biological/GB 能量）
* `E_pre_min` / `E_post_min` / `DeltaE_min`（最小化不应明显升高能量）

## cell 08b 追加（第二道门 = 数值噪声地板）
* `WT→WT ΔE`：median / MAD / P95 / max
* `WT→WT RMSD`：median / P95 / max
* repeatability：ΔE 与 feature 差的 median / max
* **比值**（真正要回答的问题）：
  `R_E = Q95(|ΔE|_WT→WT) / median(|ΔE|_WT→Mut)`
  以及 geometry 对应比值
  > 不预设绝对阈值；`R_E` 与 mutation signal 同量级 ⇒ `NUMERICAL_NOISE_TOO_HIGH` ⇒ STOP
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


# ============================================================ cell 06：子进程探针加检查
i06 = find("# ============================================================ 06 Official-XML GBn2 probe")
s06 = "".join(cells[i06]["source"])

# 1) 力清单 + 力组
s06 = s06.replace(
    '''names = [system.getForce(i).getName() for i in range(system.getNumForces())]
out["forces"] = names
out["gb_forces"] = [n for n in names if "GB" in n.upper() or "BORN" in n.upper()]''',
    '''names = [system.getForce(i).getName() for i in range(system.getNumForces())]
out["forces"] = names
out["force_groups"] = {names[i]: system.getForce(i).getForceGroup()
                       for i in range(system.getNumForces())}
out["gb_forces"] = [n for n in names if "GB" in n.upper() or "BORN" in n.upper()]''')

# 2) CustomGBForce count + GB particles == system particles + group
s06 = s06.replace(
    '''try:
    out["customgb"] = bool(any(isinstance(system.getForce(i), mm.CustomGBForce)
                               for i in range(system.getNumForces())))
except Exception:
    out["customgb"] = False
out["particles"] = system.getNumParticles()''',
    '''gb_idx = []
try:
    for i in range(system.getNumForces()):
        if isinstance(system.getForce(i), mm.CustomGBForce):
            gb_idx.append(i)
except Exception:
    gb_idx = []
out["customgb_count"] = len(gb_idx)
out["customgb"] = bool(len(gb_idx) >= 1)
out["gb_particles"] = (system.getForce(gb_idx[0]).getNumParticles()
                       if gb_idx else None)
out["gb_force_group"] = (system.getForce(gb_idx[0]).getForceGroup()
                         if gb_idx else None)
out["particles"] = system.getNumParticles()
# 硬门控：GB 力必须覆盖全部粒子（防「只有部分粒子」）
out["gb_covers_all"] = bool(gb_idx and
                            out["gb_particles"] == out["particles"])''')

# 3) E_pre / E_post / dE
s06 = s06.replace(
    '''    mm.LocalEnergyMinimizer.minimize(ctx, tolerance=10.0, maxIterations=2000)
    out["minimized"] = True
    out["energy_post"] = float(ctx.getState(getEnergy=True).getPotentialEnergy()
                              .value_in_unit(unit.kilojoule_per_mole))''',
    '''    out["E_pre_min"] = float(out["energy"])
    mm.LocalEnergyMinimizer.minimize(ctx, tolerance=10.0, maxIterations=2000)
    out["minimized"] = True
    out["E_post_min"] = float(ctx.getState(getEnergy=True).getPotentialEnergy()
                              .value_in_unit(unit.kilojoule_per_mole))
    out["DeltaE_min"] = out["E_post_min"] - out["E_pre_min"]''')

# 4) 打印块
s06 = s06.replace(
    '''    print("  CustomGBForce present: %s" % INFO.get("customgb"))
    print("  Particle count: %s" % INFO.get("particles"))
    print("  Energy: %s kJ/mol" % INFO.get("energy"))
    print("  Force finite: %s" % INFO.get("forces_finite"))
    print("  Context creation: %s" % INFO.get("context"))
    print("  Minimization: %s" % INFO.get("minimized"))''',
    '''    print("  CustomGBForce present: %s" % INFO.get("customgb"))
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
    print("  DeltaE_min: %s kJ/mol" % INFO.get("DeltaE_min"))''')

# 5) 硬门控加 gb_covers_all 与 DeltaE 不显著上升
s06 = s06.replace(
    '''GB_XML_OK = bool(
    _rc == 0
    and INFO.get("customgb") is True
    and INFO.get("energy_finite") is True
    and INFO.get("forces_finite") is True
    and INFO.get("context") is True
    and INFO.get("minimized") is True
)''',
    '''_dE = INFO.get("DeltaE_min")
GB_XML_OK = bool(
    _rc == 0
    and INFO.get("customgb") is True                 # GB 力存在
    and INFO.get("customgb_count", 0) >= 1           # >= 1 个 CustomGBForce
    and INFO.get("gb_covers_all") is True            # GB 粒子数 == system 粒子数
    and INFO.get("energy_finite") is True
    and INFO.get("forces_finite") is True
    and INFO.get("context") is True
    and INFO.get("minimized") is True
    and (_dE is not None)                            # 有 dE
    # 最小化不应明显升高能量（不设狭窄绝对区间；只要求非显著上升）
    and (_dE <= 0.0 if (_dE is not None and INFO.get("energy", 0.0) > 1e6)
         else (_dE <= 0.05 * abs(INFO.get("energy", 1.0)) if _dE is not None
               else False))
)''')
cells[i06]["source"] = s06.splitlines(True)

# ============================================================ cell 08b：完整噪声统计
i08b = find("# ============================================================ 08b Identity")
cells[i08b]["source"] = '''
# ============================================================ 08b Identity + Repeatability + NOISE FLOOR
# §31 WT->WT identity control + §32 repeatability
# 第二道门：数值噪声地板。**不预设绝对阈值**，而是与真实 mutation 信号比。
import numpy as _np

def _stats(a):
    a = _np.asarray([x for x in a if x is not None and _np.isfinite(x)], float)
    if len(a) == 0:
        return {"n": 0}
    med = float(_np.median(a))
    return {"n": int(len(a)), "median": med,
            "MAD": float(_np.median(_np.abs(a - med))),
            "P95": float(_np.percentile(_np.abs(a), 95)),
            "max": float(_np.max(_np.abs(a)))}

THREE_ = THREE
IDENT, REPEAT = [], []
N_IDENT = min(5, len(VAL))
print("identity 对照样本数:", N_IDENT)
for k in range(N_IDENT):
    v = VAL.iloc[k]
    try:
        pdb = fetch_pdb(v.target_id)
        sk = {(str(v.chain), str(v.resseq))}
        # --- WT 体系，跑两次（identity + repeatability 同源）---
        fw, ew, sw, tw = prepare_system(pdb)
        if ew:
            IDENT.append({"row": int(v.row_id), "err": f"WT_{ew}"}); continue
        sh = residues_within(tw, fw.positions, sk, RADIUS_MOVABLE)
        p1, e1pre, e1post, _r1 = minimize(sw, tw, fw.positions, sk, sh)
        p2, e2pre, e2post, _r2 = minimize(sw, tw, fw.positions, sk, sh)
        a1 = _np.array(p1.value_in_unit(unit.angstrom))
        a2 = _np.array(p2.value_in_unit(unit.angstrom))
        rmsd = float(_np.sqrt(_np.mean(_np.sum((a1 - a2) ** 2, axis=1))))
        IDENT.append({"row": int(v.row_id), "target": str(v.target_id),
                      "dE_identity": float(e2post - e1post),
                      "rmsd_identity": rmsd,
                      "shell": len(sh)})
        REPEAT.append({"row": int(v.row_id), "dE_repeat": abs(float(e2post - e1post)),
                       "rmsd_repeat": rmsd})
        print("  row %-5d %-6s dE=%+.4f kJ/mol  RMSD=%.4f A  shell=%d"
              % (int(v.row_id), str(v.target_id), e2post - e1post, rmsd, len(sh)))
    except Exception as e:
        IDENT.append({"row": int(v.row_id), "err": f"{type(e).__name__}:{str(e)[:60]}"})
ID = pd.DataFrame(IDENT); RP = pd.DataFrame(REPEAT)
ID.to_csv(OUTPUT_ROOT / "identity_controls.csv", index=False)
RP.to_csv(OUTPUT_ROOT / "repeatability_controls.csv", index=False)

# ---------- 真实 mutation 的 |dE| 与 RMSD（配对）----------
print()
print("=== 真实 mutation 信号（配对 WT/MUT，label-free）===")
MUT = []
for k in range(min(20, len(VAL))):
    v = VAL.iloc[k]
    try:
        pdb = fetch_pdb(v.target_id)
        sk = {(str(v.chain), str(v.resseq))}
        fw, ew, sw, tw = prepare_system(pdb)
        if ew: continue
        shw = residues_within(tw, fw.positions, sk, RADIUS_MOVABLE)
        pw, wpre, wpost, _ = minimize(sw, tw, fw.positions, sk, shw)
        muts = [(str(v.chain), int(v.resseq), THREE_[str(v.wt)], THREE_[str(v.mut)])]
        fm, em, sm, tm = prepare_system(pdb, muts)
        if em: continue
        shm = residues_within(tm, fm.positions, sk, RADIUS_MOVABLE)
        pm, mpre, mpost, _ = minimize(sm, tm, fm.positions, sk, shm)
        # 把 MUT 坐标对齐前先按原子索引比较（同一拓扑顺序）
        try:
            aw = _np.array(pw.value_in_unit(unit.angstrom))
            am = _np.array(pm.value_in_unit(unit.angstrom))
            n = min(len(aw), len(am))
            rmsd_m = float(_np.sqrt(_np.mean(_np.sum((aw[:n] - am[:n]) ** 2, axis=1))))
        except Exception:
            rmsd_m = _np.nan
        MUT.append({"row": int(v.row_id), "target": str(v.target_id),
                    "wt": str(v.wt), "mut": str(v.mut),
                    "dE_total": float(mpost - wpost),
                    "dE_relax_cost": float((mpost - mpre) - (wpost - wpre)),
                    "rmsd_wt_mut": rmsd_m})
    except Exception as e:
        MUT.append({"row": int(v.row_id), "err": f"{type(e).__name__}:{str(e)[:60]}"})
MU = pd.DataFrame(MUT)
MU.to_csv(OUTPUT_ROOT / "mutation_smoke.csv", index=False)
if len(MU):
    print(MU[["row", "wt", "mut", "dE_total", "rmsd_wt_mut"]].to_string(index=False))

# ---------- 噪声地板统计 ----------
print()
print("=== NUM_NOISE_FLOOR / §31 §32 ===")
ID_ST = _stats(ID.get("dE_identity", []))
RS_ST = _stats(ID.get("rmsd_identity", []))
RP_ST = _stats(RP.get("dE_repeat", []))
MU_ST = _stats(MU.get("dE_total", [])) if len(MU) else {"n": 0}
MR_ST = _stats(MU.get("rmsd_wt_mut", [])) if len(MU) else {"n": 0}
print("WT->WT dE   :", ID_ST)
print("WT->WT RMSD :", RS_ST)
print("repeat dE   :", RP_ST)
print("WT->Mut dE  :", MU_ST)
print("WT->Mut RMSD:", MR_ST)

# 比值：真正要回答的问题
R_E = R_G = None
if ID_ST.get("n") and MU_ST.get("n"):
    den = abs(MU_ST["median"]) if abs(MU_ST["median"]) > 1e-12 else _np.nan
    R_E = ID_ST["P95"] / den if _np.isfinite(den) else None
if RS_ST.get("n") and MR_ST.get("n"):
    den2 = MR_ST["median"] if MR_ST["median"] > 1e-12 else _np.nan
    R_G = RS_ST["P95"] / den2 if _np.isfinite(den2) else None
print()
print("  R_E = Q95(|dE|_WT->WT) / median(|dE|_WT->Mut) = %s" % R_E)
print("  R_G = Q95(RMSD_WT->WT) / median(RMSD_WT->Mut) = %s" % R_G)

NUM_NOISE_FLOOR_KJMOL = ID_ST.get("P95")
print()
print("NUM_NOISE_FLOOR_KJMOL = %s" % NUM_NOISE_FLOOR_KJMOL)

# ---------- 第二道门（不预设绝对阈值，只比量级）----------
print()
print("=== NOISE GATE ===")
noise_ok = None
if R_E is None:
    print("  R_E 不可算（样本不足）-> 无法判定，STOP")
    noise_ok = False
else:
    print("  R_E = %.3f" % R_E)
    if R_E >= 1.0:
        print("  => NUMERICAL_NOISE_TOO_HIGH（identity 噪声 >= mutation 信号）")
        noise_ok = False
    elif R_E >= 0.5:
        print("  => MARGINAL（同量级；需人工判断是否继续）")
        noise_ok = False
    else:
        print("  => OK（identity 噪声明显低于 mutation 信号）")
        noise_ok = True
(OUTPUT_ROOT / "noise_gate.json").write_text(json.dumps(
    {"identity_dE": ID_ST, "identity_rmsd": RS_ST, "repeat_dE": RP_ST,
     "mutation_dE": MU_ST, "mutation_rmsd": MR_ST,
     "R_E": R_E, "R_G": R_G,
     "NUM_NOISE_FLOOR_KJMOL": NUM_NOISE_FLOOR_KJMOL,
     "noise_gate_pass": noise_ok,
     "gb_model": GB_MODEL, "platform": PLAT_NAME}, indent=2))

if not noise_ok:
    print()
    print("STOP: 数值噪声门未通过 -> 不得进入全量。")
    print("  identity 对照暴露的是**数值地板**问题，与物理模型无关。")
    raise RuntimeError("NUMERICAL_NOISE_GATE_FAIL —— 依第 4 轮要求 STOP")
print()
print("两道门均已通过：GBn2_XML (cell 06) + NOISE GATE (cell 08b)")
print("NOTE: 本单元只做可行性/噪声评估；全量在 cell 10（RUN_SMOKE_ONLY=False）。")
'''
cells[i08b]["source"] = ("".join(cells[i08b]["source"])).splitlines(True)

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
for k, v in (("CustomGBForce count", "customgb_count" in t),
             ("GB particles == system particles", "gb_covers_all" in t),
             ("GB force group", "gb_force_group" in t),
             ("每个力的 force group", "force_groups" in t),
             ("E_pre/E_post/DeltaE_min", "DeltaE_min" in t),
             ("identity dE 统计", "WT->WT dE" in t),
             ("R_E 比值", "R_E" in t),
             ("R_G 比值", "R_G" in t),
             ("NOISE GATE", "NOISE GATE" in t),
             ("噪声门失败 STOP", "NUMERICAL_NOISE_GATE_FAIL" in t),
             ("20 个 mutation smoke", "mutation_smoke.csv" in t),
             ("无 OBC2 回落", "GBSAOBC2Force" not in t)):
    print(f"  {'✅' if v else '❌'} {k}")
