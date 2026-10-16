# -*- coding: utf-8 -*-
"""scripts/v7_4/02i_fix_dupes.py -- 修掉 notebook 里的两处真实缺陷.

1. **重复定义** `residues_within`：cell 07（nm 版，正确）被 cell 08（Å 版）覆盖。
   功能上同义（都返回 8 Å 内的残基键），但**同一名字两个实现**是隐患
   —— cell 08b/10 调用时用哪一个取决于执行顺序。删掉 cell 08 的那份。
2. **cell 12** 用旧的 3 参签名调 `minimize(sw, tw, fw.positions)`，
   而新签名是 `minimize(system, topology, positions, site_keys, shell_keys)`。
   且原实现只在原子数相等时才有意义（WT/MUT 原子数常不同）。
   改为**真正**提取 §25 要求的几何增量。
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


# ---------- 1) 删除 cell 08 里重复的 residues_within ----------
i08 = find("# ============================================================ 08 Paired WT/MUT Smoke")
s08 = "".join(cells[i08]["source"])
DUP_START = "def residues_within(topology, positions, site_keys, radius):"
if DUP_START in s08:
    a = s08.index(DUP_START)
    # 找到函数结束（下一个顶层 def/print/SMOKE）
    b = s08.index("SMOKE = []", a)
    s08 = s08[:a] + "# （residues_within 已在 cell 07 定义，此处不再重复）\n\n" + s08[b:]
    cells[i08]["source"] = s08.splitlines(True)
    print("已删除 cell 08 中重复的 residues_within")
else:
    print("cell 08 无重复定义（已处理过）")

# ---------- 2) 重写 cell 12：真正的几何增量 ----------
i12 = find("# ============================================================ 12 Relaxed Geometry")
cells[i12]["source"] = '''
# ============================================================ 12 Relaxed Geometry Extraction
# §25：在**松弛后**的 WT/MUT 上重算局部几何，取 mutation-centered 增量。
# 修复：原版用旧签名调 minimize 且只在原子数相等时可用；现改为配对提取。
import numpy as _np

def _site_geom(topology, positions, site_keys, radius=RADIUS_MOVABLE):
    """以突变位点为中心，算 §25 要求的局部几何量。"""
    pos = _np.array([list(p) for p in positions.value_in_unit(unit.angstrom)])
    site_atoms = [a.index for a in topology.atoms()
                  if (a.residue.chain.id, a.residue.id) in site_keys]
    other = [a for a in topology.atoms()
             if (a.residue.chain.id, a.residue.id) not in site_keys]
    if not site_atoms or not other:
        return None
    oi = _np.array([a.index for a in other])
    D = _np.sqrt(_np.sum((pos[oi][:, None, :] - pos[site_atoms][None, :, :]) ** 2,
                         axis=2))            # (n_other, n_site)
    dmin_per_other = D.min(axis=1)
    # 伙伴残基层面的最小距离
    res_min = {}
    for k, a in enumerate(other):
        key = (a.residue.chain.id, a.residue.id)
        v = float(dmin_per_other[k])
        if key not in res_min or v < res_min[key]:
            res_min[key] = v
    return {"n_contact_4": int((dmin_per_other <= 4.0).sum()),
            "n_contact_6": int((dmin_per_other <= 6.0).sum()),
            "n_contact_8": int((dmin_per_other <= 8.0).sum()),
            "min_partner_distance": float(dmin_per_other.min()),
            "mean_partner_distance": float(dmin_per_other.mean())}

def geometry_deltas(rid, target, muts):
    """配对 WT/MUT 松弛后，返回 §25 的几何增量。"""
    pdb = fetch_pdb(target)
    sk = set()
    for (ch, rs, _w, _m) in muts:
        sk.add((str(ch), str(rs)))
    if not sk and muts:
        sk = {(str(muts[0][0]), str(muts[0][1]))}
    # WT
    fw, ew, sw, tw = prepare_system(pdb)
    if ew:
        return {"row_id": rid, "error": f"WT_{ew}"}
    shw = residues_within(tw, fw.positions, sk, RADIUS_MOVABLE)
    pw, wpre, wpost, _ = minimize(sw, tw, fw.positions, sk, shw)
    gw = _site_geom(tw, pw, sk)
    # MUT
    fm, em, sm, tm = prepare_system(pdb, muts)
    if em:
        return {"row_id": rid, "error": f"MUT_{em}"}
    shm = residues_within(tm, fm.positions, sk, RADIUS_MOVABLE)
    pm, mpre, mpost, _ = minimize(sm, tm, fm.positions, sk, shm)
    gm = _site_geom(tm, pm, sk)
    out = {"row_id": rid, "target_id": target,
           "dE_total": float(mpost - wpost),
           "dE_relax_cost": float((mpost - mpre) - (wpost - wpre)),
           "shell_atoms_wt": len(shw), "shell_atoms_mut": len(shm)}
    if gw and gm:
        for k in ("n_contact_4", "n_contact_6", "n_contact_8",
                  "min_partner_distance", "mean_partner_distance"):
            out["delta_" + k] = float(gm[k] - gw[k])
    return out

print("几何增量函数已定义（§25：delta_min_partner_distance / contact counts / shell）。")
print("SASA 相关项（delta_site_SASA / delta_interface_burial）在本单元未实现 -> 标记 NOT_RUN。")
'''.strip("\n").splitlines(True)

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
print("  residues_within 定义次数:", len([1 for c in code2
                                            if "def residues_within" in "".join(c["source"])]))
print("  minimize 旧签名调用残留:",
      len([1 for m in __import__("re").finditer(r"minimize\((sw|sm), (tw|tm), (fw|fm)\.positions\)", t)]))
print("  cell 12 有 _site_geom:", "_site_geom" in t)
print("  delta_min_partner_distance:", "delta_min_partner_distance" in t or
      '"delta_" + k' in t)
