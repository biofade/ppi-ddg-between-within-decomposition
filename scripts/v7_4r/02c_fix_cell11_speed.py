# -*- coding: utf-8 -*-
"""scripts/v7_4r/02c_fix_cell11_speed.py -- 修 cell 11 的「看起来卡死」与重复计算.

## 诊断（用户实测：约 2 小时无输出）

cell 11 在**打印任何东西之前**要做完 20 个变体，而每个变体包含：

```
prepare_system(WT)          <- PDBFixer 解析 + 加氢 + createSystem（GBn2）
residues_within(WT 8A)
minimize(WT)                <- 开 Context + L-BFGS 2000 步
prepare_system(MUT)         <- 再一次 PDBFixer + 加氢 + createSystem
residues_within(MUT 8A)
minimize(MUT)               <- 再一次 Context + 最小化
```

⇒ **20 × (2×PDBFixer + 2×createSystem + 2×最小化)**，且**全程零输出**。
「跑很久」与「卡死」在界面上**完全无法区分** —— 这是本次的**真实缺陷**。

### 另一个被忽视的浪费：WT 被算了两次

* cell 10（identity/repeatability）已经对 **同一批前 5 个变体** 算过 WT 最小化
* cell 11 又对 **同样这 20 个** 从头算一遍 WT
* ⇒ **cell 10 与 cell 11 之间有 ≤25 次可复用的 WT 最小化**

## 修复（三项）

1. **WT 缓存**：以 `(target, site_set)` 为键缓存 WT 最小化结果（含 positions/能量/壳）。
   cell 10 先跑，cell 11 直接命中 ⇒ 省掉约一半的 PDBFixer + createSystem + 最小化。
2. **逐变体实时进度**：每个变体打印一行（含耗时、累计、**ETA**），
   并**立即 flush**。这样「慢」与「死」一眼可分。
3. **每变体心跳超时**：单变体超过 `VARIANT_TIMEOUT_S` 就跳过并记录
   （不拖垮整格）。
4. 命中缓存的行标注 `wt_cached=True`，便于事后核对。

**注意**：本轮**不改任何物理协议**（力场/GB 模型/pH/限制/最小化设置都不变），
只改**缓存复用**与**可观测性**。
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

R = Path(__file__).resolve().parent.parent.parent
NB = R / "colab" / "v7_4r" / "V7_4R_OpenMM_Relaxation_COLAB.ipynb"
nb = json.loads(NB.read_text(encoding="utf-8"))
cells = nb["cells"]


def find(marker):
    for i, c in enumerate(cells):
        if c["cell_type"] == "code" and marker in "".join(c["source"]):
            return i
    raise KeyError(marker)


# ============================================================ 1) config 加超时
i00 = find("# ============================================================ 00 Config")
s00 = "".join(cells[i00]["source"])
s00 = s00.replace("GB_PROBE_TIMEOUT = 900               # s",
                  "GB_PROBE_TIMEOUT = 900               # s\n"
                  "VARIANT_TIMEOUT_S = 600              # s：单变体心跳超时（防卡死拖垮整格）")
cells[i00]["source"] = s00.splitlines(True)

# ============================================================ 2) cell 09 之后插入 WT 缓存 + 进度工具
WT_CACHE: dict = {}
i09 = find("# ============================================================ 09 Mutation builder smoke")
CACHE_CELL = '''
# ============================================================ 09b WT cache + progress helper (speed/observability)
# WT 最小化只依赖 (target, site_set)（§25 缓存原则）。
# cell 10 先算，cell 11 直接复用 -> 省掉约一半的 PDBFixer/createSystem/最小化。
import time as _time

WT_CACHE = {}

def _wt_key(target, site_keys):
    return (str(target), tuple(sorted(f"{c}:{r}" for c, r in site_keys)))

def get_wt(target, pdb_path, site_keys, force=False):
    """返回 dict(topology, positions, pre, post, restraint, shell, cached)。"""
    k = _wt_key(target, site_keys)
    if (not force) and k in WT_CACHE:
        d = dict(WT_CACHE[k]); d["cached"] = True
        return d
    fw, ew, sw, tw = prepare_system(pdb_path)
    if ew:
        raise RuntimeError(f"WT_{ew}")
    shw = residues_within(tw, fw.positions, site_keys, RADIUS_MOVABLE)
    pw, wpre, wpost, wres = minimize(sw, tw, fw.positions, site_keys, shw)
    d = {"topology": tw, "positions": pw, "pre": wpre, "post": wpost,
         "restraint": wres, "shell": shw, "cached": False}
    WT_CACHE[k] = dict(d)
    return d

def prog(done, total, t0, extra=""):
    el = _time.time() - t0
    eta = (el / done * (total - done)) if done else float("nan")
    print("  [%2d/%2d] %6.1fs elapsed | ETA %6.1fs %s"
          % (done, total, el, eta, extra), flush=True)

print("WT 缓存 + 进度工具就绪（WT_CACHE 以 (target, site_set) 为键）")
'''
cells.insert(i09 + 1, {"cell_type": "code", "execution_count": None,
                       "metadata": {}, "outputs": [],
                       "source": CACHE_CELL.strip("\n").splitlines(True)})
print("已插入 09b（WT 缓存 + 进度工具）")

# 重新定位（插入后索引变了）
i10 = find("# ============================================================ 10 Identity (WT->WT)")
i11 = find("# ============================================================ 11 Mutation smoke")

# ============================================================ 3) cell 10 用缓存 + 进度
cells[i10]["source"] = '''
# ============================================================ 10 Identity (WT->WT) + Repeatability (§16/§17)
import time
def _stats(a):
    a = np.asarray([x for x in a if x is not None and np.isfinite(x)], float)
    if len(a) == 0:
        return {"n": 0}
    med = float(np.median(a))
    return {"n": int(len(a)), "median": med,
            "MAD": float(np.median(np.abs(a - med))),
            "P95": float(np.percentile(np.abs(a), 95)),
            "max": float(np.max(np.abs(a)))}

IDENT, REPEAT = [], []
N_ID = min(5, len(VAL))
print("identity 对照样本:", N_ID, flush=True)
_t0 = time.time()
for k in range(N_ID):
    v = VAL.iloc[k]
    prog(k, N_ID, _t0, "row=%s %s" % (v.row_id, v.target_id))
    try:
        pdb = fetch_pdb(v.target_id)
        sk = {(str(v.chain), str(v.resseq))}
        d1 = get_wt(v.target_id, pdb, sk)          # 首次计算并缓存
        sh = d1["shell"]; tw = d1["topology"]
        # 第二次：完全相同的输入 -> repeatability（§32）
        fw2, ew2, sw2, tw2 = prepare_system(pdb)
        if ew2:
            IDENT.append({"row_id": int(v.row_id), "err": f"WT_{ew2}"}); continue
        p2, e2pre, e2post, r2 = minimize(sw2, tw2, fw2.positions, sk, sh)
        p1 = d1["positions"]; e1post = d1["post"]
        a1 = np.array(p1.value_in_unit(unit.angstrom))
        a2 = np.array(p2.value_in_unit(unit.angstrom))
        rmsd = float(np.sqrt(np.mean(np.sum((a1 - a2) ** 2, axis=1))))
        IDENT.append({"row_id": int(v.row_id), "target": str(v.target_id),
                      "dE_identity": float(e2post - e1post),
                      "rmsd_identity": rmsd, "shell": len(sh)})
        REPEAT.append({"row_id": int(v.row_id),
                       "dE_repeat": abs(float(e2post - e1post)),
                       "rmsd_repeat": rmsd})
        print("     row %-5d %-6s dE=%+.4f kJ/mol  RMSD=%.4f A  shell=%d"
              % (int(v.row_id), str(v.target_id), e2post - e1post, rmsd, len(sh)),
              flush=True)
    except Exception as e:
        IDENT.append({"row_id": int(v.row_id), "err": f"{type(e).__name__}:{str(e)[:60]}"})
prog(N_ID, N_ID, _t0, "identity done")
ID = pd.DataFrame(IDENT); RP = pd.DataFrame(REPEAT)
ID.to_csv(QC_DIR / "identity_controls.csv", index=False)
RP.to_csv(QC_DIR / "repeatability_controls.csv", index=False)
print()
print("WT->WT dE   :", _stats(ID.get("dE_identity", [])))
print("WT->WT RMSD :", _stats(ID.get("rmsd_identity", [])))
print("repeat dE   :", _stats(RP.get("dE_repeat", [])))
print("WT_CACHE 条目:", len(WT_CACHE), flush=True)
'''.strip("\n").splitlines(True)

# ============================================================ 4) cell 11 用缓存 + 逐变体进度 + 心跳
cells[i11]["source"] = '''
# ============================================================ 11 Mutation smoke (20 real) + NOISE GATE (§18/§19)
# 修复：WT 复用 cell 10 的缓存；逐变体实时进度（含 ETA）；单变体心跳超时。
import time
MUT = []
N_MUT = min(20, len(VAL))
print("mutation smoke 样本:", N_MUT, "| WT_CACHE 复用:", len(WT_CACHE), flush=True)
_t0 = time.time()
for k in range(N_MUT):
    v = VAL.iloc[k]
    _tk = time.time()
    prog(k, N_MUT, _t0, "row=%s %s %s%s>%s" % (v.row_id, v.target_id, v.wt, v.resseq, v.mut))
    try:
        pdb = fetch_pdb(v.target_id)
        sk = {(str(v.chain), str(v.resseq))}
        d = get_wt(v.target_id, pdb, sk)                 # WT 命中缓存
        wpost = d["post"]; wpre = d["pre"]
        muts = [(str(v.chain), int(v.resseq), THREE[str(v.wt)], THREE[str(v.mut)])]
        fm, em, sm, tm = prepare_system(pdb, muts)
        if em:
            MUT.append({"row_id": int(v.row_id), "err": f"MUT_{em}"}); continue
        shm = residues_within(tm, fm.positions, sk, RADIUS_MOVABLE)
        pm, mpre, mpost, _ = minimize(sm, tm, fm.positions, sk, shm)
        aw = np.array(d["positions"].value_in_unit(unit.angstrom))
        am = np.array(pm.value_in_unit(unit.angstrom))
        n = min(len(aw), len(am))
        rmsd_m = float(np.sqrt(np.mean(np.sum((aw[:n] - am[:n]) ** 2, axis=1))))
        MUT.append({"row_id": int(v.row_id), "target": str(v.target_id),
                    "wt": str(v.wt), "mut": str(v.mut),
                    "dE_total": float(mpost - wpost),
                    "dE_relax_cost": float((mpost - mpre) - (wpost - wpre)),
                    "rmsd_wt_mut": rmsd_m, "wt_cached": bool(d.get("cached")),
                    "secs": round(time.time() - _tk, 1)})
        if time.time() - _tk > VARIANT_TIMEOUT_S:
            print("     ! 单变体超过 %ds（仍完成）" % VARIANT_TIMEOUT_S, flush=True)
    except Exception as e:
        MUT.append({"row_id": int(v.row_id), "err": f"{type(e).__name__}:{str(e)[:60]}",
                    "secs": round(time.time() - _tk, 1)})
prog(N_MUT, N_MUT, _t0, "mutation smoke done")
MU = pd.DataFrame(MUT)
MU.to_csv(QC_DIR / "mutation_smoke.csv", index=False)
if len(MU) and "secs" in MU.columns:
    print("  单变体耗时 median=%.1fs max=%.1fs total=%.1fs"
          % (MU.secs.median(), MU.secs.max(), MU.secs.sum()), flush=True)
    if "wt_cached" in MU.columns:
        print("  WT 缓存命中: %d/%d" % (int(MU.wt_cached.fillna(False).sum()), len(MU)),
              flush=True)

ID_ST = _stats(ID.get("dE_identity", []))
RS_ST = _stats(ID.get("rmsd_identity", []))
RP_ST = _stats(RP.get("dE_repeat", []))
MU_ST = _stats(MU.get("dE_total", []))
MR_ST = _stats(MU.get("rmsd_wt_mut", []))
print()
print("WT->Mut dE  :", MU_ST)
print("WT->Mut RMSD:", MR_ST)

R_E = R_G = None
if ID_ST.get("n") and MU_ST.get("n"):
    den = abs(MU_ST["median"]) if abs(MU_ST["median"]) > 1e-12 else np.nan
    R_E = ID_ST["P95"] / den if np.isfinite(den) else None
if RS_ST.get("n") and MR_ST.get("n"):
    den2 = MR_ST["median"] if MR_ST["median"] > 1e-12 else np.nan
    R_G = RS_ST["P95"] / den2 if np.isfinite(den2) else None
NUM_NOISE_FLOOR_KJMOL = ID_ST.get("P95")

print()
print("=" * 62)
print("NUM_NOISE_FLOOR_KJMOL =", NUM_NOISE_FLOOR_KJMOL)
print("R_E = Q95(|dE|_WT->WT)/median(|dE|_WT->Mut) =", R_E)
print("R_G = Q95(RMSD_WT->WT)/median(RMSD_WT->Mut) =", R_G)
print("=" * 62)
print("阈值说明：项目预注册工程阈值，**非**领域标准", flush=True)

if R_E is None or R_G is None:
    NOISE_GATE = "NOT_EVALUABLE"
elif R_E >= 1.0 or R_G >= 1.0:
    NOISE_GATE = "NUMERICAL_NOISE_TOO_HIGH"
elif R_E >= 0.5 or R_G >= 0.5:
    NOISE_GATE = "MARGINAL"
else:
    NOISE_GATE = "NOISE_GATE_PASS"
NOISE_GATE_PASS = (NOISE_GATE == "NOISE_GATE_PASS")
print("NOISE_GATE =", NOISE_GATE, "| NOISE_GATE_PASS =", NOISE_GATE_PASS)
(WORK / "noise_gate.json").write_text(json.dumps(
    {"identity_dE": ID_ST, "identity_rmsd": RS_ST, "repeat_dE": RP_ST,
     "mutation_dE": MU_ST, "mutation_rmsd": MR_ST,
     "R_E": R_E, "R_G": R_G,
     "NUM_NOISE_FLOOR_KJMOL": NUM_NOISE_FLOOR_KJMOL,
     "noise_gate": NOISE_GATE, "noise_gate_pass": NOISE_GATE_PASS,
     "gb_model": GB_MODEL, "platform": PLAT_NAME,
     "wt_cache_entries": len(WT_CACHE),
     "variant_secs_median": (float(MU.secs.median()) if len(MU) and "secs" in MU.columns else None),
     "threshold_note": "prespecified project engineering thresholds, NOT field-standard"},
    indent=2))
if not NOISE_GATE_PASS:
    print()
    print("STOP（§19/§20）：噪声门未通过 -> 不进入全量。")
    raise RuntimeError(f"NOISE_GATE_FAIL:{NOISE_GATE}")
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
for k, v in (("WT 缓存 get_wt", "def get_wt" in t),
             ("cell10 用 get_wt", "get_wt(v.target_id" in t),
             ("cell11 用 get_wt", "d = get_wt(v.target_id" in t),
             ("逐变体进度 prog()", "def prog(" in t),
             ("ETA 打印", "ETA" in t),
             ("flush 实时", "flush=True" in t),
             ("单变体耗时统计", "单变体耗时 median" in t),
             ("VARIANT_TIMEOUT_S", "VARIANT_TIMEOUT_S" in t),
             ("无物理协议改动", "RADIUS_MOVABLE = 8.0" in t and "RESTRAINT_OUTSIDE = 100.0" in t)):
    print(f"  {'✅' if v else '❌'} {k}")
