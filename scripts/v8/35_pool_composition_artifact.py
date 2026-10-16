# -*- coding: utf-8 -*-
"""scripts/v8/35_pool_composition_artifact.py -- 受控演示：**模型池组成会人为翻转指标秩相关**。

## 现象（PHASE 13 发现）

在**同一宇宙、同一批 4473 行**上：
```
67 个模型的池（含 51 个 V6.2 随机划分键）   rho(Overall r, Macro rho) = -0.2685
排除那 51 个（16 个 family-OOD 键）        rho(Overall r, Macro rho) = +0.1861
```
符号翻转。

## 受控演示设计

把 136 个模型**按训练划分的宽松度**分成两组（用键名前缀作为 provenance 标签）：

  * `PERMISSIVE`：随机 / 结构 / seq40 划分上训练的键
    （`L1_random_*`, `L2_complex_*`, `L3_seq40_*`, `A0..A6` 中的 base 变体）
  * `STRICT`：family-OOD 划分上训练的键（`L4_family_*`, `F*`, `R*`, `B*`, `p_*`, `z_*`）

三组分别计算 ρ(Overall r, Macro ρ)，并做**逐模型散点诊断**：
看 PERMISSIVE 组的负相关是否由「高 Overall r + 低 Macro ρ」的少数点驱动。

**关键**：因为三组在**同一宇宙同一行**上评价，差异**只能**来自模型池组成。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(R / "scripts" / "v8"))
import v8metrics as M  # noqa: E402

OUT = R / "experiments" / "v8"
MIN_N = 5

# ---------- 载入 4473 行宇宙上的全部模型 ----------
H = pd.read_csv(OUT / "external" / "hermes_matched_final.csv", low_memory=False)
z = np.load(R / "experiments" / "v7_1" / "phase2_rank" / "z_predictions.npz",
            allow_pickle=True)
Y = np.asarray(z["y"], float)
lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)
H = H.sort_values("pair_row").reset_index(drop=True)
TGT = H.target_id_canonical.astype(str).to_numpy()
FAM = H.family_id_canonical.astype(str).to_numpy()
assert len(H) == len(Y) == len(lab)

NONPRED = {"y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows",
           "z_raw", "y_cal", "family", "seed"}
def is_label(k):
    kl = str(k).lower()
    return kl == "y" or kl.endswith("_y") or kl.endswith("_ddg")

# provenance 标签：由键名判定训练划分宽松度
PERM_TOK = ("L1_random_", "L2_complex_", "L3_seq40_")
STRICT_TOK = ("L4_family_", "F0_", "F1_", "F2_", "F3_", "F4_", "F5_", "F6_",
              "F7_", "F8_", "F9_", "R0_", "R1_", "R2_", "R3_", "B0_", "B1_",
              "B2_", "B3", "B4_", "A2_smallPLM", "p_", "z_R", "z_A")

def prov(key: str) -> str:
    if any(key.startswith(t) for t in PERM_TOK):
        return "PERMISSIVE"
    if any(key.startswith(t) for t in STRICT_TOK):
        return "STRICT"
    return "OTHER"

ART = [
    ("experiments/v6_2/phase2_split4/split4_predictions.npz", "mas"),
    ("experiments/v6_2_1/phase3_d_only/d_only_predictions.npz", "mas"),
    ("experiments/v7_0/phase2_ablation/ablation_predictions.npz", "align"),
    ("experiments/v7_0/phase3_losses/loss_predictions.npz", "align"),
    ("experiments/v7_1/phase2_rank/z_predictions.npz", "rowid"),
    ("experiments/v7_1/phase5_lsr/lsr_predictions.npz", "align"),
    ("experiments/v7_2r_1/phase2_full_eval/full_ranker_predictions.npz", "rowid"),
    ("experiments/v7_2r_2/phase5_eval/trackB_predictions.npz", "rowid"),
    ("experiments/v7_2r_3/phase5_eval/ppi_predictions.npz", "rowid"),
    ("experiments/v7_3/phase5_rankers/family_predictions.npz", "rowid"),
]
ROWS_ALL = lab.pair_row.to_numpy(int)
pos_all = pd.Series(np.arange(len(ROWS_ALL)), index=ROWS_ALL)

# V6.2 的 mas 制品对齐到 matched_analysis_set(7523)，需映射到 canonical 的 4473
mas = pd.read_csv(R / "experiments" / "v6_2" / "phase2_split4" /
                  "matched_analysis_set.csv", low_memory=False)
mas_row = mas.row.to_numpy(int)
pos_of_mas = pos_all.reindex(mas_row)

rows = []
for rel, kind in ART:
    p = R / rel
    if not p.exists():
        continue
    zz = np.load(p, allow_pickle=True)
    rid = None
    for k in ("row_id", "pair_row", "rows", "row"):
        if k in zz.files:
            rid = np.asarray(zz[k]).ravel().astype(int); break
    for key in zz.files:
        if key in NONPRED or is_label(key):
            continue
        a = np.asarray(zz[key])
        try:
            a = a.astype(float)
        except (ValueError, TypeError):
            continue
        if a.ndim != 1:
            continue
        if kind == "align":
            if len(a) != len(ROWS_ALL):
                continue
            av = a
        elif kind == "rowid":
            if rid is None or len(a) != len(rid):
                continue
            m = pos_all.reindex(rid)
            if m.isna().any():
                continue
            av = np.full(len(ROWS_ALL), np.nan)
            av[m.to_numpy(int)] = a
        else:  # mas
            if len(a) != len(mas_row):
                continue
            av = np.full(len(ROWS_ALL), np.nan)
            ok = pos_of_mas.notna().to_numpy()
            av[pos_of_mas[ok].to_numpy(int)] = a[ok]
        mm = np.isfinite(Y) & np.isfinite(av)
        if mm.sum() < 100:
            continue
        U = M.prep_universe(Y[mm], av[mm], TGT[mm], FAM[mm])
        r = M.compute_all(U)
        r["model_name"] = key
        r["provenance"] = prov(key)
        rows.append(r)

D = pd.DataFrame(rows)
D.to_csv(OUT / "external" / "pool_composition_models.csv", index=False,
         encoding="utf-8")
print(f"模型总数 = {len(D)}")
print(D.provenance.value_counts().to_string())

# ---------- 分组秩相关 ----------
COLS = ["overall_pearson", "macro_spearman", "efficiency@10", "C_W"]
def agree(g):
    g = g.dropna(subset=COLS)
    if len(g) < 6:
        return None
    return {"n": len(g),
            "rho_r_macro": float(spearmanr(g.overall_pearson, g.macro_spearman).correlation),
            "tau_r_macro": float(kendalltau(g.overall_pearson, g.macro_spearman).correlation),
            "rho_r_eff": float(spearmanr(g.overall_pearson, g["efficiency@10"]).correlation),
            "med_r": float(g.overall_pearson.median()),
            "med_macro": float(g.macro_spearman.median())}

print("\n" + "=" * 82)
print("受控演示：同一宇宙（4473 行）、只改模型池组成")
print("=" * 82)
res = {}
groups = [("ALL", D), ("PERMISSIVE（随机/complex/seq40 划分训练）",
                       D[D.provenance == "PERMISSIVE"]),
          ("STRICT（family-OOD 划分训练）", D[D.provenance == "STRICT"]),
          ("OTHER", D[D.provenance == "OTHER"])]
print(f"{'group':<44}{'n':>5}{'rho(r,macro)':>15}{'tau':>9}{'rho(r,eff)':>12}")
for nm, g in groups:
    a = agree(g)
    if a is None:
        print(f"{nm:<44}  (模型过少)"); continue
    res[nm] = a
    print(f"{nm:<44}{a['n']:>5}{a['rho_r_macro']:>+15.4f}"
          f"{a['tau_r_macro']:>+9.4f}{a['rho_r_eff']:>+12.4f}")

# ---------- 逐模型诊断：谁在驱动负相关 ----------
per = D[D.provenance == "PERMISSIVE"].dropna(subset=COLS)
st = D[D.provenance == "STRICT"].dropna(subset=COLS)
print("\n" + "=" * 82)
print("诊断：PERMISSIVE 组中「高 Overall r + 低 Macro ρ」的极端点")
print("=" * 82)
if len(per):
    hi_r = per.overall_pearson.quantile(0.75)
    lo_m = per.macro_spearman.quantile(0.25)
    ext = per[(per.overall_pearson >= hi_r) & (per.macro_spearman <= lo_m)]
    print(f"  PERMISSIVE: Overall r p75={hi_r:+.4f}, Macro ρ p25={lo_m:+.4f}")
    print(f"  同时满足 (r>=p75 且 macro<=p25) 的模型数 = {len(ext)} / {len(per)}"
          f"  ({100*len(ext)/max(1,len(per)):.1f}%)")
    print(f"  这些点的 Overall r 中位 = {ext.overall_pearson.median():+.4f}"
          f", Macro ρ 中位 = {ext.macro_spearman.median():+.4f}")
    print(f"  它们对**全部 136 模型**池 ρ 的影响：")
    keep = D[~D.model_name.isin(ext.model_name)].dropna(subset=COLS)
    allr = spearmanr(D.dropna(subset=COLS).overall_pearson,
                     D.dropna(subset=COLS).macro_spearman).correlation
    kee = spearmanr(keep.overall_pearson, keep.macro_spearman).correlation
    print(f"    剔除前 (n={len(D.dropna(subset=COLS))}) ρ = {allr:+.4f}")
    print(f"    剔除后 (n={len(keep)})           ρ = {kee:+.4f}")
    res["ALL_minus_extremes"] = {"n": len(keep), "rho_r_macro": float(kee)}

meta = {"groups": res,
        "finding": ("模型池组成（不同训练划分的模型混池）可人为改变指标间秩相关的"
                    "**符号**，即使评价宇宙与行完全一致。"),
        "methodological_lesson": ("报告「指标间模型排名秩相关」时必须声明模型池的"
                                  "训练划分 provenance；否则数值不可解释。"),
        "caveat": ("PERMISSIVE/STRICT 标签由键名前缀推断，非上游显式字段；"
                   "极端点剔除为诊断性（post-hoc），不作确证性推断。")}
(OUT / "external" / "pool_composition_meta.json").write_text(
    json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"\n写出 pool_composition_models.csv / pool_composition_meta.json")
