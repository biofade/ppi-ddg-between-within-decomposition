# -*- coding: utf-8 -*-
"""scripts/v8/41_expand_external_pool.py -- 扩池：HERMES 多 backbone 变体 + proteinMPNN 基线。

## 为什么这不是"再试一次"，而是新证据

`hermes_experiments/skempi/results/` 目录列举显示：
  * **15 个 HERMES 变体**（`py/bp` × `000/050` × 不同 finetune / esmfold / untrained）
  * **3 个 proteinMPNN 基线**（`v_48_002 / 020 / 030`）
每个目录内都有 `…-use_mt_structure=0.csv`（逐突变预测）。

⇒ 加上 CATH-ddG 表里的 5 个基线，外部池可达 **≥8 个可评价模型**，
首次可能跨过「可检门槛」。

## 严格性说明（防伪多样性）

多个 HERMES 变体**共享同一架构与训练数据**，只是 seed / backbone / finetune 不同。
因此它们是**同一方法族的变体**，不是独立方法。
**本脚本同时给出「按变体计」与「按方法族计」两套读数**，
不得只报有利于结论的那一套。
"""

from __future__ import annotations

import json
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(R / "scripts" / "v8"))
import v8metrics as M  # noqa: E402

OUT = R / "experiments" / "v8"
EXT = OUT / "external"
DL = Path("F:/dsh/downloads/hermes_variants")
DL.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "Mozilla/5.0"}
BASE = ("https://raw.githubusercontent.com/StatPhysBio/"
        "hermes_experiments/main/skempi/results")
MIN_N = 5

VARIANTS = {
    "HERMES_py000": ("hermes_py_000_ft_cdna117k_ddg_st",
                     "skempi_v2_cleaned_NO_1KBH-hermes_py_000_ft_cdna117k_ddg_st-use_mt_structure=0.csv"),
    "HERMES_py050": ("hermes_py_050_ft_cdna117k_ddg_st",
                     "skempi_v2_cleaned_NO_1KBH-hermes_py_050_ft_cdna117k_ddg_st-use_mt_structure=0.csv"),
    "HERMES_bp000": ("hermes_bp_000_ft_cdna117k_ddg_st",
                     "skempi_v2_cleaned_NO_1KBH-hermes_bp_000_ft_cdna117k_ddg_st-use_mt_structure=0.csv"),
    "HERMES_bp050": ("hermes_bp_050_ft_cdna117k_ddg_st",
                     "skempi_v2_cleaned_NO_1KBH-hermes_bp_050_ft_cdna117k_ddg_st-use_mt_structure=0.csv"),
    "HERMES_py000_esmfold": ("hermes_py_000_ft_esmfold_cdna117k_ddg_st",
                             "skempi_v2_cleaned_NO_1KBH-hermes_py_000_ft_esmfold_cdna117k_ddg_st-use_mt_structure=0.csv"),
    "proteinMPNN_002": ("proteinmpnn_v_48_002",
                        "skempi_v2_cleaned_NO_1KBH-num_seq_per_target=10-use_mt_structure=0.csv"),
    "proteinMPNN_020": ("proteinmpnn_v_48_020",
                        "skempi_v2_cleaned_NO_1KBH-num_seq_per_target=10-use_mt_structure=0.csv"),
    "proteinMPNN_030": ("proteinmpnn_v_48_030",
                        "skempi_v2_cleaned_NO_1KBH-num_seq_per_target=10-use_mt_structure=0.csv"),
}
FAMILY = {"HERMES_py000": "HERMES", "HERMES_py050": "HERMES",
          "HERMES_bp000": "HERMES", "HERMES_bp050": "HERMES",
          "HERMES_py000_esmfold": "HERMES",
          "proteinMPNN_002": "proteinMPNN", "proteinMPNN_020": "proteinMPNN",
          "proteinMPNN_030": "proteinMPNN"}

files = {}
for tag, (d, fn) in VARIANTS.items():
    dst = DL / f"{tag}.csv"
    if dst.exists() and dst.stat().st_size > 5000:
        files[tag] = dst; continue
    try:
        req = urllib.request.Request(f"{BASE}/{d}/{fn}", headers=UA)
        with urllib.request.urlopen(req, timeout=90) as r:
            b = r.read()
        dst.write_bytes(b)
        files[tag] = dst
        print(f"  [ok] {tag:<24} {len(b):>8,} B")
    except Exception as e:  # noqa: BLE001
        print(f"  [fail] {tag:<24} {type(e).__name__}")
print(f"\n下载成功 = {len(files)} / {len(VARIANTS)}")

# ---------- canonical ----------
lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)
ms = pd.read_csv(R / "data" / "v5_0" / "mutation_sites.csv", low_memory=False)
ms["pair_row"] = pd.to_numeric(ms.pair_row, errors="coerce")
z = np.load(R / "experiments" / "v7_1" / "phase2_rank" / "z_predictions.npz",
            allow_pickle=True)
Y = np.asarray(z["y"], float)
def s(x):
    return str(x).strip().upper()
ms["tok"] = ms.wt_aa.map(s) + ms.ann_pos1.astype("Int64").astype(str) + ms.mut_aa.map(s)
agg = (ms.sort_values("pair_row").groupby("pair_row", as_index=False)
         .agg(pdb4=("pdb_id", lambda v: str(v.iloc[0]).upper().strip()[:4]),
              n_sites=("tok", "size"),
              mut_str=("tok", lambda v: ",".join(sorted(set(v))))))
C = lab.merge(agg, on="pair_row", how="left")
assert len(C) == len(lab)
C["y_true"] = Y
TGT = C.target_id_canonical.astype(str).to_numpy()
FAM = C.family_id_canonical.astype(str).to_numpy()

# canonical key：单点 -> (PDB4|链?|pos|mut)。HERMES 用 <wt><pos><mut>（无链）
_pos = pd.to_numeric(ms.groupby("pair_row").ann_pos1.first()
                     .reindex(C.pair_row.to_numpy()).to_numpy(), errors="coerce")
_wt = ms.groupby("pair_row").wt_aa.first().reindex(C.pair_row.to_numpy()).to_numpy()
_mut = ms.groupby("pair_row").mut_aa.first().reindex(C.pair_row.to_numpy()).to_numpy()
C["key_h"] = (C.pdb4.astype(str) + "|"
              + pd.Series(_wt).astype(str).str.upper().to_numpy()
              + pd.Series(_pos).astype("float").astype("Int64").astype(str).to_numpy()
              + pd.Series(_mut).astype(str).str.upper().to_numpy())

# ---------- 评价 ----------
rows = []
print("\n" + "=" * 84)
print("扩池后的外部模型（U_HERMES_CANONICAL 口径，canonical 行级）")
print("=" * 84)
print(f"{'model':<26}{'family':<14}{'n':>6}{'tgt>=5':>8}{'Overall r':>11}"
      f"{'Macro rho':>11}{'C_W':>9}{'Eff@10':>9}")
for tag, p in files.items():
    h = pd.read_csv(p, low_memory=False)
    h["pdb4"] = h.PDB_filename.astype(str).str.upper().str.strip().str[:4]
    # HERMES: `mutant` = <wt><pos><mut>；proteinMPNN: 同样格式
    mcol = "mutant" if "mutant" in h.columns else "Mutation(s)_cleaned"
    h["key_h"] = h.pdb4 + "|" + h[mcol].astype(str).str.strip().str.upper()
    lp = next((c for c in h.columns
               if ("log_proba" in c.lower() or "log_p_mt" in c.lower())), None)
    if lp is None:
        lp = next((c for c in h.columns
                   if "ddg" in c.lower() and "neg" in c.lower()), None)
    if lp is None:
        print(f"  {tag}: 无可用分数列 ({list(h.columns)[:6]})"); continue
    h["proxy"] = -pd.to_numeric(h[lp], errors="coerce")
    ha = h.groupby("key_h", as_index=False).agg(proxy=("proxy", "mean"))
    m = C.merge(ha, on="key_h", how="left")
    ok = m.proxy.notna().to_numpy() & np.isfinite(Y)
    if ok.sum() < 100:
        print(f"  {tag:<26}{FAMILY[tag]:<14}{int(ok.sum()):>6}  —— 命中不足")
        continue
    U = M.prep_universe(Y[ok], m.proxy.to_numpy()[ok], TGT[ok], FAM[ok])
    r = M.compute_all(U)
    r.update({"model": tag, "family": FAMILY[tag], "raw_col": lp})
    rows.append(r)
    print(f"{tag:<26}{FAMILY[tag]:<14}{r['n_rows']:>6}{r['n_targets_ge5']:>8}"
          f"{r['overall_pearson']:>+11.4f}{r['macro_spearman']:>+11.4f}"
          f"{r['C_W']:>+9.4f}{r['efficiency@10']:>+9.4f}")

D = pd.DataFrame(rows)
if len(D):
    D.to_csv(EXT / "external_expanded_pool.csv", index=False, encoding="utf-8")
    # 加上 CATH-ddG 表的 5 个基线（不同宇宙，单独标注）
    print(f"\n  HERMES/proteinMPNN 变体可评价数 = {len(D)}")
    sel = D.dropna(subset=["overall_pearson", "macro_spearman", "efficiency@10"])
    print(f"  其中三轴齐全 = {len(sel)}")
    if len(sel) >= 5:
        r1 = spearmanr(sel.overall_pearson, sel.macro_spearman).correlation
        r3 = spearmanr(sel.macro_spearman, sel["efficiency@10"]).correlation
        print(f"\n  ★ 按**变体**计 (n={len(sel)}): "
              f"rho(Overall r, Macro rho) = {r1:+.4f} | rho(macro,eff) = {r3:+.4f}")
        fam = sel.groupby("family").agg(
            overall_pearson=("overall_pearson", "median"),
            macro_spearman=("macro_spearman", "median"),
            efficiency=("efficiency@10", "median")).reset_index()
        print(f"  按**方法族**计 (n={len(fam)}):")
        for r in fam.itertuples():
            print(f"    {r.family:<14} Overall r={r.overall_pearson:+.4f} "
                  f"Macro rho={r.macro_spearman:+.4f} Eff@10={r.efficiency:+.4f}")
        if len(fam) >= 4:
            r1f = spearmanr(fam.overall_pearson, fam.macro_spearman).correlation
            print(f"    rho(Overall r, Macro rho) 按族 = {r1f:+.4f}")
        print(f"\n  ⚠ 严格性：{len(sel)} 个变体来自仅 {sel.family.nunique()} 个方法族")
        print("     ⇒ 不得只报「按变体」的数值；两套都要报。")

meta = {"variants_downloaded": sorted(files.keys()),
        "n_evaluable_variants": int(len(D)),
        "n_families": int(D.family.nunique()) if len(D) else 0,
        "strictness": ("多 HERMES 变体共享架构与训练数据，属同一方法族；"
                       "必须同时报告「按变体」与「按方法族」两套读数"),
        "sign_convention": "proxy = -(log_proba_mt - log_proba_wt)（或负 ddG 列）",
        "combined_with_cathddg_baselines": 5,
        "total_external_pool_if_merged": int(len(D)) + 5,
        "caveat": "与 CATH-ddG 的 5 个基线属**不同宇宙**，合并时须分别标注"}
(EXT / "external_expanded_pool_meta.json").write_text(
    json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"\n写出 external_expanded_pool.csv / _meta.json")
