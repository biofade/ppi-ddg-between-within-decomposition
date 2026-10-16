# -*- coding: utf-8 -*-
"""scripts/v8/34_multipoint_mechanism.py -- 机制检验：能力分化是否由**多点突变**驱动？

## 动机

PHASE 13 的外部检验发现：在**仅单点突变**的 `U_HERMES_CANONICAL`（2865 行）上，
内部模型的 ρ(Overall r, Macro ρ) = **+0.474**，
而在**含多点突变**的 `U_CANONICAL` 全集（4473 行）上为 **−0.20 ~ −0.26**。

⇒ 候选机制 (1)：分化由**多点突变**驱动。

## 做法（纯内部，无需外部资源）

在 `U_CANONICAL` 内部按 `n_sites` 分层，**同一批模型、同一套度量**：
  * 全部行（4473）
  * 仅单点（n_sites == 1）
  * 仅多点（n_sites > 1）
  * 单点补集 / 多点补集 的 target 规模匹配子集

对每层计算 ρ(Overall r, Macro ρ) 与 ρ(Overall r, Eff@10)。

**注意**：单点子集与多点子集的 target 集合不同，直接比较 ρ 会混入 target 组成效应。
因此额外做一次**目标规模匹配**：把两个子集分别限制在**共同的 target** 上
（即同时有单点与多点突变的 target），使比较尽可能受控。
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

# ---------- 数据 ----------
H = pd.read_csv(OUT / "external" / "hermes_matched_final.csv", low_memory=False)
z = np.load(R / "experiments" / "v7_1" / "phase2_rank" / "z_predictions.npz",
            allow_pickle=True)
Y = np.asarray(z["y"], float)
lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)
H = H.sort_values("pair_row").reset_index(drop=True)
assert len(H) == len(lab) == len(Y)
H["y_true"] = Y
H["n_sites"] = pd.to_numeric(H.n_sites, errors="coerce")
print(f"canonical {len(H)} 行")
print(f"  n_sites 分布: {H.n_sites.value_counts().sort_index().head(8).to_dict()}")
print(f"  单点 = {int((H.n_sites==1).sum())}  多点 = {int((H.n_sites>1).sum())}")

TGT_ALL = H.target_id_canonical.astype(str).to_numpy()
FAM_ALL = H.family_id_canonical.astype(str).to_numpy()

# 共同 target：同时含单点与多点的 target
g = H.groupby(TGT_ALL).n_sites.agg(["min", "max"])
common = set(g[(g["min"] == 1) & (g["max"] > 1)].index)
print(f"  同时含单点与多点的 target = {len(common)} / {H.target_id_canonical.nunique()}")

# ---------- 预测装载（复用 33 的对齐逻辑）----------
NONPRED = {"y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows",
           "z_raw", "y_cal", "family", "seed"}
def is_label(k):
    kl = str(k).lower()
    return kl == "y" or kl.endswith("_y") or kl.endswith("_ddg")

ART = [
    ("experiments/v7_1/phase2_rank/z_predictions.npz", "rowid"),
    ("experiments/v7_2r_1/phase2_full_eval/full_ranker_predictions.npz", "rowid"),
    ("experiments/v7_2r_2/phase5_eval/trackB_predictions.npz", "rowid"),
    ("experiments/v7_2r_3/phase5_eval/ppi_predictions.npz", "rowid"),
    ("experiments/v7_3/phase5_rankers/family_predictions.npz", "rowid"),
    ("experiments/v7_0/phase2_ablation/ablation_predictions.npz", "align"),
    ("experiments/v7_0/phase3_losses/loss_predictions.npz", "align"),
    ("experiments/v7_1/phase5_lsr/lsr_predictions.npz", "align"),
]
ROWS_ALL = lab.pair_row.to_numpy(int)
pos_all = pd.Series(np.arange(len(ROWS_ALL)), index=ROWS_ALL)

P, names = [], []
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
        else:
            if rid is None or len(a) != len(rid):
                continue
            m = pos_all.reindex(rid)
            if m.isna().any():
                continue
            av = np.full(len(ROWS_ALL), np.nan)
            av[m.to_numpy(int)] = a
        P.append(av); names.append(key)
P = np.vstack(P)
print(f"预测矩阵 P = {P.shape}")


def rank_agree(y, tgt, fam, Psub):
    """给定行子集，返回模型层面的 (overall_r, macro_rho, eff10) 与秩一致性。"""
    ov, mc, ef = [], [], []
    for i in range(Psub.shape[0]):
        a = Psub[i]
        m = np.isfinite(y) & np.isfinite(a)
        if m.sum() < 50:
            ov.append(np.nan); mc.append(np.nan); ef.append(np.nan); continue
        U = M.prep_universe(y[m], a[m], tgt[m], fam[m])
        r = M.compute_all(U)
        ov.append(r["overall_pearson"]); mc.append(r["macro_spearman"])
        ef.append(r["efficiency@10"])
    ov = np.asarray(ov, float); mc = np.asarray(mc, float)
    ef = np.asarray(ef, float)
    v = np.isfinite(ov) & np.isfinite(mc) & np.isfinite(ef)
    if v.sum() < 8:
        return {"n_models": int(v.sum()), "n_rows": int(len(y)),
                "n_targets_ge5": None, "rho_overall_macro": np.nan,
                "tau_overall_macro": np.nan, "rho_overall_eff10": np.nan}
    return {
        "n_models": int(v.sum()), "n_rows": int(len(y)),
        "n_targets_ge5": int((pd.Series(tgt).value_counts() >= MIN_N).sum()),
        "rho_overall_macro": float(spearmanr(ov[v], mc[v]).correlation),
        "tau_overall_macro": float(kendalltau(ov[v], mc[v]).correlation),
        "rho_overall_eff10": float(spearmanr(ov[v], ef[v]).correlation),
        "median_overall_r": float(np.nanmedian(ov[v])),
        "median_macro": float(np.nanmedian(mc[v])),
        "median_eff10": float(np.nanmedian(ef[v])),
    }


# ---------- 分层 ----------
STRATA = {}
m1 = (H.n_sites == 1).to_numpy()
m2 = (H.n_sites > 1).to_numpy()
incom = np.isin(TGT_ALL, list(common))
STRATA["ALL_rows"] = np.ones(len(H), bool)
STRATA["single_point"] = m1
STRATA["multi_point"] = m2
STRATA["single_point__common_targets"] = m1 & incom
STRATA["multi_point__common_targets"] = m2 & incom

rows = []
print("\n" + "=" * 84)
print("机制检验：按 n_sites 分层的秩一致性")
print("=" * 84)
print(f"{'stratum':<34}{'rows':>7}{'tgt>=5':>8}{'models':>8}"
      f"{'rho(r,macro)':>14}{'tau':>9}{'rho(r,eff)':>11}")
for nm, mask in STRATA.items():
    res = rank_agree(H.y_true.to_numpy(float)[mask], TGT_ALL[mask],
                     FAM_ALL[mask], P[:, mask])
    res["stratum"] = nm
    rows.append(res)
    print(f"{nm:<34}{res['n_rows']:>7}{str(res['n_targets_ge5']):>8}"
          f"{res['n_models']:>8}{res['rho_overall_macro']:>+14.4f}"
          f"{res['tau_overall_macro']:>+9.4f}{res['rho_overall_eff10']:>+11.4f}")

D = pd.DataFrame(rows)[["stratum", "n_rows", "n_targets_ge5", "n_models",
                        "rho_overall_macro", "tau_overall_macro",
                        "rho_overall_eff10", "median_overall_r",
                        "median_macro", "median_eff10"]]
D.to_csv(OUT / "split_inflation" / "multipoint_mechanism.csv", index=False,
         encoding="utf-8")

# ---------- 判定 ----------
a = D[D.stratum == "single_point"].iloc[0]
b = D[D.stratum == "multi_point"].iloc[0]
ac = D[D.stratum == "single_point__common_targets"].iloc[0]
bc = D[D.stratum == "multi_point__common_targets"].iloc[0]
print("\n" + "=" * 84)
print("判定")
print("=" * 84)
print(f"  仅单点       ρ(r,macro) = {a.rho_overall_macro:+.4f}  (n_rows={a.n_rows})")
print(f"  仅多点       ρ(r,macro) = {b.rho_overall_macro:+.4f}  (n_rows={b.n_rows})")
print(f"  共同 target 上：单点 {ac.rho_overall_macro:+.4f}  vs  多点 {bc.rho_overall_macro:+.4f}")
print()
sig = "SIGN_FLIP" if (a.rho_overall_macro * b.rho_overall_macro < 0) else "same_sign"
print(f"  单点 vs 多点符号  = {sig}")
print(f"  单点 vs 全集      = {a.rho_overall_macro:+.4f} vs "
      f"{D[D.stratum=='ALL_rows'].iloc[0].rho_overall_macro:+.4f}")
verdict = ("MULTIPOINT_DRIVES_DISSOCIATION" if sig == "SIGN_FLIP"
           else "MULTIPOINT_NOT_THE_DRIVER")
print(f"\n  ⇒ MECHANISM_VERDICT = {verdict}")

meta = {"strata": D.to_dict(orient="records"),
        "verdict": verdict,
        "caveat": ("单点/多点子集的 target 集合不同；common_targets 变体用于"
                   "尽量受控。仍非完全受控实验。")}
(OUT / "external" / "multipoint_mechanism_meta.json").write_text(
    json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"\n写出 split_inflation/multipoint_mechanism.csv + "
      f"external/multipoint_mechanism_meta.json")
