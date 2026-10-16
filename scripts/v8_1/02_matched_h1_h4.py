# -*- coding: utf-8 -*-
"""scripts/v8_1/02_matched_h1_h4.py -- V8.1 §7/§8/§9/§10/§11：匹配后的 H1/H4/H2/AXIS-C 重算。

**规则**
* 跨模型比较**只在同一 `rowset_hash` 内**进行（§3）；
* 秩一致使用 `-RMSE`（§2）；
* 主池由**performance-blind 字典序规则**选取（§3.3）：
  (1) provenance 干净度 → (2) 独立 method-family 数 → (3) n_ge5 → (4) n_rows
  **在读取指标值之前**即固定。
* AXIS C 使用解析几何 `D(n,K)=K(n-K)/[n(n+1)]`（§11.2）。

产出 `v8_1/FROZEN_RESULTS.json` 与若干 CSV。
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

V8 = R / "experiments" / "v8"
A = R / "v8_1"
MIN_N = 5
K = 10

UR = pd.read_csv(A / "MATCHED_UNIVERSE_REGISTRY.csv")
MR = pd.read_csv(A / "MODEL_FAMILY_REGISTRY.csv")
lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)
CAN = lab[["pair_row", "target_id_canonical", "family_id_canonical"]].rename(
    columns={"pair_row": "row_id", "target_id_canonical": "target",
             "family_id_canonical": "family"})
z = np.load(R / "experiments" / "v7_1" / "phase2_rank" / "z_predictions.npz",
            allow_pickle=True)
Y_ALL = np.asarray(z["y"], float)
pos_all = pd.Series(np.arange(len(CAN)), index=CAN.row_id)
mas = pd.read_csv(R / "experiments" / "v6_2" / "phase2_split4" /
                  "matched_analysis_set.csv", low_memory=False)

ARTMAP = {
    "V7.1_z": ("experiments/v7_1/phase2_rank/z_predictions.npz", "rowid"),
    "V7.2_nested": ("experiments/v7_2/phase5_nested_cv/ranker_predictions.npz", "rowid"),
    "V7.2R1_full": ("experiments/v7_2r_1/phase2_full_eval/full_ranker_predictions.npz", "rowid"),
    "V7.2R2_trackB": ("experiments/v7_2r_2/phase5_eval/trackB_predictions.npz", "rowid"),
    "V7.2R3_ppi": ("experiments/v7_2r_3/phase5_eval/ppi_predictions.npz", "rowid"),
    "V7.2R3_noise": ("experiments/v7_2r_3/phase8_noise_control/noise_control_predictions.npz", "rowid"),
    "V7.3_physics": ("experiments/v7_3/phase5_rankers/family_predictions.npz", "rowid"),
    "V7.0_ablation": ("experiments/v7_0/phase2_ablation/ablation_predictions.npz", "align"),
    "V7.0_losses": ("experiments/v7_0/phase3_losses/loss_predictions.npz", "align"),
    "V7.1_lsr": ("experiments/v7_1/phase5_lsr/lsr_predictions.npz", "align"),
    "V6_2_split4": ("experiments/v6_2/phase2_split4/split4_predictions.npz", "mas"),
    "V6_2_1_donly": ("experiments/v6_2_1/phase3_d_only/d_only_predictions.npz", "mas"),
}
NONPRED = {"y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows",
           "z_raw", "y_cal", "family", "seed"}
def is_label(k):
    kl = str(k).lower()
    return kl == "y" or kl.endswith("_y") or kl.endswith("_ddg")

# ---------- 构建每个 artifact 的 (y, tgt, fam, 预测矩阵) ----------
STORE = {}
for tag, (rel, kind) in ARTMAP.items():
    p = R / rel
    if not p.exists():
        continue
    zz = np.load(p, allow_pickle=True)
    rid = None
    for k in ("row_id", "pair_row", "rows", "row"):
        if k in zz.files:
            rid = np.asarray(zz[k]).ravel().astype(int); break
    if kind == "mas":
        yr = pd.to_numeric(mas.ddG, errors="coerce").to_numpy(float)
        tg = mas.pdb_id.astype(str).to_numpy()
        fa = mas.family_label.astype(str).to_numpy()
        keys, P = [], []
        for k in zz.files:
            if k in NONPRED or is_label(k):
                continue
            a = np.asarray(zz[k])
            try:
                a = a.astype(float)
            except (ValueError, TypeError):
                continue
            if a.ndim == 1 and len(a) == len(yr):
                keys.append(k); P.append(a)
    else:
        ids = rid if rid is not None else CAN.row_id.to_numpy(int)
        if kind == "align":
            ok = np.ones(len(CAN), bool)
            tg = CAN.target.to_numpy(); fa = CAN.family.to_numpy(); yr = Y_ALL
        else:
            m = pos_all.reindex(ids)
            ok = m.notna().to_numpy()
            if not ok.all():
                continue
            tg = CAN.target.to_numpy()[m.to_numpy(int)]
            fa = CAN.family.to_numpy()[m.to_numpy(int)]
            yr = Y_ALL[m.to_numpy(int)]
        keys, P = [], []
        for k in zz.files:
            if k in NONPRED or is_label(k):
                continue
            a = np.asarray(zz[k])
            try:
                a = a.astype(float)
            except (ValueError, TypeError):
                continue
            if a.ndim != 1:
                continue
            if kind == "align":
                if len(a) != len(CAN):
                    continue
                v = a
            else:
                if len(a) != len(ids):
                    continue
                v = a
            keys.append(k); P.append(v)
    if keys:
        STORE[tag] = {"y": yr, "tgt": tg, "fam": fa,
                      "P": np.vstack(P), "keys": keys,
                      "rowset_hash": UR[UR.artifact == tag].rowset_hash.iloc[0]}
print(f"载入 {len(STORE)} 个 artifact")

# ---------- 逐 rowset_hash 计算模型级三轴 ----------
rows = []
for tag, d in STORE.items():
    for i, k in enumerate(d["keys"]):
        pv = d["P"][i]
        m = np.isfinite(d["y"]) & np.isfinite(pv)
        if m.sum() < 100:
            continue
        U = M.prep_universe(d["y"][m], pv[m], d["tgt"][m], d["fam"][m])
        r = M.compute_all(U)
        r.update({"artifact": tag, "key": k, "rowset_hash": d["rowset_hash"]})
        rows.append(r)
D = pd.DataFrame(rows)
D = D.merge(MR[["artifact", "key", "method_family", "split_regime",
                "is_ablation", "is_random_split_model",
                "is_primary_representative"]],
            on=["artifact", "key"], how="left")
D.to_csv(A / "all_models_by_rowset.csv", index=False, encoding="utf-8")
print(f"模型级读数 = {len(D)} 行；rowset 数 = {D.rowset_hash.nunique()}")


# ---------- §3.3 performance-blind 主池选择 ----------
def elig(g, level):
    if level == "C":
        g = g[g.is_primary_representative == True]          # noqa: E712
    return g.dropna(subset=["overall_pearson", "macro_spearman",
                            "efficiency@10", "rmse"])


CANDS = []
for rh, g in D.groupby("rowset_hash"):
    for level in ("A", "C"):
        ge = elig(g, level)
        if len(ge) < 5:
            continue
        green_frac = float((ge.split_regime != "permissive").mean())
        fams = int(ge.method_family.nunique())
        CANDS.append({"rowset_hash": rh, "level": level, "n_models": len(ge),
                      "green_fraction": green_frac, "n_method_families": fams,
                      "n_ge5": int(ge.n_targets_ge5.median()),
                      "n_rows": int(ge.n_rows.median()),
                      "n_permissive": int((ge.split_regime == "permissive").sum())})
CAND = pd.DataFrame(CANDS).sort_values(
    ["green_fraction", "n_method_families", "n_ge5", "n_rows"],
    ascending=[False, False, False, False]).reset_index(drop=True)
CAND.to_csv(A / "primary_universe_candidates.csv", index=False, encoding="utf-8")
print("\n§3.3 主池候选（performance-blind 排序，读取指标值之前固定）:")
print(CAND.to_string(index=False))

# 主池 = 第一条；但若某 level-C 池 family 数>=2 且 green=1.0 优先
prim = None
for _, r in CAND.iterrows():
    if r.green_fraction == 1.0 and r.n_method_families >= 2:
        prim = r; break
if prim is None and len(CAND):
    prim = CAND.iloc[0]
print(f"\n⇒ 主池 = rowset {prim.rowset_hash} / level {prim.level} / "
      f"n_models={prim.n_models} / families={prim.n_method_families} / "
      f"n_ge5={prim.n_ge5} / rows={prim.n_rows}")

# ---------- §7 秩一致（-RMSE）逐 stratum ----------
def agree(g):
    g = g.dropna(subset=["overall_pearson", "macro_spearman", "efficiency@10", "rmse"])
    if len(g) < 5:
        return None
    g = g.assign(neg_rmse=-g.rmse)
    out = {"n_models": len(g), "n_method_families": int(g.method_family.nunique()),
            "n_green": int((g.split_regime != "permissive").sum())}
    for a, b, nm in (("overall_pearson", "macro_spearman", "overall_vs_macro"),
                     ("neg_rmse", "macro_spearman", "negrmse_vs_macro"),
                     ("overall_pearson", "efficiency@10", "overall_vs_util"),
                     ("macro_spearman", "efficiency@10", "macro_vs_util")):
        out[f"rho_{nm}"] = float(spearmanr(g[a], g[b]).correlation)
        out[f"tau_{nm}"] = float(kendalltau(g[a], g[b]).correlation)
    return out


STRAT = []
for rh, g in D.groupby("rowset_hash"):
    for level in ("A", "C"):
        ge = g if level == "A" else g[g.is_primary_representative == True]  # noqa: E712
        a = agree(ge)
        if a:
            a.update({"rowset_hash": rh, "level": level})
            STRAT.append(a)
ST = pd.DataFrame(STRAT)
ST.to_csv(A / "rank_agreement_by_stratum.csv", index=False, encoding="utf-8")
print("\n§7 逐 stratum 秩一致（含 -RMSE 修正）:")
print(ST.round(4).to_string(index=False))

# ---------- §8 H4 受控划分膨胀 ----------
V62 = STORE.get("V6_2_split4")
h4 = {}
if V62 is not None:
    yv, tg = V62["y"], V62["tgt"]
    LV = {"L1_random": "random", "L2_complex": "complex",
          "L3_seq40": "seq40", "L4_family": "family_OOD"}
    rec = {}
    for pref, nm in LV.items():
        ks = [k for k in V62["keys"] if k.startswith(pref + "_")
              and not k.endswith("_y") and not k.endswith("_pdb")]
        ovs, cbs, cws, rms, macs, efs = [], [], [], [], [], []
        for k in ks:
            i = V62["keys"].index(k)
            m = np.isfinite(yv) & np.isfinite(V62["P"][i])
            U = M.prep_universe(yv[m], V62["P"][i][m], tg[m])
            r = M.compute_all(U)
            ovs.append(r["overall_pearson"]); cbs.append(r["C_B_row"])
            cws.append(r["C_W_row"]); rms.append(r["rmse"])
            macs.append(r["macro_spearman"]); efs.append(r["efficiency@10"])
        rec[nm] = {"n_models": len(ks),
                   "overall_pearson": float(np.median(ovs)),
                   "C_B": float(np.median(cbs)), "C_W": float(np.median(cws)),
                   "rmse": float(np.median(rms)),
                   "macro_spearman": float(np.median(macs)),
                   "efficiency@10": float(np.median(efs))}
    h4 = rec
    print("\n§8 H4 受控划分膨胀（同一批 7523 行、同一批配置，仅改划分）:")
    for nm, v in rec.items():
        print(f"   {nm:<12} r={v['overall_pearson']:+.4f} C_B={v['C_B']:+.4f} "
              f"C_W={v['C_W']:+.4f} rmse={v['rmse']:.4f} "
              f"macro={v['macro_spearman']:+.4f}")
    if "random" in rec and "family_OOD" in rec:
        dr = rec["family_OOD"]["overall_pearson"] - rec["random"]["overall_pearson"]
        dcb = rec["family_OOD"]["C_B"] - rec["random"]["C_B"]
        dcw = rec["family_OOD"]["C_W"] - rec["random"]["C_W"]
        den = abs(dcb) + abs(dcw)
        h4["delta"] = {"delta_r": dr, "delta_C_B": dcb, "delta_C_W": dcw,
                       "fraction_between": abs(dcb) / den if den else np.nan,
                       "fraction_within": abs(dcw) / den if den else np.nan}
        print(f"   Δr={dr:+.4f}  ΔC_B={dcb:+.4f}  ΔC_W={dcw:+.4f}")
        print(f"   contribution fraction: between={h4['delta']['fraction_between']:.4f}"
              f"  within={h4['delta']['fraction_within']:.4f}")

# ---------- §11.2 AXIS C 解析几何 ----------
n = np.arange(1, 402)
Erand = np.where(n > K, K / n, 1.0)
Emax = np.where(n > K, K * (2 * n - K + 1) / (n * (n + 1)), 1.0)
Dn = Emax - Erand
peak = int(n[np.argmax(Dn)])
nstar = K + np.sqrt(K * K + K)
GEO = pd.DataFrame({"n": n, "E_rand": Erand, "E_max": Emax, "D": Dn})
GEO.to_csv(A / "axis_c_geometry.csv", index=False, encoding="utf-8")
print(f"\n§11.2 AXIS C 几何: D(n,K)=K(n-K)/[n(n+1)]; 峰值 n*={peak} "
      f"(解析 K+sqrt(K^2+K)={nstar:.2f}); D(50,10)={float(Dn[49]):.4f} "
      f"< D(20,10)={float(Dn[19]):.4f}")
print(f"   ⇒ 「n>=50 区分度更高」**数学上不成立**（D 在 n*≈20.5 后单调下降）")

# 实测 n<=K 占比（独立复验）
TGT_FULL = CAN.target.to_numpy()
cnt = pd.Series(TGT_FULL).value_counts()
ge5 = cnt[cnt >= MIN_N]
n_le_k = int((ge5 <= K).sum())
print(f"   n<=K=10 的 n>=5 target = {n_le_k}/{len(ge5)} = "
      f"{100*n_le_k/len(ge5):.1f}%  （历史报 40.8% ⇒ "
      f"{'复验一致' if abs(n_le_k/len(ge5)-0.408)<0.01 else '**不一致**'}）")

FR = {
    "version": "v8.1",
    "UPSTREAM_HASH_CHANGE": 0,
    "direction_bug": {"present": True, "affected_script": "12_decision_alignment.py",
                      "fix": "use -RMSE for rank agreement",
                      "sign_flip_observed_for": list(
                          json.load(open(A / "_core_meta.json",
                                         encoding="utf-8"))["bug_report"].keys())},
    "rowsets": {"n_artifacts": int(len(UR)),
                "n_distinct_rowset_hash": int(UR.rowset_hash.nunique()),
                "detail": UR.to_dict(orient="records")},
    "primary_universe": {"rowset_hash": str(prim.rowset_hash),
                         "level": str(prim.level),
                         "n_models": int(prim.n_models),
                         "n_method_families": int(prim.n_method_families),
                         "n_ge5": int(prim.n_ge5), "n_rows": int(prim.n_rows),
                         "green_fraction": float(prim.green_fraction),
                         "selection_rule": "performance-blind lexicographic: "
                                           "(1) provenance cleanliness "
                                           "(2) n independent method families "
                                           "(3) n_ge5 (4) n_rows"},
    "rank_agreement_by_stratum": ST.to_dict(orient="records"),
    "H4_controlled_split": h4,
    "axis_c_geometry": {"formula_D": "K(n-K)/(n(n+1))",
                        "n_star_analytic": float(nstar), "n_star_numeric": peak,
                        "D_at_n20": float(Dn[19]), "D_at_n50": float(Dn[49]),
                        "n_le_K_count": n_le_k, "n_ge5_count": int(len(ge5)),
                        "n_le_K_fraction": float(n_le_k / len(ge5)),
                        "forbidden_claim": "n>=50 improves fixed-K discrimination"},
    "provenance": json.load(open(A / "_core_meta.json", encoding="utf-8")),
}
json.dump(FR, open(A / "FROZEN_RESULTS.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
print(f"\n写出 FROZEN_RESULTS.json / rank_agreement_by_stratum.csv / "
      f"axis_c_geometry.csv / primary_universe_candidates.csv / all_models_by_rowset.csv")
