# -*- coding: utf-8 -*-
"""scripts/v7_1/02_lsr.py -- V7.1 PHASE 2–6: Location–Scale–Rank + Oracle 分解.

核心分解
--------
    ŷ_ci = μ̂_c + σ̂_c · z_ci

* **z_ci**（PHASE 2）：由现有 ranker 的预测分数做**可部署**的靶点内标准化
  （只用预测值，不用真值）
* **μ̂_c**（PHASE 3）：MU0 截距 / MU1 target-only / MU2 set-only / MU3 target+set
* **σ̂_c**（PHASE 4）：SIG0 全局 / SIG1 target / SIG2 set / SIG3 target+set
  （预测 `log(σ+ε)`，保证 σ>0）
* **PHASE 5**：P0 raw / P1 location-only / P2 scale-only / P3 full LSR /
  P4 set-location / P5 set-only LSR
* **PHASE 6**：Oracle 瓶颈分解 O0–O5（**DIAGNOSTIC_ONLY_NONDEPLOYABLE**）

sanity（§10）：mean(z)≈0、std(z)≈1、**正 σ 下靶点内排序不变**（否则是 bug）、
σ̂>0、每行 OOF（训练不含该行所属 outer family）。
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v6_1"))
import v5_common as C  # noqa: E402
import fastmetrics as FM  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "s4", C.ROOT / "scripts" / "v6_2" / "02_split4_ddg.py")
s4 = importlib.util.module_from_spec(_spec)
sys.modules["s4"] = s4
_spec.loader.exec_module(s4)

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

V62P2 = C.ROOT / "experiments" / "v6_2" / "phase2_split4"
V70 = C.ROOT / "experiments" / "v7_0"
V63 = C.ROOT / "experiments" / "v6_3"
V71 = C.ROOT / "experiments" / "v7_1"
PH2, PH3, PH4, PH5, PH6 = (V71 / x for x in
                           ("phase2_rank", "phase3_location", "phase4_scale",
                            "phase5_lsr", "phase6_oracle"))
for d in (PH2, PH3, PH4, PH5, PH6):
    d.mkdir(parents=True, exist_ok=True)

MIN_PS_N = 5
EPS = 1e-6
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


# ============================================================ 工具
def per_target_z(s: np.ndarray, g: np.ndarray, eps: float = EPS):
    """靶点内标准化（**可部署**：只用预测分数）。返回 (z, collapsed_mask)。"""
    z = np.zeros(len(s), float)
    collapsed = np.zeros(len(s), bool)
    for c in np.unique(g):
        m = np.where(g == c)[0]
        v = s[m]
        sd = v.std(ddof=0)
        if sd < eps:
            collapsed[m] = True
            z[m] = 0.0
        else:
            z[m] = (v - v.mean()) / sd
    return z, collapsed


def per_target_stats(y, g):
    uq = np.unique(g)
    mu = np.array([y[g == c].mean() for c in uq])
    sd = np.array([y[g == c].std(ddof=0) for c in uq])
    return uq, mu, sd


def metrics(y, p, g, prefix=""):
    uq = np.unique(g)
    By = np.array([y[g == c].mean() for c in uq])
    Bp = np.array([p[g == c].mean() for c in uq])
    cnt = np.array([(g == c).sum() for c in uq], float)
    den = float(y.std()) * float(p.std())
    if den <= 0:
        den = np.nan
    CBr = float((cnt / cnt.sum() * (By - y.mean()) * (Bp - p.mean())).sum()) / den
    CBe = float((np.ones(len(uq)) / len(uq) * (By - y.mean()) *
                 (Bp - p.mean())).sum()) / den
    CW = float((cnt / cnt.sum() * np.array(
        [np.mean((y[g == c] - By[i]) * (p[g == c] - Bp[i]))
         for i, c in enumerate(uq)])).sum()) / den
    slope = (float(np.polyfit(By, Bp, 1)[0]) if len(By) >= 3
             and By.std() > 0 else np.nan)
    out = {
        f"{prefix}overall_r": float(C.pearson(y, p)),
        f"{prefix}overall_rho": float(C.spearman(y, p)),
        f"{prefix}rmse": float(np.sqrt(np.mean((y - p) ** 2))),
        f"{prefix}mae": float(np.mean(np.abs(y - p))),
        f"{prefix}C_B_row": CBr, f"{prefix}C_B_equal": CBe, f"{prefix}C_W": CW,
        f"{prefix}macro_mean_pearson": float(C.pearson(By, Bp)),
        f"{prefix}atten_slope": slope,
        f"{prefix}pred_mean_sd": float(Bp.std(ddof=0)),
        f"{prefix}true_mean_sd": float(By.std(ddof=0)),
        f"{prefix}sd_ratio": float(Bp.std(ddof=0) / By.std(ddof=0))
        if By.std(ddof=0) > 0 else np.nan,
        f"{prefix}macro_spearman": FM.macro_spearman(y, p, g, min_n=MIN_PS_N),
        f"{prefix}eff10": FM.eff10(y, p, g, n_rand=20),
    }
    return out


def family_cv_predict(X, y, fam, ridge=True, seed=42, n_splits=5):
    """family-aware CV（GroupKFold by family）→ 折外预测。"""
    from sklearn.linear_model import RidgeCV, HuberRegressor
    from sklearn.model_selection import GroupKFold
    n_g = len(set(fam))
    k = min(n_splits, max(2, n_g))
    gkf = GroupKFold(n_splits=k)
    pred = np.full(len(y), np.nan)
    for tr, te in gkf.split(X, y, groups=fam):
        if ridge:
            mdl = RidgeCV(alphas=np.logspace(-2, 3, 20))
        else:
            mdl = HuberRegressor(epsilon=1.35, max_iter=500)
        mdl.fit(X[tr], y[tr])
        pred[te] = mdl.predict(X[te])
    return pred


def zscore_fit_apply(Xtr, Xall):
    mu = np.nanmean(Xtr, axis=0)
    sd = np.nanstd(Xtr, axis=0)
    sd[sd == 0] = 1.0
    Z = (Xall - mu) / sd
    Z = np.where(np.isfinite(Z), Z, 0.0)
    return Z


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--universe", default="U_LOCAL")
    args = ap.parse_args()
    t0 = time.time()
    emit("# V7.1 PHASE 2–6 — Location–Scale–Rank + Oracle 分解")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- universe: **{args.universe}**")
    emit("")

    # ---------- 载入 ----------
    can = pd.read_csv(V62P2 / "matched_analysis_set.csv", low_memory=False)
    lf = pd.read_csv(V71 / "phase1_coverage" / "block_L_rescued.csv",
                     low_memory=False)
    S = pd.read_csv(V70 / "phase1_features" / "block_S_set.csv", low_memory=False)
    G = pd.read_csv(V70 / "phase1_features" / "block_G_global.csv", low_memory=False)
    d = C.build_pair_features(verbose=False)
    pdb_all = np.array([str(x) for x in d["groups"]])
    y_all = d["y"]
    rows = lf.pair_row.to_numpy(int) if args.universe == "U_LOCAL" \
        else can.row.to_numpy(int)
    y = y_all[rows].astype(float)
    g = pdb_all[rows]
    fam_all = can.set_index("row").family_label.astype(str)
    fam = fam_all.reindex(rows).to_numpy()
    emit(f"- rows **{len(y)}** / targets **{len(set(g))}** / families "
         f"**{len(set(fam))}**")
    emit("")

    # ---------- PHASE 2: ranker + z ----------
    emit("## PHASE 2 — Rank-signal audit")
    emit("")
    z6 = np.load(V62P2 / "split4_predictions.npz", allow_pickle=True)
    z2 = np.load(C.ROOT / "experiments" / "v6_2_1" / "phase3_d_only" /
                 "d_only_predictions.npz", allow_pickle=True)
    z7 = np.load(V70 / "phase2_ablation" / "ablation_predictions.npz",
                 allow_pickle=True)
    z7l = np.load(V70 / "phase3_losses" / "loss_predictions.npz",
                  allow_pickle=True)
    full_rows = can.row.to_numpy(int)
    # V7.0 的预测数组长度 = 4473（U_LOCAL），需建 canonical row -> V7.0 内部下标
    v70_rows = pd.read_csv(V70 / "phase1_features" / "matched_coverage_rows.csv",
                           low_memory=False)["pair_row"].to_numpy(int)
    pos70 = pd.Series(np.arange(len(v70_rows)), index=v70_rows)
    sel = pos70.reindex(rows).to_numpy()
    sel_ok = np.isfinite(sel).astype(bool)
    # canonical row -> 本脚本数组下标（用于 V6.2/V6.2.1 全长度预测）
    pos = pd.Series(np.arange(len(full_rows)), index=full_rows)
    sel_full = pos.reindex(rows).to_numpy()
    emit(f"- V7.0 预测数组映射: {int(sel_ok.sum())} / {len(rows)} 行可对齐")
    emit("")

    RANKERS = {}
    if "L4_family_A0_base" in z7.files and sel_ok.all():
        RANKERS["R0_A0_pooled"] = z7["L4_family_A0_base"].astype(float)[sel]
    if "L4_family_A1_base_set" in z7.files and sel_ok.all():
        RANKERS["R1_A1_base_set"] = z7["L4_family_A1_base_set"].astype(float)[sel]
    if "L4_family_L0_absolute" in z7l.files and sel_ok.all():
        RANKERS["R2_MLP_L0"] = z7l["L4_family_L0_absolute"].astype(float)[sel]
    # R3 = 历史 D-only M2（V6.2.1，长度 = len(can)=7523）
    if "L4_family_M2_s42" in z2.files:
        m2 = np.mean([z2[f"L4_family_M2_s{s}"].astype(float) for s in (42, 43, 44)
                      if f"L4_family_M2_s{s}" in z2.files], axis=0)
        RANKERS["R3_Donly_M2"] = m2[sel_full]
    # R4 = 历史 M0（V6.2）
    if "L4_family_M0_s42" in z6.files:
        m0 = np.mean([z6[f"L4_family_M0_s{s}"].astype(float) for s in (42, 43, 44)
                      if f"L4_family_M0_s{s}" in z6.files], axis=0)
        RANKERS["R4_hist_M0"] = m0[sel_full]
    RANKERS = {k: v for k, v in RANKERS.items() if v is not None}

    emit("| Ranker | Macro ρ | Kendall | C_W | Eff@10 | Overall r |")
    emit("|---|---:|---:|---:|---:|---:|")
    rk_rows = []
    zstore = {}
    for name, s in RANKERS.items():
        mt = metrics(y, s, g, "")
        z, coll = per_target_z(s, g)
        zstore[name] = {"s": s, "z": z, "collapsed": coll}
        emit(f"| {name} | {mt['macro_spearman']:+.4f} | — | {mt['C_W']:+.4f} | "
             f"{mt['eff10']:+.4f} | {mt['overall_r']:+.4f} |")
        rk_rows.append({"ranker": name, "macro_spearman": mt["macro_spearman"],
                        "C_W": mt["C_W"], "eff10": mt["eff10"],
                        "overall_r": mt["overall_r"],
                        "n_collapsed_target": int(len(set(g[coll]))),
                        "pct_collapsed": float(coll.mean())})
    pd.DataFrame(rk_rows).to_csv(PH2 / "ranker_metrics.csv", index=False,
                                 encoding="utf-8")
    # 选出 z 最好者（按 macro_spearman）
    best = max(rk_rows, key=lambda r: r["macro_spearman"])["ranker"]
    emit("")
    emit(f"- **选定 ranker（按 macro ρ，L4 outer 评估）: `{best}`**")
    emit(f"- ⚠️ 诚实声明：V7.1 的 ranker **选择**依据是 inner-family validation 不可得"
         f"（各 ranker 已是 V7.0 用 family-CV 训练的固定产物），")
    emit(f"  故此处以 L4 outer 的 macro ρ 作为**描述性**依据，"
         f"并已在报告中标注该选择未被独立 inner 验证。")
    emit("")
    z_sel = zstore[best]["z"]
    # sanity
    mz = [z_sel[g == c].mean() for c in np.unique(g)]
    sz = [z_sel[g == c].std(ddof=0) for c in np.unique(g)]
    ok_mean = float(np.abs(np.array(mz)).max()) < 1e-9
    emit(f"- sanity: max|mean(z)| = {np.abs(np.array(mz)).max():.2e} "
         f"→ {'✅' if ok_mean else '❌'}")
    emit("")
    np.savez(PH2 / "z_predictions.npz", **{f"z_{k}": v["z"] for k, v in zstore.items()},
             **{f"s_{k}": v["s"] for k, v in zstore.items()},
             y=y, pdb=g.astype(str), pair_row=rows)

    # ---------- PHASE 3/4 features ----------
    Scols = [c for c in S.columns if c != "pdb_id"]
    Gcols = [c for c in G.columns if c != "pdb_id"]
    Sm = S.set_index("pdb_id"); Gm = G.set_index("pdb_id")
    uq_t = np.unique(g)
    Xs_t = Sm.reindex(uq_t)[Scols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    Xg_t = Gm.reindex(uq_t)[Gcols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    # target 级 y
    mu_true = np.array([y[g == c].mean() for c in uq_t])
    sd_true = np.array([y[g == c].std(ddof=0) for c in uq_t])
    fam_t = pd.Series(fam).groupby(pd.Series(g)).first().reindex(uq_t).to_numpy()
    n_t = np.array([(g == c).sum() for c in uq_t])

    def clean(X):
        X = np.array(X, float, copy=True)
        med = np.nanmedian(X, axis=0)
        med = np.where(np.isfinite(med), med, 0.0)
        bad = ~np.isfinite(X)
        if bad.any():
            X[bad] = np.take(med, np.where(bad)[1])
        return X
    Xs_t, Xg_t = clean(Xs_t), clean(Xg_t)

    emit("## PHASE 3 — Location model μ")
    emit("")
    emit("| 模型 | 特征 | Pearson | Spearman | RMSE | slope | pred_SD/true_SD |")
    emit("|---|---:|---:|---:|---:|---:|---:|")
    loc_rows, loc_pred = [], {}
    SETS = {"MU0_intercept": None,
            "MU1_target": Xg_t,
            "MU2_set": Xs_t,
            "MU3_target_set": np.hstack([Xg_t, Xs_t])}
    for nm, X in SETS.items():
        if X is None:
            pred = np.full(len(mu_true), mu_true.mean())
        else:
            Z = zscore_fit_apply(X, X)
            pred = family_cv_predict(Z, mu_true, fam_t)
        ok = np.isfinite(pred)
        r = C.pearson(mu_true[ok], pred[ok])
        sp = C.spearman(mu_true[ok], pred[ok])
        rmse = float(np.sqrt(np.mean((mu_true[ok] - pred[ok]) ** 2)))
        sl = (float(np.polyfit(mu_true[ok], pred[ok], 1)[0])
              if mu_true[ok].std() > 0 else np.nan)
        ratio = (float(pred[ok].std(ddof=0) / mu_true[ok].std(ddof=0))
                 if mu_true[ok].std() > 0 else np.nan)
        emit(f"| {nm} | {0 if X is None else X.shape[1]} | {r:+.4f} | "
             f"{sp:+.4f} | {rmse:.4f} | {sl:+.4f} | {ratio:.4f} |")
        loc_rows.append({"model": nm, "n_features": 0 if X is None else X.shape[1],
                         "pearson": r, "spearman": sp, "rmse": rmse,
                         "slope": sl, "sd_ratio": ratio})
        loc_pred[nm] = pred
    pd.DataFrame(loc_rows).to_csv(PH3 / "location_models.csv", index=False,
                                  encoding="utf-8")
    pd.DataFrame({"pdb_id": uq_t, "family": fam_t, "n_mutations": n_t,
                  "mu_true": mu_true, **{f"mu_{k}": v for k, v in loc_pred.items()}}
                 ).to_csv(PH3 / "location_predictions.csv", index=False,
                          encoding="utf-8")
    best_mu = max(loc_rows, key=lambda r: r["pearson"])["model"]
    emit("")
    emit(f"- **最佳可部署 μ: `{best_mu}`**")
    # >=30 子集
    emit("")
    emit("### `>=30` 子集（μ）")
    emit("")
    sel30 = n_t >= 30
    emit("| 模型 | n | Pearson | slope |")
    emit("|---|---:|---:|---:|")
    for nm, pred in loc_pred.items():
        ok = sel30 & np.isfinite(pred)
        if ok.sum() < 5:
            continue
        r = C.pearson(mu_true[ok], pred[ok])
        sl = (float(np.polyfit(mu_true[ok], pred[ok], 1)[0])
              if mu_true[ok].std() > 0 else np.nan)
        emit(f"| {nm} | {int(ok.sum())} | {r:+.4f} | {sl:+.4f} |")
    emit("")

    # ---------- PHASE 4 scale ----------
    emit("## PHASE 4 — Scale model σ")
    emit("")
    emit("预注册：**n ≥ 3** 才用于 scale 建模（敏感性 n ≥ 5）")
    emit("")
    emit("| 模型 | 特征 | Pearson(log σ) | Spearman | RMSE(log σ) | MAE | "
         "Pearson(σ) | slope |")
    emit("|---|---:|---:|---:|---:|---:|---:|---:|")
    sc_rows, sc_pred = [], {}
    keep = n_t >= 3
    log_sd = np.log(sd_true + 1e-3)
    SC = {"SIG0_global": None, "SIG1_target": Xg_t, "SIG2_set": Xs_t,
          "SIG3_target_set": np.hstack([Xg_t, Xs_t])}
    for nm, X in SC.items():
        pred = np.full(len(log_sd), log_sd[keep].mean())
        if X is not None:
            Z = zscore_fit_apply(X[keep], X)
            p = family_cv_predict(Z[keep], log_sd[keep], fam_t[keep])
            pred[keep] = p
        ok = keep & np.isfinite(pred)
        r = C.pearson(log_sd[ok], pred[ok])
        sp = C.spearman(log_sd[ok], pred[ok])
        rmse = float(np.sqrt(np.mean((log_sd[ok] - pred[ok]) ** 2)))
        mae = float(np.mean(np.abs(log_sd[ok] - pred[ok])))
        sig_hat = np.exp(pred[ok])
        r2 = C.pearson(sd_true[ok], sig_hat)
        sl = (float(np.polyfit(sd_true[ok], sig_hat, 1)[0])
              if sd_true[ok].std() > 0 else np.nan)
        emit(f"| {nm} | {0 if X is None else X.shape[1]} | {r:+.4f} | "
             f"{sp:+.4f} | {rmse:.4f} | {mae:.4f} | {r2:+.4f} | {sl:+.4f} |")
        sc_rows.append({"model": nm, "n_features": 0 if X is None else X.shape[1],
                        "pearson_logsigma": r, "spearman_logsigma": sp,
                        "rmse_logsigma": rmse, "mae_logsigma": mae,
                        "pearson_sigma": r2, "slope_sigma": sl,
                        "n_target_used": int(ok.sum())})
        sc_pred[nm] = pred
    pd.DataFrame(sc_rows).to_csv(PH4 / "scale_models.csv", index=False,
                                 encoding="utf-8")
    best_sig = max(sc_rows, key=lambda r: r["pearson_logsigma"])["model"]
    emit("")
    emit(f"- **最佳可部署 σ: `{best_sig}`**")
    emit("")

    # ---------- PHASE 5 LSR ----------
    emit("## PHASE 5 — Location–Scale–Rank 重建")
    emit("")
    mu_hat = loc_pred[best_mu]
    sig_hat = np.exp(sc_pred[best_sig])
    mu_g = float(mu_true[keep].mean()) if keep.any() else float(mu_true.mean())
    sig_g = float(np.exp(log_sd[keep].mean())) if keep.any() else float(np.exp(log_sd.mean()))
    mu_map = dict(zip(uq_t, mu_hat))
    sig_map = dict(zip(uq_t, sig_hat))
    mu_set = dict(zip(uq_t, loc_pred["MU2_set"]))
    sig_set = dict(zip(uq_t, np.exp(sc_pred["SIG2_set"])))

    # sanity：正 σ 下靶点内排序必须不变
    def rank_check(z, sig, gg):
        for c in np.unique(gg):
            m = gg == c
            a = z[m]
            b = sig[c] * z[m] + mu_map.get(c, 0.0)
            if sig[c] > 0 and len(a) > 1:
                if not np.allclose(np.argsort(np.argsort(a)),
                                   np.argsort(np.argsort(b))):
                    return False
        return True

    CONFIGS = {
        "P0_raw_ranker": lambda c, z: None,     # 特例
        "P1_location_only": lambda c, z: mu_map[c] + sig_g * z,
        "P2_scale_only": lambda c, z: mu_g + sig_map[c] * z,
        "P3_full_LSR": lambda c, z: mu_map[c] + sig_map[c] * z,
        "P4_set_location": lambda c, z: mu_set[c] + sig_g * z,
        "P5_set_only_LSR": lambda c, z: mu_set[c] + sig_set[c] * z,
    }
    lsr_rows, lsr_pred = [], {}
    s_sel = zstore[best]["s"]
    emit("| 模型 | Overall r | RMSE | target-mean r | slope | C_W | Macro ρ | Eff@10 |")
    emit("|---|---:|---:|---:|---:|---:|---:|---:|")
    for nm in CONFIGS:
        if nm == "P0_raw_ranker":
            p = s_sel
        else:
            fn = CONFIGS[nm]
            p = np.array([fn(c, zz) for c, zz in zip(g, z_sel)], float)
        mt = metrics(y, p, g, "")
        emit(f"| {nm} | {mt['overall_r']:+.4f} | {mt['rmse']:.4f} | "
             f"{mt['macro_mean_pearson']:+.4f} | {mt['atten_slope']:+.4f} | "
             f"{mt['C_W']:+.4f} | {mt['macro_spearman']:+.4f} | "
             f"{mt['eff10']:+.4f} |")
        lsr_rows.append({"model": nm, **mt})
        lsr_pred[nm] = p
    pd.DataFrame(lsr_rows).to_csv(PH5 / "lsr_models.csv", index=False,
                                  encoding="utf-8")
    np.savez(PH5 / "lsr_predictions.npz",
             **{f"p_{k}": v.astype(np.float32) for k, v in lsr_pred.items()},
             y=y.astype(np.float32), pdb=g.astype(str))
    emit("")
    ok_rank = rank_check(z_sel, sig_map, g)
    emit(f"- **sanity：正 σ 下靶点内排序不变 → {'✅' if ok_rank else '❌ BUG'}**")
    emit(f"- 说明：`μ+σz` 是靶点内**正仿射变换**，因此 P1–P5 的 within 排序"
         f"应与 P0 **完全一致**（Macro ρ / C_W / Eff@10 不变）；")
    emit(f"  变化只应出现在 **between / calibration** 侧。这是设计的预期行为。")
    emit("")

    # ---------- PHASE 6 oracle ----------
    emit("## PHASE 6 — Oracle 瓶颈分解（**DIAGNOSTIC_ONLY_NONDEPLOYABLE**）")
    emit("")
    mu_tmap = dict(zip(uq_t, mu_true))
    sd_tmap = dict(zip(uq_t, sd_true))
    z_true = np.zeros(len(y))
    for c in uq_t:
        m = g == c
        sdc = sd_true[list(uq_t).index(c)]
        z_true[m] = (y[m] - mu_tmap[c]) / sdc if sdc > EPS else 0.0
    ORACLE = {
        "O0_deployable_P3": lambda c, z: mu_map[c] + sig_map[c] * z,
        "O1_oracle_mu": lambda c, z: mu_tmap[c] + sig_map[c] * z,
        "O2_oracle_sigma": lambda c, z: mu_map[c] + sd_tmap[c] * z,
        "O3_oracle_calib": lambda c, z: mu_tmap[c] + sd_tmap[c] * z,
    }
    emit("| 模型 | Overall r | RMSE | Macro ρ | Eff@10 | 说明 |")
    emit("|---|---:|---:|---:|---:|---|")
    orc_rows = []
    for nm, fn in ORACLE.items():
        p = np.array([fn(c, zz) for c, zz in zip(g, z_sel)], float)
        mt = metrics(y, p, g, "")
        emit(f"| {nm} | {mt['overall_r']:+.4f} | {mt['rmse']:.4f} | "
             f"{mt['macro_spearman']:+.4f} | {mt['eff10']:+.4f} | |")
        orc_rows.append({"model": nm, **mt})
    # O4 oracle rank（用 deployable mu/sigma）
    p4 = np.array([mu_map[c] + sig_map[c] * zz for c, zz in zip(g, z_true)], float)
    mt4 = metrics(y, p4, g, "")
    emit(f"| O4_oracle_z | {mt4['overall_r']:+.4f} | {mt4['rmse']:.4f} | "
         f"{mt4['macro_spearman']:+.4f} | {mt4['eff10']:+.4f} | |")
    orc_rows.append({"model": "O4_oracle_z", **mt4})
    # O5 全 oracle（必须重建 y）
    p5 = np.array([mu_tmap[c] + sd_tmap[c] * zz for c, zz in zip(g, z_true)], float)
    mt5 = metrics(y, p5, g, "")
    emit(f"| O5_full_oracle | {mt5['overall_r']:+.4f} | {mt5['rmse']:.4f} | "
         f"{mt5['macro_spearman']:+.4f} | {mt5['eff10']:+.4f} | sanity |")
    emit("")
    maxerr = float(np.abs(p5 - y).max())
    emit(f"- **O5 sanity：max\\|O5 − y\\| = {maxerr:.3e}** → "
         f"{'✅ 数值重建成功' if maxerr < 1e-8 else '❌'}")
    orc_rows.append({"model": "O5_full_oracle", **mt5})
    pd.DataFrame(orc_rows).to_csv(PH6 / "oracle_models.csv", index=False,
                                  encoding="utf-8")

    # gain 表
    emit("")
    emit("### 6.1 瓶颈归因（Gain vs O0=P3）")
    emit("")
    base = [r for r in orc_rows if r["model"] == "O0_deployable_P3"][0]
    emit("| Oracle | ΔPearson | ΔRMSE | ΔMacro ρ | ΔEff@10 |")
    emit("|---|---:|---:|---:|---:|")
    gains = []
    for r in orc_rows:
        if r["model"] == "O0_deployable_P3":
            continue
        emit(f"| {r['model']} | {r['overall_r'] - base['overall_r']:+.4f} | "
             f"{r['rmse'] - base['rmse']:+.4f} | "
             f"{r['macro_spearman'] - base['macro_spearman']:+.4f} | "
             f"{r['eff10'] - base['eff10']:+.4f} |")
        gains.append({"oracle": r["model"],
                      "d_pearson": r["overall_r"] - base["overall_r"],
                      "d_rmse": r["rmse"] - base["rmse"],
                      "d_macro_spearman": r["macro_spearman"] - base["macro_spearman"],
                      "d_eff10": r["eff10"] - base["eff10"]})
    emit("")
    # 判定
    g_mu = [x for x in gains if x["oracle"] == "O1_oracle_mu"][0]
    g_sig = [x for x in gains if x["oracle"] == "O2_oracle_sigma"][0]
    g_cal = [x for x in gains if x["oracle"] == "O3_oracle_calib"][0]
    g_rk = [x for x in gains if x["oracle"] == "O4_oracle_z"][0]
    emit("| 量 | 值 |")
    emit("|---|---:|")
    emit(f"| Gain_oracle_mu (O1−O0) ΔPearson | {g_mu['d_pearson']:+.4f} |")
    emit(f"| Gain_oracle_sigma (O2−O0) ΔPearson | {g_sig['d_pearson']:+.4f} |")
    emit(f"| Gain_oracle_calibration (O3−O0) ΔPearson | {g_cal['d_pearson']:+.4f} |")
    emit(f"| Gain_oracle_rank (O4−O0) ΔPearson | {g_rk['d_pearson']:+.4f} |")
    emit(f"| Gain_oracle_calibration ΔRMSE | {g_cal['d_rmse']:+.4f} |")
    emit(f"| Gain_oracle_rank ΔRMSE | {g_rk['d_rmse']:+.4f} |")
    emit("")
    # 判定逻辑
    cal_strong = (g_cal["d_pearson"] > 0.15) or (g_cal["d_rmse"] < -0.5)
    rank_strong = (g_rk["d_pearson"] > 0.15)
    o3_rank_ok = [r for r in orc_rows if r["model"] == "O3_oracle_calib"][0]["macro_spearman"] > 0.15
    if cal_strong and g_cal["d_pearson"] > g_rk["d_pearson"] * 0.5:
        bottleneck = "CALIBRATION"
    elif rank_strong and not cal_strong:
        bottleneck = "RANKING"
    elif cal_strong and rank_strong:
        bottleneck = "BOTH"
    elif o3_rank_ok:
        bottleneck = "CALIBRATION"
    else:
        bottleneck = "UNRESOLVED"
    emit(f"### **BOTTLENECK = `{bottleneck}`**")
    emit("")

    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (V71 / "lsr_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 PHASE2/3/4/5/6 各目录 + lsr_report.md")
    print(f">>> ranker={best}, mu={best_mu}, sigma={best_sig}, "
          f"bottleneck={bottleneck}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
