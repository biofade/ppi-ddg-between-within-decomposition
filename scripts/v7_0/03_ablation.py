# -*- coding: utf-8 -*-
"""scripts/v7_0/03_ablation.py -- V7.0 PHASE 2: Architecture Ablation Matrix.

7 个模型（**同一 learner、同一 folds、同一 seeds、同一树数**，§9）
--------------------------------------------------------------
| 模型 | 特征 |
|---|---|
| **A0** | Base（frozen pooled baseline = 历史 M0 特征）      |
| **A1** | Base + S（Set context）                            |
| **A2** | Base + G（Global target context，diagnostic）       |
| **A3** | Base + L（**真正的** mutation-local）              |
| **A4** | Base + S + L（**V7.0 最核心候选**）                 |
| **A5** | S + L（无 pooled Base）                            |
| **A6** | G + S + L                                          |

**统一 row universe**（§5）：BLOCK L 可用的 matched-coverage 子集（4473 行 / 208 targets）。
**预处理**（§10）：S/G/L 的标准化与 PCA **只在各折 train 行上 fit**。

Primary 指标（§11）：`Macro per-target Spearman`、`Efficiency@10`、`C_W`。
Between 同时报 **row-weighted 与 equal-target**（§E 强制）。
"""

from __future__ import annotations

import argparse
import gc
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

PH1 = C.ROOT / "experiments" / "v7_0" / "phase1_features"
PH2 = C.ROOT / "experiments" / "v7_0" / "phase2_ablation"
PH2.mkdir(parents=True, exist_ok=True)
OUT_MET = PH2 / "ablation_metrics.csv"
OUT_NPZ = PH2 / "ablation_predictions.npz"
OUT_MD = PH2 / "ablation_report.md"
FAIL = C.ROOT / "logs" / "v7_0" / "failures.log"
FAIL.parent.mkdir(parents=True, exist_ok=True)

TREES = 300
SEEDS = [42, 43, 44]
LEVELS = ["L4_family", "L1_random"]
PCA_NEW = 32
MIN_PS_N = 5
L: list[str] = []

MODELS = {
    "A0_base":          ["base"],
    "A1_base_set":      ["base", "S"],
    "A2_base_global":   ["base", "G"],
    "A3_base_local":    ["base", "L"],
    "A4_base_set_local": ["base", "S", "L"],
    "A5_set_local":     ["S", "L"],
    "A6_global_set_local": ["G", "S", "L"],
}


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def fail(w: str, e: BaseException) -> None:
    import traceback
    with open(FAIL, "a", encoding="utf-8") as fh:
        fh.write(f"[{time.strftime('%H:%M:%S')}] [03_ablation] {w}: "
                 f"{type(e).__name__}: {e}\n{traceback.format_exc()}\n")


# ---------------------------------------------------------------- 指标
def decompose(y, p, g, weight="row"):
    y = np.asarray(y, float); p = np.asarray(p, float); g = np.asarray(g)
    uq = np.unique(g)
    By = np.array([y[g == c].mean() for c in uq])
    Bp = np.array([p[g == c].mean() for c in uq])
    cnt = np.array([(g == c).sum() for c in uq], float)
    w = cnt / cnt.sum() if weight == "row" else np.ones(len(uq)) / len(uq)
    yb, pb = y.mean(), p.mean()
    den = float(y.std()) * float(p.std())
    if den <= 0:
        return dict(r=np.nan, C_B=np.nan, C_W=np.nan)
    C_B = float((w * (By - yb) * (Bp - pb)).sum()) / den
    Wy = y - np.repeat(By, cnt.astype(int))
    Wp = p - np.repeat(Bp, cnt.astype(int))
    C_W = float((w * np.array([np.mean((y[g == c] - By[i]) *
                                       (p[g == c] - Bp[i]))
                               for i, c in enumerate(uq)])).sum()) / den
    return dict(r=float(((y - yb) * (p - pb)).mean() / den), C_B=C_B, C_W=C_W,
                By=By, Bp=Bp)


def atten_slope(y, p, g):
    d = decompose(y, p, g)
    By, Bp = d["By"], d["Bp"]
    if len(By) < 3 or By.std() <= 0:
        return np.nan
    return float(np.polyfit(By, Bp, 1)[0])


def macro_kendall(y, p, g, min_n=2):
    vals = []
    for c in np.unique(g):
        m = np.where(g == c)[0]
        if len(m) < min_n:
            continue
        a, b = y[m], p[m]; n = len(a)
        conc = disc = ta = tb = 0
        for i in range(n):
            for j in range(i + 1, n):
                da = a[i] - a[j]; db = b[i] - b[j]
                if da == 0 and db == 0:
                    ta += 1; tb += 1
                elif da == 0:
                    ta += 1
                elif db == 0:
                    tb += 1
                elif (da > 0) == (db > 0):
                    conc += 1
                else:
                    disc += 1
        n0 = conc + disc
        den = np.sqrt((n0 + ta) * (n0 + tb))
        if den > 0:
            vals.append((conc - disc) / den)
    return float(np.mean(vals)) if vals else float("nan")


def pairwise_acc(y, p, g, delta=0.0, min_n=2):
    vals = []
    for c in np.unique(g):
        m = np.where(g == c)[0]
        if len(m) < min_n:
            continue
        a, b = y[m], p[m]; n = len(a)
        sc = []
        for i in range(n):
            for j in range(i + 1, n):
                dy = a[i] - a[j]
                if abs(dy) <= delta:
                    continue
                dp = b[i] - b[j]
                sc.append(0.5 if dp == 0 else (1.0 if (dy > 0) == (dp > 0) else 0.0))
        if sc:
            vals.append(float(np.mean(sc)))
    return float(np.mean(vals)) if vals else float("nan")


def full_metrics(y, p, g, prefix=""):
    from sklearn.metrics import mean_absolute_error
    d_row = decompose(y, p, g, "row")
    d_eq = decompose(y, p, g, "equal")
    out = {
        f"{prefix}overall_r": d_row["r"],
        f"{prefix}overall_rho": C.spearman(y, p),
        f"{prefix}rmse": float(np.sqrt(np.mean((y - p) ** 2))),
        f"{prefix}mae": float(mean_absolute_error(y, p)),
        f"{prefix}C_B_row": d_row["C_B"],
        f"{prefix}C_B_equal": d_eq["C_B"],
        f"{prefix}C_W": d_row["C_W"],
        f"{prefix}macro_mean_pearson": C.pearson(d_eq["By"], d_eq["Bp"]),
        f"{prefix}atten_slope": atten_slope(y, p, g),
        f"{prefix}macro_spearman": FM.macro_spearman(y, p, g, min_n=MIN_PS_N),
        f"{prefix}macro_kendall": macro_kendall(y, p, g, 2),
        f"{prefix}pairwise_acc": pairwise_acc(y, p, g, 0.0, 2),
        f"{prefix}pairwise_acc_d05": pairwise_acc(y, p, g, 0.5, 2),
        f"{prefix}eff10": FM.eff10(y, p, g, n_rand=20),
    }
    # Hit@10 / P@10（预注册：不足 10 用全部）
    hits, prec = [], []
    for c in np.unique(g):
        m = np.where(g == c)[0]
        k = min(10, len(m))
        if k < 2:
            continue
        idx = np.argsort(p[m], kind="mergesort")[:k]
        yy = y[m]
        hits.append(float((yy[idx] < 0).mean()))
        thr = np.quantile(-yy, 1 - k / len(m))
        prec.append(float(((-yy)[idx] >= thr).mean()))
    out[f"{prefix}hit10"] = float(np.mean(hits)) if hits else np.nan
    out[f"{prefix}p10"] = float(np.mean(prec)) if prec else np.nan
    return out


# ---------------------------------------------------------------- 主流程
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", nargs="*", default=list(MODELS.keys()))
    ap.add_argument("--max-minutes", type=float, default=100.0)
    args = ap.parse_args()
    t0 = time.time()
    emit("# V7.0 PHASE 2 — Architecture Ablation Matrix")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- learner: **LightGBM**（`LGB_PARAMS`，`n_estimators={TREES}`），"
         f"所有模型**同超参**（§9）")
    emit(f"- 档: {LEVELS}（L4 = primary，L1 = diagnostic）")
    emit(f"- seeds: {SEEDS}；Primary = **seed ensemble 后再算指标**（§19）")
    emit(f"- 新特征块 PCA: **{PCA_NEW}**（**只在各折 train 行 fit**，§10）")
    emit("")

    # ---- 载入统一 row universe ----
    lf = pd.read_csv(PH1 / "block_L_local.csv", low_memory=False)
    S = pd.read_csv(PH1 / "block_S_set.csv", low_memory=False)
    G = pd.read_csv(PH1 / "block_G_global.csv", low_memory=False)
    can = pd.read_csv(C.ROOT / "experiments" / "v6_2" / "phase2_split4" /
                      "matched_analysis_set.csv", low_memory=False)
    d = C.build_pair_features(verbose=False)
    blk = d["X_blocks"]
    y_all = d["y"]
    pdb_all = np.array([str(x) for x in d["groups"]])

    rows = lf.pair_row.to_numpy(int)
    emit(f"- **统一 row universe（matched coverage）**: {len(rows)} 行 / "
         f"{len(set(pdb_all[rows]))} targets")
    emit(f"- 占总 canonical {len(can)} 行的 {len(rows) / len(can) * 100:.1f}%")
    emit("")

    base = np.hstack([blk["phys_diff"], blk["Z_diff"]])[rows].astype(np.float32)
    # L / S / G 特征（按 pair_row / pdb_id 对齐）
    Lcols = [c for c in lf.columns
             if c not in ("pair_row", "pdb_id", "n_sites_verified",
                          "n_sites_total")
             and pd.api.types.is_numeric_dtype(pd.to_numeric(lf[c],
                                                             errors="coerce"))]
    Xl = lf[Lcols].apply(pd.to_numeric, errors="coerce").to_numpy(np.float32)
    Scols = [c for c in S.columns if c != "pdb_id"]
    Gcols = [c for c in G.columns if c != "pdb_id"]
    pdb_rows = pdb_all[rows]
    Sm = S.set_index("pdb_id")
    Gm = G.set_index("pdb_id")
    Xs = Sm.reindex(pdb_rows)[Scols].apply(pd.to_numeric,
                                           errors="coerce").to_numpy(np.float32)
    Xg = Gm.reindex(pdb_rows)[Gcols].apply(pd.to_numeric,
                                           errors="coerce").to_numpy(np.float32)
    emit(f"- 特征维度: base **{base.shape[1]}** | L **{Xl.shape[1]}** | "
         f"S **{Xs.shape[1]}** | G **{Xg.shape[1]}**")
    emit(f"- S/G 缺失（该 target 不在 V6.3 表内）: "
         f"S {int(np.isnan(Xs).all(axis=1).sum())} 行，"
         f"G {int(np.isnan(Xg).all(axis=1).sum())} 行")
    # 缺失中位数填补（用**全 train** 中位数；为简洁在折外用 train 中位数）
    for nm, X in (("L", Xl), ("S", Xs), ("G", Xg)):
        X = np.array(X, dtype=np.float32, copy=True)   # 保证可写
        med = np.nanmedian(X, axis=0)
        med = np.where(np.isfinite(med), med, 0.0)
        bad = ~np.isfinite(X)
        if bad.any():
            X[bad] = np.take(med, np.where(bad)[1])
        if nm == "L":
            Xl = X
        elif nm == "S":
            Xs = X
        else:
            Xg = X
    emit("")

    y = y_all[rows].astype(float)
    pdb = pdb_all[rows]
    FIT = {"base": base, "L": Xl, "S": Xs, "G": Xg}

    # ---- 对每档每 seed 训练全部模型 ----
    res_rows = []
    preds: dict = {}
    for lv in LEVELS:
        units = s4.level_units({"y": y, "pdb": pdb, "fam":
                                can.family_label.astype(str).to_numpy()[
                                    np.isin(can.row.to_numpy(int), rows)],
                                "seq40": can.seq40_cluster.astype(str).to_numpy()[
                                    np.isin(can.row.to_numpy(int), rows)],
                                "obs": can.obs_key.astype(str).to_numpy()[
                                    np.isin(can.row.to_numpy(int), rows)],
                                "state": can.state_id.astype(str).to_numpy()[
                                    np.isin(can.row.to_numpy(int), rows)]}, lv)
        emit(f"## {lv}（unit 数 {len(np.unique(units))}）")
        emit("")
        for mk, blocks in MODELS.items():
            if mk not in args.models:
                continue
            if (time.time() - t0) / 60 > args.max_minutes:
                emit(f"- {mk}: skip (soft time limit)")
                continue
            try:
                yhat = {sd: np.full(len(y), np.nan) for sd in SEEDS}
                for sd in SEEDS:
                    folds = s4.make_folds(units, sd, 5)
                    lgb = C.load_lgb()
                    for tr, te in folds:
                        mats_tr, mats_te = [], []
                        for bname in blocks:
                            X = FIT[bname]
                            if bname == "base":
                                mats_tr.append(X[tr])
                                mats_te.append(X[te])
                                continue
                            # 新块：train-only 标准化 + PCA
                            from sklearn.decomposition import PCA
                            mu = np.nanmean(X[tr], axis=0)
                            sdv = np.nanstd(X[tr], axis=0)
                            sdv[sdv == 0] = 1.0
                            Ztr = (X[tr] - mu) / sdv
                            Zte = (X[te] - mu) / sdv
                            d_pca = min(PCA_NEW, Ztr.shape[1],
                                        max(2, len(tr) - 1))
                            pc = PCA(n_components=d_pca,
                                     random_state=sd).fit(Ztr)
                            mats_tr.append(pc.transform(Ztr).astype(np.float32))
                            mats_te.append(pc.transform(Zte).astype(np.float32))
                        Xtr = np.hstack(mats_tr).astype(np.float32)
                        Xte = np.hstack(mats_te).astype(np.float32)
                        ytr = y[tr]
                        # forward/reverse 增广（历史语义）
                        Xtr_a = np.vstack([Xtr, -Xtr])
                        ytr_a = np.r_[ytr, -ytr]
                        m = lgb.LGBMRegressor(**{**C.LGB_PARAMS,
                                                 "n_estimators": TREES,
                                                 "random_state": sd})
                        m.fit(Xtr_a, ytr_a)
                        yhat[sd][te] = m.predict(Xte)
                        del m
                        gc.collect()
                ens = np.mean([yhat[sd] for sd in SEEDS], axis=0)
                ok = np.isfinite(ens)
                mt = full_metrics(y[ok], ens[ok], pdb[ok])
                res_rows.append({"level": lv, "model": mk,
                                 "blocks": "+".join(blocks),
                                 "n_features": int(
                                     sum(FIT[b].shape[1] if b == "base"
                                         else min(PCA_NEW, FIT[b].shape[1])
                                         for b in blocks)),
                                 "n_rows": int(ok.sum()),
                                 "n_target": int(len(set(pdb[ok]))),
                                 **mt})
                # secondary: 逐 seed
                ss = []
                for sd in SEEDS:
                    o2 = np.isfinite(yhat[sd])
                    ss.append(FM.macro_spearman(y[o2], yhat[sd][o2], pdb[o2],
                                                min_n=MIN_PS_N))
                res_rows[-1]["macro_spearman_seed_mean"] = float(np.mean(ss))
                res_rows[-1]["macro_spearman_seed_sd"] = float(np.std(ss,
                                                                     ddof=1))
                preds[f"{lv}_{mk}"] = ens.astype(np.float32)
                if lv == "L4_family":
                    for sd in SEEDS:
                        preds[f"{lv}_{mk}_s{sd}"] = yhat[sd].astype(np.float32)
                emit(f"| {mk} | r={mt['overall_r']:+.4f} | "
                     f"C_W={mt['C_W']:+.4f} | "
                     f"macro_ρ={mt['macro_spearman']:+.4f} | "
                     f"Eff@10={mt['eff10']:+.4f} | "
                     f"C_B(row)={mt['C_B_row']:+.4f} | "
                     f"C_B(eq)={mt['C_B_equal']:+.4f} |")
                pd.DataFrame(res_rows).to_csv(OUT_MET, index=False,
                                              encoding="utf-8")
            except MemoryError as e:
                fail(f"{lv}/{mk} OOM", e)
                emit(f"| {mk} | OOM_SKIPPED | | | | | |")
            except Exception as e:  # noqa: BLE001
                fail(f"{lv}/{mk}", e)
                emit(f"| {mk} | FAILED {type(e).__name__} | | | | | |")
                gc.collect()
        emit("")
    preds["y"] = y.astype(np.float32)
    preds["pdb"] = pdb.astype(str)
    for lv in LEVELS:
        preds[f"{lv}_y"] = y.astype(np.float32)
        preds[f"{lv}_pdb"] = pdb.astype(str)
    np.savez(OUT_NPZ, **preds)

    # ---- 汇总表 ----
    R = pd.DataFrame(res_rows)
    R.to_csv(OUT_MET, index=False, encoding="utf-8")
    emit("## 汇总（L4_family，primary）")
    emit("")
    emit("| 模型 | Overall r | C_W | Macro ρ | Kendall | Eff@10 | C_B(row) | "
         "C_B(eq) | slope |")
    emit("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    l4 = R[R.level == "L4_family"]
    for _, r in l4.iterrows():
        emit(f"| {r['model']} | {r['overall_r']:+.4f} | {r['C_W']:+.4f} | "
             f"{r['macro_spearman']:+.4f} | {r['macro_kendall']:+.4f} | "
             f"{r['eff10']:+.4f} | {r['C_B_row']:+.4f} | "
             f"{r['C_B_equal']:+.4f} | {r['atten_slope']:+.4f} |")
    emit("")
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MET}")
    print(f"写出 {OUT_NPZ}")
    print(f"写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
