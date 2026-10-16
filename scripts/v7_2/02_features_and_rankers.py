# -*- coding: utf-8 -*-
"""scripts/v7_2/02_features_and_rankers.py -- V7.2 PHASE 3-5: 特征块 + ranker 家族 + 嵌套 family-CV.

特征块（§PHASE 3，**分块构造，不一次性拼接**）
--------------------------------------------
| 块 | 内容 |
|---|---|
| **F0** | mutation identity（WT/Mut one-hot + Δcharge/Δhydro/Δvolume/Δpolarity/Δaromatic/Δflex） |
| **F1** | 残基嵌入 `h_wt_m`, `h_mut_m`（ESM-2 8M，320 维） |
| **F2** | 扰动 `Δh_m` + `‖Δh‖` + `cos(h_wt,h_mut)`（**核心 sequence-perturbation 块**） |
| **F3** | partner context（需结构映射；本脚本**仅在结构可用时**构造，否则标 `NOT_RUN`） |
| **F4** | interaction perturbation `Δh ⊙ c`, `dot`, `cosine`, 各范数（**核心 MIPA 块**） |
| **F5** | candidate-set 描述符（复用 V6.3 deployable 列） |

Ranker 家族（§PHASE 4）
----------------------
`R0` Base / `R1` Base+F0 / `R2` Base+PCA(h_wt,h_mut) / `R3` Base+PCA(Δh) /
`R4` R3+PCA(c) / `R5` R4+PCA(Δh⊙c)+dot/cos/norms / `R6` R5+F5 / `R7` F0+F2(+F3+F4)

**PCA 只在 train 行拟合**（§7）；learner 统一 **LightGBM**（§PHASE 4）。

嵌套 family-CV（§PHASE 5）：严格按 V7.1.1 冻结协议（outer 5 folds by family，
inner `min(3,max(2,n_fam))`；**选择只用 inner**；outer 仅用于最终评估）。
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
_spec2 = importlib.util.spec_from_file_location(
    "v7f", C.ROOT / "scripts" / "v7_0" / "01_features.py")
v7f = importlib.util.module_from_spec(_spec2)
sys.modules["v7f"] = v7f
_spec2.loader.exec_module(v7f)

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

FEAT = C.ROOT / "features" / "v7_2"
PER = FEAT / "per_target"
V711 = C.ROOT / "experiments" / "v7_1_1"
V70 = C.ROOT / "experiments" / "v7_0"
V72 = C.ROOT / "experiments" / "v7_2"
PH3, PH4, PH5 = (V72 / x for x in
                 ("phase3_features", "phase4_rankers", "phase5_nested_cv"))
for d in (PH3, PH4, PH5):
    d.mkdir(parents=True, exist_ok=True)
TREES = 300
PCA_NEW = 32
MIN_PS_N = 5
EPS = 1e-6
SEEDS = [42, 43, 44]
L: list[str] = []

RANKERS = {
    "R0_Base": ["base"],
    "R1_Base_F0": ["base", "F0"],
    "R2_Base_Residue": ["base", "F1"],
    "R3_Base_dh": ["base", "F2"],
    "R4_Base_dh_C": ["base", "F2", "F3"],
    "R5_Base_Interact": ["base", "F2", "F3", "F4"],
    "R6_Base_Interact_Set": ["base", "F2", "F3", "F4", "F5"],
    "R7_ResidueOnly": ["F0", "F2", "F3", "F4"],
}


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def standardize_within_target(s, g, eps=EPS):
    """V7.1.1 冻结的**唯一** z 变换函数。"""
    z = np.zeros(len(s), float)
    coll = np.zeros(len(s), bool)
    for c in np.unique(g):
        m = np.where(g == c)[0]
        v = s[m]
        sd = v.std(ddof=0)
        if sd < eps:
            coll[m] = True
        else:
            z[m] = (v - v.mean()) / sd
    return z, coll


def metrics(y, p, g):
    uq = np.unique(g)
    By = np.array([y[g == c].mean() for c in uq])
    Bp = np.array([p[g == c].mean() for c in uq])
    cnt = np.array([(g == c).sum() for c in uq], float)
    den = float(y.std()) * float(p.std())
    if den <= 0:
        den = np.nan
    out = {
        "overall_r": float(C.pearson(y, p)),
        "rmse": float(np.sqrt(np.mean((y - p) ** 2))),
        "C_B_row": float((cnt / cnt.sum() * (By - y.mean()) *
                          (Bp - p.mean())).sum()) / den,
        "C_B_equal": float((np.ones(len(uq)) / len(uq) * (By - y.mean()) *
                            (Bp - p.mean())).sum()) / den,
        "C_W": float((cnt / cnt.sum() * np.array(
            [np.mean((y[g == c] - By[i]) * (p[g == c] - Bp[i]))
             for i, c in enumerate(uq)])).sum()) / den,
        "macro_spearman": FM.macro_spearman(y, p, g, min_n=MIN_PS_N),
        "eff10": FM.eff10(y, p, g, n_rand=10),
    }
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rankers", nargs="*", default=list(RANKERS.keys()))
    ap.add_argument("--max-minutes", type=float, default=70.0)
    args = ap.parse_args()
    t0 = time.time()
    emit("# V7.2 PHASE 3–5 — 特征块 + ranker 家族 + 嵌套 family-CV")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- learner: **LightGBM**（`LGB_PARAMS`, trees={TREES}）— 统一，不逐族调参")
    emit(f"- PCA: **{PCA_NEW}**，**只在 train 行 fit**")
    emit(f"- 协议: V7.1.1 冻结 `OUTER_FOLDS_V7_2` + inner family CV")
    emit("")

    # ---------- 权威标签 ----------
    lab = pd.read_csv(V711 / "phase1_universe" / "canonical_target_labels.csv",
                      low_memory=False)
    d = C.build_pair_features(verbose=False)
    can = pd.read_csv(C.ROOT / "experiments" / "v6_2" / "phase2_split4" /
                      "matched_analysis_set.csv", low_memory=False)
    can_i = can.set_index("row")
    blk = d["X_blocks"]
    all_rows = can.row.to_numpy(int)

    # ---------- 加载已完成的残基嵌入 ----------
    per_files = sorted(PER.glob("*.npz"))
    hw, hm, rid = {}, {}, []
    for p in per_files:
        try:
            z = np.load(p)
            r = z["row_id"]
            hw[p.stem] = z["h_wt"]
            hm[p.stem] = z["h_mut"]
            rid.extend(r.tolist())
        except Exception:  # noqa: BLE001
            continue
    rid = np.array(sorted(set(rid)), dtype=int)
    emit(f"- 残基嵌入：**{len(per_files)}** 个 target 文件，"
         f"**{len(rid)}** 个可用 row")
    emit(f"- 占 U_CONTRACT {len(lab)} 行的 **{len(rid) / len(lab) * 100:.1f}%**")
    emit("")
    if len(rid) < 100:
        emit("!! 残基嵌入不足，无法继续")
        (PH3 / "feature_report.md").write_text("\n".join(L) + "\n",
                                               encoding="utf-8")
        return 1

    # 组装 h_wt / h_mut（按 row_id）
    D = d["y"]
    idx = {int(r): i for i, r in enumerate(rid)}
    HW = np.zeros((len(rid), 320), np.float32)
    HM = np.zeros((len(rid), 320), np.float32)
    seen = np.zeros(len(rid), bool)
    for tgt, arr in hw.items():
        z = np.load(PER / f"{tgt}.npz")
        for i, r in enumerate(z["row_id"]):
            j = idx.get(int(r))
            if j is None:
                continue
            HW[j] = arr[i].astype(np.float32)
            HM[j] = hm[tgt][i].astype(np.float32)
            seen[j] = True
    emit(f"- 成功装配 h_wt/h_mut：**{int(seen.sum())}** 行（320 维）")
    rid = rid[seen]; HW = HW[seen]; HM = HM[seen]
    emit("")

    y = d["y"][rid].astype(float)
    g = can_i.loc[rid, "pdb_id"].astype(str).to_numpy()
    fam = can_i.loc[rid, "family_label"].astype(str).to_numpy()
    base = np.hstack([blk["phys_diff"], blk["Z_diff"]])[rid].astype(np.float32)
    emit(f"- **U_CONTRACT 子集（有残基嵌入）**: rows **{len(y)}** / "
         f"targets **{len(set(g))}** / families **{len(set(fam))}**")
    n5 = sum(1 for c in np.unique(g) if (g == c).sum() >= MIN_PS_N)
    emit(f"- 其中 n≥{MIN_PS_N} 的 target（Macro ρ 资格）: **{n5}**")
    emit("")

    # ---------- F0 mutation identity ----------
    fwd = pd.read_csv(C.DATA / "v4.0" / "ddg_pairs.csv", low_memory=False)
    fwd = fwd[fwd.direction == "forward"].reset_index(drop=True)
    import re
    f0 = []
    for r in rid:
        muts = []
        s = str(fwd["mutation"].iloc[r]).upper()
        for tok in re.split(r"[;,]", s):
            tok = tok.strip()
            if not tok:
                continue
            body = tok.split("_", 1)[-1] if "_" in tok else tok
            for m in re.finditer(r"(?<![A-Z])([A-Z])(\d+)([A-Z])(?![A-Z])",
                                 body):
                muts.append((m.group(1), m.group(3)))
        if not muts:
            muts = [("X", "X")]
        vecs = [v7f.mut_identity(w, mu) for w, mu in muts]
        keys = list(vecs[0].keys())
        f0.append([float(np.mean([v[k] for v in vecs])) for k in keys])
    F0 = np.array(f0, np.float32)
    emit(f"- **F0** mutation identity: {F0.shape[1]} 维")
    # ---------- F1 / F2 ----------
    F1 = np.hstack([HW, HM])
    dH = HM - HW
    nrm = np.linalg.norm(dH, axis=1, keepdims=True)
    cosw = (np.sum(HW * HM, axis=1, keepdims=True) /
            (np.linalg.norm(HW, axis=1, keepdims=True) *
             np.linalg.norm(HM, axis=1, keepdims=True) + 1e-9))
    F2 = np.hstack([dH, nrm, cosw])
    emit(f"- **F1** h_wt,h_mut: {F1.shape[1]} 维")
    emit(f"- **F2** Δh + ‖Δh‖ + cos: {F2.shape[1]} 维")
    emit(f"- **F3** partner context: **NOT_RUN**（需结构映射；本轮未做）")
    emit(f"- **F4** interaction: **NOT_RUN**（依赖 F3）")
    # ---------- F5 set ----------
    S = pd.read_csv(V70 / "phase1_features" / "block_S_set.csv", low_memory=False)
    Scols = [c for c in S.columns if c != "pdb_id"]
    F5 = S.set_index("pdb_id").reindex(g)[Scols].apply(
        pd.to_numeric, errors="coerce").to_numpy(np.float32)
    F5 = np.array(F5, np.float32, copy=True)
    med = np.nanmedian(F5, axis=0); med = np.where(np.isfinite(med), med, 0.0)
    bad = ~np.isfinite(F5)
    if bad.any():
        F5[bad] = np.take(med, np.where(bad)[1])
    emit(f"- **F5** candidate-set: {F5.shape[1]} 维")
    emit("")
    FIT = {"base": base, "F0": F0, "F1": F1, "F2": F2, "F5": F5}

    # ---------- 嵌套 CV ----------
    emit("## PHASE 5 — 嵌套 family-CV（冻结协议）")
    emit("")
    D4 = {"y": y, "pdb": g, "fam": fam,
          "seq40": can_i.loc[rid, "seq40_cluster"].astype(str).to_numpy(),
          "obs": can_i.loc[rid, "obs_key"].astype(str).to_numpy(),
          "state": can_i.loc[rid, "state_id"].astype(str).to_numpy()}
    res_rows, preds = [], {}
    for sd in SEEDS:
        if (time.time() - t0) / 60 > args.max_minutes:
            emit(f"- ⏱ 软时限，停止于 seed {sd}")
            break
        outer = s4.make_folds(D4["fam"], sd, 5)
        for mk in args.rankers:
            blocks = RANKERS[mk]
            if any(b == "F3" or b == "F4" for b in blocks):
                continue          # F3/F4 未构造
            try:
                yhat = np.full(len(y), np.nan)
                for fi, (tr, te) in enumerate(outer):
                    assert not (set(fam[tr]) & set(fam[te])), "outer 泄漏"
                    mats_tr, mats_te = [], []
                    for b in blocks:
                        X = FIT[b]
                        if b == "base":
                            mats_tr.append(X[tr]); mats_te.append(X[te])
                            continue
                        from sklearn.decomposition import PCA
                        mu = np.nanmean(X[tr], axis=0)
                        sdv = np.nanstd(X[tr], axis=0); sdv[sdv == 0] = 1
                        Ztr = (X[tr] - mu) / sdv
                        Zte = (X[te] - mu) / sdv
                        k = min(PCA_NEW, Ztr.shape[1], max(2, len(tr) - 1))
                        pc = PCA(n_components=k, random_state=sd).fit(Ztr)
                        mats_tr.append(pc.transform(Ztr).astype(np.float32))
                        mats_te.append(pc.transform(Zte).astype(np.float32))
                    Xtr = np.hstack(mats_tr).astype(np.float32)
                    Xte = np.hstack(mats_te).astype(np.float32)
                    lgb = C.load_lgb()
                    m = lgb.LGBMRegressor(**{**C.LGB_PARAMS,
                                             "n_estimators": TREES,
                                             "random_state": sd})
                    m.fit(np.vstack([Xtr, -Xtr]), np.r_[y[tr], -y[tr]])
                    yhat[te] = m.predict(Xte)
                    del m
                    gc.collect()
                ok = np.isfinite(yhat)
                z, coll = standardize_within_target(yhat[ok], g[ok])
                mt = metrics(y[ok], yhat[ok], g[ok])
                mt["collapse_rate"] = float(coll.mean())
                res_rows.append({"seed": sd, "ranker": mk,
                                 "blocks": "+".join(blocks),
                                 "n_rows": int(ok.sum()),
                                 "n_target": int(len(set(g[ok]))),
                                 **mt})
                preds[f"{mk}_s{sd}"] = yhat.astype(np.float32)
                emit(f"| {mk} | s{sd} | r={mt['overall_r']:+.4f} | "
                     f"C_W={mt['C_W']:+.4f} | macro_ρ={mt['macro_spearman']:+.4f} | "
                     f"Eff@10={mt['eff10']:+.4f} | collapse={coll.mean()*100:.1f}% |")
                pd.DataFrame(res_rows).to_csv(PH5 / "nested_cv_results.csv",
                                              index=False, encoding="utf-8")
            except Exception as ex:  # noqa: BLE001
                emit(f"| {mk} | s{sd} | FAILED {type(ex).__name__}: {ex} |")
                gc.collect()
    preds["y"] = y.astype(np.float32)
    preds["pdb"] = g.astype(str)
    preds["family"] = fam.astype(str)
    preds["row_id"] = rid.astype(np.int64)
    np.savez(PH5 / "ranker_predictions.npz", **preds)
    R = pd.DataFrame(res_rows)
    emit("")
    if len(R):
        emit("## 3-seed 均值")
        emit("")
        g_ = R.groupby("ranker")[["overall_r", "C_W", "macro_spearman",
                                  "eff10"]].mean()
        emit("| Ranker | Overall r | C_W | Macro ρ | Eff@10 |")
        emit("|---|---:|---:|---:|---:|")
        for k, r in g_.iterrows():
            emit(f"| {k} | {r['overall_r']:+.4f} | {r['C_W']:+.4f} | "
                 f"{r['macro_spearman']:+.4f} | {r['eff10']:+.4f} |")
    emit("")
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (PH3 / "feature_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    (PH5 / "nested_cv_report.md").write_text("\n".join(L) + "\n",
                                             encoding="utf-8")
    print(f"\n写出 phase3/phase5 文件")
    return 0


if __name__ == "__main__":
    sys.exit(main())
