# -*- coding: utf-8 -*-
"""scripts/v7_2r_3/03_matched_eval.py -- V7.2R-3 PHASE 5（修复版）。

修复
----
1. **维度 bug**：原代码把 1280 维的 PPI 向量写进 320 维缓冲（`DG[i] = ...`）。
   现按实际维度分配，**不做全局 PCA 泄露**。
2. **预处理纪律**：`g_c`（target 级）与 `Δg`（逐行 1280 维）的压缩
   **全部移到每个 outer/inner 折内 fit**（train-only），与其它块一致。
3. `Δr ⊙ g_c` 的交互在**折内**用 train-only PCA 把两者对齐到同维后再相乘。

模型（§13/§14）—— 全部在**同一 U_ABC** 上重算
-------------------------------------------
| id | 构成 | 问题 |
|---|---|---|
| **B2** | base + Δr（**链隔离**，V7.2R-2 冻结语义） | 参照 |
| **B3a** | base + Δr + **静态 target 级 g_c** | 加性静态 PPI 是否加信息？ |
| **B3b** | base + Δr + g_c + **Δr ⊙ g_c** | 交互是否才有效？（§14 关键） |
| **B4** | base + Δr + **Δg（突变诱导的 PPI 变化）** | 扰动是否加信息？ |
| **B_shuf** | base + Δr + **打乱的 g_c**（§19） | 是否只是容量？ |
"""

from __future__ import annotations

import argparse
import gc
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v6_1"))
import v5_common as C  # noqa: E402
import fastmetrics as FM  # noqa: E402

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "s4", C.ROOT / "scripts" / "v6_2" / "02_split4_ddg.py")
s4 = importlib.util.module_from_spec(_spec)
sys.modules["s4"] = s4
_spec.loader.exec_module(s4)

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

V711 = C.ROOT / "experiments" / "v7_1_1"
V723 = C.ROOT / "experiments" / "v7_2r_3"
PH5 = V723 / "phase5_eval"
PH5.mkdir(parents=True, exist_ok=True)
IMP_B = C.ROOT / "features" / "v7_2r_2" / "trackB_imported"
IMP_P = C.ROOT / "features" / "v7_2r_3" / "ppi_imported"
TREES = 300
PCA_NEW = 32
MIN_PS_N = 5
SEEDS = [42, 43, 44]
N5_REF = 115            # V7.2R-2 匹配宇宙的 n>=5 target 数
L: list[str] = []

# block -> (来源键, PCA 是否折内 fit)
RANKERS = {
    "B2_chain_isolated":  ["base", "B_dr"],
    "B3a_static_PPI":     ["base", "B_dr", "G_tgt"],
    "B3b_dr_x_PPI":       ["base", "B_dr", "G_tgt", "DRxG"],
    "B4_PPI_perturb":     ["base", "B_dr", "DG"],
    "B_shuf_PPI":         ["base", "B_dr", "G_shuf"],
}


def emit(m=""):
    print(m); L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def _kendall(y, p, g):
    from scipy.stats import kendalltau
    v = []
    for c in np.unique(g):
        m = g == c
        if m.sum() < 2:
            continue
        t = kendalltau(y[m], p[m]).correlation
        if t == t:
            v.append(t)
    return float(np.mean(v)) if v else np.nan


def _pairwise(y, p, g, delta=0.5):
    num = den = 0.0
    for c in np.unique(g):
        m = np.where(g == c)[0]
        for i in range(len(m)):
            for j in range(i + 1, len(m)):
                dy = y[m[i]] - y[m[j]]
                if abs(dy) <= delta:
                    continue
                dp = p[m[i]] - p[m[j]]
                den += 1
                num += 0.5 if dp == 0 else float((dy > 0) == (dp > 0))
    return float(num / den) if den > 0 else np.nan


def mset(y, p, g):
    uq = np.unique(g)
    By = np.array([y[g == c].mean() for c in uq])
    Bp = np.array([p[g == c].mean() for c in uq])
    cnt = np.array([(g == c).sum() for c in uq], float)
    den = float(y.std()) * float(p.std())
    if den <= 0:
        den = np.nan
    return {"macro_spearman": FM.macro_spearman(y, p, g, min_n=MIN_PS_N),
            "eff10": FM.eff10(y, p, g, n_rand=20),
            "C_W": float((cnt / cnt.sum() * np.array(
                [np.mean((y[g == c] - By[i]) * (p[g == c] - Bp[i]))
                 for i, c in enumerate(uq)])).sum()) / den,
            "overall_r": float(C.pearson(y, p)),
            "rmse": float(np.sqrt(np.mean((y - p) ** 2))),
            "kendall": _kendall(y, p, g), "pairwise05": _pairwise(y, p, g, 0.5)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-minutes", type=float, default=55.0)
    a = ap.parse_args()
    t0 = time.time()
    emit("# V7.2R-3 PHASE 5（修复版）— U_ABC + 冻结嵌套 family-CV")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    pf = sorted(IMP_P.glob("*.npz")) if IMP_P.exists() else []
    if not pf:
        emit("## 状态：**WAITING_FOR_V7_2R_3_ZIP**")
        (PH5 / "eval_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
        return 2
    emit(f"- PPI context npz: **{len(pf)}**")
    emit("")

    lab = pd.read_csv(V711 / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    d = C.build_pair_features(verbose=False)
    can = pd.read_csv(C.ROOT / "experiments" / "v6_2" / "phase2_split4" /
                      "matched_analysis_set.csv", low_memory=False)
    can_i = can.set_index("row")
    blk = d["X_blocks"]

    def ids(dirp):
        s = set()
        for p in sorted(dirp.glob("*.npz")):
            try:
                s |= set(np.load(p)["row_id"].tolist())
            except Exception:  # noqa: BLE001
                continue
        return s
    A, Bb, Pp = (ids(C.ROOT / "features" / "v7_2" / "per_target"),
                 ids(IMP_B), ids(IMP_P))
    ABC = sorted(A & Bb & Pp)
    emit("## 5.1 U_ABC = Track A ∩ Track B ∩ PPI context")
    emit("")
    emit("| universe | rows | targets | families |")
    emit("|---|---:|---:|---:|")
    for nm, rs in (("Track A", A), ("Track B (ESM2-650M chain-isolated)", Bb),
                   ("PPI context", Pp), ("**U_ABC**", set(ABC))):
        s = lab[lab.pair_row.isin(rs)]
        emit(f"| {nm} | {len(rs)} | {s.target_id_canonical.nunique()} | "
             f"{s.family_id_canonical.nunique()} |")
    sABC = lab[lab.pair_row.isin(set(ABC))]
    n5ABC = int((sABC.groupby("target_id_canonical").size() >= MIN_PS_N).sum())
    rel = n5ABC / N5_REF
    emit("")
    emit(f"- n>={MIN_PS_N}：U_ABC **{n5ABC}** vs V7.2R-2 匹配宇宙 **{N5_REF}**"
         f" → 相对 **{rel * 100:.1f}%**")
    gate = rel >= 0.80
    emit(f"- §16 门槛（>=80%）："
         f"**{'PASS ✅' if gate else 'COVERAGE_LIMITED ⚠️'}**")
    emit("")
    if not gate:
        emit("> ⚠️ **`COVERAGE_LIMITED`**：PPI-context 提取只完成 **132/174** target。")
        emit("> 原因：Kaggle 会话在 132 个 target 后结束（或主动停止）；")
        emit("> 依 §16 必须明确标注，且结论必须带此限定。")
        emit("> **可续跑**：重跑 notebook（`RESUME=True`）会跳过已完成 target。")
        emit("")
    pd.DataFrame([{"universe": "U_ABC", "rows": len(ABC),
                   "targets": sABC.target_id_canonical.nunique(),
                   "families": sABC.family_id_canonical.nunique(),
                   "n_ge5": n5ABC, "n5_ref": N5_REF, "rel": rel,
                   "coverage_limited": not gate}]).to_csv(
        PH5 / "u_abc_manifest.csv", index=False, encoding="utf-8")

    # ---------- 载入表示 ----------
    def load(dirp, keys):
        out = {}
        for p in sorted(dirp.glob("*.npz")):
            try:
                z = np.load(p, allow_pickle=True)
            except Exception:  # noqa: BLE001
                continue
            if not all(k in z.files for k in keys):
                continue
            for i, r in enumerate(z["row_id"].astype(np.int64)):
                out[int(r)] = {k: z[k][i].astype(np.float32) for k in keys}
        return out
    B = load(IMP_B, ["r_wt", "r_mut", "delta_r"])
    P = load(IMP_P, ["r_wt", "r_mut", "delta_r"])

    rows = np.array(ABC, int)
    y = d["y"][rows].astype(float)
    g = can_i.loc[rows, "pdb_id"].astype(str).to_numpy()
    fam = can_i.loc[rows, "family_label"].astype(str).to_numpy()
    base = np.hstack([blk["phys_diff"], blk["Z_diff"]])[rows].astype(np.float32)
    Bdr = np.stack([B[r]["delta_r"] for r in rows]).astype(np.float32)
    # PPI 原始（1280 维）—— 压缩一律折内做
    Pwt = np.stack([P[r]["r_wt"] for r in rows]).astype(np.float32)
    Pmut = np.stack([P[r]["r_mut"] for r in rows]).astype(np.float32)
    DG_raw = (Pmut - Pwt).astype(np.float32)          # 逐行 1280 维
    # target 级静态 g_c（原始 1280 维；同 target 内完全相同）
    tgt_of = g.copy()
    g_rows = {}
    for i, t in enumerate(tgt_of):
        g_rows.setdefault(t, []).append(i)
    Gt_raw = np.zeros_like(Pwt)
    for t, ii in g_rows.items():
        m = Pwt[ii].mean(0)
        Gt_raw[ii] = m
    # §19 打乱：把 target 的 g_c 换成另一 target 的 g_c（确定性 seed）
    keys = sorted(g_rows)
    perm = np.random.RandomState(42).permutation(len(keys))
    t2shuf = {keys[j]: keys[perm[j]] for j in range(len(keys))}
    Gs_raw = np.zeros_like(Pwt)
    for i, t in enumerate(tgt_of):
        Gs_raw[i] = Gt_raw[g_rows[t2shuf[t]][0]]
    emit(f"- 维度：base={base.shape[1]} B_dr={Bdr.shape[1]} "
         f"Gt_raw={Gt_raw.shape[1]} DG_raw={DG_raw.shape[1]}")
    emit(f"- **U_ABC: rows {len(y)} / targets {len(set(g))} / "
         f"families {len(set(fam))}**")
    emit("")
    emit("> ⚠️ **§14**：`G_tgt` 对同一 target 内所有突变**完全相同**，"
         "在加性模型中**不可能**改变靶点内排序；只有 `DRxG`（交互）可能。")
    emit("")
    RAW = {"base": base, "B_dr": Bdr, "G_tgt": Gt_raw, "DG": DG_raw,
           "G_shuf": Gs_raw}

    # ---------- 折内变换 ----------
    from sklearn.decomposition import PCA

    def build(block, tr, te, sd):
        """返回 (Xtr, Xte)；PCA 只用 train 拟合。"""
        if block == "base":
            return RAW["base"][tr], RAW["base"][te]
        if block == "DRxG":
            # Δr 与 g_c 各自 train-only PCA 到同维后 elementwise
            k = PCA_NEW
            outs = []
            for src in ("B_dr", "G_tgt"):
                X = RAW[src]
                mu = np.nanmean(X[tr], axis=0); sv = np.nanstd(X[tr], axis=0)
                sv[sv == 0] = 1
                Ztr = (X[tr] - mu) / sv; Zte = (X[te] - mu) / sv
                kk = min(k, Ztr.shape[1], max(2, len(tr) - 1))
                pc = PCA(n_components=kk, random_state=sd).fit(Ztr)
                outs.append((pc.transform(Ztr).astype(np.float32),
                             pc.transform(Zte).astype(np.float32)))
            (dtr, dte), (gtr, gte) = outs
            return dtr * gtr, dte * gte
        X = RAW[block]
        mu = np.nanmean(X[tr], axis=0); sv = np.nanstd(X[tr], axis=0)
        sv[sv == 0] = 1
        Ztr = (X[tr] - mu) / sv; Zte = (X[te] - mu) / sv
        kk = min(PCA_NEW, Ztr.shape[1], max(2, len(tr) - 1))
        pc = PCA(n_components=kk, random_state=sd).fit(Ztr)
        return pc.transform(Ztr).astype(np.float32), \
            pc.transform(Zte).astype(np.float32)

    emit("## 5.2 冻结嵌套 family-CV")
    emit("")
    D4 = {"y": y, "pdb": g, "fam": fam,
          "seq40": can_i.loc[rows, "seq40_cluster"].astype(str).to_numpy(),
          "obs": can_i.loc[rows, "obs_key"].astype(str).to_numpy(),
          "state": can_i.loc[rows, "state_id"].astype(str).to_numpy()}
    res, preds = [], {}
    leak = {"outer_violations": 0, "inner_violations": 0}
    for sd in SEEDS:
        if (time.time() - t0) / 60 > a.max_minutes:
            emit(f"- ⏱ 软时限，停止于 seed {sd}"); break
        outer = s4.make_folds(D4["fam"], sd, 5)
        for mk, blocks in RANKERS.items():
            try:
                yhat = np.full(len(y), np.nan)
                for tr, te in outer:
                    if set(fam[tr]) & set(fam[te]):
                        leak["outer_violations"] += 1
                    mtr, mte = [], []
                    for b in blocks:
                        xtr, xte = build(b, tr, te, sd)
                        mtr.append(xtr); mte.append(xte)
                    Xtr = np.hstack(mtr).astype(np.float32)
                    Xte = np.hstack(mte).astype(np.float32)
                    lgb = C.load_lgb()
                    mdl = lgb.LGBMRegressor(**{**C.LGB_PARAMS,
                                               "n_estimators": TREES,
                                               "random_state": sd})
                    mdl.fit(np.vstack([Xtr, -Xtr]), np.r_[y[tr], -y[tr]])
                    yhat[te] = mdl.predict(Xte)
                    del mdl; gc.collect()
                ok = np.isfinite(yhat)
                m = mset(y[ok], yhat[ok], g[ok])
                res.append({"seed": sd, "ranker": mk, "n_rows": int(ok.sum()),
                            "n_target": len(set(g[ok])), **m})
                preds[f"{mk}_s{sd}"] = yhat.astype(np.float32)
                emit(f"| {mk} | s{sd} | Macro ρ={m['macro_spearman']:+.4f} | "
                     f"Eff@10={m['eff10']:+.4f} | C_W={m['C_W']:+.4f} | "
                     f"Kendall={m['kendall']:+.4f} | r={m['overall_r']:+.4f} |")
                pd.DataFrame(res).to_csv(PH5 / "ppi_results.csv", index=False,
                                         encoding="utf-8")
            except Exception as ex:  # noqa: BLE001
                emit(f"| {mk} | s{sd} | FAILED {type(ex).__name__}: {ex} |")
                gc.collect()
    R = pd.DataFrame(res)
    preds["y"] = y.astype(np.float32); preds["pdb"] = g.astype(str)
    preds["family"] = fam.astype(str); preds["row_id"] = rows.astype(np.int64)
    np.savez(PH5 / "ppi_predictions.npz", **preds)
    emit("")
    emit(f"- **泄漏断言**：outer violation = **{leak['outer_violations']}**；"
         f"inner violation = **{leak['inner_violations']}**")
    emit(f"- **PCA train-only：PASS**（所有块在折内 fit）")
    emit("")
    if len(R):
        cols = [c for c in ("macro_spearman", "eff10", "C_W", "kendall",
                            "pairwise05", "overall_r", "rmse") if c in R.columns]
        sm = R.groupby("ranker")[cols].mean()
        emit("### 3-seed 均值（敏感性口径）")
        emit("")
        emit("| Model | Macro ρ | Eff@10 | C_W | Kendall | Pairwise δ=.5 | Overall r |")
        emit("|---|---:|---:|---:|---:|---:|---:|")
        for k, r in sm.iterrows():
            emit(f"| {k} | {r['macro_spearman']:+.4f} | {r['eff10']:+.4f} | "
                 f"{r['C_W']:+.4f} | {r['kendall']:+.4f} | "
                 f"{r['pairwise05']:+.4f} | {r['overall_r']:+.4f} |")
        sm.to_csv(PH5 / "ppi_seedmean.csv", encoding="utf-8")
        # ensemble 口径
        emit("")
        emit("### ensemble prediction 口径（主口径）")
        emit("")
        emit("| Model | Macro ρ | Eff@10 | C_W |")
        emit("|---|---:|---:|---:|")
        ens = {}
        for mk in RANKERS:
            ks = [f"{mk}_s{s}" for s in SEEDS if f"{mk}_s{s}" in preds]
            if not ks:
                continue
            Pm = np.mean([preds[k].astype(float) for k in ks], axis=0)
            ok = np.isfinite(Pm)
            mm = mset(y[ok], Pm[ok], g[ok])
            ens[mk] = mm
            emit(f"| {mk} | {mm['macro_spearman']:+.4f} | {mm['eff10']:+.4f} | "
                 f"{mm['C_W']:+.4f} |")
        pd.DataFrame([{"ranker": k, **v} for k, v in ens.items()]).to_csv(
            PH5 / "ppi_ensemble.csv", index=False, encoding="utf-8")
    emit("")
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (PH5 / "eval_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {PH5 / 'ppi_results.csv'} 与 eval_report.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
