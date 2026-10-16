# -*- coding: utf-8 -*-
"""scripts/v8/38b_boundary_bootstrap_fast.py -- 快且**经等价性验证**的边界 bootstrap。

## 为什么需要它

`38_boundary_bootstrap_clean.py` 用规范实现 `v8metrics.compute_all`，**正确但慢**
（14.5 s/迭代 ⇒ 2000 次约 8 小时）。

## 本脚本的做法（**先证等价，再提速**）

1. **向量化** `overall_pearson`：用逐 target 预存的一/二阶量
   （`YS/PS/Y2/P2/YP`），只在 `n>=MIN_N` 的 target 上累加。
2. **向量化** `efficiency@10`：逐 target 预存 `EFF`，宏平均。
3. **强制等价性门（gate）**：在全量 4473 行上，把向量化实现与
   `v8metrics.compute_all` 的逐模型读数**逐位比对**（容差 1e-9）。
   **不通过则直接退出，不做 bootstrap。**
4. 通过后才跑 target-cluster bootstrap。

> 上一轮 `36` 的教训：内联公式**必须先过等价性门**才能用于推断。
"""

from __future__ import annotations

import argparse
import json
import sys
import time
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
SEED = 42
MIN_N = 5
K = 10
TOL = 1e-9


def load_all():
    H = pd.read_csv(EXT / "hermes_matched_final.csv", low_memory=False)
    z = np.load(R / "experiments" / "v7_1" / "phase2_rank" / "z_predictions.npz",
                allow_pickle=True)
    Y = np.asarray(z["y"], float)
    lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    H = H.sort_values("pair_row").reset_index(drop=True)
    TGT = H.target_id_canonical.astype(str).to_numpy()
    FAM = H.family_id_canonical.astype(str).to_numpy()
    assert len(H) == len(Y) == len(lab)

    POOL = pd.read_csv(EXT / "pool_composition_models.csv", low_memory=False)
    POOL = POOL[POOL.provenance.isin(["PERMISSIVE", "STRICT"])].dropna(
        subset=["overall_pearson", "efficiency@10"]).reset_index(drop=True)

    NONPRED = {"y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows",
               "z_raw", "y_cal", "family", "seed"}
    def is_label(k):
        kl = str(k).lower()
        return kl == "y" or kl.endswith("_y") or kl.endswith("_ddg")

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
    mas = pd.read_csv(R / "experiments" / "v6_2" / "phase2_split4" /
                      "matched_analysis_set.csv", low_memory=False)
    mas_row = mas.row.to_numpy(int)
    pos_of_mas = pos_all.reindex(mas_row)

    want = set(POOL.model_name)
    P, NM, GR = [], [], []
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
            if key not in want or key in NONPRED or is_label(key):
                continue
            av = np.asarray(zz[key])
            try:
                av = av.astype(float)
            except (ValueError, TypeError):
                continue
            if av.ndim != 1:
                continue
            if kind == "align":
                if len(av) != len(ROWS_ALL):
                    continue
                vec = av
            elif kind == "rowid":
                if rid is None or len(av) != len(rid):
                    continue
                m = pos_all.reindex(rid)
                if m.isna().any():
                    continue
                vec = np.full(len(ROWS_ALL), np.nan)
                vec[m.to_numpy(int)] = av
            else:
                if len(av) != len(mas_row):
                    continue
                vec = np.full(len(ROWS_ALL), np.nan)
                ok = pos_of_mas.notna().to_numpy()
                vec[pos_of_mas[ok].to_numpy(int)] = av[ok]
            P.append(vec); NM.append(key)
            GR.append(0 if POOL.loc[POOL.model_name == key,
                                    "provenance"].iloc[0] == "PERMISSIVE" else 1)
    return Y, TGT, FAM, np.vstack(P), np.asarray(NM), np.asarray(GR), POOL


def precompute(Y, TGT, P, tcode, n_t, cnt, starts, elig):
    """逐 eligible target 预存：EFF, YS, PS, Y2, P2, YP。"""
    M_ = P.shape[0]
    EFF = np.full((M_, n_t), np.nan)
    YS = np.zeros(n_t); PS = np.zeros((M_, n_t))
    Y2 = np.zeros(n_t); P2 = np.zeros((M_, n_t)); YP = np.zeros((M_, n_t))
    for g in elig:
        s = starts[g]; e = s + cnt[g]; n = cnt[g]
        yg = Y[s:e]; pg = P[:, s:e]
        ar = np.arange(n)
        true_order = np.lexsort((ar, yg))
        rank = np.empty(n, int); rank[true_order] = ar
        gain = 1.0 + (n - rank)
        ideal = np.lexsort((ar, -gain))
        kk = max(1, min(K, n))
        gid = gain[ideal[:kk]].sum()
        pred = np.argsort(pg, axis=1, kind="stable")[:, :kk]
        if gid > 0:
            EFF[:, g] = gain[pred].sum(axis=1) / gid
        YS[g] = yg.sum(); PS[:, g] = pg.sum(axis=1)
        Y2[g] = (yg ** 2).sum()
        P2[:, g] = (pg ** 2).sum(axis=1)
        YP[:, g] = (pg * yg).sum(axis=1)
    return EFF, YS, PS, Y2, P2, YP


def fast_metrics(c, elig, EFF, YS, PS, Y2, P2, YP, M_, tcode, Y):
    """c: 长度 n_t 的目标计数向量（已限制在 elig 上）。

    `overall_pearson` 用**组均值精确分解**：
        Cov    = Σ_t c_t (ȳ_t − ȳ)(p̄_t − p̄)      （跨项 ≡ 0，精确）
        Var_y  = Σ_t c_t (ȳ_t − ȳ)²
        Var_p  = Σ_t c_t (p̄_t − p̄)²
    其中 ȳ_t / p̄_t 为**目标均值**，ȳ / p̄ 为**行级整体均值**（按 c 加权）。
    该式按构造等于 `v8metrics._pearson`，由等价性门（容差 1e-9）验证。

    ⚠️ 前一版错误：对 y、p 行级中心化后点乘，未按 target 正确加权。
    等价性门正是为此设置。
    """
    n_tot = c.sum()
    if n_tot < 10:
        return np.full(M_, np.nan), np.full(M_, np.nan)
    use = c > 0
    cc = np.where(use, c, 0.0)
    y_row = np.bincount(tcode, weights=Y, minlength=len(c))
    ybar = (y_row * c).sum() / n_tot
    yt = np.where(use, y_row / np.maximum(c, 1e-12), 0.0)      # (n_t,)
    pt = np.where(use, PS / np.maximum(c, 1e-12), 0.0)         # (M_, n_t)
    pt_bar = (pt * cc).sum(axis=1) / n_tot                     # (M_,)
    dy = cc * (yt - ybar)                                      # (n_t,)
    dp = (pt - pt_bar[:, None]) * cc                           # (M_, n_t)
    cov = (dp * dy[None, :]).sum(axis=1)
    vy = (dy * (yt - ybar)).sum()
    vp = (dp * (pt - pt_bar[:, None])).sum(axis=1)
    den = np.sqrt(np.maximum(vy, 1e-12) * np.maximum(vp, 1e-12))
    ov = cov / den
    ne = c[elig].sum()
    ef = (EFF[:, elig] * c[elig]).sum(axis=1) / max(ne, 1.0)
    return ov, ef


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--max-minutes", type=float, default=30.0)
    a = ap.parse_args()

    Y, TGT, FAM, P, NM, GR, POOL = load_all()
    M_ = P.shape[0]
    print(f"模型池 = {M_} (PERMISSIVE {int((GR==0).sum())} / STRICT {int((GR==1).sum())})"
          f" | 行 = {len(Y)}")
    uq, tcode = np.unique(TGT, return_inverse=True)
    n_t = len(uq)
    cnt = np.bincount(tcode, minlength=n_t)
    starts = np.concatenate([[0], np.cumsum(cnt)[:-1]])
    elig = np.where(cnt >= MIN_N)[0]
    print(f"n>={MIN_N} 的 target = {len(elig)}")

    EFF, YS, PS, Y2, P2, YP = precompute(Y, TGT, P, tcode, n_t, cnt, starts, elig)
    print("预计算完成")

    # ---------- 等价性门 ----------
    cfull = np.zeros(n_t); cfull[elig] = cnt[elig]
    ov_f, ef_f = fast_metrics(cfull, elig, EFF, YS, PS, Y2, P2, YP, M_, tcode, Y)
    print("\n" + "=" * 78)
    print("等价性门：向量化 vs v8metrics.compute_all（全量 4473 行）")
    print("=" * 78)
    d_ov, d_ef = [], []
    for i in range(M_):
        m = np.isfinite(Y) & np.isfinite(P[i])
        if m.sum() < 50:
            continue
        U = M.prep_universe(Y[m], P[i][m], TGT[m], FAM[m])
        r = M.compute_all(U)
        d_ov.append(abs(r["overall_pearson"] - ov_f[i]))
        d_ef.append(abs(r["efficiency@10"] - ef_f[i]))
    d_ov = np.asarray(d_ov); d_ef = np.asarray(d_ef)
    print(f"  比对模型数 = {len(d_ov)}")
    print(f"  max |Δ overall_pearson| = {d_ov.max():.3e}")
    print(f"  max |Δ efficiency@10|   = {d_ef.max():.3e}")
    if d_ov.max() > TOL or d_ef.max() > TOL:
        print(f"\n  ❌ 等价性门**未通过**（容差 {TOL:.0e}）⇒ 不做 bootstrap。")
        print("  这保护我们不再犯 36 号脚本的错误。")
        return 1
    print(f"  ✅ 等价性门通过（容差 {TOL:.0e}）⇒ 允许 bootstrap")

    def rho_of(ov, ef, mask):
        v = np.isfinite(ov) & np.isfinite(ef) & mask
        if v.sum() < 6:
            return np.nan
        return float(spearmanr(ov[v], ef[v]).correlation)

    mperm = GR == 0; mstr = GR == 1
    p0 = rho_of(ov_f, ef_f, mperm); s0 = rho_of(ov_f, ef_f, mstr)
    print(f"\n点估计: PERMISSIVE rho={p0:+.4f} | STRICT rho={s0:+.4f} | Δ={p0-s0:+.4f}")
    print(f"  对照 pool: PERMISSIVE +0.8088 / STRICT -0.0289  "
          f"⇒ {'✅ 一致' if abs(p0-0.8088)<0.02 and abs(s0+0.0289)<0.02 else '⚠ 不一致'}")

    tidx = [np.where(tcode == u)[0] for u in range(n_t)]
    rs = np.random.RandomState(SEED)
    ps, ss, ds = [], [], []
    t0 = time.time(); valid = 0
    for it in range(a.n_boot):
        pick = rs.randint(0, len(elig), len(elig))
        c = np.bincount(elig[pick], minlength=n_t).astype(float)
        ov, ef = fast_metrics(c, elig, EFF, YS, PS, Y2, P2, YP, M_, tcode, Y)
        rp = rho_of(ov, ef, mperm); rq = rho_of(ov, ef, mstr)
        if not (np.isfinite(rp) and np.isfinite(rq)):
            continue
        valid += 1; ps.append(rp); ss.append(rq); ds.append(rp - rq)
        if (it + 1) % 250 == 0:
            el = (time.time() - t0) / 60
            print(f"  {it+1}/{a.n_boot} valid={valid} {el:.2f}min "
                  f"ETA {el/(it+1)*(a.n_boot-it-1):.2f}min", flush=True)
            if el > a.max_minutes:
                print("  软时限，停止"); break
    ps = np.asarray(ps); ss = np.asarray(ss); ds = np.asarray(ds)
    print(f"\n有效迭代 = {valid}/{a.n_boot} ({(time.time()-t0)/60:.2f} min)")
    out = []
    for nm, v, pt in (("PERMISSIVE", ps, p0), ("STRICT", ss, s0),
                      ("DELTA", ds, p0 - s0)):
        if len(v) < 20:
            continue
        lo, hi = np.percentile(v, [2.5, 97.5])
        ex = bool(lo > 0 or hi < 0)
        print(f"  {nm:<11} point={pt:+.4f}  target-cluster 95%CI=[{lo:+.4f},{hi:+.4f}]  "
              f"{'排除0' if ex else '含0'}")
        out.append({"quantity": nm, "point": float(pt), "ci_low": float(lo),
                    "ci_high": float(hi), "excludes_zero": ex})
    pd.DataFrame(out).to_csv(EXT / "boundary_bootstrap_final.csv", index=False,
                             encoding="utf-8")
    pl, ph = np.percentile(ps, [2.5, 97.5]); sl, sh = np.percentile(ss, [2.5, 97.5])
    overlap = not (ph < sl or sh < pl)
    meta = {"implementation": "vectorized, PASSED equivalence gate vs v8metrics",
            "equivalence_max_abs_diff_overall_pearson": float(d_ov.max()),
            "equivalence_max_abs_diff_eff10": float(d_ef.max()),
            "n_boot_requested": a.n_boot, "n_valid": valid, "seed": SEED,
            "cluster_unit": "target", "K": K, "min_n": MIN_N,
            "point_PERMISSIVE": float(p0), "point_STRICT": float(s0),
            "ci_overlap": bool(overlap),
            "verdict": ("BOUNDARY_ROBUST" if not overlap else "BOUNDARY_NOT_ROBUST")}
    (EXT / "boundary_bootstrap_final_meta.json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n两组 target-cluster CI {'重叠' if overlap else '**不重叠** ✅'}")
    print("写出 boundary_bootstrap_final.csv / _meta.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
