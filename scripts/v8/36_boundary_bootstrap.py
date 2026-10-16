# -*- coding: utf-8 -*-
"""scripts/v8/36_boundary_bootstrap.py -- §6.2.2 的**区间**检验（PHASE 15 扩展）。

## 被检验的量

同一宇宙（canonical 4473 行），按训练划分 provenance 分两组：
    PERMISSIVE (n=37) : 随机 / complex / seq40 划分训练
    STRICT     (n=133): family-OOD 划分训练

点估计：ρ(Overall r, Eff@10)  = **+0.8088** (PERMISSIVE) vs **−0.0289** (STRICT)

**问题**：这个 0.84 的差距是稳健的吗？

## 方法（target-cluster bootstrap，§55/§56）

重采样单位 = **target**，N=5000，seed=42。
每次迭代：
  1. 有放回抽 target，拼出新的行集合
  2. 在该行集合上**重算**每个模型的 (Overall r, Eff@10)
  3. 在两组内分别求 ρ(Overall r, Eff@10)
  4. 记录 Δ = ρ_PERMISSIVE − ρ_STRICT

**关键**：两组模型**互不重叠**，故二者可视为独立，Δ 的区间可由
两个边际分布直接相减得到（也同时记录直接用 Δ 分布的结果做交叉校验）。

CPU-first：预计算每模型逐 target 的中间量不可行（Eff@10 需 target 内排序），
故本脚本用**向量化的逐 target 贡献**实现：对每个模型预存每个 target 的
`gain[np.argsort(pg)[:k]].sum()`，重采样时只做查表求和。
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
OUT = R / "experiments" / "v8"
EXT = OUT / "external"
SEED = 42
MIN_N = 5
K = 10


def build():
    H = pd.read_csv(EXT / "hermes_matched_final.csv", low_memory=False)
    z = np.load(R / "experiments" / "v7_1" / "phase2_rank" / "z_predictions.npz",
                allow_pickle=True)
    Y = np.asarray(z["y"], float)
    lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    H = H.sort_values("pair_row").reset_index(drop=True)
    TGT = H.target_id_canonical.astype(str).to_numpy()
    assert len(H) == len(Y) == len(lab)

    NONPRED = {"y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows",
               "z_raw", "y_cal", "family", "seed"}
    def is_label(k):
        kl = str(k).lower()
        return kl == "y" or kl.endswith("_y") or kl.endswith("_ddg")

    PERM_TOK = ("L1_random_", "L2_complex_", "L3_seq40_")
    STRICT_TOK = ("L4_family_", "F0_", "F1_", "F2_", "F3_", "F4_", "F5_",
                  "F6_", "F7_", "F8_", "F9_", "R0_", "R1_", "R2_", "R3_",
                  "B0_", "B1_", "B2_", "B3", "B4_", "A2_smallPLM", "p_",
                  "z_R", "z_A")
    def grp(key):
        if any(key.startswith(t) for t in PERM_TOK):
            return 0
        if any(key.startswith(t) for t in STRICT_TOK):
            return 1
        return -1

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
            else:
                if len(a) != len(mas_row):
                    continue
                av = np.full(len(ROWS_ALL), np.nan)
                ok = pos_of_mas.notna().to_numpy()
                av[pos_of_mas[ok].to_numpy(int)] = a[ok]
            if (np.isfinite(Y) & np.isfinite(av)).sum() < 100:
                continue
            g = grp(key)
            if g < 0:
                continue
            P.append(av); NM.append(key); GR.append(g)
    return Y, TGT, np.vstack(P), np.asarray(NM), np.asarray(GR)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-boot", type=int, default=5000)
    ap.add_argument("--max-minutes", type=float, default=45.0)
    a = ap.parse_args()

    Y, TGT, P, NM, GR = build()
    M_, N_ = P.shape
    print(f"模型 = {M_} (PERMISSIVE {int((GR==0).sum())} / STRICT {int((GR==1).sum())})"
          f" | 行 = {N_} | targets = {len(set(TGT))}")

    uq, tcode = np.unique(TGT, return_inverse=True)
    n_t = len(uq)
    order = np.argsort(tcode, kind="stable")
    cnt = np.bincount(tcode, minlength=n_t)
    starts = np.concatenate([[0], np.cumsum(cnt)[:-1]])
    elig = np.where(cnt >= MIN_N)[0]
    print(f"n>={MIN_N} 的 target = {len(elig)}")

    # 预计算：每个模型在每个 target 的 (gain 和, Eff@10, overall 贡献)
    # gain 与 target 无关的排序无关部分只依赖 y
    EFF = np.full((M_, n_t), np.nan)     # Eff@10 逐 target
    GSUM = np.zeros((M_, n_t))           # Σ(y - ȳ_target) 需要整体均值 -> 用原始 y 和
    YS = np.zeros((M_, n_t)); YC = np.zeros(n_t); PS = np.zeros((M_, n_t))
    for g in elig:
        s = starts[g]; e = s + cnt[g]; n = cnt[g]
        yg = Y[s:e]
        ar = np.arange(n)
        true_order = np.lexsort((ar, yg))
        rank = np.empty(n, int); rank[true_order] = ar
        gain = 1.0 + (n - rank)
        ideal = np.lexsort((ar, -gain))
        kk = max(1, min(K, n))
        gid = gain[ideal[:kk]].sum()
        pg = P[:, s:e]
        pred = np.argsort(pg, axis=1, kind="stable")[:, :kk]
        if gid > 0:
            EFF[:, g] = gain[pred].sum(axis=1) / gid
        YS[:, g] = yg.sum()
        PS[:, g] = pg.sum(axis=1)
    YC[:] = cnt
    print("预计算完成")

    def metrics(sel_t):
        """sel_t: 被抽中的 target 索引（可重复）。返回 (overall_r, eff) 逐模型。

        **只在 n>=MIN_N 的 target 上计算**（与 `v8metrics.compute_all` 口径一致），
        否则会混入小 target 而改变 overall Pearson。
        """
        c = np.bincount(sel_t, minlength=n_t).astype(float)
        c = c * np.isin(np.arange(n_t), elig)     # ★ 只保留 eligible target
        n_tot = c.sum()
        if n_tot < 10:
            return (np.full(M_, np.nan), np.full(M_, np.nan))
        sy = (YS * c).sum(axis=1)                 # Σ y      (M,)
        sp = (PS * c).sum(axis=1)                 # Σ p      (M,)
        syy = (Y2 * c).sum()                      # Σ y^2    scalar
        spp = (P2 * c).sum(axis=1)                # Σ p^2    (M,)
        syp = (YP * c).sum(axis=1)                # Σ y·p    (M,)
        ybar = (np.bincount(tcode, weights=Y, minlength=n_t) * c).sum() / n_tot
        pbar = sp / n_tot
        vy = syy / n_tot - ybar ** 2
        vp = spp / n_tot - pbar ** 2
        cv = syp / n_tot - ybar * pbar
        den = np.sqrt(np.maximum(vy, 1e-12) * np.maximum(vp, 1e-12))
        ov = cv / den
        # Eff@10：只在**被抽中的 eligible target** 上取宏平均
        ne = c[elig].sum()
        ef = (EFF[:, elig] * c[elig]).sum(axis=1) / max(ne, 1.0)
        return ov, ef

    print("二阶量预计算完成")
    Y2 = np.zeros(n_t); P2 = np.zeros((M_, n_t)); YP = np.zeros((M_, n_t))
    for g in elig:
        s = starts[g]; e = s + cnt[g]
        yg = Y[s:e]; pg = P[:, s:e]
        Y2[g] = (yg ** 2).sum()
        P2[:, g] = (pg ** 2).sum(axis=1)
        YP[:, g] = (pg * yg).sum(axis=1)
    print("二阶量预计算完成")

    def metrics2_old(sel_t):
        c = np.bincount(sel_t, minlength=n_t).astype(float)
        n_tot = c.sum()
        sy = (YS * c).sum(axis=1)          # (M,)
        sp = (PS * c).sum(axis=1)
        syy = (Y2 * c).sum()               # scalar
        spp = (P2 * c).sum(axis=1)
        syp = (YP * c).sum(axis=1)
        ybar = syy * 0 + (np.bincount(tcode, weights=Y,
                                     minlength=n_t) * c).sum() / n_tot
        pbar = sp / n_tot
        vy = syy / n_tot - ybar ** 2
        vp = spp / n_tot - pbar ** 2
        cv = syp / n_tot - ybar * pbar
        den = np.sqrt(np.maximum(vy, 1e-12) * np.maximum(vp, 1e-12))
        ov = cv / den
        ef = (EFF[:, :] * c).sum(axis=1) / c.sum()
        return ov, ef

    ov0, ef0 = metrics(np.arange(n_t))
    def rho(ov, ef, mask):
        v = np.isfinite(ov) & np.isfinite(ef) & mask
        if v.sum() < 6:
            return np.nan
        return float(spearmanr(ov[v], ef[v]).correlation)
    mperm = GR == 0; mstr = GR == 1
    p0 = rho(ov0, ef0, mperm); s0 = rho(ov0, ef0, mstr)
    print(f"\n点估计: PERMISSIVE rho={p0:+.4f} (n={int(mperm.sum())}) | "
          f"STRICT rho={s0:+.4f} (n={int(mstr.sum())}) | Δ={p0-s0:+.4f}")

    rs = np.random.RandomState(SEED)
    ps, ss, ds = [], [], []
    t0 = time.time(); valid = 0
    for it in range(a.n_boot):
        pick = rs.randint(0, len(elig), len(elig))
        sel_t = elig[pick]
        ov, ef = metrics(sel_t)
        rp = rho(ov, ef, mperm); rsx = rho(ov, ef, mstr)
        if not (np.isfinite(rp) and np.isfinite(rsx)):
            continue
        valid += 1
        ps.append(rp); ss.append(rsx); ds.append(rp - rsx)
        if (it + 1) % 500 == 0:
            el = (time.time() - t0) / 60
            print(f"  {it+1}/{a.n_boot} valid={valid} {el:.1f}min "
                  f"ETA {el/(it+1)*(a.n_boot-it-1):.1f}min", flush=True)
            if el > a.max_minutes:
                print("  软时限，停止"); break

    ps = np.asarray(ps); ss = np.asarray(ss); ds = np.asarray(ds)
    print(f"\n有效迭代 = {valid}/{a.n_boot} ({(time.time()-t0)/60:.1f} min)")
    out = []
    for nm, v, pt in (("PERMISSIVE", ps, p0), ("STRICT", ss, s0),
                      ("DELTA", ds, p0 - s0)):
        lo, hi = np.percentile(v, [2.5, 97.5])
        ex = bool(lo > 0 or hi < 0)
        print(f"  {nm:<11} point={pt:+.4f}  95%CI=[{lo:+.4f},{hi:+.4f}]  "
              f"{'排除0' if ex else '含0'}")
        out.append({"quantity": nm, "point": pt, "ci_low": float(lo),
                    "ci_high": float(hi), "excludes_zero": ex,
                    "n_valid": valid})
    # 两组区间是否重叠
    pl, ph = np.percentile(ps, [2.5, 97.5])
    sl, sh = np.percentile(ss, [2.5, 97.5])
    overlap = not (ph < sl or sh < pl)
    print(f"\n  两组 95% CI {'**重叠**' if overlap else '**不重叠** ✅'}")
    D = pd.DataFrame(out)
    D.to_csv(EXT / "boundary_bootstrap.csv", index=False, encoding="utf-8")
    meta = {"n_boot_requested": a.n_boot, "n_valid": valid, "seed": SEED,
            "cluster_unit": "target", "K": K, "min_n": MIN_N,
            "point_PERMISSIVE": p0, "point_STRICT": s0,
            "ci_overlap": bool(overlap),
            "verdict": ("BOUNDARY_ROBUST" if not overlap else "BOUNDARY_NOT_ROBUST"),
            "caveat": ("两组模型不重叠但非随机分配（provenance 由键名推断）；"
                       "bootstrap 只反映 target 抽样不确定性，不反映模型池选择不确定性。")}
    (EXT / "boundary_bootstrap_meta.json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n写出 boundary_bootstrap.csv / boundary_bootstrap_meta.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
