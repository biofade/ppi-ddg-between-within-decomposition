# -*- coding: utf-8 -*-
"""scripts/v8/09_weighting_audit.py -- V8 PHASE 9 / 核心实验 C：目标加权敏感性（§32/§33/§34）。

预测不变，**只改聚合权重**：
    ROW    : 每个突变等权（= 行权重，n 大的 target 主导）
    EQUAL  : 每个 target 等权
    SQRT   : 权重 ∝ sqrt(n)
    LOG    : 权重 ∝ log(1+n)
    CAP30  : 权重 ∝ min(n,30)

对每个模型报告 C_B / C_W / 加权 Pearson 在各方案下的取值、极差与符号变化，
以及模型**排名**是否随加权方案改变。

另做 §34 目标规模集中度审计：剔除最大 5/10/20 个 target 后 Overall r / C_B / C_W 的变化。
"""

from __future__ import annotations

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
OUT = R / "experiments" / "v8"
W = OUT / "weighting"
W.mkdir(parents=True, exist_ok=True)

SCHEMES = ["ROW", "EQUAL", "SQRT", "LOG", "CAP30"]


def scheme_weights(cnt, scheme):
    c = np.asarray(cnt, float)
    if scheme == "ROW":
        return c
    if scheme == "EQUAL":
        return np.ones_like(c)
    if scheme == "SQRT":
        return np.sqrt(c)
    if scheme == "LOG":
        return np.log1p(c)
    if scheme == "CAP30":
        return np.minimum(c, 30.0)
    raise ValueError(scheme)


def weighted_metrics(y, p, tgt, scheme):
    """返回 (w_pearson, C_B, C_W, C_B_equal, C_W_equal)。

    加权 Pearson：对行加权后求相关（等价于对中心和尺度也重新定义）。
    C_B / C_W：按 §28 分解，但在**加权**测度下计算。
    """
    y = np.asarray(y, float); p = np.asarray(p, float)
    t = np.asarray(tgt).astype(str)
    uq, inv = np.unique(t, return_inverse=True)
    n_t = np.bincount(inv)
    w_t = scheme_weights(n_t, scheme)
    w = w_t[inv]
    W = w.sum()
    my = np.sum(w * y) / W; mp = np.sum(w * p) / W
    vy = np.sum(w * (y - my) ** 2) / W; vp = np.sum(w * (p - mp) ** 2) / W
    if vy <= 0 or vp <= 0:
        return (np.nan,) * 5
    cov = np.sum(w * (y - my) * (p - mp)) / W
    r = cov / np.sqrt(vy * vp)
    # 目标均值（加权）
    By = np.array([np.average(y[inv == k], weights=w[inv == k]) for k in range(len(uq))])
    Bp = np.array([np.average(p[inv == k], weights=w[inv == k]) for k in range(len(uq))])
    sd = np.sqrt(vy) * np.sqrt(vp)
    C_B = np.sum(w_t * (By - my) * (Bp - mp)) / (W * sd)
    Wy = y - By[inv]; Wp = p - Bp[inv]
    C_W = np.sum(w * Wy * Wp) / (W * sd)
    # 等权变体
    Bye = np.array([y[inv == k].mean() for k in range(len(uq))])
    Bpe = np.array([p[inv == k].mean() for k in range(len(uq))])
    C_Be = np.mean((Bye - Bye.mean()) * (Bpe - Bpe.mean())) / sd
    C_We = np.mean(Wy * Wp) / sd
    return float(r), float(C_B), float(C_W), float(C_Be), float(C_We)


def drop_largest(y, p, tgt, k):
    t = np.asarray(tgt).astype(str)
    uq, cnt = np.unique(t, return_counts=True)
    order = uq[np.argsort(-cnt)]
    drop = set(order[:k].tolist())
    m = ~np.isin(t, list(drop))
    return y[m], p[m], t[m]


def main() -> int:
    S = pd.read_csv(OUT / "BENCHMARK_SCOREBOARD.csv", low_memory=False)
    print("scoreboard rows =", len(S))

    # 重建每个模型的 (y, p, tgt) —— 复用 06 的装载逻辑
    sys.path.insert(0, str(R / "scripts" / "v8"))
    lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    mas = pd.read_csv(R / "experiments" / "v6_2" / "phase2_split4" /
                      "matched_analysis_set.csv", low_memory=False)
    CAN = lab[["pair_row", "target_id_canonical", "family_id_canonical"]].rename(
        columns={"pair_row": "row_id", "target_id_canonical": "target_id",
                 "family_id_canonical": "family_id"})
    ROWS = CAN.row_id.to_numpy(int)
    TGT = CAN.target_id.to_numpy(str)
    z = np.load(R / "experiments" / "v7_1" / "phase2_rank" / "z_predictions.npz",
               allow_pickle=True)
    Y = np.asarray(z["y"], float)
    MAS_ROWS = mas.row.to_numpy(int); MAS_TGT = mas.pdb_id.astype(str).to_numpy()
    Y_MAS = pd.to_numeric(mas["ddG"], errors="coerce").to_numpy(float)

    NONPRED = {"y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows",
               "z_raw", "y_cal", "family", "seed", "masks"}
    def is_label_key(k):
        kl = str(k).lower()
        return kl == "y" or kl.endswith("_y") or kl.endswith("_ddg")

    _cache = {}
    rows, rank_rows = [], []
    for r in S.itertuples():
        if r.universe == "U_CANONICAL":
            y, tgt, base = Y, TGT, ROWS
        else:
            y, tgt, base = Y_MAS, MAS_TGT, MAS_ROWS
        rel = r.prediction_file
        if rel not in _cache:
            _cache[rel] = np.load(R / rel, allow_pickle=True)
        zz = _cache[rel]
        if r.model_name not in zz.files or is_label_key(r.model_name):
            continue
        a = np.asarray(zz[r.model_name])
        if a.ndim != 1:
            continue
        try:
            a = a.astype(float)
        except (ValueError, TypeError):
            continue
        if len(a) == len(Y) and np.isfinite(a).all():
            if abs(np.corrcoef(a, Y)[0, 1]) > 0.98:
                continue
        if r.universe == "U_CANONICAL" and len(a) != len(base):
            continue
        m = np.isfinite(y) & np.isfinite(a)
        yy, pp, tt = y[m], a[m], tgt[m]
        rec = {"model_id": r.model_id, "model_name": r.model_name,
               "model_family": r.model_family, "universe": r.universe}
        for s in SCHEMES:
            wr, cb, cw, cbe, cwe = weighted_metrics(yy, pp, tt, s)
            rec[f"pearson_{s}"] = wr
            rec[f"C_B_{s}"] = cb
            rec[f"C_W_{s}"] = cw
        rows.append(rec)
        # 排名一致性：各方案下按加权 Pearson 排名
        rank_rows.append(rec)
    D = pd.DataFrame(rows)
    D.to_csv(W / "weighting_sensitivity.csv", index=False, encoding="utf-8")

    print(f"\n评估 {len(D)} 个模型 × {len(SCHEMES)} 加权方案")
    print("\n=== C_B 在各方案下的极差（按极差排序前 12）===")
    cbcols = [f"C_B_{s}" for s in SCHEMES]
    cwcols = [f"C_W_{s}" for s in SCHEMES]
    D["C_B_range"] = D[cbcols].max(axis=1) - D[cbcols].min(axis=1)
    D["C_W_range"] = D[cwcols].max(axis=1) - D[cwcols].min(axis=1)
    D["C_B_sign_change"] = (D[cbcols] > 0).sum(axis=1).between(1, len(SCHEMES) - 1)
    D["C_W_sign_change"] = (D[cwcols] > 0).sum(axis=1).between(1, len(SCHEMES) - 1)
    show = ["model_name", "model_family", "universe"] + cbcols + ["C_B_range", "C_B_sign_change"]
    print(D.nlargest(12, "C_B_range")[show].to_string(index=False))
    print()
    print("=== 汇总 ===")
    print(f"  C_B 极差中位数 = {D.C_B_range.median():.4f} | 最大 = {D.C_B_range.max():.4f}")
    print(f"  C_W 极差中位数 = {D.C_W_range.median():.4f} | 最大 = {D.C_W_range.max():.4f}")
    print(f"  C_B 出现符号变化的模型数 = {int(D.C_B_sign_change.sum())} / {len(D)}")
    print(f"  C_W 出现符号变化的模型数 = {int(D.C_W_sign_change.sum())} / {len(D)}")

    # 模型**排名**是否随方案改变
    print("\n=== 模型排名随加权方案的变化（每个 universe）===")
    rank_out = []
    for uni, g in D.groupby("universe"):
        g = g.dropna(subset=[f"pearson_{s}" for s in SCHEMES])
        if len(g) < 8:
            continue
        rk = {}
        for s in SCHEMES:
            v = g[f"pearson_{s}"].to_numpy(float)
            o = np.argsort(-v, kind="stable")
            rr = np.empty(len(v), int); rr[o] = np.arange(1, len(v) + 1)
            rk[s] = rr
        RK = pd.DataFrame(rk)
        sp = spearmanr(RK["ROW"], RK["EQUAL"]).correlation
        ke = kendalltau(RK["ROW"], RK["EQUAL"]).correlation
        mx = int(RK["ROW"].sub(RK["EQUAL"]).abs().max())
        rank_out.append({"universe": uni, "n_models": len(g),
                         "spearman_ROW_vs_EQUAL": sp, "kendall_ROW_vs_EQUAL": ke,
                         "max_rank_shift": mx})
        print(f"  {uni}: n={len(g)} Spearman(ROW,EQUAL)={sp:+.3f} "
              f"Kendall={ke:+.3f} 最大秩位移={mx}")
    pd.DataFrame(rank_out).to_csv(W / "weighting_rank_shift.csv", index=False,
                                  encoding="utf-8")

    # §34 目标规模集中度：剔除最大 target
    print("\n=== §34 目标规模集中度审计（U_CANONICAL，取 macro_spearman 前 3 的模型）===")
    uc = S[(S.universe == "U_CANONICAL") & (S.n_rows == len(Y))].nlargest(
        3, "macro_spearman")
    con = []
    for r in uc.itertuples():
        zz = _cache.get(r.prediction_file)
        if zz is None:
            zz = np.load(R / r.prediction_file, allow_pickle=True)
            _cache[r.prediction_file] = zz
        if r.model_name not in zz.files:
            continue
        a = np.asarray(zz[r.model_name]).astype(float)
        if len(a) != len(Y):
            continue
        m = np.isfinite(Y) & np.isfinite(a)
        yy, pp, tt = Y[m], a[m], TGT[m]
        for k in (0, 5, 10, 20):
            y2, p2, t2 = drop_largest(yy, pp, tt, k) if k else (yy, pp, tt)
            wr, cb, cw, _, _ = weighted_metrics(y2, p2, t2, "ROW")
            con.append({"model": r.model_name, "drop_top_n": k,
                        "n_rows": len(y2), "n_targets": len(set(t2)),
                        "w_pearson": wr, "C_B": cb, "C_W": cw})
    C = pd.DataFrame(con)
    C.to_csv(W / "target_size_concentration.csv", index=False, encoding="utf-8")
    if len(C):
        print(C.to_string(index=False))

    # 报告
    lines = ["# V8 Target-Weighting Sensitivity Report\n",
             "> 预测**完全不变**；只改聚合权重。\n",
             f"> 方案：{' / '.join(SCHEMES)}\n",
             "## 1. C_B / C_W 对加权的敏感性\n",
             f"- C_B 极差中位 {D.C_B_range.median():.4f}，最大 {D.C_B_range.max():.4f}",
             f"- C_W 极差中位 {D.C_W_range.median():.4f}，最大 {D.C_W_range.max():.4f}",
             f"- C_B 出现符号变化的模型：{int(D.C_B_sign_change.sum())} / {len(D)}",
             f"- C_W 出现符号变化的模型：{int(D.C_W_sign_change.sum())} / {len(D)}",
             "\n## 2. 模型排名是否随加权改变\n"]
    for x in rank_out:
        lines.append(f"- `{x['universe']}`: n={x['n_models']}, "
                     f"Spearman(ROW,EQUAL)={x['spearman_ROW_vs_EQUAL']:+.3f}, "
                     f"最大秩位移={x['max_rank_shift']}")
    lines.append("\n## 3. 目标规模集中度（剔除最大 target）\n")
    if len(C):
        lines.append("```")
        lines.append(C.to_string(index=False))
        lines.append("```")
    (OUT / "WEIGHTING_SENSITIVITY_REPORT.md").write_text("\n".join(lines) + "\n",
                                                         encoding="utf-8")
    print("\n写出 weighting_sensitivity.csv / weighting_rank_shift.csv / "
          "target_size_concentration.csv / WEIGHTING_SENSITIVITY_REPORT.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
