# -*- coding: utf-8 -*-
"""scripts/v8_1/102_weighting_semantics.py -- 判定五种权重方案的真实数学语义。

两个候选解释（用户提出的核心问题）：

  (P) PER-ROW MULTIPLIER —— 代码当前实现
        w_ci = a_c ,  a_c in {n_c, 1, sqrt(n_c), log1p(n_c), min(n_c,30)}
        => target 的总质量 = n_c * a_c
        => 名为 EQUAL 的方案实际是「每行等权」，大 target 仍主导
        => 名为 ROW 的方案每行权重 n_c，是纯行权重

  (T) TARGET TOTAL MASS —— 与 scheme 名称一致的解释
        A_c = a_c（target 总权重），行权重 w_ci = A_c / n_c
        => target 的总质量 = A_c
        => EQUAL(a=1) 才真正是「每 target 等权」
        => ROW(a=n_c) 等价于每行等权（每 target 总质量 n_c）

对两种解释分别计算 109 个配置 × 5 schemes 的 r / C_B / C_W，并比较：
  * 恒等式是否精确
  * C_B / C_W 的跨 scheme 极差与符号翻转计数
  * 模型排名在不同 scheme 间是否一致
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(".")
SCHEMES = ["ROW", "EQUAL", "SQRT", "LOG", "CAP30"]


def a_of(n, scheme):
    c = np.asarray(n, float)
    return {"ROW": c, "EQUAL": np.ones_like(c), "SQRT": np.sqrt(c),
            "LOG": np.log1p(c), "CAP30": np.minimum(c, 30.0)}[scheme]


def weighted(y, p, tgt, scheme, mode):
    """mode='P' -> w=a_c ; mode='T' -> w=a_c/n_c。返回 (r, C_B, C_W)。"""
    y = np.asarray(y, float); p = np.asarray(p, float)
    t = np.asarray(tgt).astype(str)
    uq, inv = np.unique(t, return_inverse=True)
    n_c = np.bincount(inv).astype(float)
    a = a_of(n_c, scheme)
    w_row = a if mode == "P" else a / n_c        # w_ci
    w = w_row[inv]
    W = w.sum()
    my = np.sum(w * y) / W; mp = np.sum(w * p) / W
    vy = np.sum(w * (y - my) ** 2) / W; vp = np.sum(w * (p - mp) ** 2) / W
    if vy <= 0 or vp <= 0:
        return (np.nan,) * 3
    sd = np.sqrt(vy * vp)
    By = np.empty(len(uq)); Bp = np.empty(len(uq)); wbar = np.empty(len(uq))
    for k in range(len(uq)):
        m = inv == k
        wbar[k] = w[m].sum()
        By[k] = np.sum(w[m] * y[m]) / wbar[k]
        Bp[k] = np.sum(w[m] * p[m]) / wbar[k]
    C_B = np.sum(wbar * (By - my) * (Bp - mp)) / (W * sd)
    C_W = np.sum(w * (y - By[inv]) * (p - Bp[inv])) / (W * sd)
    return float(C_B + C_W), float(C_B), float(C_W)


lab = pd.read_csv(R / "experiments/v7_1_1/phase1_universe/canonical_target_labels.csv",
                  low_memory=False)
TGT = lab.target_id_canonical.to_numpy(str)
zcan = np.load(R / "experiments/v7_1/phase2_rank/z_predictions.npz", allow_pickle=True)
Y_CAN = np.asarray(zcan["y"], float)
mas = pd.read_csv(R / "experiments/v6_2/phase2_split4/matched_analysis_set.csv",
                  low_memory=False)
Y_MAS = pd.to_numeric(mas["ddG"], errors="coerce").to_numpy(float)
TGT_MAS = mas.pdb_id.astype(str).to_numpy()
S = pd.read_csv(R / "experiments/v8/BENCHMARK_SCOREBOARD.csv", low_memory=False)

cache, recs = {}, []
LAB = {"y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows",
       "z_raw", "y_cal", "family", "seed", "masks"}
for r in S.itertuples():
    rel = r.prediction_file
    if rel not in cache:
        try:
            cache[rel] = np.load(R / rel, allow_pickle=True)
        except Exception:  # noqa: BLE001
            cache[rel] = None
    zz = cache[rel]
    if zz is None or r.model_name not in zz.files or r.model_name in LAB \
            or r.model_name.lower().endswith("_y"):
        continue
    a = np.asarray(zz[r.model_name])
    if a.ndim != 1:
        continue
    try:
        a = a.astype(float)
    except (ValueError, TypeError):
        continue
    y = Y_CAN if r.join_mode.startswith("canonical") else Y_MAS
    tgt = TGT if r.join_mode.startswith("canonical") else TGT_MAS
    if len(a) != len(y):
        continue
    m = np.isfinite(y) & np.isfinite(a)
    yy, pp, tt = y[m], a[m], tgt[m]
    if len(set(tt)) < 3:
        continue
    rec = {"model_name": r.model_name, "universe": r.universe}
    for mode in ("P", "T"):
        for s in SCHEMES:
            rr, cb, cw = weighted(yy, pp, tt, s, mode)
            rec[f"{mode}_r_{s}"] = rr
            rec[f"{mode}_C_B_{s}"] = cb
            rec[f"{mode}_C_W_{s}"] = cw
    recs.append(rec)

D = pd.DataFrame(recs)
print(f"配置数 = {len(D)}")
D.to_csv(R / "v8_1/WEIGHTING_SEMANTICS.csv", index=False, encoding="utf-8")

out = {}
for mode, label in (("P", "PER-ROW MULTIPLIER  w=a_c  (代码当前实现)"),
                    ("T", "TARGET TOTAL MASS    w=a_c/n_c (与名称一致)")):
    print()
    print("=" * 86)
    print(f"解释 {mode}: {label}")
    print("=" * 86)
    cbc = [f"{mode}_C_B_{s}" for s in SCHEMES]
    cwc = [f"{mode}_C_W_{s}" for s in SCHEMES]
    rc = [f"{mode}_r_{s}" for s in SCHEMES]
    # 恒等式
    res = max(float((D[f"{mode}_r_{s}"] - D[f"{mode}_C_B_{s}"]
                     - D[f"{mode}_C_W_{s}"]).abs().max()) for s in SCHEMES)
    print(f"  恒等式 max|r-C_B-C_W| = {res:.3e}   全部<=1e-12: {res<=1e-12}")
    rngB = (D[cbc].max(axis=1) - D[cbc].min(axis=1))
    rngW = (D[cwc].max(axis=1) - D[cwc].min(axis=1))
    print(f"  逐配置跨 scheme 极差中位数: C_B = {rngB.median():.5f}   C_W = {rngW.median():.5f}")
    sfB = {str(lo): int((((D[cbc] > lo).sum(axis=1) > 0) &
                         ((D[cbc] < -lo).sum(axis=1) > 0)).sum())
           for lo in (0, 0.001, 0.005, 0.01)}
    sfW = {str(lo): int((((D[cwc] > lo).sum(axis=1) > 0) &
                         ((D[cwc] < -lo).sum(axis=1) > 0)).sum())
           for lo in (0, 0.001, 0.005, 0.01)}
    print(f"  C_B 符号翻转: {sfB}")
    print(f"  C_W 符号翻转: {sfW}")
    for s in SCHEMES:
        v = D[f"{mode}_C_W_{s}"].dropna()
        print(f"    {s:6s} C_W 中位数 = {v.median():+.5f}")
    out[mode] = dict(identity_residual=res, range_C_B_median=float(rngB.median()),
                     range_C_W_median=float(rngW.median()),
                     signflip_C_B=sfB, signflip_C_W=sfW)

# ---- 排名一致性：两种解释下，模型在不同 scheme 间的排名是否一致 ----
print()
print("=" * 86)
print("排名对 scheme 的敏感性（排名一致的方案对）")
print("=" * 86)
for mode in ("P", "T"):
    rc = [f"{mode}_r_{s}" for s in SCHEMES]
    sub = D.dropna(subset=rc)
    print(f"\n  解释 {mode}:")
    for i in range(len(SCHEMES)):
        for j in range(i + 1, len(SCHEMES)):
            rho = spearmanr(sub[rc[i]], sub[rc[j]]).correlation
            print(f"    {SCHEMES[i]:6s} vs {SCHEMES[j]:6s}  rho = {rho:+.4f}")

out["note"] = ("P = 代码实现 w=a_c；T = 与名称一致 w=a_c/n_c。"
               "两者均满足恒等式（因为分子分母同测度）。"
               "差异在于 target 总质量是 n_c*a_c 还是 a_c。")
json.dump(out, open(R / "v8_1/WEIGHTING_SEMANTICS.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
print("\nwrote v8_1/WEIGHTING_SEMANTICS.{csv,json}")
