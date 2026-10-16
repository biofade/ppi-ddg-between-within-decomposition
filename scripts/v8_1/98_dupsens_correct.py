# -*- coding: utf-8 -*-
"""scripts/v8_1/98_dupsens_correct.py -- 重复键敏感性（用与冻结源一致的分解）。

先验证：本脚本的分解能否逐位复现 all_models_by_rowset.csv 的
C_B_row / C_W_row（V6_2_split4, 7,523 行）。只有验证通过，才把它
用于 FIRST / MEAN 两种去重口径。

分解（未加权，逐配置精确）：
    N = 行数,  n_c = target c 的行数
    m_y = mean(y), m_p = mean(p),  sigma = std(y)*std(p)   (除以 N)
    B_y,c = mean(y in c), B_p,c = mean(p in c)
    C_B = sum_c (n_c/N) (B_y,c - m_y)(B_p,c - m_p) / sigma
    C_W = sum_c (n_c/N) mean_i[(y-B_y,c)(p-B_p,c)] / sigma
    r   = C_B + C_W
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(".")


def decomp(y, p, tgt, min_n_macro=5):
    y = np.asarray(y, float); p = np.asarray(p, float)
    t = np.asarray(tgt).astype(str)
    uq, inv = np.unique(t, return_inverse=True)
    n_c = np.bincount(inv).astype(float)
    N = len(y)
    my, mp = y.mean(), p.mean()
    sd = y.std() * p.std()          # population sd (ddof=0)
    By = np.array([y[inv == k].mean() for k in range(len(uq))])
    Bp = np.array([p[inv == k].mean() for k in range(len(uq))])
    C_B = np.sum((n_c / N) * (By - my) * (Bp - mp)) / sd
    C_W = np.sum((n_c / N) * np.array(
        [(y[inv == k] - By[k]).dot(p[inv == k] - Bp[k]) / n_c[k]
         for k in range(len(uq))])) / sd
    sp = [spearmanr(y[inv == k], p[inv == k]).correlation
          for k in range(len(uq)) if n_c[k] >= min_n_macro]
    sp = [x for x in sp if np.isfinite(x)]
    return float(C_B + C_W), float(C_B), float(C_W), (float(np.mean(sp)) if sp else np.nan)


# ---------------- 验证阶段 ----------------
mas = pd.read_csv(R / "experiments/v6_2/phase2_split4/matched_analysis_set.csv",
                  low_memory=False)
mas["key"] = mas.pdb_id.astype(str) + "|" + mas.mutation.astype(str)
y_all = pd.to_numeric(mas["ddG"], errors="coerce").to_numpy(float)
tgt_all = mas.pdb_id.astype(str).to_numpy()
groups = mas.groupby("key").indices
keys = list(groups.keys())
first_idx = np.array([groups[k][0] for k in keys])
agg_rows = [groups[k] for k in keys]
y_key_mean = np.array([y_all[r].mean() for r in agg_rows])
tgt_key = np.array([tgt_all[r[0]] for r in agg_rows])

z = np.load(R / "experiments/v6_2/phase2_split4/split4_predictions.npz",
            allow_pickle=True)
disk = pd.read_csv(R / "v8_1/all_models_by_rowset.csv", low_memory=False)
disk = disk[(disk.n_rows == 7523) & (disk.artifact == "V6_2_split4")] \
    .drop_duplicates("key").set_index("key")

print("=" * 86)
print("验证：本脚本分解 vs all_models_by_rowset.csv")
print("=" * 86)
mx = 0.0
nchk = 0
for k in sorted(disk.index):
    if k not in z.files:
        continue
    p_all = np.asarray(z[k], float)
    rr, cb, cw, _ = decomp(y_all, p_all, tgt_all)
    ref = disk.loc[k]
    d = max(abs(cb - ref.C_B_row), abs(cw - ref.C_W_row), abs(rr - ref.overall_pearson))
    mx = max(mx, d); nchk += 1
print(f"  比对配置 = {nchk}   最大绝对差 = {mx:.3e}")
print("  => 分解与冻结源一致" if mx < 1e-9 else "  => 不一致，需继续核对")

# ---------------- 敏感性阶段 ----------------
def fold_mean(p_all):
    return np.array([p_all[r].mean() for r in agg_rows])


rows = []
for level, regime in (("L1", "Random"), ("L4", "Family-OOD")):
    names = sorted(k for k in z.files
                   if k.startswith(level + "_") and not k.endswith(("_y", "_pdb")))
    for nm in names:
        p_all = np.asarray(z[nm], float)
        for mode in ("ALL", "FIRST", "MEAN"):
            if mode == "ALL":
                rr, cb, cw, ms = decomp(y_all, p_all, tgt_all); n = len(y_all)
            elif mode == "FIRST":
                rr, cb, cw, ms = decomp(y_all[first_idx], p_all[first_idx],
                                        tgt_all[first_idx]); n = len(first_idx)
            else:
                rr, cb, cw, ms = decomp(y_key_mean, fold_mean(p_all), tgt_key)
                n = len(keys)
            rows.append(dict(regime=regime, model=nm,
                             pairkey="".join(nm.split("_")[-2:]),
                             mode=mode, n=n, r=rr, C_B=cb, C_W=cw, macro=ms))
D = pd.DataFrame(rows)

print()
print("=" * 86)
print("三种口径的中位数")
print("=" * 86)
for regime in ("Random", "Family-OOD"):
    print(f"\n--- {regime} ---")
    for mode in ("ALL", "FIRST", "MEAN"):
        s = D[(D.regime == regime) & (D["mode"] == mode)]
        print(f"  {mode:6s} n={s.n.iloc[0]:5d}  median r={s.r.median():+.4f}  "
              f"C_B={s.C_B.median():+.4f}  C_W={s.C_W.median():+.4f}  "
              f"macro={s.macro.median():.4f}")

print()
print("=" * 86)
print("Random -> Family-OOD 配对归因")
print("=" * 86)
summ = {}
for mode in ("ALL", "FIRST", "MEAN"):
    a = D[(D.regime == "Random") & (D["mode"] == mode)].set_index("pairkey")
    b = D[(D.regime == "Family-OOD") & (D["mode"] == mode)].set_index("pairkey")
    com = a.index.intersection(b.index)
    dr = (b.loc[com, "r"] - a.loc[com, "r"]).median()
    dcb = (b.loc[com, "C_B"] - a.loc[com, "C_B"]).median()
    dcw = (b.loc[com, "C_W"] - a.loc[com, "C_W"]).median()
    ratio = abs(dcb) / (abs(dcb) + abs(dcw))
    summ[mode] = dict(n_pairs=int(len(com)), d_r=float(dr), d_C_B=float(dcb),
                      d_C_W=float(dcw), attribution=float(ratio),
                      median_r_random=float(a.loc[com, "r"].median()),
                      median_r_ood=float(b.loc[com, "r"].median()))
    print(f"  {mode:6s} n={len(com):2d}  Δr={dr:+.4f}  ΔC_B={dcb:+.4f}  ΔC_W={dcw:+.4f}"
          f"  归因={100*ratio:.1f}%")

ok = all(v["attribution"] > 0.5 and v["d_C_B"] < 0 for v in summ.values())
print()
print(f"  三口径下 ΔC_B 均为负 且 归因均 > 50%: {ok}")
print(f"  归因比例范围: {min(v['attribution'] for v in summ.values())*100:.1f}%"
      f" .. {max(v['attribution'] for v in summ.values())*100:.1f}%")

json.dump(dict(decomp_verified_max_abs_diff=float(mx), n_checked=int(nchk),
               by_mode=summ, direction_preserved=bool(ok),
               caveat=("测试集内部敏感性分析；预测 artifact 未保存 fold 归属，"
                       "无法重跑 split 层面的键级分组")),
          open(R / "v8_1/DUPLICATE_KEY_SENSITIVITY.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
D.to_csv(R / "v8_1/DUPLICATE_KEY_SENSITIVITY.csv", index=False, encoding="utf-8")
print("\nwrote v8_1/DUPLICATE_KEY_SENSITIVITY.{json,csv}")
