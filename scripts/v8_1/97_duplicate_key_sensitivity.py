# -*- coding: utf-8 -*-
"""scripts/v8_1/97_duplicate_key_sensitivity.py -- 重复 (pdb, mutation) 键的敏感性分析。

背景：7,523 行只对应 4,812 个不同 (pdb, mutation) 键；2,011 个键有多条实验记录，
其中 1,405 个键的 ΔΔG 不一致。预测目录没有保存 fold 归属，因此**无法**重跑 split
来检验 duplicate-key leakage。但可以检验一个等价的关键问题：

    结论是否依赖「同一键的多条记录在测试集里被重复计入」？

做法：在存储的测试预测上（7,523 行）比较三种口径
    ALL   : 全部 7,523 行
    FIRST : 每个键只保留第一条记录（4,812 行）
    MEAN  : 每个键的 y 与 p 都取组内均值（4,812 行，键级完全聚合）
对 L1(random) 与 L4(family-OOD) 两档，计算 r / C_B / C_W / macro Spearman，
并检查 81% 归因方向是否保持。

注意：这是**测试集内部**的敏感性分析，不能替代重新生成 split；
若审稿需要严格结论，必须重跑带键级分组的训练/评估（见报告）。
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(".")

mas = pd.read_csv(R / "experiments/v6_2/phase2_split4/matched_analysis_set.csv",
                  low_memory=False)
mas["key"] = mas.pdb_id.astype(str) + "|" + mas.mutation.astype(str)
y_all = pd.to_numeric(mas["ddG"], errors="coerce").to_numpy(float)
tgt_all = mas.pdb_id.astype(str).to_numpy()
key_all = mas.key.to_numpy(str)

# 键 -> 该键在 7523 行中的行号列表
groups = mas.groupby("key").indices
keys = list(groups.keys())
first_idx = np.array([groups[k][0] for k in keys])
# 键级代表行（用于把 (7523,) 的预测折叠为 (4812,)）
agg_rows = [groups[k] for k in keys]
y_key_mean = np.array([y_all[r].mean() for r in agg_rows])
tgt_key = np.array([tgt_all[r[0]] for r in agg_rows])

z = np.load(R / "experiments/v6_2/phase2_split4/split4_predictions.npz",
            allow_pickle=True)


def decomp(y, p, tgt):
    """未加权的 r = C_B + C_W（逐配置精确）。

    注意：target 索引必须先在**全量行**上建立，再按有限性掩码截取，
    否则掩码会平移 target 标签（与 09_weighting_audit 的做法保持一致）。
    """
    t = np.asarray(tgt).astype(str)
    uq, inv = np.unique(t, return_inverse=True)
    y = np.asarray(y, float); p = np.asarray(p, float)
    m = np.isfinite(y) & np.isfinite(p)
    y, p, inv = y[m], p[m], inv[m]
    my, mp = y.mean(), p.mean()
    vy, vp = y.var(), p.var()
    sd = np.sqrt(vy * vp)
    r = np.cov(y, p, bias=True)[0, 1] / sd
    By = np.array([y[inv == k].mean() for k in range(len(uq))])
    Bp = np.array([p[inv == k].mean() for k in range(len(uq))])
    C_B = np.cov(By, Bp, bias=True)[0, 1] / sd
    C_W = r - C_B
    sp = [spearmanr(y[inv == k], p[inv == k]).correlation
          for k in range(len(uq)) if (inv == k).sum() >= 5]
    sp = [x for x in sp if np.isfinite(x)]
    return float(r), float(C_B), float(C_W), float(np.mean(sp)) if sp else np.nan


def fold_mean(p_all):
    return np.array([p_all[r].mean() for r in agg_rows])


rows = []
for level, regime in (("L1", "Random"), ("L4", "Family-OOD")):
    names = [k for k in z.files if k.startswith(level) and k.endswith(tuple(
        f"_s{s}" for s in (42, 43, 44)))]
    for nm in sorted(names):
        p_all = np.asarray(z[nm], float)
        for mode in ("ALL", "FIRST", "MEAN"):
            if mode == "ALL":
                rr, cb, cw, ms = decomp(y_all, p_all, tgt_all)
                n = len(y_all)
            elif mode == "FIRST":
                rr, cb, cw, ms = decomp(y_all[first_idx], p_all[first_idx],
                                        tgt_all[first_idx])
                n = len(first_idx)
            else:
                rr, cb, cw, ms = decomp(y_key_mean, fold_mean(p_all), tgt_key)
                n = len(keys)
            rows.append(dict(regime=regime, model=nm,
                             pairkey=nm.split("_", 1)[1],
                             mode=mode, n=n, r=rr, C_B=cb, C_W=cw, macro=ms))

D = pd.DataFrame(rows)
print("=" * 88)
print("重复键敏感性：三种口径的对比")
print("=" * 88)
for regime in ("Random", "Family-OOD"):
    g = D[D.regime == regime]
    print(f"\n--- {regime} ---")
    piv = g.pivot_table(index="model", columns="mode", values=["r", "C_B", "C_W"])
    for mode in ("ALL", "FIRST", "MEAN"):
        s = g[g["mode"] == mode]
        print(f"  {mode:6s} n={s.n.iloc[0]:5d}  median r={s.r.median():+.4f}  "
              f"median C_B={s.C_B.median():+.4f}  median C_W={s.C_W.median():+.4f}  "
              f"median macro={s.macro.median():.4f}")

print()
print("=" * 88)
print("81% 归因方向是否保持（Random -> Family-OOD，同 model 配对取中位数）")
print("=" * 88)
summ = {}
for mode in ("ALL", "FIRST", "MEAN"):
    a = D[(D.regime == "Random") & (D["mode"] == mode)].set_index("pairkey")
    b = D[(D.regime == "Family-OOD") & (D["mode"] == mode)].set_index("pairkey")
    com = a.index.intersection(b.index)
    dcb = (b.loc[com, "C_B"] - a.loc[com, "C_B"]).median()
    dcw = (b.loc[com, "C_W"] - a.loc[com, "C_W"]).median()
    dr = (b.loc[com, "r"] - a.loc[com, "r"]).median()
    ratio = abs(dcb) / (abs(dcb) + abs(dcw))
    summ[mode] = dict(n_pairs=len(com), d_r=float(dr), d_C_B=float(dcb),
                      d_C_W=float(dcw), attribution=float(ratio))
    print(f"  {mode:6s} n_pairs={len(com):2d}  Δr={dr:+.4f}  ΔC_B={dcb:+.4f}  "
          f"ΔC_W={dcw:+.4f}  归因比例={ratio:.4f} ({100*ratio:.1f}%)")
print()
ok = all(v["attribution"] > 0.5 for v in summ.values())
print(f"  三种口径下 between-target 归因均 > 50%: {ok}")
print(f"  归因比例范围: {min(v['attribution'] for v in summ.values()):.4f}"
      f" .. {max(v['attribution'] for v in summ.values()):.4f}")

json.dump(dict(by_mode=summ, all_gt_half=bool(ok),
               note=("测试集内部敏感性分析；预测 artifact 未保存 fold 归属，"
                     "无法重跑 split 层面的键级分组")),
          open(R / "v8_1/DUPLICATE_KEY_SENSITIVITY.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
D.to_csv(R / "v8_1/DUPLICATE_KEY_SENSITIVITY.csv", index=False, encoding="utf-8")
print("\nwrote v8_1/DUPLICATE_KEY_SENSITIVITY.{json,csv}")
