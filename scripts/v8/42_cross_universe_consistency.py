# -*- coding: utf-8 -*-
"""scripts/v8/42_cross_universe_consistency.py -- 跨宇宙一致性：外部方法在 3 个宇宙上的横向对照。

## 为什么这比"再下载一个模型"更有价值

V8 已建立 **3 个外部宇宙**（不同行集合、不同 target 定义、不同来源）：
  * `U_HERMES_CANONICAL`  —— HERMES 团队复算，target = canonical 家族
  * `U_EXT_CANONICAL`     —— CATH-ddG 表，target = canonical 家族（单点匹配）
  * `U_EXT_HELDOUT_CATH`  —— CATH-ddG 官方 held-out，target = `complex_PPI`

同一个外部方法在多个宇宙上被**独立评价** ⇒ 可以检验
**「方法的相对强弱排序是否跨宇宙稳定」**。

**这是一个真正的稳健性检验**（§30），且完全不依赖新资源。

**纪律**：不跨宇宙比较**绝对数值**（§21）；只比较**排序**（秩一致）。
"""

from __future__ import annotations

import json
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
EXT = OUT / "external"

# ---------- 收集 3 个宇宙的读数 ----------
uni = {}

p1 = EXT / "external_expanded_pool.csv"
if p1.exists():
    d = pd.read_csv(p1)
    d["universe"] = "U_HERMES_CANONICAL"
    d = d.rename(columns={"model": "method"})
    # 把变体归并到方法族（按方法族更保守）
    d["method"] = d["family"]
    d = d.groupby(["method", "universe"], as_index=False).agg(
        overall_pearson=("overall_pearson", "median"),
        macro_spearman=("macro_spearman", "median"),
        efficiency=("efficiency@10", "median"),
        n=("n_rows", "max"))
    uni["U_HERMES_CANONICAL"] = d[["method", "universe", "overall_pearson",
                                   "macro_spearman", "efficiency", "n"]]

p2 = EXT / "external_multimodel_scoreboard.csv"
if p2.exists():
    d = pd.read_csv(p2).rename(columns={"efficiency@10": "efficiency",
                                       "model": "method"})
    d["universe"] = "U_EXT_CANONICAL"
    d = d[["method", "universe", "overall_pearson", "macro_spearman",
           "efficiency", "n_rows"]].rename(columns={"n_rows": "n"})
    uni["U_EXT_CANONICAL"] = d

p3 = EXT / "external_heldout_cath.csv"
if p3.exists():
    d = pd.read_csv(p3).rename(columns={"efficiency@10": "efficiency",
                                       "model": "method"})
    d["universe"] = "U_EXT_HELDOUT_CATH"
    d = d.rename(columns={"model": "method", "n_rows": "n"})
    d = d[["method", "universe", "overall_pearson", "macro_spearman",
           "efficiency", "n"]]
    uni["U_EXT_HELDOUT_CATH"] = d

print(f"收集到 {len(uni)} 个外部宇宙")
for k, v in uni.items():
    print(f"  {k:<24} 方法数 = {len(v)}  方法: {sorted(v.method.unique())}")

A = pd.concat(uni.values(), ignore_index=True)
A.to_csv(EXT / "cross_universe_external.csv", index=False, encoding="utf-8")

# ---------- 方法 × 宇宙 矩阵 ----------
print("\n" + "=" * 78)
print("Overall Pearson（方法 × 宇宙）—— 只读排序，不比绝对值（§21）")
print("=" * 78)
P = A.pivot_table(index="method", columns="universe",
                  values="overall_pearson")
M = A.pivot_table(index="method", columns="universe",
                  values="macro_spearman")
print(P.round(4).to_string())
print("\nMacro Spearman（方法 × 宇宙）")
print(M.round(4).to_string())

# ---------- 跨宇宙秩一致 ----------
print("\n" + "=" * 78)
print("跨宇宙方法排序一致性（仅在两个宇宙都有读数的方法上）")
print("=" * 78)
cols = list(P.columns)
rows = []
for i in range(len(cols)):
    for j in range(i + 1, len(cols)):
        a, b = cols[i], cols[j]
        for nm, T in (("Overall r", P), ("Macro rho", M)):
            sub = T[[a, b]].dropna()
            if len(sub) < 3:
                print(f"  {nm:<10} {a} vs {b}: 共同方法仅 {len(sub)} ⇒ 无法计算")
                continue
            sp = spearmanr(sub[a], sub[b]).correlation
            kt = kendalltau(sub[a], sub[b]).correlation
            print(f"  {nm:<10} {a:<22} vs {b:<22} n={len(sub)}  "
                  f"rho={sp:+.4f}  tau={kt:+.4f}")
            rows.append({"metric": nm, "universe_a": a, "universe_b": b,
                         "n_common": len(sub), "spearman": float(sp),
                         "kendall": float(kt)})
pd.DataFrame(rows).to_csv(EXT / "cross_universe_agreement.csv", index=False,
                          encoding="utf-8")

# ---------- 逐方法跨宇宙画像 ----------
print("\n" + "=" * 78)
print("逐方法跨宇宙画像（同一方法在几个宇宙上被评价）")
print("=" * 78)
cnt = A.groupby("method").universe.nunique().sort_values(ascending=False)
for m, k in cnt.items():
    sub = A[A.method == m]
    if k >= 2:
        line = f"  {m:<14} 出现在 {k} 个宇宙: "
        line += "  ".join(f"{r.universe[:14]}: r={r.overall_pearson:+.3f}"
                          for r in sub.itertuples())
        print(line)
print()
print("  ⚠ 只有 1 个宇宙读数的方法无法做跨宇宙检验：")
print("   ", sorted(cnt[cnt == 1].index.tolist()))

meta = {"universes": list(uni.keys()),
        "methods_per_universe": {k: sorted(v.method.unique())
                                 for k, v in uni.items()},
        "cross_universe_pairs": rows,
        "discipline": ("只比较方法**排序**，不比较绝对数值（§21 不同宇宙不可比）"),
        "caveat": ("三个宇宙的 target 定义不同（canonical 家族 vs complex_PPI）、"
                   "行数不同（4473 / 3822 / 813）、模型数不同 ⇒ "
                   "跨宇宙秩一致只能作为**弱**证据")}
(EXT / "cross_universe_meta.json").write_text(
    json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"\n写出 cross_universe_external.csv / cross_universe_agreement.csv / _meta.json")
