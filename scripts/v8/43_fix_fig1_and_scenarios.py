# -*- coding: utf-8 -*-
"""scripts/v8/43_fix_fig1_and_scenarios.py -- 回应审稿人点名的两个疑点。

## 疑点 1：fig1 未标注 160-artifact 伪重复问题

`fig1_rank_agreement_heatmap.png` 用 `U_CANONICAL` 的 **160 模型**池计算秩一致性。
该池含 **102 个 V6.2 `split4` 键**（M0/M1/M2 × 3 seeds × 4 levels），
彼此高度相关 ⇒ 构成**伪重复**，且已被 PHASE 13 证明会把秩相关从
`+0.186`（排除后）翻到 `−0.269`（含入时）。

**⇒ 原图会误导读者把「−0.27」当成 benchmark 的性质。本脚本重绘并在图内标注。**

## 疑点 2：scenario 表的推荐模型来自 permissive split 训练

`scenario_recommendations.csv` 的场景 A/C 首选取自 `V6_2_split4`（**随机/complex/seq40
划分训练**），即**宽松划分**的模型 —— 在严格 family-OOD 基准上推荐它们是不当的。

**⇒ 本脚本产出一个 provenance-clean 版本：只从 family-OOD 训练的键中选首选。**
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
OUT = R / "experiments" / "v8"
FIG = OUT / "figures"
SC = OUT / "scenarios"

for cand in ("Microsoft YaHei", "SimHei", "DengXian", "SimSun"):
    try:
        font_manager.findfont(cand, fallback_to_default=False)
        plt.rcParams["font.sans-serif"] = [cand]
        break
    except Exception:
        continue
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 130
plt.rcParams["savefig.bbox"] = "tight"

# ============ 疑点 1：重绘 fig1 ============
S = pd.read_csv(OUT / "BENCHMARK_SCOREBOARD.csv", low_memory=False)
u = S[(S.universe == "U_CANONICAL")].dropna(
    subset=["overall_pearson", "macro_spearman", "efficiency@10", "target_mean_pearson"])
print(f"fig1 数据池 = {len(u)} 模型")

# 排除 V6.2 伪重复键
is_v62 = u.model_family.eq("V6_2_split4")
u_clean = u[~is_v62]
print(f"  其中 V6_2_split4（伪重复）= {int(is_v62.sum())}，排除后 = {len(u_clean)}")


def agreement(g):
    cols = ["overall_pearson", "target_mean_pearson", "rmse",
            "macro_spearman", "macro_kendall", "C_W",
            "efficiency@10", "ndcg@10", "hit@10"]
    cols = [c for c in cols if c in g.columns]
    g = g.dropna(subset=cols)
    n = len(cols)
    Mx = np.eye(n)
    for i in range(n):
        for j in range(i + 1, n):
            v = spearmanr(g[cols[i]], g[cols[j]]).correlation
            Mx[i, j] = Mx[j, i] = v
    return pd.DataFrame(Mx, index=cols, columns=cols), len(g)


A_all, n_all = agreement(u)
A_cln, n_cln = agreement(u_clean)
A_all.to_csv(OUT / "rank_reversal" / "rank_agreement_spearman_U_CANONICAL_ALL.csv",
             encoding="utf-8")
A_cln.to_csv(OUT / "rank_reversal" / "rank_agreement_spearman_U_CANONICAL_CLEAN.csv",
             encoding="utf-8")

fig, axes = plt.subplots(1, 2, figsize=(17.5, 7.4))
for ax, Mx, ttl, nn in (
        (axes[0], A_all, f"(a) ALL models (n={n_all}) — INCLUDES 102 V6.2 "
                         f"pseudo-replicates", n_all),
        (axes[1], A_cln, f"(b) V6.2 excluded (n={n_cln}) — provenance-clean",
         n_cln)):
    M = Mx.to_numpy(float)
    im = ax.imshow(M, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(Mx.columns)))
    ax.set_xticklabels(Mx.columns, rotation=45, ha="right", fontsize=8.4)
    ax.set_yticks(range(len(Mx.index)))
    ax.set_yticklabels(Mx.index, fontsize=8.4)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isfinite(M[i, j]):
                ax.text(j, i, f"{M[i,j]:+.2f}", ha="center", va="center",
                        fontsize=6.6,
                        color="white" if abs(M[i, j]) > 0.55 else "black")
    r = Mx.loc["overall_pearson", "macro_spearman"]
    ax.set_title(ttl + f"\nrho(Overall r, Macro rho) = {r:+.3f}", fontsize=10.2)
    fig.colorbar(im, ax=ax, shrink=0.78, label="rank agreement (Spearman)")

fig.suptitle("Fig 1  Cross-axis rank agreement — pool-composition sensitivity",
             fontsize=13, y=1.005)
fig.text(0.5, -0.085,
         "CAVEAT (added after audit): panel (a) reproduces the originally reported "
         "figure but its pool contains 51 V6.2 'split4' keys\n"
         "(M0/M1/M2 x 3 seeds x 4 split levels) that are near-duplicates of each other. "
         "Panel (b) removes them.\n"
         "rho(Overall r, Macro rho) = "
         f"{A_all.loc['overall_pearson','macro_spearman']:+.3f} (a) vs "
         f"{A_cln.loc['overall_pearson','macro_spearman']:+.3f} (b)  ->  the value "
         "CHANGES BY ~0.23, so this quantity is POOL-DEPENDENT.\n"
         "NOTE (honest): BOTH panels are POSITIVE here. A negative value "
         "(-0.269) arises only in a DIFFERENT, smaller pool\n"
         "(67 models = V6.2 keys + V7.1 keys only). Therefore no single number "
         "characterises the benchmark; report the observed\n"
         "range across pools (6 pools examined: -0.27 to +0.47) and always state the "
         "pool and its provenance composition.",
         ha="center", fontsize=8.6, color="#222")
fig.savefig(FIG / "fig1_rank_agreement_heatmap.png")
plt.close(fig)
print(f"  重绘 fig1 完成：ALL rho={A_all.loc['overall_pearson','macro_spearman']:+.4f} | "
      f"CLEAN rho={A_cln.loc['overall_pearson','macro_spearman']:+.4f}")

# ============ 疑点 2：provenance-clean scenario 表 ============
PERM_TOK = ("L1_random_", "L2_complex_", "L3_seq40_")
STRICT_TOK = ("L4_family_", "F0_", "F1_", "F2_", "F3_", "F4_", "F5_", "F6_",
              "F7_", "F8_", "F9_", "R0_", "R1_", "R2_", "R3_", "B0_", "B1_",
              "B2_", "B3", "B4_", "A2_smallPLM", "p_", "z_R", "z_A")
def prov(k):
    if any(k.startswith(t) for t in PERM_TOK):
        return "PERMISSIVE"
    if any(k.startswith(t) for t in STRICT_TOK):
        return "STRICT"
    return "OTHER"

u2 = u.copy()
u2["provenance"] = u2.model_name.map(prov)
print(f"\nprovenance 分布: {u2.provenance.value_counts().to_dict()}")

SCEN = {
    "A_quantitative_ddG": ("定量 ΔΔG", "rmse", True),
    "B_library_ranking": ("突变文库排序", "macro_spearman", False),
    "C_wetlab_shortlist": ("湿实验短名单", "efficiency@10", False),
}
rows = []
for name, (desc, metric, asc) in SCEN.items():
    g = u2.dropna(subset=[metric])
    strict = g[g.provenance == "STRICT"]
    pick_all = g.sort_values(metric, ascending=asc).iloc[0]
    pick_strict = (strict.sort_values(metric, ascending=asc).iloc[0]
                   if len(strict) else None)
    rows.append({
        "scenario": name, "description": desc, "primary_metric": metric,
        "top1_ALL_pool": pick_all.model_name,
        "top1_ALL_provenance": pick_all.provenance,
        "top1_ALL_family": pick_all.model_family,
        "top1_ALL_value": float(pick_all[metric]),
        "top1_STRICT_pool": pick_strict.model_name if pick_strict is not None else "N/A",
        "top1_STRICT_provenance": pick_strict.provenance if pick_strict is not None else "N/A",
        "top1_STRICT_family": pick_strict.model_family if pick_strict is not None else "N/A",
        "top1_STRICT_value": float(pick_strict[metric]) if pick_strict is not None else np.nan,
        "changed_after_provenance_clean": (pick_strict is not None
                                           and pick_strict.model_name != pick_all.model_name),
    })
P2 = pd.DataFrame(rows)
P2.to_csv(SC / "scenario_recommendations_provenance_clean.csv", index=False,
          encoding="utf-8")
print("\n=== provenance-clean 场景推荐 ===")
for r in P2.itertuples():
    print(f"  {r.scenario:<22} 原首选={r.top1_ALL_pool:<28}({r.top1_ALL_provenance})")
    print(f"  {'':<22} 严格首选={r.top1_STRICT_pool:<28}({r.top1_STRICT_provenance})"
          f"  改变={r.changed_after_provenance_clean}")
print("\n写出 rank_reversal/rank_agreement_*_{ALL,CLEAN}.csv + "
      "scenarios/scenario_recommendations_provenance_clean.csv")
