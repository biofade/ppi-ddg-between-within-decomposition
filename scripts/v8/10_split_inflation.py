# -*- coding: utf-8 -*-
"""scripts/v8/10_split_inflation.py -- V8 PHASE 10 / 核心实验 D：划分膨胀（§35）。

**同一批模型族**在**同一批行**上，只改**划分方式**：
    L1_random   : 随机划分（最宽松）
    L2_complex  : 复合物/结构划分
    L3_seq40    : 序列 40% 同一性簇划分
    L4_family   : family-OOD（最严格）

比较 Overall Pearson / RMSE / Macro ρ / C_W / Efficiency@10 / NDCG@10。
问题（H4）：宽松划分是否**更多**抬高表观全局性能，而**更少**改善真实靶内泛化？

同时给出「能力轴一致性」随划分严格度的变化（呼应 checkpoint 的发现）。
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
sys.path.insert(0, str(R / "scripts" / "v8"))
import v8metrics as M  # noqa: E402

OUT = R / "experiments" / "v8"
SI = OUT / "split_inflation"
SI.mkdir(parents=True, exist_ok=True)

mas = pd.read_csv(R / "experiments" / "v6_2" / "phase2_split4" /
                  "matched_analysis_set.csv", low_memory=False)
TGT = mas.pdb_id.astype(str).to_numpy()
FAM = mas.family_label.astype(str).to_numpy()
CLU = mas.seq40_cluster.astype(str).to_numpy()
Y = pd.to_numeric(mas["ddG"], errors="coerce").to_numpy(float)
print(f"U_MAS rows={len(Y)} targets={len(set(TGT))} families={len(set(FAM))} "
      f"seq40_clusters={len(set(CLU))}")

z = np.load(R / "experiments" / "v6_2" / "phase2_split4" /
            "split4_predictions.npz", allow_pickle=True)
keys = list(z.files)
LEVELS = ["L1_random", "L2_complex", "L3_seq40", "L4_family"]
LEVEL_LABEL = {"L1_random": "random (最宽松)",
               "L2_complex": "complex/structure",
               "L3_seq40": "seq40 cluster",
               "L4_family": "family-OOD (最严格)"}
# 每个 level 的 y / pdb 键（用于确认是同一批行）
print("\n=== 每个 level 的键 ===")
for lv in LEVELS:
    lk = [k for k in keys if k.startswith(lv + "_")]
    print(f"  {lv}: {len(lk)} 键 -> {lk[:8]}")
    yk = f"{lv}_y"
    if yk in keys:
        yy = np.asarray(z[yk], float)
        same = np.allclose(yy, Y, equal_nan=True)
        print(f"      {yk}: 与 matched_analysis_set 的 ddG 一致 = {same}")

rows = []
for lv in LEVELS:
    lk = [k for k in keys if k.startswith(lv + "_")]
    models = [k for k in lk if not k.endswith("_y") and not k.endswith("_pdb")]
    for mk in models:
        a = np.asarray(z[mk])
        try:
            a = a.astype(float)
        except (ValueError, TypeError):
            continue
        if a.ndim != 1 or len(a) != len(Y):
            continue
        m = np.isfinite(Y) & np.isfinite(a)
        if m.sum() < 100:
            continue
        U = M.prep_universe(Y[m], a[m], TGT[m], FAM[m])
        res = M.compute_all(U)
        res.update({"level": lv, "level_label": LEVEL_LABEL[lv], "model": mk})
        rows.append(res)

S = pd.DataFrame(rows)
S.to_csv(SI / "split_inflation_raw.csv", index=False, encoding="utf-8")
print(f"\n评价 {len(S)} 个 (level, model) 组合")

COLS = ["overall_pearson", "rmse", "macro_spearman", "C_W",
        "efficiency@10", "ndcg@10", "target_mean_pearson", "C_B_row"]
print("\n=== 每个 level 的中位数 ===")
agg = S.groupby("level")[COLS].median()
print(agg.round(4).to_string())

# 以 L1_random 为基准的膨胀幅度
base = agg.loc["L1_random"]
print("\n=== 相对 random 的变化（中位数差）===")
delta = agg - base
print(delta.round(4).to_string())

# 能力轴一致性随 level 变化
print("\n=== 能力轴秩一致性（Overall r vs Macro ρ / vs Eff@10）===")
cons = []
for lv, g in S.groupby("level"):
    g2 = g.dropna(subset=["overall_pearson", "macro_spearman", "efficiency@10"])
    if len(g2) < 5:
        continue
    s1 = spearmanr(g2.overall_pearson, g2.macro_spearman).correlation
    k1 = kendalltau(g2.overall_pearson, g2.macro_spearman).correlation
    s2 = spearmanr(g2.overall_pearson, g2["efficiency@10"]).correlation
    k2 = kendalltau(g2.overall_pearson, g2["efficiency@10"]).correlation
    cons.append({"level": lv, "n_models": len(g2),
                 "spearman_overall_vs_macro": s1, "kendall_overall_vs_macro": k1,
                 "spearman_overall_vs_eff10": s2, "kendall_overall_vs_eff10": k2})
    print(f"  {lv:12s} n={len(g2):3d}  r~macro: ρ={s1:+.3f} τ={k1:+.3f}   "
          f"r~eff10: ρ={s2:+.3f} τ={k2:+.3f}")
C = pd.DataFrame(cons)
C.to_csv(SI / "capability_consistency_by_level.csv", index=False, encoding="utf-8")

# 同族模型跨 level 对比（相同 base 名）
fam_ct = [c for c in S.columns if c == "level"]
S["model_base"] = S.model.str.replace(r"^(L\d_[a-z0-9]+_)(.*)_s\d+$", r"\2", regex=True)
piv = S.pivot_table(index="model_base", columns="level",
                    values=["overall_pearson", "macro_spearman"], aggfunc="median")
piv.to_csv(SI / "model_level_pivot.csv", encoding="utf-8")
print("\n=== 同族模型：Overall r 与 Macro ρ 跨 level（前 12 行）===")
print(piv.head(12).round(4).to_string())

# 报告
lines = ["# V8 Split-Inflation Report（§35 / H4）\n",
         "> **同一批行（U_MAS 7523）、同一批模型族**，只改划分方式。\n",
         "> 划分：`L1_random`（宽松）→ `L2_complex` → `L3_seq40` → `L4_family`（严格）\n",
         "## 1. 各划分下的中位指标\n", "```", agg.round(4).to_string(), "```\n",
         "## 2. 相对 random 的变化（中位数差）\n", "```", delta.round(4).to_string(), "```\n",
         "## 3. 能力轴一致性随划分严格度的变化\n",
         "| level | n_models | ρ(r, macro ρ) | τ(r, macro ρ) | ρ(r, Eff@10) | τ(r, Eff@10) |",
         "|---|---:|---:|---:|---:|---:|"]
for r in C.itertuples():
    lines.append(f"| `{r.level}` | {r.n_models} | {r.spearman_overall_vs_macro:+.3f} | "
                 f"{r.kendall_overall_vs_macro:+.3f} | {r.spearman_overall_vs_eff10:+.3f} | "
                 f"{r.kendall_overall_vs_eff10:+.3f} |")
(SI.parent / "SPLIT_INFLATION_REPORT.md").write_text("\n".join(lines) + "\n",
                                                     encoding="utf-8")
print("\n写出 split_inflation_raw.csv / capability_consistency_by_level.csv / "
      "model_level_pivot.csv / SPLIT_INFLATION_REPORT.md")
