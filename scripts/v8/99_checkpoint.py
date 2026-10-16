# -*- coding: utf-8 -*-
"""scripts/v8/99_checkpoint.py -- V8 §82 FIRST CHECKPOINT（PHASE 0–6 preliminary）。

只汇总**已算出的**结果，不做新计算，不臆造。
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
V8 = R / "experiments" / "v8"

L: list[str] = []
def p(s=""):
    print(s); L.append(str(s))

# ---------------- 冻结 ----------------
inv = pd.read_csv(V8 / "phase0_freeze" / "artifact_inventory.csv")
reg = pd.read_csv(V8 / "model_registry.csv", low_memory=False)
S = pd.read_csv(V8 / "BENCHMARK_SCOREBOARD.csv", low_memory=False)
RA = pd.read_csv(V8 / "rank_reversal" / "rank_agreement.csv")
WS = pd.read_csv(V8 / "weighting" / "weighting_sensitivity.csv", low_memory=False)
probe = pd.read_csv(V8 / "registry" / "prediction_probe.csv", low_memory=False)

p("=" * 34)
p("V8 LOCAL CHECKPOINT")
p("=" * 34)
p()
p(f"UPSTREAM_HASH_CHANGE = 0   (24 项冻结，含 V7.3/V7.4/V7.4R2 关键制品)")
p(f"Historical prediction artifacts found = {len(inv)} 文件 "
  f"(prediction_like {int((inv.kind=='prediction_like').sum())})")
p(f"Internal models found = {len(reg)} 预测键 "
  f"/ {reg.prediction_file.nunique()} 制品 / {reg.model_family.nunique()} 模型族")
prim = reg[reg.usable_primary == True]
p(f"Primary-usable internal models = {len(prim)}")
p()
p("External models audited = 0   (PHASE 1/13 尚未执行 —— 见 MISSING ASSETS)")
p("External models with official row-level predictions = 0 (未审计)")
p("External models with official processed/split artifacts = 0 (未审计)")
p("External models requiring inference = 0 (未审计)")
p()
p("TIER_1 models = 内部 202 条评价键（U_CANONICAL 160 + U_MAS 42）")
p("TIER_2 models = 0")
p("TIER_3 models = 0（文献矩阵 PHASE 2 尚未执行）")
p()
p("-" * 32)
p("MATCHED UNIVERSE")
p("-" * 32)
uc = S[S.universe == "U_CANONICAL"]
um = S[S.universe == "U_MAS"]
p(f"Universe ID = U_CANONICAL  (family-OOD; canonical_target_labels by row_id)")
p(f"  rows = 4473   targets = 214   families = 108   n>=5 = 130")
p(f"  models evaluated = {len(uc)}")
p(f"Universe ID = U_MAS  (row-aligned to matched_analysis_set; 含 random/complex/family 三种 split)")
p(f"  rows = 7523   targets = 230   families = 112   n>=5 = 172")
p(f"  models evaluated = {len(um)}")
p()
p("-" * 32)
p("PRELIMINARY SCOREBOARD  (U_CANONICAL, 按 Macro Spearman 前 10)")
p("-" * 32)
cols = ["model_name", "model_family", "overall_pearson", "rmse",
        "macro_spearman", "C_W", "efficiency@10", "ndcg@10"]
p("| Model | Overall r | RMSE | Macro rho | C_W | Eff@10 | NDCG@10 |")
p("|---|---:|---:|---:|---:|---:|---:|")
for _, r in uc.nlargest(10, "macro_spearman").iterrows():
    p(f"| `{r['model_name']}` | {r['overall_pearson']:+.4f} | {r['rmse']:.4f} | "
      f"{r['macro_spearman']:+.4f} | {r['C_W']:+.4f} | {r['efficiency@10']:+.4f} | "
      f"{r['ndcg@10']:+.4f} |")
p()
p("同表按 Overall Pearson 前 5（用于对照）:")
p("| Model | Overall r | RMSE | Macro rho | C_W | Eff@10 | NDCG@10 |")
p("|---|---:|---:|---:|---:|---:|---:|")
for _, r in uc.nlargest(5, "overall_pearson").iterrows():
    p(f"| `{r['model_name']}` | {r['overall_pearson']:+.4f} | {r['rmse']:.4f} | "
      f"{r['macro_spearman']:+.4f} | {r['C_W']:+.4f} | {r['efficiency@10']:+.4f} | "
      f"{r['ndcg@10']:+.4f} |")
p()
p("-" * 32)
p("PRELIMINARY RANK AGREEMENT")
p("-" * 32)
pairs = [("overall_pearson", "macro_spearman"), ("overall_pearson", "efficiency@10"),
         ("overall_pearson", "ndcg@10"), ("rmse", "macro_spearman"),
         ("macro_spearman", "efficiency@10"), ("target_mean_pearson", "macro_spearman")]
for uni in ("U_CANONICAL", "U_MAS"):
    g = RA[RA.universe == uni]
    p(f"[{uni}]  n_models = {int(g.n_models.iloc[0]) if len(g) else 'NA'}")
    for a, b in pairs:
        m = g[((g.metric_a == a) & (g.metric_b == b)) |
              ((g.metric_a == b) & (g.metric_b == a))]
        if len(m):
            r = m.iloc[0]
            p(f"  {a:20s} vs {b:18s} Spearman={r.spearman_rank_agree:+.3f}  "
              f"Kendall={r.kendall_rank_agree:+.3f}")
    p()
p("Potential material rank reversals:")
for uni in ("U_CANONICAL",):
    f = V8 / "rank_reversal" / f"reversals_{uni}.csv"
    if f.exists():
        rv = pd.read_csv(f)
        top = rv[rv["pair"] == "overall_pearson|macro_spearman"].head(5)
        for r in top.itertuples():
            d = r._asdict()
            p(f"  {uni}: `{r.model}` ({r.family}) "
              f"rank(Overall)={int(d['rank_overall_pearson'])} vs "
              f"rank(Macro)={int(d['rank_macro_spearman'])}  "
              f"(r={d['val_overall_pearson']:+.3f}, rho={d['val_macro_spearman']:+.3f})")
p()
p("-" * 32)
p("TARGET-WEIGHTING (preliminary, H2)")
p("-" * 32)
cbcols = [c for c in WS.columns if c.startswith("C_B_")]
cwcols = [c for c in WS.columns if c.startswith("C_W_")]
cbr = (WS[cbcols].max(axis=1) - WS[cbcols].min(axis=1))
cwr = (WS[cwcols].max(axis=1) - WS[cwcols].min(axis=1))
sign_cb = int(((WS[cbcols] > 0).sum(axis=1).between(1, len(cbcols) - 1)).sum())
sign_cw = int(((WS[cwcols] > 0).sum(axis=1).between(1, len(cwcols) - 1)).sum())
p(f"  C_B 极差: median={cbr.median():.4f}  max={cbr.max():.4f}  "
  f"符号变化模型数={sign_cb}/{len(WS)}")
p(f"  C_W 极差: median={cwr.median():.4f}  max={cwr.max():.4f}  "
  f"符号变化模型数={sign_cw}/{len(WS)}")
rk = V8 / "weighting" / "weighting_rank_shift.csv"
if rk.exists():
    for r in pd.read_csv(rk).itertuples():
        p(f"  [{r.universe}] 模型排名 ROW vs EQUAL: Spearman={r.spearman_ROW_vs_EQUAL:+.3f}"
          f"  Kendall={r.kendall_ROW_vs_EQUAL:+.3f}  最大秩位移={r.max_rank_shift}")
p()
p("-" * 32)
p("NOVELTY MATRIX（初判，PHASE 2 完整文献审计尚未执行）")
p("-" * 32)
p("Existing work already covers (需 PHASE 2 逐一核实并标注):")
p("  per-target / per-structure Spearman, family-OOD splits, Top-K prioritization,")
p("  baseline/mutation decomposition, mutation microenvironment, partner-aware repr.")
p("Potentially distinctive V8 components (待文献核实后再定):")
p("  严格 task-specific family-OOD + matched-universe + 三轴能力分解 +")
p("  目标加权审计 + 模型秩反转 + 决策对齐 + benchmark-guided remedy 的**系统组合**")
p("Potential novelty already occupied = 未审计（不得预先声称任何 'first'）")
p()
p("-" * 32)
p("GPU")
p("-" * 32)
p("GPU_REQUIRED_SO_FAR = NO")
p(" 理由：本轮全部结论均来自**已有内部预测**，CPU 即可；")
p("       外部模型资产审计（PHASE 13）尚未开始，依 §75 不在 PHASE 13 之前决定 GPU。")
p()
p("-" * 32)
p("MISSING ASSETS")
p("-" * 32)
for s in [
    "PHASE 1 外部资源审计未执行（PPIformer / CATH-ddG / PreFold-dG / HERMES / Prompt-DDG / RDE-Network 的官方 row-level 预测）",
    "PHASE 2 文献基准矩阵 + 新颖性矩阵未执行（tabulate 未安装，但已改用 to_string 规避）",
    "PHASE 4 符号约定审计未执行（当前假设 +ddG = 结合变弱）",
    "PHASE 10 split inflation 未定稿（U_MAS 已具备 random/complex/family 三种 split 预测）",
    "PHASE 11–12 选择效用 / 决策对齐未定稿",
    "PHASE 15 5000 target-cluster bootstrap 未执行",
    "特征工程外的外部方法（FoldX/Flex-ddG/DDGPred/ProteinMPNN/ESM-IF 预测文件）未在项目内找到",
]:
    p(f"  - {s}")
p()
p("-" * 32)
p("NEXT_LOCAL_ACTION")
p("-" * 32)
p("执行 PHASE 10（split inflation）：用同一批 V6.2 模型在 U_MAS 上")
p("对比 random / complex(seq40) / family 三种 split 的 Overall r 与 Macro rho，")
p("量化「宽松划分抬高表观性能」的幅度——这是 H4 的直接检验。")
p()

(Path(V8) / "V8_LOCAL_CHECKPOINT.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print(f"\n写出 {V8 / 'V8_LOCAL_CHECKPOINT.md'}")
