# -*- coding: utf-8 -*-
"""scripts/v8/44_pool_provenance_audit.py -- 模型池 provenance 审计（审稿人可复核）。

## 为什么必须有这份文件

`rho(Overall r, Macro rho)` 在本次工作中被报出过 **+0.18 / +0.19 / −0.26 / +0.28 / +0.41 / +0.47 / +0.63 / +0.68**
等多个值。审稿人无法判断哪个对。根因是**模型池组成不同**，且此前未把池定义固定下来。

**本脚本把每个池的定义、成员数、provenance 构成、来源脚本全部列出**，
使每个数字都可被独立复核。

## 纪律

* 池成员**必须**由脚本明确写出（不靠记忆）；
* 每个池报告 `n_models` / `permissive` / `strict` / `other` / 来源脚本 / 点估计；
* 同时给出 `rho(Overall r, Macro rho)` 与 `rho(Overall r, Eff@10)`。
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

PERM_TOK = ("L1_random_", "L2_complex_", "L3_seq40_")
STRICT_TOK = ("L4_family_", "F0_", "F1_", "F2_", "F3_", "F4_", "F5_", "F6_",
              "F7_", "F8_", "F9_", "R0_", "R1_", "R2_", "R3_", "B0_", "B1_",
              "B2_", "B3", "B4_", "A2_smallPLM", "p_", "z_R", "z_A")


def prov(k: str) -> str:
    if any(k.startswith(t) for t in PERM_TOK):
        return "PERMISSIVE"
    if any(k.startswith(t) for t in STRICT_TOK):
        return "STRICT"
    return "OTHER"


S = pd.read_csv(OUT / "BENCHMARK_SCOREBOARD.csv", low_memory=False)
POOL = pd.read_csv(EXT / "pool_composition_models.csv", low_memory=False)


def stats(g, pool_id, source, definition):
    g = g.dropna(subset=["overall_pearson", "macro_spearman", "efficiency@10"])
    if len(g) < 5:
        return None
    pv = g.model_name.map(prov) if "provenance" not in g.columns else g.provenance
    rec = {
        "pool_id": pool_id, "source_script": source, "definition": definition,
        "n_models": int(len(g)),
        "n_PERMISSIVE": int((pv == "PERMISSIVE").sum()),
        "n_STRICT": int((pv == "STRICT").sum()),
        "n_OTHER": int((pv == "OTHER").sum()),
        "rho_overall_macro": float(spearmanr(g.overall_pearson, g.macro_spearman).correlation),
        "tau_overall_macro": float(kendalltau(g.overall_pearson, g.macro_spearman).correlation),
        "rho_overall_eff10": float(spearmanr(g.overall_pearson, g["efficiency@10"]).correlation),
        "rho_macro_eff10": float(spearmanr(g.macro_spearman, g["efficiency@10"]).correlation),
        "med_overall_r": float(g.overall_pearson.median()),
        "med_macro": float(g.macro_spearman.median()),
    }
    return rec


rows = []

u = S[S.universe == "U_CANONICAL"]
rows.append(stats(u, "P1_U_CANONICAL_all",
                  "03_model_registry / 06_capability_metrics",
                  "all U_CANONICAL scoreboard rows"))
u67 = u[u.n_rows == 4473]
rows.append(stats(u67, "P2_U_CANONICAL_nrows4473",
                  "08_rank_reversal",
                  "U_CANONICAL AND n_rows==4473 (V6.2 + V7.1 families only)"))
u67c = u67[~u67.model_family.eq("V6_2_split4")]
rows.append(stats(u67c, "P3_P2_minus_V6_2",
                  "decomposition of P2",
                  "P2 with V6_2_split4 removed"))
v7 = u[~u.model_family.eq("V6_2_split4")]
rows.append(stats(v7, "P4_U_CANONICAL_minus_V6_2",
                  "43_fix_fig1_and_scenarios (fig1 panel b)",
                  "U_CANONICAL with V6_2_split4 removed"))
rows.append(stats(POOL, "P5_pool_composition_ALL",
                  "35_pool_composition_artifact",
                  "all scoreboards incl. mas-align artifacts, 10 artifacts"))
rows.append(stats(POOL[POOL.provenance == "PERMISSIVE"],
                  "P5a_pool_composition_PERMISSIVE",
                  "35_pool_composition_artifact", "PERMISSIVE subset"))
rows.append(stats(POOL[POOL.provenance == "STRICT"],
                  "P5b_pool_composition_STRICT",
                  "35_pool_composition_artifact", "STRICT subset"))

for f, pid, defn in (
        ("external_expanded_pool.csv", "P6_HERMES_variants",
         "7 HERMES/proteinMPNN variants (2 families)"),
        ("external_multimodel_scoreboard.csv", "P7_EXT_CANONICAL",
         "5 CATH-ddG baselines on canonical rows"),
        ("external_heldout_cath.csv", "P8_EXT_HELDOUT_CATH",
         "5 CATH-ddG baselines on official held-out CATH")):
    p = EXT / f
    if not p.exists():
        continue
    d = pd.read_csv(p)
    d["model_name"] = d.get("model", pd.Series(dtype=str))
    if "model_name" not in d or d.model_name.isna().all():
        continue
    r = stats(d, pid, "39/40/41 external scripts", defn)
    if r:
        rows.append(r)

D = pd.DataFrame(rows)
D.to_csv(EXT / "pool_provenance_audit.csv", index=False, encoding="utf-8")

print("=" * 108)
print("模型池 provenance 审计 —— 每个 rho 数字的可复核来源")
print("=" * 108)
print(f"{'pool_id':<32}{'n':>5}{'perm':>6}{'strict':>8}{'other':>7}"
      f"{'rho(r,macro)':>14}{'rho(r,eff)':>12}{'rho(macro,eff)':>15}")
for r in D.itertuples():
    print(f"{r.pool_id:<32}{r.n_models:>5}{r.n_PERMISSIVE:>6}{r.n_STRICT:>8}"
          f"{r.n_OTHER:>7}{r.rho_overall_macro:>+14.4f}"
          f"{r.rho_overall_eff10:>+12.4f}{r.rho_macro_eff10:>+15.4f}")
print()
print("  可复核范围:", f"{D.rho_overall_macro.min():+.3f} .. {D.rho_overall_macro.max():+.3f}")
print("  ⇒ 任何单一数字都必须在附带 pool_id 的情况下引用。")

md = ["# V8 模型池 Provenance 审计（审稿人可复核）\n",
      "> `rho(Overall r, Macro rho)` 在本次工作中被报出过多个不同值。",
      "> 根因是**模型池组成不同**。本表把每个池的定义、成员构成与来源脚本全部固定下来。\n",
      "| pool_id | 定义 | 来源脚本 | n | PERMISSIVE | STRICT | OTHER | **rho(r,macro)** | rho(r,eff) | rho(macro,eff) |",
      "|---|---|---|---:|---:|---:|---:|---:|---:|---:|"]
for r in D.itertuples():
    md.append(f"| `{r.pool_id}` | {r.definition} | `{r.source_script}` | {r.n_models} | "
              f"{r.n_PERMISSIVE} | {r.n_STRICT} | {r.n_OTHER} | "
              f"**{r.rho_overall_macro:+.4f}** | {r.rho_overall_eff10:+.4f} | "
              f"{r.rho_macro_eff10:+.4f} |")
md += ["", f"**观测范围：{D.rho_overall_macro.min():+.4f} .. "
           f"{D.rho_overall_macro.max():+.4f}**", "",
       "## 使用规则（写入论文方法部分）", "",
       "1. 引用任何秩相关时**必须**同时给出 `pool_id` 与其 provenance 构成；",
       "2. **不得**把任一单值表述为 benchmark 的性质；",
       "3. 报告时给出**池间范围**而非点估计；",
       "4. 最小池（n<20）的秩相关**不应**作为结论依据（有效样本量受伪重复削弱）。"]
(OUT / "POOL_PROVENANCE_AUDIT.md").write_text("\n".join(md) + "\n", encoding="utf-8")
meta = {"pools": rows,
        "observed_range_rho_overall_macro":
            [float(D.rho_overall_macro.min()), float(D.rho_overall_macro.max())]}
(EXT / "pool_provenance_audit.json").write_text(
    json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"\n写出 POOL_PROVENANCE_AUDIT.md / external/pool_provenance_audit.csv")
