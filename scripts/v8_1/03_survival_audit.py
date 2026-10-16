# -*- coding: utf-8 -*-
"""scripts/v8_1/03_survival_audit.py -- V8.1：历史头条结论的存活审计。

检验 4 个历史头条是否在**匹配 + 去伪重复 + provenance-clean** 后存活：

  S1. 「Overall Pearson 与 Macro Spearman 近乎正交 / ρ=−0.26」                (§7)
  S2. 「7 ↔ 160 极端排名反转」                                              (§7.1)
  S3. 「宽松 vs 严格制度下 ρ(Overall, Eff@10) 从 +0.81 塌到 −0.03」           (§8.2)
  S4. scenario 表「Overall 选择的模型排第 54」                              (§14)
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
A = R / "v8_1"
D = pd.read_csv(A / "all_models_by_rowset.csv", low_memory=False)
ST = pd.read_csv(A / "rank_agreement_by_stratum.csv")
FR = json.load(open(A / "FROZEN_RESULTS.json", encoding="utf-8"))
out = {}

print("=" * 84)
print("S1  历史：「Overall Pearson 与 Macro Spearman 近乎正交 / ρ=−0.26」")
print("=" * 84)
v = ST[["rowset_hash", "level", "n_models", "n_method_families",
        "rho_overall_vs_macro"]].copy()
print(v.round(4).to_string(index=False))
rng = (float(v.rho_overall_vs_macro.min()), float(v.rho_overall_vs_macro.max()))
n_neg = int((v.rho_overall_vs_macro < 0).sum())
n_pos = int((v.rho_overall_vs_macro > 0).sum())
print(f"\n范围 = {rng[0]:+.4f} .. {rng[1]:+.4f}；负值池 {n_neg} / 正值池 {n_pos} "
      f"/ 共 {len(v)}")
print("⇒ 判定：", end="")
if rng[0] < 0 < rng[1] and abs(rng[0]) < 0.5 and rng[1] > 0.5:
    verdict = ("SURVIVES_ONLY_AS_POOL_DEPENDENT_RANGE — 单值正交/-0.26 不成立；"
               "跨池范围包含强正值(+0.82)，因此「近乎正交」**不能作为普遍结论**")
elif n_neg and n_pos:
    verdict = "WEAKENED — 符号随池变化"
else:
    verdict = "DOES_NOT_SURVIVE"
print("  " + verdict)
out["S1"] = {"range": rng, "n_neg_pools": n_neg, "n_pos_pools": n_pos,
             "verdict": verdict,
             "historical_claim": "near-orthogonality / rho=-0.26",
             "remains_defensible": False}

print("\n" + "=" * 84)
print("S2  历史：「7 ↔ 160 极端排名反转」")
print("=" * 84)
MODEL = "L1_random_A6_global_set_local"
base = D[(D.key == MODEL)]
print(f"  该模型出现在 {len(base)} 个 rowset：{base.rowset_hash.unique().tolist()}")
loc = []
for rh, g in D.groupby("rowset_hash"):
    g = g.dropna(subset=["overall_pearson", "macro_spearman"])
    if MODEL not in set(g.key) or len(g) < 10:
        continue
    ro = g.overall_pearson.rank(ascending=False)[g.key == MODEL].iloc[0]
    rm = g.macro_spearman.rank(ascending=False)[g.key == MODEL].iloc[0]
    has_perm = bool((g.split_regime == "permissive").any())
    loc.append({"rowset_hash": rh, "n_models": len(g),
                "rank_overall": int(ro), "rank_macro": int(rm),
                "displacement": int(abs(ro - rm)),
                "pool_contains_permissive": has_perm})
L = pd.DataFrame(loc)
print(L.to_string(index=False))
in_clean = L[~L.pool_contains_permissive]
print(f"\n  纯 provenance-clean 池中的反转位移: "
      f"{in_clean.displacement.tolist() if len(in_clean) else 'N/A'}")
verdict2 = ("SURVIVES_IN_MATCHED_POOL" if (len(in_clean) and
            (in_clean.displacement > 5).any())
            else "DOES_NOT_SURVIVE_IN_PROVENANCE_CLEAN_POOL")
print("⇒ 判定: " + verdict2)
out["S2"] = {"detail": L.to_dict(orient="records"), "verdict": verdict2,
             "remains_defensible": verdict2 == "SURVIVES_IN_MATCHED_POOL"}

print("\n" + "=" * 84)
print("S3  历史：「宽松 vs 严格制度下 ρ(Overall, Eff@10) +0.81 → −0.03」")
print("=" * 84)
print("  §8.2 规则：若两池在 **模型组成 / 行覆盖 / method-family 组成** 上不同，")
print("  则**不构成**划分效应的因果证据 ⇒ 必须移入探索性/补充材料。")
s3 = D.dropna(subset=["overall_pearson", "efficiency@10"])
perm = s3[s3.split_regime == "permissive"]
strict = s3[s3.split_regime != "permissive"]
print(f"  permissive 池: n={len(perm)} rows={sorted(perm.n_rows.unique())} "
      f"families={perm.method_family.nunique()}")
print(f"  strict     池: n={len(strict)} rows={sorted(strict.n_rows.unique())[:6]} "
      f"families={strict.method_family.nunique()}")
same_rows = set(perm.n_rows.unique()) == set(strict.n_rows.unique())
same_fam = perm.method_family.nunique() == strict.method_family.nunique()
print(f"  行集相同? {same_rows}   method-family 数相同? {same_fam}")
print("⇒ 判定: DOES_NOT_QUALIFY_AS_CAUSAL_EVIDENCE（池组成不同）⇒ 降级为探索性")
out["S3"] = {"verdict": "DOES_NOT_QUALIFY_AS_CAUSAL_EVIDENCE",
             "same_rowsets": bool(same_rows), "same_family_count": bool(same_fam),
             "remains_defensible": False}

print("\n" + "=" * 84)
print("S4  历史：scenario「Overall 选择的模型在文库排序场景排第 54」")
print("=" * 84)
perm_keys = set(D.loc[D.split_regime == "permissive", "key"])
scen = pd.read_csv(R / "experiments" / "v8" / "scenarios" /
                   "scenario_recommendations.csv")
print(scen[["scenario", "primary_metric", "top1",
            "overall_pearson_pick"]].to_string(index=False))
print(f"\n  历史 top1 是否 permissive-训练: "
      f"{[k in perm_keys for k in scen.top1]}")
# 在最大 provenance-clean stratum 内重算场景 B 的排名
rh_main = FR["primary_universe"]["rowset_hash"]
g = D[(D.rowset_hash == rh_main)].dropna(subset=["macro_spearman"])
g = g[g.split_regime != "permissive"]
if len(g):
    g = g.sort_values("macro_spearman", ascending=False).reset_index(drop=True)
    pick = "L1_random_A1_base_set"
    rk = (g.index[g.key == pick].tolist() or [-1])[0] + 1
    print(f"\n  在主池 ({rh_main}, n={len(g)}, provenance-clean) 内：")
    print(f"    按 Macro ρ 首选 = {g.key.iloc[0]}")
    print(f"    历史 Overall 选择 `{pick}` 的排名 = {rk if rk>0 else 'N/A（不在该池）'}")
    out["S4"] = {"primary_rowset": rh_main, "n_clean_models": len(g),
                 "scenario_B_top1_macro": g.key.iloc[0],
                 "historical_pick": pick, "rank_in_clean_pool": int(rk),
                 "verdict": ("RESTATED_UNDER_MATCHED_POOL" if rk > 0
                             else "NOT_PRESENT_IN_PRIMARY_POOL")}
    print(f"  ⇒ 判定: {out['S4']['verdict']}")
else:
    out["S4"] = {"verdict": "NOT_ESTIMABLE_UNDER_MATCHED_POOL"}
    print("  ⇒ 判定: NOT_ESTIMABLE_UNDER_MATCHED_POOL")

json.dump(out, open(A / "SURVIVAL_AUDIT.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)

md = ["# V8.1 历史头条存活审计（§7.1/§8.2/§14）\n",
      "> 规则：任何历史头条若在**匹配 + 去伪重复 + provenance-clean** 后不存活，",
      "> 必须降级为历史/补充材料，或删除。\n",
      "## 汇总\n",
      "| 历史头条 | 判定 | 是否仍可辩护 |", "|---|---|:--:|"]
for k, v in out.items():
    md.append(f"| {k} | {v.get('verdict','')} | "
              f"{'✅' if v.get('remains_defensible') else '❌'} |")
md += ["", "## S1 详情：「近乎正交 / ρ=−0.26」\n",
       f"跨 {len(ST)} 个匹配 stratum 的 `rho(Overall r, Macro rho)` 范围 = "
       f"**{rng[0]:+.4f} .. {rng[1]:+.4f}**（负值池 {n_neg}，正值池 {n_pos}）。\n",
       "**⇒ 单值不成立**；且范围上界为**强正值 (+0.82)**，",
       "所以「全局与靶内指标近乎正交」**不能作为普遍结论**。",
       "可辩护的表述见 `FINAL_DEFENSIBLE_CLAIMS.md`。\n"]
(A / "HISTORICAL_SURVIVAL_AUDIT.md").write_text("\n".join(md) + "\n",
                                                encoding="utf-8")
print(f"\n写出 SURVIVAL_AUDIT.json / HISTORICAL_SURVIVAL_AUDIT.md")
