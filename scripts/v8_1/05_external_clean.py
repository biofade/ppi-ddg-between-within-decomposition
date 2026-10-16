# -*- coding: utf-8 -*-
"""scripts/v8_1/05_external_clean.py -- V8.1 §13.3/§13.5：在 **EXACT_UNIQUE only** 上重算外部能力画像。

## 为什么必须重算

`04_external_audits.py` 证明此前的 CATH-ddG 映射（约 85% 匹配）**主要由歧义键构成**：
canonical 4473 行中只有 **1781 EXACT_UNIQUE**，另有 **2048 AMBIGUOUS_KEY**
（链与 WT 被丢弃后一个键对应多行）。**歧义键不得进入主匹配评估**（§13.3）。

⇒ 旧的外部读数（`U_EXT_CANONICAL` 上 Overall r +0.26 ~ +0.34 等）**作废**，本脚本重算。

## 输出

* `v8_1/external_clean_profiles.csv` —— EXACT_UNIQUE 上的三轴读数
* `v8_1/EXTERNAL_CLEAN.md`
* `FROZEN_RESULTS.json` 增加 `external_clean` 段
"""

from __future__ import annotations

import json
import re
import sys
from itertools import combinations
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

A = R / "v8_1"
V8 = R / "experiments" / "v8"
MIN_N = 5

lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)
ms = pd.read_csv(R / "data" / "v5_0" / "mutation_sites.csv", low_memory=False)
ms["pair_row"] = pd.to_numeric(ms.pair_row, errors="coerce")
z = np.load(R / "experiments" / "v7_1" / "phase2_rank" / "z_predictions.npz",
            allow_pickle=True)
Y_ALL = np.asarray(z["y"], float)

C = lab.merge(ms.sort_values("pair_row").groupby("pair_row", as_index=False).agg(
    pdb4=("pdb_id", lambda v: str(v.iloc[0]).upper().strip()[:4]),
    n_sites=("wt_aa", "size")), on="pair_row", how="left")
_pos = pd.to_numeric(ms.groupby("pair_row").ann_pos1.first()
                     .reindex(C.pair_row.to_numpy()).to_numpy(), errors="coerce")
_mut = ms.groupby("pair_row").mut_aa.first().reindex(
    C.pair_row.to_numpy()).to_numpy()
C["key_norm"] = (C.pdb4.astype(str) + "|"
                 + pd.Series(_pos).astype("float").astype("Int64").astype(str).to_numpy()
                 + "|" + pd.Series(_mut).astype(str).str.upper().to_numpy())
kc = C.key_norm.value_counts()

E = pd.read_csv("F:/dsh/downloads/cathddg/skempi_v2.csv", low_memory=False)


def parse_ext(txt):
    t = str(txt).strip().upper().replace(" ", "")
    m = re.match(r"^([A-Z]?)([A-Z])(\d+)([A-Z])$", t)
    return (None, None) if not m else (int(m.group(3)), m.group(4))


E[["pos", "mutaa"]] = pd.DataFrame([parse_ext(x) for x in E["Mutation(s)_cleaned"]],
                                   index=E.index)
E["pdb4"] = E.complex.astype(str).str.upper().str.strip().str[:4]
E["key_norm"] = (E.pdb4 + "|" + E.pos.astype("float").astype("Int64").astype(str)
                 + "|" + E.mutaa)

MODELS = ["Flex-ddG", "DDGPred", "ddG_FoldX", "ESM-IF", "proteinMPNN"]
for c in MODELS + ["ddG"]:
    if c in E.columns:
        E[c] = pd.to_numeric(E[c], errors="coerce")
gagg = {c: "mean" for c in MODELS + ["ddG"]}
EXT = E.groupby("key_norm", as_index=False).agg(gagg)
C2 = C.merge(EXT, on="key_norm", how="left", suffixes=("", "_ext"))
C2["y"] = Y_ALL
C2["status"] = np.where(~C2.key_norm.isin(set(EXT.key_norm)), "UNMATCHED",
                        np.where(C2.key_norm.map(kc) > 1, "AMBIGUOUS_KEY",
                                 "EXACT_UNIQUE"))
EU = C2[C2.status == "EXACT_UNIQUE"].copy()
print(f"canonical {len(C2)} 行 -> EXACT_UNIQUE {len(EU)}")
TGT = EU.target_id_canonical.astype(str).to_numpy()
FAM = EU.family_id_canonical.astype(str).to_numpy()
print(f"  EXACT_UNIQUE: targets={pd.Series(TGT).nunique()} "
      f"ge5={int((pd.Series(TGT).value_counts()>=MIN_N).sum())} "
      f"families={pd.Series(FAM).nunique()}")

rows = []
print("\n" + "=" * 86)
print("EXACT_UNIQUE 上的外部方法三轴读数（§13.3 主匹配评估）")
print("=" * 86)
print(f"{'model':<13}{'n':>6}{'tgt>=5':>8}{'Overall r':>11}{'Macro rho':>11}"
      f"{'C_B':>9}{'C_W':>9}{'Eff@10':>9}{'RMSE':>9}")
for c in MODELS:
    if c not in EU.columns:
        continue
    m = EU[c].notna().to_numpy() & np.isfinite(EU.y.to_numpy(float))
    if m.sum() < 50:
        print(f"{c:<13}{int(m.sum()):>6}  —— 行数不足")
        continue
    U = M.prep_universe(EU.y.to_numpy(float)[m], EU[c].to_numpy(float)[m],
                        TGT[m], FAM[m])
    r = M.compute_all(U)
    r["model"] = c
    r["rmse"] = r.get("rmse", np.nan)
    rows.append(r)
    print(f"{c:<13}{r['n_rows']:>6}{r['n_targets_ge5']:>8}"
          f"{r['overall_pearson']:>+11.4f}{r['macro_spearman']:>+11.4f}"
          f"{r['C_B_row']:>+9.4f}{r['C_W']:>+9.4f}"
          f"{r['efficiency@10']:>+9.4f}{r['rmse']:>9.4f}")

D = pd.DataFrame(rows)
D.to_csv(A / "external_clean_profiles.csv", index=False, encoding="utf-8")

# 外部方法族秩一致（小 n，含置换 p）
ext = D.dropna(subset=["overall_pearson", "macro_spearman", "efficiency@10"])
ext_out = {"n_methods": int(len(ext))}
print(f"\n  可评价外部方法 n = {len(ext)}")
if len(ext) >= 4:
    for a, b in (("overall_pearson", "macro_spearman"),
                 ("overall_pearson", "efficiency@10"),
                 ("macro_spearman", "efficiency@10")):
        rho = float(spearmanr(ext[a], ext[b]).correlation)
        tau = float(kendalltau(ext[a], ext[b]).correlation)
        xa = ext[a].to_numpy(float); xb = ext[b].to_numpy(float)
        try:
            from itertools import permutations as _perm
            obs = abs(spearmanr(xa, xb).correlation)
            cnt = tot = 0
            for pm in _perm(xb):
                tot += 1
                if abs(spearmanr(xa, np.asarray(pm)).correlation) >= obs - 1e-12:
                    cnt += 1
            p = cnt / tot
        except Exception:
            p = np.nan
        print(f"  rho({a}, {b}) = {rho:+.4f}  tau={tau:+.4f}  perm_p={p:.4f}")
        ext_out[f"rho_{a}_vs_{b}"] = rho
        ext_out[f"tau_{a}_vs_{b}"] = tau
        ext_out[f"perm_p_{a}_vs_{b}"] = p
    ext_out["low_power_warning"] = (
        f"n_methods={len(ext)} < 8 ⇒ rank correlation is UNDERPOWERED; "
        "report method-level profiles rather than a correlation headline")

# 逐方法画像（§13.5 首选呈现方式）
prof = []
for r in D.itertuples():
    prof.append({"model": r.model, "n_rows": r.n_rows,
                 "n_targets_ge5": r.n_targets_ge5,
                 "overall_pearson": r.overall_pearson, "rmse": r.rmse,
                 "macro_spearman": r.macro_spearman, "C_B": r.C_B_row,
                 "C_W": r.C_W, "efficiency@10": r.efficiency_at_10
                 if hasattr(r, "efficiency_at_10") else r._asdict().get("efficiency@10")})
P = pd.DataFrame(prof)
if len(P):
    ow = P.loc[P.overall_pearson.idxmax()]
    rw = P.loc[P.macro_spearman.idxmax()]
    print(f"\n  Overall 赢家 = {ow.model} ({ow.overall_pearson:+.4f})")
    print(f"  靶内排序赢家 = {rw.model} ({rw.macro_spearman:+.4f})")
    print(f"  同一方法? {ow.model == rw.model}")
    ext_out["overall_winner"] = str(ow.model)
    ext_out["ranking_winner"] = str(rw.model)
    ext_out["same_winner"] = bool(ow.model == rw.model)
    P.to_csv(A / "external_clean_profiles_flat.csv", index=False, encoding="utf-8")

# 与旧（含歧义键）读数对照
old = pd.read_csv(V8 / "external" / "external_multimodel_scoreboard.csv") \
    if (V8 / "external" / "external_multimodel_scoreboard.csv").exists() else None
if old is not None and len(D):
    cmp = D[["model", "overall_pearson", "macro_spearman"]].merge(
        old[["model", "overall_pearson", "macro_spearman"]],
        on="model", suffixes=("_cleanEU", "_inclAmbiguous"))
    cmp["delta_overall"] = cmp.overall_pearson_cleanEU - cmp.overall_pearson_inclAmbiguous
    cmp["delta_macro"] = cmp.macro_spearman_cleanEU - cmp.macro_spearman_inclAmbiguous
    cmp.to_csv(A / "external_clean_vs_old.csv", index=False, encoding="utf-8")
    print("\n  与旧读数对照（含歧义键 vs 仅 EXACT_UNIQUE）:")
    print(cmp.round(4).to_string(index=False))
    ext_out["comparison_to_ambiguous_mapping"] = cmp.to_dict(orient="records")

FR = json.load(open(A / "FROZEN_RESULTS.json", encoding="utf-8"))
FR["external_clean"] = {
    "universe": "U_EXT_EXACT_UNIQUE",
    "n_rows": int(len(EU)), "n_targets": int(pd.Series(TGT).nunique()),
    "n_ge5": int((pd.Series(TGT).value_counts() >= MIN_N).sum()),
    "n_families": int(pd.Series(FAM).nunique()),
    "excluded_ambiguous_rows": int((C2.status == "AMBIGUOUS_KEY").sum()),
    "excluded_unmatched_rows": int((C2.status == "UNMATCHED").sum()),
    "profiles": D.to_dict(orient="records"),
    "rank_agreement_external": ext_out,
}
FR.setdefault("forbidden_claims", []).append(
    "any external capability claim computed on chain-ambiguous mappings")
json.dump(FR, open(A / "FROZEN_RESULTS.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)

md = ["# 外部证据（EXACT_UNIQUE only，§13.3/§13.5）\n",
      f"**宇宙**：`U_EXT_EXACT_UNIQUE` — {len(EU)} 行 / "
      f"{pd.Series(TGT).nunique()} targets / "
      f"{(pd.Series(TGT).value_counts()>=MIN_N).sum()} 个 n≥5 / "
      f"{pd.Series(FAM).nunique()} families\n",
      f"**排除**：`AMBIGUOUS_KEY` {int((C2.status=='AMBIGUOUS_KEY').sum())} 行、"
      f"`UNMATCHED` {int((C2.status=='UNMATCHED').sum())} 行\n",
      "> 此前基于完整键匹配（约 85%）的读数**主要由歧义键构成**，已作废。\n",
      "## 方法级能力画像（首选呈现）\n",
      "| model | n | n≥5 | Overall r | Macro ρ | C_B | C_W | Eff@10 | RMSE |",
      "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
for r in D.itertuples():
    d = r._asdict()
    md.append(f"| `{r.model}` | {r.n_rows} | {r.n_targets_ge5} | "
              f"{r.overall_pearson:+.4f} | {r.macro_spearman:+.4f} | "
              f"{r.C_B_row:+.4f} | {r.C_W:+.4f} | "
              f"{d.get('efficiency@10', float('nan')):+.4f} | {r.rmse:.4f} |")
md += ["", "## 外部秩一致（**低功效**）\n", "```", json.dumps(
    ext_out, indent=2, ensure_ascii=False), "```",
    "\n> n 过小 ⇒ 秩相关**不宜作为头条**；以上方法级画像为准（§13.5）。\n"]
(A / "EXTERNAL_CLEAN.md").write_text("\n".join(md) + "\n", encoding="utf-8")
print(f"\n写出 external_clean_profiles.csv / EXTERNAL_CLEAN.md / "
      f"external_clean_profiles_flat.csv（+ FROZEN_RESULTS.external_clean）")
