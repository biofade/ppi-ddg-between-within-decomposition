# -*- coding: utf-8 -*-
"""scripts/v8/32b_hermes_match.py -- HERMES × canonical 精确匹配（用 mutation_sites 构造 key）。

canonical 只有 (pair_row, target, family)；突变身份从
`data/v5_0/mutation_sites.csv` 取（key = `pair_row`）。

构造 HERMES 风格的 key：`<aa><ann_pos1><mut_aa>`（如 `L38G`），
其含义为 **mutation_sites.wt_aa 应为 HERMES 的 wt**。
→ 用两种归一化做匹配并报告命中率：
   (a) `wt_aa + ann_pos1 + mut_aa`   —— 假定编号一致
   (b) 只用 `mut_aa` + 按 mutation_sites 顺序对到 HERMES 同 PDB 的行
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
DL = Path("F:/dsh/downloads/hermes")
OUT = R / "experiments" / "v8" / "external"

lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)
ms = pd.read_csv(R / "data" / "v5_0" / "mutation_sites.csv", low_memory=False)
print(f"canonical={len(lab)}  mutation_sites={len(ms)}")

ms["pair_row"] = pd.to_numeric(ms.pair_row, errors="coerce")
lab = lab.merge(ms[["pair_row", "pdb_id", "ann_chain", "wt_aa", "mut_aa",
                    "ann_pos1", "seq_index0"]],
                left_on="pair_row", right_on="pair_row", how="left")
print(f"合并后: {len(lab)}, mutation 列缺失 = {lab.wt_aa.isna().sum()}")

def s(x):
    return str(x).strip().upper()
lab["pdb4"] = lab.pdb_id.astype(str).str.upper().str.strip().str[:4]
lab["mut_h"] = (lab.wt_aa.map(s) + lab.ann_pos1.astype("Int64").astype(str)
                + lab.mut_aa.map(s))
lab["key_h"] = lab.pdb4 + "|" + lab.mut_h
print(f"canonical key 样例: {lab.key_h.head(4).tolist()}")

# HERMES
hs = {}
for tag in ("hermes_py000_cdna117k", "hermes_py050_cdna117k"):
    p = DL / f"{tag}.csv"
    if p.exists():
        hs[tag] = pd.read_csv(p, low_memory=False)
h = list(hs.values())[0].copy()
LP = "log_proba_mt__minus__log_proba_wt"
h["pdb4"] = h.PDB_filename.astype(str).str.upper().str.strip().str[:4]
h["mut_h"] = h["mutant"].astype(str).str.strip().str.upper()
h["key_h"] = h.pdb4 + "|" + h.mut_h
h["proxy"] = -pd.to_numeric(h[LP], errors="coerce")   # 方向对齐 V8
print(f"HERMES key 样例: {h.key_h.head(4).tolist()}")

hagg = h.groupby("key_h", as_index=False).agg(proxy=("proxy", "mean"),
                                              n=("proxy", "size"))
M = lab.merge(hagg, on="key_h", how="left")
hit = M.proxy.notna()
print("\n" + "=" * 70)
print("精确匹配（wt_aa + ann_pos1 + mut_aa）")
print("=" * 70)
print(f"  canonical 行 = {len(M)}")
print(f"  命中 = {int(hit.sum())} ({hit.mean():.4f})")
print(f"  覆盖 target = {M.loc[hit,'target_id_canonical'].nunique()} / "
      f"{M.target_id_canonical.nunique()}")

# 回退：若命中率低，试 seq_index0 无 chain 只用 (pdb4, mut_aa, pos)
if hit.mean() < 0.5:
    print("\n  命中率偏低 ⇒ 试 HERMES 的 `Mutation(s)_cleaned` 前缀链格式")
    h2 = list(hs.values())[0].copy()
    h2["pdb4"] = h2.PDB_filename.astype(str).str.upper().str.strip().str[:4]
    h2["clean"] = h2["Mutation(s)_cleaned"].astype(str).str.strip().str.upper()
    h2["key_c"] = h2.pdb4 + "|" + h2.clean
    h2["proxy"] = -pd.to_numeric(h2[LP], errors="coerce")
    # canonical: 用 ann_chain + ann_pos1 + mut_aa  (HERMES cleaned = <chain><pos><mut>)
    lab["clean_c"] = (lab.ann_chain.map(s)
                      + lab.ann_pos1.astype("Int64").astype(str)
                      + lab.mut_aa.map(s))
    lab["key_c"] = lab.pdb4 + "|" + lab.clean_c
    h2a = h2.groupby("key_c", as_index=False).agg(proxy=("proxy", "mean"))
    M2 = lab.merge(h2a, on="key_c", how="left")
    hit2 = M2.proxy.notna()
    print(f"  变体(b) chain+pos+mut: 命中 = {int(hit2.sum())} ({hit2.mean():.4f})")
    if hit2.mean() > hit.mean():
        M, hit = M2, hit2
        print("  ⇒ 采用变体(b)")

# 相关（方向校验）
yc = "ddg" if "ddg" in M.columns else None
print(f"\ncanonical 真值列 = {yc!r}")
if yc and hit.sum() > 10:
    yy = pd.to_numeric(M.loc[hit, yc], errors="coerce")
    pp = M.loc[hit, "proxy"]
    mm = yy.notna() & pp.notna()
    r = np.corrcoef(yy[mm], pp[mm])[0, 1]
    print(f"  r(canonical ddG, HERMES proxy) = {r:+.4f}  (n={int(mm.sum())})")
    print("  ⇒ 期望为正（proxy 已取负对齐 ΔΔG>0=变弱）")

M.to_csv(OUT / "hermes_matched.csv", index=False, encoding="utf-8")
meta = {
    "match_rule": "canonical(mutation_sites) x HERMES on (PDB4, wt+pos+mut or chain+pos+mut)",
    "canonical_rows": int(len(M)),
    "matched_rows": int(hit.sum()),
    "match_rate": float(hit.mean()),
    "matched_targets": int(M.loc[hit, "target_id_canonical"].nunique()),
    "total_targets": int(M.target_id_canonical.nunique()),
    "sign_convention": "hermes proxy = -(log_proba_mt - log_proba_wt); V8 ddG>0=weaker",
    "evidence_tier": "TIER 2 (SKEMPI v2.0 finetuned; §50 requires family-overlap audit)",
}
(OUT / "hermes_matched_meta.json").write_text(
    json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"\n写出 hermes_matched.csv ({len(M)} 行) / hermes_matched_meta.json")
