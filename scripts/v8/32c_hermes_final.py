# -*- coding: utf-8 -*-
"""scripts/v8/32c_hermes_final.py -- HERMES × canonical **严格 1:1** 匹配（修正 32b 的行膨胀）。

32b 的 `pair_row` 合并使 4473 行膨胀到 8705（mutation_sites 对多突变体有多行）。
本脚本按 `pair_row` 聚合后**强制回到 4473 行**，再报告匹配率。

同时给出**每条 canonical 行的完整突变身份**用于复核。
"""

from __future__ import annotations

import json
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
ms["pair_row"] = pd.to_numeric(ms.pair_row, errors="coerce")
print(f"canonical={len(lab)}  mutation_sites={len(ms)}  "
      f"pair_row 唯一={ms.pair_row.nunique()}")

# 按 pair_row 聚合突变身份（多突变体 -> 逗号连接）
def s(x):
    return str(x).strip().upper()

ms["tok"] = ms.wt_aa.map(s) + ms.ann_pos1.astype("Int64").astype(str) + ms.mut_aa.map(s)
agg = (ms.sort_values("pair_row")
         .groupby("pair_row", as_index=False)
         .agg(pdb4=("pdb_id", lambda v: str(v.iloc[0]).upper().strip()[:4]),
              n_sites=("tok", "size"),
              mut_str=("tok", lambda v: ",".join(sorted(set(v))))))
print(f"聚合后 pair_row 数 = {len(agg)}")

L = lab.merge(agg, left_on="pair_row", right_on="pair_row", how="left")
assert len(L) == len(lab), f"行数膨胀! {len(L)} != {len(lab)}"
print(f"✅ 严格 1:1: {len(L)} 行")
print(f"   突变身份缺失 = {int(L.mut_str.isna().sum())}")
print(f"   单点突变行 = {int((L.n_sites == 1).sum())}, "
      f"多点 = {int((L.n_sites > 1).sum())}")
L["key_h"] = L.pdb4 + "|" + L.mut_str

h = pd.read_csv(DL / "hermes_py000_cdna117k.csv", low_memory=False)
LP = "log_proba_mt__minus__log_proba_wt"
h["pdb4"] = h.PDB_filename.astype(str).str.upper().str.strip().str[:4]
h["key_h"] = h.pdb4 + "|" + h["mutant"].astype(str).str.strip().str.upper()
h["proxy"] = -pd.to_numeric(h[LP], errors="coerce")
hagg = h.groupby("key_h", as_index=False).agg(proxy=("proxy", "mean"),
                                              n_hermes=("proxy", "size"))
M = L.merge(hagg, on="key_h", how="left")
hit = M.proxy.notna()

print("\n" + "=" * 72)
print("HERMES × canonical 严格匹配结果")
print("=" * 72)
print(f"  canonical 行           = {len(M)}")
print(f"  命中 HERMES            = {int(hit.sum())}  ({hit.mean():.4f})")
print(f"  覆盖 target            = {M.loc[hit,'target_id_canonical'].nunique()}"
      f" / {M.target_id_canonical.nunique()}")
print(f"  覆盖 family            = {M.loc[hit,'family_id_canonical'].nunique()}"
      f" / {M.family_id_canonical.nunique()}")
single = M.n_sites == 1
print(f"\n  仅单点突变子集: n={int(single.sum())}, "
      f"命中={int(M.loc[single,'proxy'].notna().sum())} "
      f"({M.loc[single,'proxy'].notna().mean():.4f})")
n5 = M[hit].target_id_canonical.value_counts()
print(f"  命中子集中 n>=5 的 target = {int((n5>=5).sum())}")

# 方向校验：用 y（来自 row-aligned z_predictions）
z = np.load(R / "experiments" / "v7_1" / "phase2_rank" / "z_predictions.npz",
            allow_pickle=True)
Y = np.asarray(z["y"], float)
assert len(Y) == len(lab)
M = M.sort_values("pair_row").reset_index(drop=True)
M["y_true"] = Y
mm = hit.to_numpy() & np.isfinite(Y)
if mm.sum() > 10:
    r = np.corrcoef(Y[mm], M.proxy.to_numpy()[mm])[0, 1]
    print(f"\n  方向校验: r(y_true ΔΔG, HERMES proxy) = {r:+.4f} (n={int(mm.sum())})")
    print(f"  ⇒ 期望**为正**（proxy 已取负以对齐 ΔΔG>0=结合变弱）")

M.to_csv(OUT / "hermes_matched_final.csv", index=False, encoding="utf-8")
meta = {
    "hermes_source": ("https://github.com/StatPhysBio/hermes_experiments/"
                      "skempi/results/hermes_py_000_ft_cdna117k_ddg_st/"),
    "hermes_raw_column": LP,
    "sign_convention": ("hermes proxy = -(log_proba_mt - log_proba_wt); "
                        "V8: ddG>0 = weaker binding"),
    "canonical_rows": int(len(M)),
    "matched_rows": int(hit.sum()),
    "match_rate": float(hit.mean()),
    "matched_targets": int(M.loc[hit, "target_id_canonical"].nunique()),
    "total_targets": int(M.target_id_canonical.nunique()),
    "matched_families": int(M.loc[hit, "family_id_canonical"].nunique()),
    "total_families": int(M.family_id_canonical.nunique()),
    "matched_targets_ge5": int((n5 >= 5).sum()),
    "evidence_tier": ("TIER 2 — HERMES 在 SKEMPI v2.0 上微调(§50)，"
                      "不能自动进入 Tier 1"),
    "official_predictions": True,
    "inference_needed": False,
    "gpu_needed": False,
}
(OUT / "hermes_matched_meta.json").write_text(
    json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"\n写出 hermes_matched_final.csv / hermes_matched_meta.json")
print("\n★ 关键: HERMES 是**官方发布的 row-level 预测** ⇒ 无需推理、无需 GPU。")
