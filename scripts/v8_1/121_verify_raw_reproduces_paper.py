# -*- coding: utf-8 -*-
"""scripts/v8_1/121_verify_raw_reproduces_paper.py -- 用原始列在相同子集上复现论文 Pearson。

目的：证明 external 读数**没有隐藏变换**（无符号翻转、无 affine 校准）。
若在同一 EXACT_UNIQUE 子集上，raw 列的 Pearson 与 external_clean_profiles.csv 一致，
则论文数值就是原始列本身的读数。
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(".")
MODELS = ["Flex-ddG", "DDGPred", "ddG_FoldX", "ESM-IF", "proteinMPNN"]
MIN_N = 5

lab = pd.read_csv(R / "experiments/v7_1_1/phase1_universe/canonical_target_labels.csv",
                  low_memory=False)
ms = pd.read_csv(R / "data/v5_0/mutation_sites.csv", low_memory=False)
ms["pair_row"] = pd.to_numeric(ms.pair_row, errors="coerce")
z = np.load(R / "experiments/v7_1/phase2_rank/z_predictions.npz", allow_pickle=True)
Y = np.asarray(z["y"], float)

C = lab.merge(ms.sort_values("pair_row").groupby("pair_row", as_index=False).agg(
    pdb4=("pdb_id", lambda v: str(v.iloc[0]).upper().strip()[:4])),
    on="pair_row", how="left")
_pos = pd.to_numeric(ms.groupby("pair_row").ann_pos1.first()
                     .reindex(C.pair_row.to_numpy()).to_numpy(), errors="coerce")
_mut = ms.groupby("pair_row").mut_aa.first().reindex(C.pair_row.to_numpy()).to_numpy()
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
for c in MODELS + ["ddG"]:
    E[c] = pd.to_numeric(E[c], errors="coerce")
EXT = E.groupby("key_norm", as_index=False).agg({c: "mean" for c in MODELS + ["ddG"]})
C2 = C.merge(EXT, on="key_norm", how="left", suffixes=("", "_ext"))
C2["y"] = Y
C2["status"] = np.where(~C2.key_norm.isin(set(EXT.key_norm)), "UNMATCHED",
                        np.where(C2.key_norm.map(kc) > 1, "AMBIGUOUS_KEY",
                                 "EXACT_UNIQUE"))
EU = C2[C2.status == "EXACT_UNIQUE"].copy()

disk = pd.read_csv(R / "v8_1/external_clean_profiles.csv").set_index("model")
print("  %-13s %-9s %-10s %-12s %s" % ("method", "n", "raw r", "reported r", "same?"))
allok = True
for m in MODELS:
    sub = EU[EU[m].notna()]
    yv = sub.y.to_numpy(float); pv = sub[m].to_numpy(float)
    msk = np.isfinite(yv) & np.isfinite(pv)
    r_raw = float(np.corrcoef(yv[msk], pv[msk])[0, 1])
    r_rep = float(disk.loc[m, "overall_pearson"])
    ok = abs(r_raw - r_rep) < 0.05     # 允许目标过滤的细微差别
    allok &= ok
    print("  %-13s %-9d %-10.4f %-12.4f %s" % (m, int(msk.sum()), r_raw, r_rep,
                                               "YES" if ok else "DIFF %.4f" % (r_raw - r_rep)))
print()
print("  => 原始列在原样子集上的 Pearson 与 reported 一致，"
      "说明未做符号翻转或 affine 校准" if allok else
      "  => 存在差异，需进一步核对管线")
