# -*- coding: utf-8 -*-
"""scripts/v8/00c_probe_alignment.py -- 确定「无 row_id」制品能否按**行序**对齐（只读）。

多个制品有 4473 / 7523 行，与 canonical(4473) 或 matched_analysis_set(7523)
行数一致。需要判定它们是否**行对齐**（同序），否则不能使用。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)
mas = pd.read_csv(R / "experiments" / "v6_2" / "phase2_split4" /
                  "matched_analysis_set.csv", low_memory=False)
CANON_ROWS = lab.pair_row.to_numpy(int)
MAS_ROWS = mas.row.to_numpy(int) if "row" in mas.columns else None
print("canonical rows =", len(CANON_ROWS), "| matched_analysis_set rows =",
      (len(MAS_ROWS) if MAS_ROWS is not None else None))

TARGETS = [
    "experiments/v7_0/phase2_ablation/ablation_predictions.npz",
    "experiments/v7_0/phase3_losses/loss_predictions.npz",
    "experiments/v7_1/phase5_lsr/lsr_predictions.npz",
    "experiments/v7_1/phase2_rank/z_predictions.npz",
    "experiments/v7_3/phase5_rankers/family_predictions.npz",
    "experiments/v7_3/phase9_calibration/calibrated_predictions.npz",
    "experiments/v6_2/phase2_split4/split4_predictions.npz",
    "experiments/v6_2_1/phase3_d_only/d_only_predictions.npz",
]
for rel in TARGETS:
    p = R / rel
    if not p.exists():
        print(f"\n--- {rel}\n    MISSING"); continue
    z = np.load(p, allow_pickle=True)
    keys = list(z.files)
    print(f"\n--- {rel}")
    print(f"    keys({len(keys)}): {keys[:18]}")
    for k in keys:
        a = np.asarray(z[k])
        if a.ndim == 2:
            print(f"      {k}: shape={a.shape} dtype={a.dtype}")
    # 找行标识
    for cand in ("row_id", "row", "rows", "pair_row", "pdb", "pdb_id", "y"):
        if cand in keys:
            a = np.asarray(z[cand]).ravel()
            print(f"    [{cand}] n={a.size} sample={a[:6]}")
            if cand in ("row_id", "row", "rows", "pair_row"):
                ri = a.astype(int)
                if len(ri) == len(CANON_ROWS):
                    print(f"      vs canonical 逐位相同: {bool(np.array_equal(ri, CANON_ROWS))}")
                if MAS_ROWS is not None and len(ri) == len(MAS_ROWS):
                    print(f"      vs matched_set 逐位相同: {bool(np.array_equal(ri, MAS_ROWS))}")
            if cand in ("pdb", "pdb_id"):
                print(f"      n_unique={len(set(a.astype(str)))}")
                # 与 canonical target 序列比较
                if len(a) == len(CANON_ROWS):
                    eq = (a.astype(str) == lab.target_id_canonical.to_numpy(str))
                    print(f"      与 canonical target 序列一致率 = {eq.mean():.4f}")
