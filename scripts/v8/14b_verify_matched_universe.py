# -*- coding: utf-8 -*-
"""scripts/v8/14b_verify_matched_universe.py -- 验证 remedy 序列的**宇宙是否真的匹配**。

`14_remedy_case_study.py` 显示 V7.3 全部 F 族都在 `n_targets_ge5 = 102` 层，
但这**不等于**它们用的是同一批行。本脚本直接比对**行身份**：

  * V7.3 family_predictions.npz 有 `row_id`（3351 行）-> 与 canonical 匹配
  * V7.2R1 full_ranker_predictions.npz 有 `row_id`（4233 行）
  * V7.2R2 trackB_predictions.npz 有 `row_id`（4219 行）
  * V7.2R3 ppi_predictions.npz 有 `row_id`（3452 行）

输出两两 Jaccard / 交集，确认哪些 remedy 比较是**严格匹配**的。
"""

from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
OUT = R / "experiments" / "v8"
lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)
CANON = set(lab.pair_row.astype(int))

ART = {
    "V7.2R1_full": "experiments/v7_2r_1/phase2_full_eval/full_ranker_predictions.npz",
    "V7.2R2_trackB": "experiments/v7_2r_2/phase5_eval/trackB_predictions.npz",
    "V7.2R3_ppi": "experiments/v7_2r_3/phase5_eval/ppi_predictions.npz",
    "V7.3_physics": "experiments/v7_3/phase5_rankers/family_predictions.npz",
    "V7.2_nested": "experiments/v7_2/phase5_nested_cv/ranker_predictions.npz",
    "V7.1_z": "experiments/v7_1/phase2_rank/z_predictions.npz",
}
rows = []
for name, rel in ART.items():
    p = R / rel
    if not p.exists():
        print(f"  {name}: MISSING"); continue
    z = np.load(p, allow_pickle=True)
    rid = None
    for k in ("row_id", "pair_row", "rows", "row"):
        if k in z.files:
            rid = set(np.asarray(z[k]).ravel().astype(int)); break
    if rid is None:
        print(f"  {name}: 无 row_id 键 -> 需按行序对齐，本脚本跳过")
        continue
    rows.append({"name": name, "n": len(rid),
                 "in_canonical": len(rid & CANON),
                 "not_in_canonical": len(rid - CANON)})

D = pd.DataFrame(rows)
print("=== 各制品的行集合 ===")
print(D.to_string(index=False))

print("\n=== 两两交集 / Jaccard ===")
sets = {}
for name, rel in ART.items():
    p = R / rel
    if not p.exists():
        continue
    z = np.load(p, allow_pickle=True)
    for k in ("row_id", "pair_row", "rows", "row"):
        if k in z.files:
            sets[name] = set(np.asarray(z[k]).ravel().astype(int)); break

out = []
names = sorted(sets)
print(f"{'A':<16}{'B':<16}{'|A|':>7}{'|B|':>7}{'|A∩B|':>8}{'Jaccard':>9}")
for a, b in combinations(names, 2):
    A, B = sets[a], sets[b]
    inter = len(A & B); uni = len(A | B)
    j = inter / uni if uni else float("nan")
    print(f"{a:<16}{b:<16}{len(A):>7}{len(B):>7}{inter:>8}{j:>9.4f}")
    out.append({"A": a, "B": b, "nA": len(A), "nB": len(B),
                "n_inter": inter, "jaccard": j,
                "exact_match": bool(A == B)})
pd.DataFrame(out).to_csv(OUT / "remedy" / "universe_overlap.csv", index=False,
                         encoding="utf-8")

print("\n=== 判读 ===")
ex = [o for o in out if o["exact_match"]]
print(f"  完全相同的行集合对：{len(ex)}")
for o in ex:
    print(f"    {o['A']} == {o['B']}")
print()
print("  ⇒ 只有**完全相同**的行集合才能直接比较 remedy 的增益；")
print("    否则必须声明为**跨宇宙比较**（§21）。")
