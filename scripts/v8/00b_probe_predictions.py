# -*- coding: utf-8 -*-
"""scripts/v8/00b_probe_predictions.py -- 探查每个预测制品的结构与可连接性（只读）。

回答三个问题：
  1. 里面有哪些 key / 列？
  2. 是否含 `row_id`（权威连接键）？
  3. 能否连到 canonical_target_labels（row_id）？

这是 §9 reuse-first 的**前置调查** —— 先知道能借什么，再决定借什么。
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
OUT = R / "experiments" / "v8" / "registry"
OUT.mkdir(parents=True, exist_ok=True)

lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)
ROW = set(lab.pair_row.astype(int))
print(f"canonical rows = {len(ROW)} | targets = {lab.target_id_canonical.nunique()}")

inv = pd.read_csv(R / "experiments" / "v8" / "phase0_freeze" /
                  "artifact_inventory.csv")
cands = inv[(inv.kind == "prediction_like") &
            (inv.ext.isin([".npz", ".csv"]))].sort_values("size_bytes",
                                                          ascending=False)
print(f"候选制品 = {len(cands)}")

rows = []
for r in cands.itertuples():
    p = R / r.path
    if not p.exists():
        continue
    rec = {"path": r.path, "ext": r.ext, "size_kb": round(r.size_bytes / 1024, 1),
           "keys_or_cols": "", "n_rows": None, "has_row_id": False,
           "row_id_match_canonical": None, "note": ""}
    try:
        if r.ext == ".npz":
            z = np.load(p, allow_pickle=True)
            keys = list(z.files)
            rec["keys_or_cols"] = "|".join(keys[:14])
            # 找 row_id 与可能的预测向量
            rid = None
            for k in ("row_id", "row_ids", "rows", "pair_row"):
                if k in keys:
                    rid = np.asarray(z[k]).ravel(); break
            if rid is not None:
                rec["n_rows"] = int(rid.size)
                ri = set(int(x) for x in rid)
                rec["has_row_id"] = True
                rec["row_id_match_canonical"] = len(ri & ROW)
            else:
                # 退而求其次：用任意 1-D 整型数组长度
                for k in keys:
                    a = np.asarray(z[k])
                    if a.ndim == 1:
                        rec["n_rows"] = int(a.size); break
                rec["note"] = "no row_id key"
        else:
            d = pd.read_csv(p, low_memory=False)
            rec["keys_or_cols"] = "|".join(map(str, d.columns[:16]))
            rec["n_rows"] = len(d)
            c_low = {c.lower(): c for c in d.columns}
            hit = None
            for k in ("row_id", "pair_row", "row"):
                if k in c_low:
                    hit = c_low[k]; break
            if hit:
                ri = set(pd.to_numeric(d[hit], errors="coerce").dropna().astype(int))
                rec["has_row_id"] = True
                rec["row_id_match_canonical"] = len(ri & ROW)
            else:
                rec["note"] = "no row_id column"
    except Exception as e:  # noqa: BLE001
        rec["note"] = f"ERR {type(e).__name__}:{str(e)[:60]}"
    rows.append(rec)

P = pd.DataFrame(rows)
P.to_csv(OUT / "prediction_probe.csv", index=False, encoding="utf-8")

print()
print("=== 含 row_id 且能连上 canonical 的制品 ===")
ok = P[(P.has_row_id == True) & (P.row_id_match_canonical > 0)].sort_values(
    "row_id_match_canonical", ascending=False)
for r in ok.itertuples():
    print("  %-72s rows=%-6s match=%-6s  %s" %
          (r.path[-72:], r.n_rows, r.row_id_match_canonical, r.keys_or_cols[:60]))
print()
print("=== 不含 row_id 的制品（需另行映射）===")
bad = P[P.has_row_id != True]
for r in bad.itertuples():
    print("  %-72s n=%-6s %s  [%s]" %
          (r.path[-72:], r.n_rows, r.keys_or_cols[:52], r.note[:40]))
print()
print(f"可连 canonical 的制品数 = {len(ok)} / {len(P)}")
