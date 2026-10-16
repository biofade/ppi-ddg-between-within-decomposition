# -*- coding: utf-8 -*-
"""scripts/v8/_test_metrics.py -- v8metrics 正确性自检（合成数据，确定性）。

构造三个已知答案的模型：
  perfect   : p == y                     -> 各项应为最优
  reversed  : p == -y                    -> 排序应为最差、C_W 应为负
  noise     : p 为确定性伪随机            -> 应接近 0
  const     : p 恒定                      -> 应恰当返回 NaN（不崩溃）
并校验 §28 的分解恒等式 r_full ≈ C_B + C_W。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v8metrics as M  # noqa: E402

rng = np.random.RandomState(0)
T = 8
NS = [6, 8, 10, 12, 6, 20, 7, 9]
y = np.concatenate([np.round(rng.randn(n) * 1.5, 3) for n in NS])
tgt = np.concatenate([[f"T{i:02d}"] * n for i, n in enumerate(NS)])
fam = np.concatenate([[f"F{i//2}"] * n for i, n in enumerate(NS)])
print(f"synthetic: rows={len(y)} targets={T} families={len(set(fam))}")

cases = {}
for name, p in (("perfect", y.copy()),
                ("reversed", -y.copy()),
                ("noise", rng.randn(len(y))),
                ("const", np.zeros(len(y)))):
    U = M.prep_universe(y, p, tgt, fam)
    cases[name] = M.compute_all(U)

cols = ["overall_pearson", "rmse", "macro_spearman", "C_W_row", "C_B_row",
        "efficiency@10", "ndcg@10", "hit@10", "ef@5"]
print()
hdr = "%-9s" % "case" + "".join("%14s" % c for c in cols)
print(hdr); print("-" * len(hdr))
for name, r in cases.items():
    line = "%-9s" % name
    for c in cols:
        v = r.get(c)
        line += "%14s" % ("nan" if v is None or (isinstance(v, float) and np.isnan(v))
                          else f"{v:+.4f}")
    print(line)

print()
ok = True
def chk(cond, msg):
    global ok
    print(("  ✅ " if cond else "  ❌ ") + msg)
    if not cond: ok = False

pf, rv, nz = cases["perfect"], cases["reversed"], cases["noise"]
chk(abs(pf["overall_pearson"] - 1) < 1e-9, "perfect: overall_pearson == 1")
chk(abs(pf["rmse"]) < 1e-9, "perfect: rmse == 0")
chk(abs(pf["macro_spearman"] - 1) < 1e-9, "perfect: macro_spearman == 1")
chk(abs(pf["ndcg@10"] - 1) < 1e-9, "perfect: ndcg@10 == 1")
chk(abs(pf["efficiency@10"] - 1) < 1e-9, "perfect: efficiency@10 == 1")
chk(abs(pf["hit@10"] - 1) < 1e-9, "perfect: hit@10 == 1")
chk(abs(rv["overall_pearson"] + 1) < 1e-9, "reversed: overall_pearson == -1")
chk(rv["macro_spearman"] < 0, "reversed: macro_spearman < 0")
chk(rv["C_W_row"] < 0, "reversed: C_W < 0")
chk(abs(nz["overall_pearson"]) < 0.35, "noise: |overall_pearson| small")
chk(np.isfinite(cases["const"]["rmse"]), "const: rmse 有限（不崩溃）")
for name, r in cases.items():
    d = r.get("decomp_residual")
    if d is not None and np.isfinite(d):
        chk(abs(d) < 1e-9, f"{name}: §28 分解 r_full ≈ C_B + C_W (resid={d:.2e})")
print()
print("METRICS_SELFTEST =", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
