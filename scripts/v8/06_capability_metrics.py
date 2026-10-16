# -*- coding: utf-8 -*-
"""scripts/v8/06_capability_metrics.py -- V8 PHASE 6：三轴能力记分板。

对注册表中每个预测键，在其**原生连接宇宙**上计算三轴指标。
prediction 已固定，只做评估。

输出：
    experiments/v8/BENCHMARK_SCOREBOARD.csv
    experiments/v8/metrics/scoreboard_<universe>.csv
    experiments/v8/metrics/decomposition_qc.csv
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(R / "scripts" / "v8"))
import v8metrics as M  # noqa: E402

OUT = R / "experiments" / "v8"
MET = OUT / "metrics"
MET.mkdir(parents=True, exist_ok=True)

lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)
mas = pd.read_csv(R / "experiments" / "v6_2" / "phase2_split4" /
                  "matched_analysis_set.csv", low_memory=False)
CANON = lab[["pair_row", "target_id_canonical", "family_id_canonical"]].rename(
    columns={"pair_row": "row_id", "target_id_canonical": "target_id",
             "family_id_canonical": "family_id"})
CANON_ROWS = CANON.row_id.to_numpy(int)
CANON_TGT = CANON.target_id.to_numpy(str)
MAS_ROWS = mas.row.to_numpy(int)
MAS_TGT = mas.pdb_id.astype(str).to_numpy()
MAS_FAM = mas.family_label.astype(str).to_numpy()

# y（真值 ΔΔG）取自已验证制品
_z = np.load(R / "experiments" / "v7_1" / "phase2_rank" / "z_predictions.npz",
             allow_pickle=True)
Y_CANON = np.asarray(_z["y"], float)
assert len(Y_CANON) == len(CANON_ROWS), "y 与 canonical 行数不一致"
print(f"canonical y: n={len(Y_CANON)} mean={Y_CANON.mean():.4f} sd={Y_CANON.std():.4f}")

# matched_analysis_set 的 y（mas 自带 ddG）
Y_MAS = pd.to_numeric(mas["ddG"], errors="coerce").to_numpy(float)
print(f"mas y: n={len(Y_MAS)} finite={np.isfinite(Y_MAS).sum()}")

NONPRED = {"y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows",
           "z_raw", "y_cal", "family", "seed", "masks"}


def is_label_key(k: str) -> bool:
    """判定某键是否**就是真值标签**（必须排除，否则污染记分板）。

    本次实测泄漏键全部形如 `L1_random_y` / `L2_complex_y` /
    `L4_family_y` / `L3_seq40_y` —— 即真值数组被以 `<level>_<split>_y` 命名。
    """
    kl = str(k).lower()
    return kl == "y" or kl.endswith("_y") or kl.endswith("_ddg") \
        or kl in ("true", "label", "labels", "target")


# 不可能的相关性上限：任何真实模型都不该超过此值（否则视为泄漏）
IMPOSSIBLE_R = 0.98

reg = pd.read_csv(OUT / "model_registry.csv", low_memory=False)
print(f"注册条目 = {len(reg)}")
_cache: dict = {}
t0 = time.time()
rows, qc, leaks = [], [], []
for i, r in enumerate(reg.itertuples()):
    rel = r.prediction_file
    key = r.prediction_key
    if key in NONPRED or is_label_key(key):
        continue
    if rel not in _cache:
        _cache[rel] = np.load(R / rel, allow_pickle=True)
    z = _cache[rel]
    if key not in z.files:
        continue
    a = np.asarray(z[key])
    if a.ndim != 1:
        continue
    try:
        a = a.astype(np.float64)
    except (ValueError, TypeError):
        continue
    # 硬性防泄漏：预测与真值近乎完全相同 -> 视为标签混入，排除
    if len(a) == len(Y_CANON) and np.isfinite(a).all():
        _r = np.corrcoef(a, Y_CANON)[0, 1]
        if np.isfinite(_r) and abs(_r) > IMPOSSIBLE_R:
            leaks.append({"model_id": r.model_id, "key": key, "r_vs_y": float(_r),
                          "action": "EXCLUDED_AS_LABEL_LEAK"})
            continue

    if r.join_mode == "canonical_rowid":
        rid = None
        for kk in ("row_id", "pair_row", "rows", "row"):
            if kk in z.files:
                rid = np.asarray(z[kk]).ravel().astype(int); break
        if rid is None:
            y = Y_CANON; tgt = CANON_TGT; fam = CANON.family_id.to_numpy()
            universe = "U_CANONICAL_ALIGN"
        else:
            pos = pd.Series(np.arange(len(CANON_ROWS)), index=CANON_ROWS).reindex(rid)
            ok = pos.notna().to_numpy()
            if ok.sum() != len(rid):
                y = Y_CANON[pos[ok].to_numpy(int)]
                tgt = CANON_TGT[pos[ok].to_numpy(int)]
                fam = CANON.family_id.to_numpy()[pos[ok].to_numpy(int)]
                a = a[ok]; universe = "U_CANONICAL_SUBSET"
            else:
                y = Y_CANON[pos.to_numpy(int)]
                tgt = CANON_TGT[pos.to_numpy(int)]
                fam = CANON.family_id.to_numpy()[pos.to_numpy(int)]
                universe = "U_CANONICAL"
    elif r.join_mode == "canonical_align":
        y = Y_CANON; tgt = CANON_TGT; fam = CANON.family_id.to_numpy()
        universe = "U_CANONICAL"
    else:  # mas_align
        y = Y_MAS; tgt = MAS_TGT; fam = MAS_FAM
        universe = "U_MAS"

    m = np.isfinite(y) & np.isfinite(a)
    if m.sum() < 50:
        continue
    U = M.prep_universe(y[m], a[m], tgt[m], fam[m])
    res = M.compute_all(U)
    res.update({"model_id": r.model_id, "model_name": key,
                "model_family": r.model_family,
                "resource_tier": r.resource_tier,
                "provenance_status": r.provenance_status,
                "prediction_file": rel, "universe": universe,
                "join_mode": r.join_mode})
    rows.append(res)
    qc.append({"model_id": r.model_id, "decomp_residual": res.get("decomp_residual"),
               "n_rows": res.get("n_rows"), "n_targets_ge5": res.get("n_targets_ge5")})

S = pd.DataFrame(rows)
front = ["model_id", "model_name", "model_family", "resource_tier",
         "provenance_status", "prediction_file", "join_mode", "universe",
         "n_rows", "n_targets", "n_targets_ge5", "min_n"]
cols = front + [c for c in M.METRIC_COLS if c in S.columns]
S = S[[c for c in cols if c in S.columns]]
S.to_csv(OUT / "BENCHMARK_SCOREBOARD.csv", index=False, encoding="utf-8")
for u, g in S.groupby("universe"):
    g.to_csv(MET / f"scoreboard_{u}.csv", index=False, encoding="utf-8")
pd.DataFrame(qc).to_csv(MET / "decomposition_qc.csv", index=False, encoding="utf-8")
pd.DataFrame(leaks).to_csv(MET / "label_leak_guard.csv", index=False, encoding="utf-8")
print(f"标签泄漏守卫排除 {len(leaks)} 个键")
for L in leaks[:10]:
    print("   EXCLUDED %-28s r_vs_y=%.6f" % (L["key"], L["r_vs_y"]))

print(f"\n完成 {len(S)} 个模型评估 ({time.time() - t0:.1f}s)")
print("\n=== 按 universe ===")
print(S.groupby("universe").agg(models=("model_id", "size"),
                                max_rows=("n_rows", "max"),
                                max_ge5=("n_targets_ge5", "max")).to_string())
print("\n=== 分解恒等式 QC（|resid| 最大值应为 ~1e-15）===")
q = pd.DataFrame(qc)
print(f"  max |resid| = {q.decomp_residual.abs().max():.3e}")
print("\n=== U_CANONICAL 上按 macro_spearman 排名前 20 ===")
u = S[S.universe == "U_CANONICAL"].copy()
show = ["model_name", "model_family", "overall_pearson", "rmse",
        "macro_spearman", "C_W", "efficiency@10", "ndcg@10"]
print(u.nlargest(20, "macro_spearman")[show].to_string(index=False))
print("\n=== U_CANONICAL 上按 overall_pearson 排名前 12 ===")
print(u.nlargest(12, "overall_pearson")[show].to_string(index=False))
