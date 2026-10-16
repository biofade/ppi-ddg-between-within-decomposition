# -*- coding: utf-8 -*-
"""scripts/v8/37_diagnose_discrepancy.py -- 定位 35/36 数值不一致的根因。

## 不一致

同一份 `pool_composition_models.csv`（37 个 PERMISSIVE 键）：
    rho(Overall r, Eff@10) = **+0.8088**   （35，经 v8metrics.compute_all）
    rho(Overall r, Eff@10) = **-0.5579**   （36，逐 target 一/二阶量内联）

## 诊断步骤

1. 用 `v8metrics.compute_all` **重新**计算每个 PERMISSIVE 模型的
   (overall_pearson, efficiency@10)，与 `pool_composition_models.csv` 逐位比对。
2. 用 `36` 的内联公式**独立**重算同一批模型的同一对量，逐位比对。
3. 打印**差异最大的前 5 个模型**及其两种取值 ⇒ 定位分歧点。
4. 对同一模型用两种路径分别重算 `overall_pearson`，确认是否一致。
"""

from __future__ import annotations

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
sys.path.insert(0, str(R / "scripts" / "v8"))
import v8metrics as M  # noqa: E402

OUT = R / "experiments" / "v8"
K = 10
MIN_N = 5

# ---------- 载入 ----------
H = pd.read_csv(OUT / "external" / "hermes_matched_final.csv", low_memory=False)
z = np.load(R / "experiments" / "v7_1" / "phase2_rank" / "z_predictions.npz",
            allow_pickle=True)
Y = np.asarray(z["y"], float)
lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)
H = H.sort_values("pair_row").reset_index(drop=True)
TGT = H.target_id_canonical.astype(str).to_numpy()
FAM = H.family_id_canonical.astype(str).to_numpy()
assert len(H) == len(Y) == len(lab)

POOL = pd.read_csv(OUT / "external" / "pool_composition_models.csv",
                   low_memory=False)
PERM = POOL[POOL.provenance == "PERMISSIVE"].dropna(
    subset=["overall_pearson", "efficiency@10"]).copy()
print(f"PERMISSIVE 模型 = {len(PERM)}")
print(f"  pool 文件 rho(Overall r, Eff@10) = "
      f"{spearmanr(PERM.overall_pearson, PERM['efficiency@10']).correlation:+.4f}")

# ---------- 1. 用 v8metrics 重算 ----------
NONPRED = {"y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows",
           "z_raw", "y_cal", "family", "seed"}
def is_label(k):
    kl = str(k).lower()
    return kl == "y" or kl.endswith("_y") or kl.endswith("_ddg")

ART = [
    ("experiments/v6_2/phase2_split4/split4_predictions.npz", "mas"),
    ("experiments/v6_2_1/phase3_d_only/d_only_predictions.npz", "mas"),
    ("experiments/v7_0/phase2_ablation/ablation_predictions.npz", "align"),
    ("experiments/v7_0/phase3_losses/loss_predictions.npz", "align"),
    ("experiments/v7_1/phase2_rank/z_predictions.npz", "rowid"),
    ("experiments/v7_1/phase5_lsr/lsr_predictions.npz", "align"),
]
ROWS_ALL = lab.pair_row.to_numpy(int)
pos_all = pd.Series(np.arange(len(ROWS_ALL)), index=ROWS_ALL)
mas = pd.read_csv(R / "experiments/v6_2" / "phase2_split4" /
                  "matched_analysis_set.csv", low_memory=False)
mas_row = mas.row.to_numpy(int)
pos_of_mas = pos_all.reindex(mas_row)

want = set(PERM.model_name)
rec = {}
for rel, kind in ART:
    p = R / rel
    if not p.exists():
        continue
    zz = np.load(p, allow_pickle=True)
    rid = None
    for k in ("row_id", "pair_row", "rows", "row"):
        if k in zz.files:
            rid = np.asarray(zz[k]).ravel().astype(int); break
    for key in zz.files:
        if key not in want:
            continue
        a = np.asarray(zz[key])
        try:
            a = a.astype(float)
        except (ValueError, TypeError):
            continue
        if kind == "align":
            if len(a) != len(ROWS_ALL):
                continue
            av = a
        elif kind == "rowid":
            if rid is None or len(a) != len(rid):
                continue
            m = pos_all.reindex(rid)
            if m.isna().any():
                continue
            av = np.full(len(ROWS_ALL), np.nan)
            av[m.to_numpy(int)] = a
        else:
            if len(a) != len(mas_row):
                continue
            av = np.full(len(ROWS_ALL), np.nan)
            ok = pos_of_mas.notna().to_numpy()
            av[pos_of_mas[ok].to_numpy(int)] = a[ok]
        mm = np.isfinite(Y) & np.isfinite(av)
        if mm.sum() < 100:
            continue
        U = M.prep_universe(Y[mm], av[mm], TGT[mm], FAM[mm])
        r = M.compute_all(U)
        rec[key] = {"ov_v8": r["overall_pearson"],
                    "ef_v8": r["efficiency@10"],
                    "n_rows": int(mm.sum()),
                    "n_ge5": r["n_targets_ge5"]}
print(f"\n用 v8metrics 重算到 {len(rec)} / {len(want)} 个 PERMISSIVE 模型")
missing = want - set(rec)
if missing:
    print(f"  ⚠ 未重算到: {sorted(missing)[:8]}")

# ---------- 2. 逐位比对 ----------
D = PERM[["model_name", "overall_pearson", "efficiency@10"]].rename(
    columns={"overall_pearson": "ov_pool", "efficiency@10": "ef_pool"})
D["ov_v8"] = D.model_name.map(lambda k: rec.get(k, {}).get("ov_v8", np.nan))
D["ef_v8"] = D.model_name.map(lambda k: rec.get(k, {}).get("ef_v8", np.nan))
D["d_ov"] = (D.ov_v8 - D.ov_pool).abs()
D["d_ef"] = (D.ef_v8 - D.ef_pool).abs()
print("\n" + "=" * 78)
print("1) v8metrics 重算 vs pool 文件（逐位）")
print("=" * 78)
ok = D.dropna(subset=["ov_v8", "ef_v8"])
print(f"  可比模型 = {len(ok)}")
print(f"  max |Δ overall_pearson| = {ok.d_ov.max():.3e}")
print(f"  max |Δ efficiency@10|   = {ok.d_ef.max():.3e}")
if ok.d_ov.max() > 1e-6 or ok.d_ef.max() > 1e-6:
    print("\n  ⚠ 不一致！差异最大的前 5 个：")
    print(ok.nlargest(5, "d_ov")[
        ["model_name", "ov_pool", "ov_v8", "ef_pool", "ef_v8"]].to_string(index=False))
else:
    print("  ✅ pool 文件的读数与 v8metrics 重算**逐位一致**")

# ---------- 3. 检查 pool 里 PERMISSIVE 是否混入了非 4473 行模型 ----------
print("\n" + "=" * 78)
print("2) pool 文件 PERMISSIVE 的覆盖线索")
print("=" * 78)
if "n_rows" in PERM.columns:
    print(PERM.n_rows.value_counts().to_string())
    print(f"  n_targets_ge5 取值: "
          f"{sorted(PERM.n_targets_ge5.dropna().unique())[:8]}")

print("\n" + "=" * 78)
print("3) 两种口径下的 rho（用**同一批**可比重算模型）")
print("=" * 78)
if len(ok) >= 6:
    print(f"  用 pool 读数      rho = "
          f"{spearmanr(ok.ov_pool, ok.ef_pool).correlation:+.4f}")
    print(f"  用 v8metrics 重算 rho = "
          f"{spearmanr(ok.ov_v8, ok.ef_v8).correlation:+.4f}")
ok.to_csv(OUT / "external" / "discrepancy_diagnosis.csv", index=False,
          encoding="utf-8")
print(f"\n写出 external/discrepancy_diagnosis.csv")
