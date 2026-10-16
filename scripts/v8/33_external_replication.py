# -*- coding: utf-8 -*-
"""scripts/v8/33_external_replication.py -- PHASE 13 定稿：**外部模型是否复现能力分化**（§68）。

## 问题

V8 的核心发现是：在严格 family-OOD 下，**校准轴（Overall Pearson）与排序轴
（Macro per-target Spearman）对模型的排序几乎不相关**（内部模型上 ρ ≈ −0.20 ~ −0.26）。

§68 的 `STRONGLY_SUPPORTED` 要求这一分化在**外部模型**上也可复现。

## 做法（严格匹配宇宙，§21）

1. 用 HERMES 官方行级预测构造 `U_HERMES_CANONICAL`（单点突变子集，
   2865 行 / 195 target / 94 个 n≥5）。
2. 在**同一批行**上重算：HERMES + 全部可映射的内部模型。
3. 比较 HERMES 与内部模型群在三轴上的位置，并检验秩一致性
   ρ(Overall r, Macro ρ) 在**加入 HERMES 后**是否仍为低/负值。

**符号约定（§19）**：`hermes_proxy = −(log_proba_mt − log_proba_wt)`，已由
`32c` 的方向校验（r=+0.1688 > 0）确认与 V8 的 `ΔΔG>0=结合变弱` 一致。
"""

from __future__ import annotations

import json
import sys
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

OUT = R / "experiments" / "v8"
EXT = OUT / "external"
MIN_N = 5

# ---------- 1. 构造 U_HERMES_CANONICAL ----------
H = pd.read_csv(EXT / "hermes_matched_final.csv", low_memory=False)
print(f"hermes_matched_final: {len(H)} 行, 列={list(H.columns)[:14]}")
z = np.load(R / "experiments" / "v7_1" / "phase2_rank" / "z_predictions.npz",
            allow_pickle=True)
Y_ALL = np.asarray(z["y"], float)
lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)

H = H.sort_values("pair_row").reset_index(drop=True)
assert len(H) == len(lab) == len(Y_ALL), (len(H), len(lab), len(Y_ALL))
H["y_true"] = Y_ALL

sub = H[H.proxy.notna() & H.n_sites.eq(1)].copy()
print(f"\nU_HERMES_CANONICAL (单点且命中): n={len(sub)}")
TGT = sub.target_id_canonical.astype(str).to_numpy()
FAM = sub.family_id_canonical.astype(str).to_numpy()
YT = sub.y_true.to_numpy(float)
HP = sub.proxy.to_numpy(float)
cnt = pd.Series(TGT).value_counts()
elig = cnt[cnt >= MIN_N]
print(f"  targets={len(cnt)}, n>={MIN_N} = {len(elig)}")
print(f"  y_true: mean={YT.mean():.4f} sd={YT.std():.4f}")
print(f"  hermes : mean={HP.mean():.4f} sd={HP.std():.4f}")

# ---------- 2. HERMES 的三轴 ----------
UH = M.prep_universe(YT, HP, TGT, FAM)
Rh = M.compute_all(UH)
print("\n" + "=" * 74)
print("HERMES 三轴读数（U_HERMES_CANONICAL）")
print("=" * 74)
for k in ("n_rows", "n_targets_ge5", "overall_pearson", "rmse", "mae",
          "target_mean_pearson", "macro_spearman", "macro_kendall",
          "pairwise_d05", "C_B_row", "C_W", "efficiency@10", "ndcg@10",
          "decomp_residual"):
    v = Rh.get(k)
    print(f"  {k:<22} {v:+.4f}" if isinstance(v, float) else f"  {k:<22} {v}")

# ---------- 3. 同一批行上重算内部模型 ----------
NONPRED = {"y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows",
           "z_raw", "y_cal", "family", "seed"}
def is_label(k):
    kl = str(k).lower()
    return kl == "y" or kl.endswith("_y") or kl.endswith("_ddg")

ART = [
    ("experiments/v7_1/phase2_rank/z_predictions.npz", "rowid"),
    ("experiments/v7_2r_1/phase2_full_eval/full_ranker_predictions.npz", "rowid"),
    ("experiments/v7_2r_2/phase5_eval/trackB_predictions.npz", "rowid"),
    ("experiments/v7_2r_3/phase5_eval/ppi_predictions.npz", "rowid"),
    ("experiments/v7_3/phase5_rankers/family_predictions.npz", "rowid"),
    ("experiments/v7_0/phase2_ablation/ablation_predictions.npz", "align"),
    ("experiments/v7_0/phase3_losses/loss_predictions.npz", "align"),
    ("experiments/v7_1/phase5_lsr/lsr_predictions.npz", "align"),
]
ROWS_ALL = lab.pair_row.to_numpy(int)
pos_all = pd.Series(np.arange(len(ROWS_ALL)), index=ROWS_ALL)
sel_idx = pos_all.reindex(sub.pair_row.to_numpy(int)).to_numpy()

cache, rows = {}, []
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
        if key in NONPRED or is_label(key):
            continue
        a = np.asarray(zz[key])
        try:
            a = a.astype(float)
        except (ValueError, TypeError):
            continue
        if a.ndim != 1:
            continue
        if kind == "align":
            if len(a) != len(ROWS_ALL):
                continue
            av = a[sel_idx]
        else:
            if rid is None or len(a) != len(rid):
                continue
            m = pos_all.reindex(rid)
            if m.isna().any():
                continue
            # 先把制品对齐到 canonical 位置，再取子集
            canon_pos = m.to_numpy(int)             # 长度 = len(a)
            av_full = np.full(len(ROWS_ALL), np.nan)
            av_full[canon_pos] = a
            av = av_full[sel_idx]
        mm = np.isfinite(YT) & np.isfinite(av)
        if mm.sum() < 100:
            continue
        U = M.prep_universe(YT[mm], av[mm], TGT[mm], FAM[mm])
        r = M.compute_all(U)
        r.update({"model_name": key, "source": rel, "origin": "INTERNAL"})
        rows.append(r)

# HERMES 也作为一个"模型"进入同一张表
Rh2 = dict(Rh)
Rh2.update({"model_name": "HERMES_official_predictions",
            "source": "hermes_experiments (official row-level)",
            "origin": "EXTERNAL"})
rows.append(Rh2)

D = pd.DataFrame(rows)
D.to_csv(EXT / "external_replication_scoreboard.csv", index=False,
         encoding="utf-8")
print(f"\n同一批行上评价的模型数 = {len(D)} "
      f"(内部 {int((D.origin=='INTERNAL').sum())} + 外部 1)")

# ---------- 4. 秩一致性：加入 HERMES 前后 ----------
sel = D.dropna(subset=["overall_pearson", "macro_spearman", "efficiency@10"])
print("\n" + "=" * 74)
print("秩一致性 ρ(Overall r, Macro ρ) —— 加入外部模型前后")
print("=" * 74)
for tag, g in (("仅内部模型", sel[sel.origin == "INTERNAL"]),
               ("内部 + HERMES", sel)):
    if len(g) < 6:
        continue
    s = spearmanr(g.overall_pearson, g.macro_spearman).correlation
    k = kendalltau(g.overall_pearson, g.macro_spearman).correlation
    s2 = spearmanr(g.overall_pearson, g["efficiency@10"]).correlation
    print(f"  {tag:<16} n={len(g):3d}  ρ(r,Macro ρ)={s:+.4f} τ={k:+.4f}   "
          f"ρ(r,Eff@10)={s2:+.4f}")

# ---------- 5. HERMES 在各轴上的排名 ----------
print("\n" + "=" * 74)
print("HERMES 在三轴上的位置（百分位，越大=越好）")
print("=" * 74)
for col, asc in (("overall_pearson", False), ("macro_spearman", False),
                 ("C_W", False), ("efficiency@10", False), ("rmse", True)):
    v = sel[col].astype(float)
    hv = float(sel.loc[sel.origin == "EXTERNAL", col].iloc[0])
    pct = float((v <= hv).mean()) if asc else float((v <= hv).mean())
    print(f"  {col:<20} HERMES={hv:+.4f}  内部中位={v[sel.origin!='EXTERNAL'].median():+.4f}  "
          f"百分位={100*pct:.1f}%")

# ---------- 6. 判定 ----------
print("\n" + "=" * 74)
print("§68 外部复现判定")
print("=" * 74)
gint = sel[sel.origin == "INTERNAL"]
r_int = spearmanr(gint.overall_pearson, gint.macro_spearman).correlation
g_all = sel
r_all = spearmanr(g_all.overall_pearson, g_all.macro_spearman).correlation
h = sel[sel.origin == "EXTERNAL"].iloc[0]
macro_med = gint.macro_spearman.median()
over_med = gint.overall_pearson.median()
same_quadrant = ((h.overall_pearson - over_med) *
                 (h.macro_spearman - macro_med))
print(f"  内部模型 ρ(r, Macro ρ)          = {r_int:+.4f}")
print(f"  加入 HERMES 后 ρ(r, Macro ρ)     = {r_all:+.4f}")
print(f"  HERMES: Overall r={h.overall_pearson:+.4f} (内部中位 {over_med:+.4f})")
print(f"  HERMES: Macro ρ  ={h.macro_spearman:+.4f} (内部中位 {macro_med:+.4f})")
print(f"  HERMES 相对内部中位的象限      = "
      f"{'同象限' if same_quadrant >= 0 else '**跨象限**'}")
verdict = ("EXTERNAL_REPLICATION_CONSISTENT"
           if (r_all < 0.35 and same_quadrant >= 0)
           else "EXTERNAL_REPLICATION_DISCORDANT")
print(f"\n  ⇒ EXTERNAL_REPLICATION_VERDICT = {verdict}")
print(f"     判据: ρ(r,Macro ρ) < 0.35（低）**且** HERMES 不跨象限")
meta = {
    "universe": "U_HERMES_CANONICAL",
    "n_rows": int(len(sub)), "n_targets": int(len(cnt)),
    "n_targets_ge5": int(len(elig)),
    "hermes_source": "StatPhysBio/hermes_experiments official row-level CSV",
    "sign_convention": "hermes_proxy = -(log_proba_mt - log_proba_wt)",
    "rho_overall_macro_internal": float(r_int),
    "rho_overall_macro_with_hermes": float(r_all),
    "hermes_overall_pearson": float(h.overall_pearson),
    "hermes_macro_spearman": float(h.macro_spearman),
    "hermes_percentile_overall": float((gint.overall_pearson <= h.overall_pearson).mean()),
    "hermes_percentile_macro": float((gint.macro_spearman <= h.macro_spearman).mean()),
    "verdict": verdict,
}
(EXT / "external_replication_meta.json").write_text(
    json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"\n写出 external_replication_scoreboard.csv / external_replication_meta.json")
