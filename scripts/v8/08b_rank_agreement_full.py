# -*- coding: utf-8 -*-
"""scripts/v8/08b_rank_agreement_full.py -- 解决**模型集不一致**导致的关键矛盾。

## 矛盾

* `08_rank_reversal.py` 在 **160 个模型**（含 canonical_align 键）上得
  ρ(Overall r, Macro ρ) = **+0.181**
* `15_bootstrap_fast.py` 在 **67 个模型**（仅 canonical_rowid 键）上得
  ρ(Overall r, Macro ρ) = **−0.264**

符号相反。必须查清是**模型集效应**还是**实现差异**。

## 做法

1. 用**同一套解析口径**（本脚本自带，不依赖 scoreboard 的来源过滤），
   枚举**全部**可用模型键：canonical_rowid + canonical_align。
2. 分别报告：全集 / 仅 rowid / 仅 align 的秩一致性。
3. 找出驱动差异的具体模型。
"""

from __future__ import annotations

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
lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)
CAN = lab[["pair_row", "target_id_canonical", "family_id_canonical"]].rename(
    columns={"pair_row": "row_id", "target_id_canonical": "target_id",
             "family_id_canonical": "family_id"})
ROWS = CAN.row_id.to_numpy(int)
TGT = CAN.target_id.to_numpy(str)
FAM = CAN.family_id.to_numpy(str)
z = np.load(R / "experiments" / "v7_1" / "phase2_rank" / "z_predictions.npz",
           allow_pickle=True)
Y = np.asarray(z["y"], float)

ART = [
    ("experiments/v7_1/phase2_rank/z_predictions.npz", "rowid"),
    ("experiments/v7_2r_1/phase2_full_eval/full_ranker_predictions.npz", "rowid"),
    ("experiments/v7_2r_2/phase5_eval/trackB_predictions.npz", "rowid"),
    ("experiments/v7_2r_3/phase5_eval/ppi_predictions.npz", "rowid"),
    ("experiments/v7_2r_3/phase8_noise_control/noise_control_predictions.npz", "rowid"),
    ("experiments/v7_3/phase5_rankers/family_predictions.npz", "rowid"),
    ("experiments/v7_2/phase5_nested_cv/ranker_predictions.npz", "rowid"),
    ("experiments/v7_0/phase2_ablation/ablation_predictions.npz", "align"),
    ("experiments/v7_0/phase3_losses/loss_predictions.npz", "align"),
    ("experiments/v7_1/phase5_lsr/lsr_predictions.npz", "align"),
]
NONPRED = {"y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows",
           "z_raw", "y_cal", "family", "seed"}


def is_label(k):
    kl = str(k).lower()
    return kl == "y" or kl.endswith("_y") or kl.endswith("_ddg")


def metrics_for(P, y, tgt):
    """P: (M, N)。返回 (M, n) 指标矩阵 [overall_r, macro_sp, eff10]。"""
    M_, N_ = P.shape
    uq, inv = np.unique(tgt, return_inverse=True)
    cnt = np.bincount(inv)
    starts = np.concatenate([[0], np.cumsum(cnt)[:-1]])
    elig = np.where(cnt >= 5)[0]
    Yc = y - y.mean()
    Pc = P - P.mean(axis=1, keepdims=True)
    den = np.sqrt((Yc ** 2).sum()) * np.sqrt((Pc ** 2).sum(axis=1))
    den[den == 0] = np.nan
    overall = (Pc @ Yc) / den
    macro = np.zeros(M_); eff = np.zeros(M_)
    for g in elig:
        s = starts[g]; e = s + cnt[g]; n = cnt[g]
        yg = y[s:e]; pg = P[:, s:e]
        true_order = np.argsort(yg, kind="stable")
        rank = np.empty(n, int); rank[true_order] = np.arange(n)
        gain = 1.0 + (n - rank)
        ideal = np.argsort(-gain, kind="stable")
        kk = max(1, min(10, n))
        disc = np.log2(np.arange(2, kk + 2))
        pred = np.argsort(pg, axis=1, kind="stable")[:, :kk]
        gid = (gain[ideal[:kk]] / disc).sum()          # 理想 DCG（带折扣）
        gsel_disc = (gain[pred] / disc).sum(axis=1)    # 选中 DCG（同折扣）
        if gid > 0:
            eff += gsel_disc / gid
        ry = np.argsort(np.argsort(yg, kind="stable"), kind="stable").astype(float)
        rp = np.argsort(np.argsort(pg, axis=1, kind="stable"),
                        axis=1, kind="stable").astype(float)
        ryc = ry - ry.mean(); rpc = rp - rp.mean(axis=1, keepdims=True)
        d = np.sqrt((ryc ** 2).sum()) * np.sqrt((rpc ** 2).sum(axis=1))
        d[d == 0] = np.nan
        macro += (rpc @ ryc) / d
    ne = len(elig)
    if ne:
        macro /= ne; eff /= ne
    return np.vstack([overall, macro, eff]).T


def collect(which):
    P, names = [], []
    for rel, kind in ART:
        if which != "all" and kind != which:
            continue
        p = R / rel
        if not p.exists():
            continue
        zz = np.load(p, allow_pickle=True)
        for k in zz.files:
            if k in NONPRED or is_label(k):
                continue
            a = np.asarray(zz[k])
            try:
                a = a.astype(float)
            except (ValueError, TypeError):
                continue
            if a.ndim != 1 or len(a) != len(Y):
                continue
            if abs(np.corrcoef(a, Y)[0, 1]) > 0.98:
                continue
            P.append(a); names.append(k)
    return np.vstack(P), names


for which, label in (("all", "rowid + align 全集"),
                     ("rowid", "仅 canonical_rowid"),
                     ("align", "仅 canonical_align")):
    P, names = collect(which)
    mt = metrics_for(P, Y, TGT)
    v = np.isfinite(mt[:, 0]) & np.isfinite(mt[:, 1]) & np.isfinite(mt[:, 2])
    r = spearmanr(mt[v, 0], mt[v, 1]).correlation
    k = kendalltau(mt[v, 0], mt[v, 1]).correlation
    re = spearmanr(mt[v, 0], mt[v, 2]).correlation
    ke = kendalltau(mt[v, 0], mt[v, 2]).correlation
    print(f"\n[{label}] M={v.sum()}")
    print(f"   ρ(Overall r, Macro ρ) = {r:+.4f}   τ = {k:+.4f}")
    print(f"   ρ(Overall r, Eff@10)  = {re:+.4f}   τ = {ke:+.4f}")
    if which == "align":
        # 打印这些键的读数，看它们是否驱动了符号
        idx = np.where(v)[0]
        for i in idx[:40]:
            print(f"     {names[i][:34]:36s} r={mt[i,0]:+.3f} rho={mt[i,1]:+.3f} "
                  f"eff={mt[i,2]:+.4f}")

print("\n=== 判读 ===")
print("若 align 子集的 ρ(Overall r, Macro ρ) 明显为负、而 rowid 子集为正（或反之），")
print("则该秩一致性是**模型集组成效应**，不是稳健的能力关系；")
print("⇒ V8 必须报告「模型集条件下的区间」，不得给出单一数字。")
