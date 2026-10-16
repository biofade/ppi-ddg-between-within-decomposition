# -*- coding: utf-8 -*-
"""scripts/v7_0/05_losses.py -- V7.0 PHASE 3: Loss Ablation.

只对 PHASE 2 的 top 架构做（§14）。PHASE 2 结果显示 **A1(+Set)** 是唯一在三项
primary 指标上同向改善者，故 PHASE 3 在 **A1 的特征集**上做 loss 消融。

Loss 变体
---------
| 变体 | 定义 |
|---|---|
| **L0** | 绝对 MSE（对照） |
| **L1** | MSE + `λ_center ·` **target-centered** 项（只用**本折 train** 的 target 均值） |
| **L2** | L1 + **pairwise ranking**（同 target 内，δ=0.5 kcal/mol，每 target 上限 200 对） |
| **L3** | **target-balanced** 训练（每 target 等总权重，修正 V6.3 weighting artifact） |

实现：小型 MLP（CPU，`torch` 已装）。**同一 folds、同一 seeds、同一 row universe**
（与 PHASE 2 一致），以便公平对比。

**禁止**任何 label-derived 变量进入特征；target 均值只用**本折 train**。
"""

from __future__ import annotations

import argparse
import gc
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v6_1"))
import v5_common as C  # noqa: E402
import fastmetrics as FM  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "s4", C.ROOT / "scripts" / "v6_2" / "02_split4_ddg.py")
s4 = importlib.util.module_from_spec(_spec)
sys.modules["s4"] = s4
_spec.loader.exec_module(s4)

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

PH1 = C.ROOT / "experiments" / "v7_0" / "phase1_features"
PH3 = C.ROOT / "experiments" / "v7_0" / "phase3_losses"
PH3.mkdir(parents=True, exist_ok=True)
OUT_MET = PH3 / "loss_metrics.csv"
OUT_NPZ = PH3 / "loss_predictions.npz"
OUT_MD = PH3 / "loss_report.md"
FAIL = C.ROOT / "logs" / "v7_0" / "failures.log"
FAIL.parent.mkdir(parents=True, exist_ok=True)

SEEDS = [42, 43, 44]
LEVEL = "L4_family"
PCA_NEW = 32
MIN_PS_N = 5
DELTA_PAIR = 0.5
MAX_PAIRS = 200
EPOCHS = 60
LR = 1e-3
HID = 128
L: list[str] = []

VARIANTS = ["L0_absolute", "L1_centered", "L2_pairwise", "L3_target_balanced"]


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def fail(w, e):
    import traceback
    with open(FAIL, "a", encoding="utf-8") as fh:
        fh.write(f"[{time.strftime('%H:%M:%S')}] [05_losses] {w}: "
                 f"{type(e).__name__}: {e}\n{traceback.format_exc()}\n")


class MLP(nn.Module):
    def __init__(self, d_in, hid=HID):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(d_in, hid), nn.ReLU(), nn.Dropout(0.1),
            nn.Linear(hid, hid // 2), nn.ReLU(),
            nn.Linear(hid // 2, 1))

    def forward(self, x):
        return self.net(x).squeeze(-1)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-minutes", type=float, default=60.0)
    args = ap.parse_args()
    t0 = time.time()
    emit("# V7.0 PHASE 3 — Loss Ablation（A1 特征集）")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- 架构: **A1 = Base + Set**（PHASE 2 唯一三指标同向改善者）")
    emit(f"- 实现: MLP({HID}) CPU torch；epochs={EPOCHS}, lr={LR}")
    emit(f"- δ_pair = {DELTA_PAIR} kcal/mol；每 target 上限 {MAX_PAIRS} 对")
    emit(f"- 强制：target 均值**只用本折 train**；无 label-derived 特征")
    emit("")
    torch.set_num_threads(2)

    lf = pd.read_csv(PH1 / "block_L_local.csv", low_memory=False)
    S = pd.read_csv(PH1 / "block_S_set.csv", low_memory=False)
    can = pd.read_csv(C.ROOT / "experiments" / "v6_2" / "phase2_split4" /
                      "matched_analysis_set.csv", low_memory=False)
    d = C.build_pair_features(verbose=False)
    blk = d["X_blocks"]
    pdb_all = np.array([str(x) for x in d["groups"]])
    rows = lf.pair_row.to_numpy(int)
    y = d["y"][rows].astype(np.float32)
    pdb = pdb_all[rows]
    base = np.hstack([blk["phys_diff"], blk["Z_diff"]])[rows].astype(np.float32)
    # base 也有约 0.03% NaN（phys_diff 缺失）→ 与 Xs 同样做中位数填补，
    # 否则标准化后 NaN 会扩散到整个 MLP（修复 PHASE 3 全 NaN 的根因）
    base = np.array(base, dtype=np.float32, copy=True)
    bmed = np.nanmedian(base, axis=0); bmed = np.where(np.isfinite(bmed), bmed, 0.0)
    bbad = ~np.isfinite(base)
    if bbad.any():
        base[bbad] = np.take(bmed, np.where(bbad)[1])
    Scols = [c for c in S.columns if c != "pdb_id"]
    Xs = S.set_index("pdb_id").reindex(pdb)[Scols].apply(
        pd.to_numeric, errors="coerce").to_numpy(np.float32)
    Xs = np.array(Xs, dtype=np.float32, copy=True)
    med = np.nanmedian(Xs, axis=0); med = np.where(np.isfinite(med), med, 0.0)
    bad = ~np.isfinite(Xs)
    if bad.any():
        Xs[bad] = np.take(med, np.where(bad)[1])
    emit(f"- row universe: **{len(y)}** 行 / **{len(set(pdb))}** targets")
    emit(f"- 特征: base {base.shape[1]} + S {Xs.shape[1]}")
    emit("")

    fam_all = can.family_label.astype(str).to_numpy()
    # 与 PHASE 2 相同的 unit 定义（复用 s4）
    D4 = {"y": y, "pdb": pdb, "fam": fam_all[
              np.isin(can.row.to_numpy(int), rows)],
          "seq40": can.seq40_cluster.astype(str).to_numpy()[
              np.isin(can.row.to_numpy(int), rows)],
          "obs": can.obs_key.astype(str).to_numpy()[
              np.isin(can.row.to_numpy(int), rows)],
          "state": can.state_id.astype(str).to_numpy()[
              np.isin(can.row.to_numpy(int), rows)]}
    units = s4.level_units(D4, LEVEL)
    emit(f"- {LEVEL} unit 数: {len(np.unique(units))}")
    emit("")

    res, preds = [], {}
    for vname in VARIANTS:
        if (time.time() - t0) / 60 > args.max_minutes:
            emit(f"- {vname}: skip (soft time limit)")
            continue
        try:
            yhat = {sd: np.full(len(y), np.nan) for sd in SEEDS}
            for sd in SEEDS:
                folds = s4.make_folds(units, sd, 5)
                torch.manual_seed(sd)
                np.random.seed(sd)
                for tr, te in folds:
                    # ---- train-only 标准化 + PCA ----
                    mu = base[tr].mean(0); sg = base[tr].std(0); sg[sg == 0] = 1
                    Btr = (base[tr] - mu) / sg; Bte = (base[te] - mu) / sg
                    mu2 = Xs[tr].mean(0); sg2 = Xs[tr].std(0); sg2[sg2 == 0] = 1
                    Ztr = (Xs[tr] - mu2) / sg2; Zte = (Xs[te] - mu2) / sg2
                    from sklearn.decomposition import PCA
                    k = min(PCA_NEW, Ztr.shape[1], max(2, len(tr) - 1))
                    pc = PCA(n_components=k, random_state=sd).fit(Ztr)
                    Xtr = np.hstack([Btr, pc.transform(Ztr)]).astype(np.float32)
                    Xte = np.hstack([Bte, pc.transform(Zte)]).astype(np.float32)
                    # ---- 增广（fwd/rev）----
                    Xa = np.vstack([Xtr, -Xtr]); ya = np.r_[y[tr], -y[tr]]
                    ga = np.r_[pdb[tr], pdb[tr]]
                    xt = torch.tensor(Xa); yt = torch.tensor(ya)
                    xv = torch.tensor(Xte)
                    model = MLP(Xa.shape[1])
                    opt = torch.optim.Adam(model.parameters(), lr=LR)
                    # 训练集内 target 均值（**只用 train**）
                    tmean = pd.Series(ya).groupby(pd.Series(ga)).transform("mean")
                    tmean = tmean.to_numpy(np.float32)
                    w = np.ones(len(ya), np.float32)
                    if vname == "L3_target_balanced":
                        cnt = pd.Series(ga).map(
                            pd.Series(ga).value_counts()).to_numpy(float)
                        w = (1.0 / cnt)
                        w = w / w.mean()
                    wt = torch.tensor(w)
                    # 预构 pairwise 对（仅 L2）
                    pairs = []
                    if vname == "L2_pairwise":
                        rs = np.random.RandomState(sd)
                        for c in np.unique(ga):
                            idx = np.where(ga == c)[0]
                            n = len(idx)
                            allp = [(i, j) for i in range(n)
                                    for j in range(i + 1, n)]
                            if not allp:
                                continue
                            if len(allp) > MAX_PAIRS:
                                sel = rs.choice(len(allp), MAX_PAIRS,
                                                replace=False)
                                allp = [allp[i] for i in sel]
                            for i, j in allp:
                                if abs(ya[idx[i]] - ya[idx[j]]) > DELTA_PAIR:
                                    pairs.append((idx[i], idx[j]))
                    for ep in range(EPOCHS):
                        model.train(); opt.zero_grad()
                        out = model(xt)
                        # L_abs（target-balanced 用样本权重）
                        loss_abs = (wt * (out - yt) ** 2).mean()
                        loss = loss_abs
                        if vname in ("L1_centered", "L2_pairwise"):
                            ga_t = torch.tensor(pd.factorize(ga)[0])
                            om = torch.zeros_like(out)
                            ym = torch.zeros_like(yt)
                            for gv in torch.unique(ga_t):
                                m = ga_t == gv
                                om[m] = out[m].mean()
                                ym[m] = yt[m].mean()
                            loss = loss + 0.5 * (((out - om) - (yt - ym)) ** 2).mean()
                        if vname == "L2_pairwise" and pairs:
                            pi = torch.tensor([p[0] for p in pairs])
                            pj = torch.tensor([p[1] for p in pairs])
                            s = torch.sign(yt[pi] - yt[pj])
                            loss = loss + 0.5 * torch.log1p(
                                torch.exp(-s * (out[pi] - out[pj]))).mean()
                        loss.backward(); opt.step()
                    model.eval()
                    with torch.no_grad():
                        yhat[sd][te] = model(xv).numpy()
                    del model
                    gc.collect()
            ens = np.mean([yhat[sd] for sd in SEEDS], axis=0)
            ok = np.isfinite(ens)
            yy, gg, pp = y[ok], pdb[ok], ens[ok]
            mt = _metrics(yy, pp, gg)
            # self-centered 版本（用模型自己预测的 target 均值）
            pm = pd.Series(pp).groupby(pd.Series(gg)).transform("mean").to_numpy()
            mt_c = _metrics(yy, pp - pm, gg)
            res.append({"variant": vname, "level": LEVEL, "n_rows": int(ok.sum()),
                        "n_target": int(len(set(gg))),
                        "macro_spearman_seed_sd": float(np.std(
                            [FM.macro_spearman(y[np.isfinite(yhat[s])],
                                               yhat[s][np.isfinite(yhat[s])],
                                               pdb[np.isfinite(yhat[s])],
                                               min_n=MIN_PS_N) for s in SEEDS],
                            ddof=1)),
                        **mt, **{f"selfcentered_{k}": v for k, v in mt_c.items()}})
            emit(f"| {vname} | r={mt['overall_r']:+.4f} | C_W={mt['C_W']:+.4f} | "
                 f"macro_ρ={mt['macro_spearman']:+.4f} | Eff@10={mt['eff10']:+.4f} | "
                 f"| self-centered macro_ρ={mt_c['macro_spearman']:+.4f} |")
            preds[f"{LEVEL}_{vname}"] = ens.astype(np.float32)
            for s in SEEDS:
                preds[f"{LEVEL}_{vname}_s{s}"] = yhat[s].astype(np.float32)
            pd.DataFrame(res).to_csv(OUT_MET, index=False, encoding="utf-8")
        except Exception as e:  # noqa: BLE001
            fail(vname, e)
            emit(f"| {vname} | FAILED {type(e).__name__} | | | | |")
            gc.collect()
    preds[f"{LEVEL}_y"] = y.astype(np.float32)
    preds[f"{LEVEL}_pdb"] = pdb.astype(str)
    np.savez(OUT_NPZ, **preds)
    R = pd.DataFrame(res)
    R.to_csv(OUT_MET, index=False, encoding="utf-8")
    emit("")
    emit("## 汇总")
    emit("")
    if len(R):
        emit("| Loss | Overall r | C_W | Macro ρ | Kendall | Eff@10 | C_B(row) | "
             "C_B(eq) |")
        emit("|---|---:|---:|---:|---:|---:|---:|---:|")
        for _, r in R.iterrows():
            emit(f"| {r['variant']} | {r['overall_r']:+.4f} | {r['C_W']:+.4f} | "
                 f"{r['macro_spearman']:+.4f} | {r['macro_kendall']:+.4f} | "
                 f"{r['eff10']:+.4f} | {r['C_B_row']:+.4f} | "
                 f"{r['C_B_equal']:+.4f} |")
    emit("")
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MET}\n写出 {OUT_NPZ}\n写出 {OUT_MD}")
    return 0


def _metrics(y, p, g):
    uq = np.unique(g)
    By = np.array([y[g == c].mean() for c in uq])
    Bp = np.array([p[g == c].mean() for c in uq])
    cnt = np.array([(g == c).sum() for c in uq], float)
    den = float(y.std()) * float(p.std())
    if den <= 0:
        return {}
    CBr = float(((cnt / cnt.sum()) * (By - y.mean()) *
                 (Bp - p.mean())).sum()) / den
    CBe = float(((np.ones(len(uq)) / len(uq)) * (By - y.mean()) *
                 (Bp - p.mean())).sum()) / den
    CW = float(((cnt / cnt.sum()) * np.array(
        [np.mean((y[g == c] - By[i]) * (p[g == c] - Bp[i]))
         for i, c in enumerate(uq)])).sum()) / den
    slope = (float(np.polyfit(By, Bp, 1)[0]) if len(By) >= 3
             and By.std() > 0 else np.nan)
    return {"overall_r": float(C.pearson(y, p)),
            "overall_rho": float(C.spearman(y, p)),
            "rmse": float(np.sqrt(np.mean((y - p) ** 2))),
            "mae": float(np.mean(np.abs(y - p))),
            "C_B_row": CBr, "C_B_equal": CBe, "C_W": CW,
            "atten_slope": slope,
            "macro_mean_pearson": float(C.pearson(By, Bp)),
            "macro_spearman": FM.macro_spearman(y, p, g, min_n=MIN_PS_N),
            "macro_kendall": _kendall(y, p, g),
            "eff10": FM.eff10(y, p, g, n_rand=20)}


def _kendall(y, p, g, min_n=2):
    vals = []
    for c in np.unique(g):
        m = np.where(g == c)[0]
        if len(m) < min_n:
            continue
        a, b = y[m], p[m]; n = len(a)
        conc = disc = ta = tb = 0
        for i in range(n):
            for j in range(i + 1, n):
                da = a[i] - a[j]; db = b[i] - b[j]
                if da == 0 and db == 0:
                    ta += 1; tb += 1
                elif da == 0:
                    ta += 1
                elif db == 0:
                    tb += 1
                elif (da > 0) == (db > 0):
                    conc += 1
                else:
                    disc += 1
        n0 = conc + disc
        den = np.sqrt((n0 + ta) * (n0 + tb))
        if den > 0:
            vals.append((conc - disc) / den)
    return float(np.mean(vals)) if vals else float("nan")


if __name__ == "__main__":
    sys.exit(main())
