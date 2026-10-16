# -*- coding: utf-8 -*-
"""scripts/v7_2r_2/03_matched_trackB_eval.py -- V7.2R-2 PHASE 5-8.

若 Kaggle Track B 特征**尚未导入** → 明确输出 `WAITING_FOR_KAGGLE_ZIP`
并**停止**（§54 末条：本地科学评估必须等该 ZIP）。

导入后：
1. 构造 **U_AB = Track A valid ∩ Track B valid**（§40）
2. 在 U_AB 上跑 **B0 / A2 / B1 / B2 / B3** 的冻结嵌套 family-CV
   （learner 统一 LightGBM；PCA=32 **train-only**）
3. 主口径 ensemble prediction，敏感性 seed-mean
"""

from __future__ import annotations

import argparse
import gc
import importlib.util
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v6_1"))
import v5_common as C  # noqa: E402
import fastmetrics as FM  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "s4", C.ROOT / "scripts" / "v6_2" / "02_split4_ddg.py")
s4 = importlib.util.module_from_spec(_spec)
sys.modules["s4"] = s4
_spec.loader.exec_module(s4)
_spec2 = importlib.util.spec_from_file_location(
    "v7f", C.ROOT / "scripts" / "v7_0" / "01_features.py")
v7f = importlib.util.module_from_spec(_spec2)
sys.modules["v7f"] = v7f
_spec2.loader.exec_module(v7f)

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

V711 = C.ROOT / "experiments" / "v7_1_1"
V722 = C.ROOT / "experiments" / "v7_2r_2"
PH5 = V722 / "phase5_eval"
PH5.mkdir(parents=True, exist_ok=True)
IMP = C.ROOT / "features" / "v7_2r_2" / "trackB_imported"
PER_A = C.ROOT / "features" / "v7_2" / "per_target"

TREES = 300
PCA_NEW = 32
MIN_PS_N = 5
EPS = 1e-6
SEEDS = [42, 43, 44]
L: list[str] = []

# §8 模型定义。
# ⚠️ B3/B4 需要 **伙伴/MINT 表示**；本 ZIP **未提供**（见 §3 语义审计：
#    `partner_repr` 字段不存在，`partner_sanity.csv` 为空，
#    且表示是 CHAIN_ISOLATED）→ **NOT_RUN**，绝不用零填充伪造。
RANKERS = {
    "B0_Base": ["base"],
    "A2_smallPLM_dh": ["base", "F2_small"],
    "B1_strong_repr": ["base", "B_wt", "B_mut"],
    "B2_strong_dr": ["base", "B_dr"],
}
RANKERS_NOT_RUN = {
    "B3_strong_dr_partner": "无 partner/MINT 表示（partner_repr 缺失，partner_sanity 为空）",
    "B4_interaction": "无伙伴向量 p，无法构造 Δr⊙p（§8 明确禁止在维度不兼容时强行做 elementwise）",
}


def emit(m=""):
    print(m); L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def mset(y, p, g):
    uq = np.unique(g)
    By = np.array([y[g == c].mean() for c in uq])
    Bp = np.array([p[g == c].mean() for c in uq])
    cnt = np.array([(g == c).sum() for c in uq], float)
    den = float(y.std()) * float(p.std())
    if den <= 0:
        den = np.nan
    return {
        "macro_spearman": FM.macro_spearman(y, p, g, min_n=MIN_PS_N),
        "eff10": FM.eff10(y, p, g, n_rand=10),
        "C_W": float((cnt / cnt.sum() * np.array(
            [np.mean((y[g == c] - By[i]) * (p[g == c] - Bp[i]))
             for i, c in enumerate(uq)])).sum()) / den,
        "overall_r": float(C.pearson(y, p)),
        "rmse": float(np.sqrt(np.mean((y - p) ** 2))),
        # §13 要求的附加排名指标
        "kendall": float(np.mean([C.kendall(y[g == c], p[g == c])
                                  for c in uq if (g == c).sum() >= 2]))
        if hasattr(C, "kendall") else _kendall_macro(y, p, g),
        "pairwise05": _pairwise(y, p, g, 0.5),
        "hit10": _hit10(y, p, g),
        "prec10": _prec10(y, p, g),
    }


def _kendall_macro(y, p, g):
    from scipy.stats import kendalltau
    vals = []
    for c in np.unique(g):
        m = g == c
        if m.sum() < 2:
            continue
        t = kendalltau(y[m], p[m]).correlation
        if t == t:
            vals.append(t)
    return float(np.mean(vals)) if vals else np.nan


def _pairwise(y, p, g, delta=0.5):
    """靶点内配对准确率：|Δy| > delta 的配对中，预测方向正确（含并列记 0.5）的比例。"""
    num = den = 0.0
    for c in np.unique(g):
        m = np.where(g == c)[0]
        if len(m) < 2:
            continue
        for i in range(len(m)):
            for j in range(i + 1, len(m)):
                dy = y[m[i]] - y[m[j]]
                if abs(dy) <= delta:
                    continue
                dp = p[m[i]] - p[m[j]]
                den += 1
                if dp == 0:
                    num += 0.5
                elif (dy > 0) == (dp > 0):
                    num += 1
    return float(num / den) if den > 0 else np.nan


def _hit10(y, p, g):
    """每 target 真实 ΔΔG 最小的 k 个中，被预测 top-k 命中的比例（k=min(10,len)）。"""
    vals = []
    for c in np.unique(g):
        m = np.where(g == c)[0]
        k = min(10, len(m))
        if k < 2:
            continue
        true_top = set(np.argsort(y[m])[:k].tolist())
        pred_top = set(np.argsort(p[m])[:k].tolist())
        vals.append(len(true_top & pred_top) / k)
    return float(np.mean(vals)) if vals else np.nan


def _prec10(y, p, g):
    """P@10：预测 top-k 中真实为「最优 k 个」的比例（与 Hit@10 同义口径，k=min(10,len)）。"""
    return _hit10(y, p, g)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-minutes", type=float, default=45.0)
    args = ap.parse_args()
    t0 = time.time()
    emit("# V7.2R-2 PHASE 5 — U_AB 匹配宇宙 + Track B 评估")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    # ---------- Track B 是否已导入 ----------
    bfiles = sorted(IMP.glob("*.npz")) if IMP.exists() else []
    if not bfiles:
        emit("")
        emit("## 状态：**WAITING_FOR_KAGGLE_ZIP**")
        emit("")
        emit(f"`{IMP}` 中无 Track B 特征文件。")
        emit("")
        emit("**本地科学评估依 §54 必须等待 Kaggle ZIP。** 需要完成的步骤：")
        emit("")
        emit("1. 上传 `kaggle/v7_2r_2/kaggle_input_package.zip` → 建 Kaggle Dataset")
        emit("2. 上传 `kaggle/v7_2r_2/V7_2R_2_TrackB_Kaggle.ipynb` → "
             "Accelerator `GPU T4 x2` + `Internet ON` → **Run All**")
        emit("3. 下载 `V7_2R_2_TRACKB_FEATURES.zip`")
        emit("4. 放到 `F:\\dsh\\downloads\\v7_2r_2\\`")
        emit("5. `python scripts/v7_2r_2/02_import_kaggle_features.py`")
        emit("6. 重跑本脚本")
        emit("")
        emit("> 详见 `kaggle/v7_2r_2/README_KAGGLE.md`。")
        emit("")
        emit("### 已就绪的本地准备（可核验）")
        emit("")
        for p in ("kaggle/v7_2r_2/kaggle_input_package.zip",
                  "kaggle/v7_2r_2/V7_2R_2_TrackB_Kaggle.ipynb",
                  "kaggle/v7_2r_2/README_KAGGLE.md",
                  "kaggle/v7_2r_2/kaggle_input_manifest.csv"):
            q = C.ROOT / p
            emit(f"- {'✅' if q.exists() else '❌'} `{p}` "
                 f"({q.stat().st_size:,} B)" if q.exists() else f"- ❌ `{p}`")
        emit("")
        (PH5 / "eval_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
        return 2

    emit(f"- Track B npz: **{len(bfiles)}**")
    emit("")

    # ---------- 载入 ----------
    lab = pd.read_csv(V711 / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    d = C.build_pair_features(verbose=False)
    can = pd.read_csv(C.ROOT / "experiments" / "v6_2" / "phase2_split4" /
                      "matched_analysis_set.csv", low_memory=False)
    can_i = can.set_index("row")
    blk = d["X_blocks"]

    # Track A 行
    a_rows = set()
    for p in sorted(PER_A.glob("*.npz")):
        try:
            a_rows |= set(np.load(p)["row_id"].tolist())
        except Exception:  # noqa: BLE001
            continue
    # Track B 行 + 表示
    B = {}
    b_rows = set()
    for p in bfiles:
        try:
            z = np.load(p, allow_pickle=True)
            rid = z["row_id"].astype(np.int64)
            for i, r in enumerate(rid):
                B[int(r)] = {
                    "wt": z["r_wt"][i].astype(np.float32),
                    "mut": z["r_mut"][i].astype(np.float32),
                    "dr": z["delta_r"][i].astype(np.float32)}
            b_rows |= set(rid.tolist())
        except Exception:  # noqa: BLE001
            continue
    AB = sorted(a_rows & b_rows)
    emit("## 5.1 U_AB（Track A valid ∩ Track B valid）")
    emit("")
    emit("| universe | rows | targets | families |")
    emit("|---|---:|---:|---:|")
    for nm, rs in (("Track A (full)", a_rows), ("Track B (imported)", b_rows),
                   ("**U_AB (matched)**", set(AB))):
        sub = lab[lab.pair_row.isin(rs)]
        emit(f"| {nm} | {len(rs)} | {sub.target_id_canonical.nunique()} | "
             f"{sub.family_id_canonical.nunique()} |")
    subA = lab[lab.pair_row.isin(a_rows)]
    subAB = lab[lab.pair_row.isin(set(AB))]
    n5A = (subA.groupby("target_id_canonical").size() >= MIN_PS_N).sum()
    n5AB = (subAB.groupby("target_id_canonical").size() >= MIN_PS_N).sum()
    rel = n5AB / max(1, n5A)
    emit("")
    emit(f"- n≥{MIN_PS_N} target 覆盖：**{n5AB}/{n5A} = {rel * 100:.1f}%**")
    emit(f"- §40 门槛（≥80% of Track A reference）："
         f"**{'PASS ✅' if rel >= 0.80 else 'COVERAGE_LIMITED ⚠️'}**")
    emit("")
    if len(AB) < 200:
        emit(f"!! U_AB 仅 {len(AB)} 行，样本不足 → 停止")
        (PH5 / "eval_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
        return 1

    # ---------- 特征 ----------
    rows = np.array(AB, int)
    y = d["y"][rows].astype(float)
    g = can_i.loc[rows, "pdb_id"].astype(str).to_numpy()
    fam = can_i.loc[rows, "family_label"].astype(str).to_numpy()
    base = np.hstack([blk["phys_diff"], blk["Z_diff"]])[rows].astype(np.float32)
    # A2: 小 PLM Δh
    HW = np.zeros((len(rows), 320), np.float32)
    HM = np.zeros((len(rows), 320), np.float32)
    pos = {int(r): i for i, r in enumerate(rows)}
    for p in sorted(PER_A.glob("*.npz")):
        z = np.load(p)
        arr_w, arr_m = z["h_wt"], z["h_mut"]
        for i, r in enumerate(z["row_id"]):
            j = pos.get(int(r))
            if j is not None:
                HW[j] = arr_w[i].astype(np.float32)
                HM[j] = arr_m[i].astype(np.float32)
    dH = HM - HW
    F2_small = np.hstack([dH, np.linalg.norm(dH, axis=1, keepdims=True)])
    # B 表示
    Bwt = np.stack([B[r]["wt"] for r in rows]).astype(np.float32)
    Bmut = np.stack([B[r]["mut"] for r in rows]).astype(np.float32)
    Bdr = np.stack([B[r]["dr"] for r in rows]).astype(np.float32)
    FIT = {"base": base, "F2_small": F2_small, "B_wt": Bwt, "B_mut": Bmut,
           "B_dr": Bdr}
    emit(f"- 维度：base={base.shape[1]} F2_small={F2_small.shape[1]} "
         f"B_wt={Bwt.shape[1]} B_dr={Bdr.shape[1]}")
    emit(f"- **U_AB universe: rows {len(y)} / targets {len(set(g))} / "
         f"families {len(set(fam))}**")
    emit("")

    # ---------- 嵌套 CV ----------
    emit("### 5.1b NOT_RUN 模型（如实标注，不伪造）")
    emit("")
    for k, why in RANKERS_NOT_RUN.items():
        emit(f"- **{k}** → NOT_RUN：{why}")
    emit("")
    emit("## 5.2 冻结嵌套 family-CV（B0/A2/B1/B2）")
    emit("")
    D4 = {"y": y, "pdb": g, "fam": fam,
          "seq40": can_i.loc[rows, "seq40_cluster"].astype(str).to_numpy(),
          "obs": can_i.loc[rows, "obs_key"].astype(str).to_numpy(),
          "state": can_i.loc[rows, "state_id"].astype(str).to_numpy()}
    res, preds = [], {}
    for sd in SEEDS:
        if (time.time() - t0) / 60 > args.max_minutes:
            emit(f"- ⏱ 软时限，停止于 seed {sd}")
            break
        outer = s4.make_folds(D4["fam"], sd, 5)
        for mk, blocks in RANKERS.items():
            try:
                yhat = np.full(len(y), np.nan)
                for tr, te in outer:
                    assert not (set(fam[tr]) & set(fam[te])), "outer 泄漏"
                    mtr, mte = [], []
                    for b in blocks:
                        X = FIT[b]
                        if b == "base":
                            mtr.append(X[tr]); mte.append(X[te]); continue
                        from sklearn.decomposition import PCA
                        mu = np.nanmean(X[tr], axis=0)
                        sv = np.nanstd(X[tr], axis=0); sv[sv == 0] = 1
                        Ztr = (X[tr] - mu) / sv; Zte = (X[te] - mu) / sv
                        k = min(PCA_NEW, Ztr.shape[1], max(2, len(tr) - 1))
                        pc = PCA(n_components=k, random_state=sd).fit(Ztr)
                        mtr.append(pc.transform(Ztr).astype(np.float32))
                        mte.append(pc.transform(Zte).astype(np.float32))
                    Xtr = np.hstack(mtr).astype(np.float32)
                    Xte = np.hstack(mte).astype(np.float32)
                    lgb = C.load_lgb()
                    mdl = lgb.LGBMRegressor(**{**C.LGB_PARAMS,
                                               "n_estimators": TREES,
                                               "random_state": sd})
                    mdl.fit(np.vstack([Xtr, -Xtr]), np.r_[y[tr], -y[tr]])
                    yhat[te] = mdl.predict(Xte)
                    del mdl; gc.collect()
                ok = np.isfinite(yhat)
                m = mset(y[ok], yhat[ok], g[ok])
                res.append({"seed": sd, "ranker": mk, "n_rows": int(ok.sum()),
                            "n_target": len(set(g[ok])), **m})
                preds[f"{mk}_s{sd}"] = yhat.astype(np.float32)
                emit(f"| {mk} | s{sd} | Macro ρ={m['macro_spearman']:+.4f} | "
                     f"Eff@10={m['eff10']:+.4f} | C_W={m['C_W']:+.4f} | "
                     f"r={m['overall_r']:+.4f} |")
                pd.DataFrame(res).to_csv(PH5 / "trackB_results.csv", index=False,
                                         encoding="utf-8")
            except Exception as ex:  # noqa: BLE001
                emit(f"| {mk} | s{sd} | FAILED {type(ex).__name__}: {ex} |")
                gc.collect()
    R = pd.DataFrame(res)
    preds["y"] = y.astype(np.float32); preds["pdb"] = g.astype(str)
    preds["family"] = fam.astype(str); preds["row_id"] = rows.astype(np.int64)
    np.savez(PH5 / "trackB_predictions.npz", **preds)
    emit("")
    if len(R):
        emit("### 3-seed 均值")
        emit("")
        cols = ["macro_spearman", "eff10", "C_W", "kendall", "pairwise05",
                "hit10", "overall_r"]
        cols = [c for c in cols if c in R.columns]
        sm = R.groupby("ranker")[cols].mean()
        emit("| Ranker | Macro ρ | Eff@10 | C_W | Kendall | Pairwise δ=.5 | "
             "Hit@10 | Overall r |")
        emit("|---|---:|---:|---:|---:|---:|---:|---:|")
        for k, r in sm.iterrows():
            emit(f"| {k} | {r['macro_spearman']:+.4f} | {r['eff10']:+.4f} | "
                 f"{r['C_W']:+.4f} | {r.get('kendall', float('nan')):+.4f} | "
                 f"{r.get('pairwise05', float('nan')):+.4f} | "
                 f"{r.get('hit10', float('nan')):+.4f} | {r['overall_r']:+.4f} |")
        sm.to_csv(PH5 / "trackB_seedmean.csv", encoding="utf-8")
    emit("")
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (PH5 / "eval_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {PH5 / 'trackB_results.csv'} 与 eval_report.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
