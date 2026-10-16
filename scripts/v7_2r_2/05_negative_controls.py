# -*- coding: utf-8 -*-
"""scripts/v7_2r_2/05_negative_controls.py -- V7.2R-2 PHASE：负对照 NC1–NC4（§46）.

| 对照 | 内容 |
|---|---|
| **NC1** | `r_wt` only vs `delta_r`（判定是否真需要**扰动**） |
| **NC2** | mutation identity / 理化向量（判定是否只是编码氨基酸替换） |
| **NC3** | `‖Δr‖` only（判定是否只是**幅度**信息） |
| **NC4** | partner shuffle（**仅当 B3 声称 partner-specific 信号时**） |

全部在**同一 U_AB** 上、同一冻结嵌套 family-CV 下评估。
"""

from __future__ import annotations

import gc
import importlib.util
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
PH7 = V722 / "phase7_controls"
PH7.mkdir(parents=True, exist_ok=True)
IMP = C.ROOT / "features" / "v7_2r_2" / "trackB_imported"
PER_A = C.ROOT / "features" / "v7_2" / "per_target"
TREES = 300
PCA_NEW = 32
MIN_PS_N = 5
SEEDS = [42, 43, 44]
L: list[str] = []

CONTROLS = {
    "NC1_B2_dr":        ["base", "B_dr"],
    "NC1_B1_wt_only":   ["base", "B_wt"],
    "NC2_mut_identity": ["base", "F0"],
    "NC3_dr_norm_only": ["base", "B_drnorm"],
    "NC4_random_partner": ["base", "B_dr", "B_randpartner"],
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
    }


def main() -> int:
    t0 = time.time()
    emit("# V7.2R-2 — 负对照（NC1–NC4，§46）")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    pred_p = V722 / "phase5_eval" / "trackB_predictions.npz"
    if not pred_p.exists():
        emit("")
        emit("## 状态：**WAITING_FOR_KAGGLE_ZIP**")
        emit("")
        emit("负对照依赖 Track B 预测（`trackB_predictions.npz`）。")
        emit("请先完成 Kaggle 侧 → `02_import_kaggle_features.py` → "
             "`03_matched_trackB_eval.py`。")
        (PH7 / "negative_controls.md").write_text("\n".join(L) + "\n",
                                                  encoding="utf-8")
        return 2
    z = np.load(pred_p, allow_pickle=True)
    rows = z["row_id"].astype(np.int64)
    y = z["y"].astype(float)
    g = np.array([str(x) for x in z["pdb"]])
    fam = np.array([str(x) for x in z["family"]])

    lab = pd.read_csv(V711 / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    d = C.build_pair_features(verbose=False)
    can = pd.read_csv(C.ROOT / "experiments" / "v6_2" / "phase2_split4" /
                      "matched_analysis_set.csv", low_memory=False)
    can_i = can.set_index("row")
    blk = d["X_blocks"]
    base = np.hstack([blk["phys_diff"], blk["Z_diff"]])[rows].astype(np.float32)

    # Track B 表示
    B = {}
    for p in sorted(IMP.glob("*.npz")):
        try:
            q = np.load(p, allow_pickle=True)
            for i, r in enumerate(q["row_id"].astype(np.int64)):
                B[int(r)] = {"wt": q["r_wt"][i].astype(np.float32),
                             "mut": q["r_mut"][i].astype(np.float32),
                             "dr": q["delta_r"][i].astype(np.float32)}
        except Exception:  # noqa: BLE001
            continue
    Bwt = np.stack([B[int(r)]["wt"] for r in rows]).astype(np.float32)
    Bdr = np.stack([B[int(r)]["dr"] for r in rows]).astype(np.float32)
    Bdrnorm = np.linalg.norm(Bdr, axis=1, keepdims=True)
    # NC4：确定性随机 partner（同维）
    rng = np.random.RandomState(42)
    perm = rng.permutation(len(rows))
    Brand = Bdr[perm]
    # NC2：mutation identity
    fwd = pd.read_csv(C.DATA / "v4.0" / "ddg_pairs.csv", low_memory=False)
    fwd = fwd[fwd.direction == "forward"].reset_index(drop=True)
    f0 = []
    for r in rows:
        s = str(fwd["mutation"].iloc[int(r)]).upper(); ms = []
        for tok in re.split(r"[;,]", s):
            tok = tok.strip()
            if not tok:
                continue
            body = tok.split("_", 1)[-1] if "_" in tok else tok
            for m in re.finditer(r"(?<![A-Z])([A-Z])(\d+)([A-Z])(?![A-Z])", body):
                ms.append((m.group(1), m.group(3)))
        if not ms:
            ms = [("X", "X")]
        vs = [v7f.mut_identity(a, b) for a, b in ms]
        ks = list(vs[0].keys())
        f0.append([float(np.mean([v[k] for v in vs])) for k in ks])
    F0 = np.array(f0, np.float32)
    FIT = {"base": base, "B_dr": Bdr, "B_wt": Bwt, "F0": F0,
           "B_drnorm": Bdrnorm, "B_randpartner": Brand}

    emit(f"- U_AB: rows **{len(y)}** / targets **{len(set(g))}** / "
         f"families **{len(set(fam))}**")
    emit("")
    emit("| 对照 | 构成 | Macro ρ | Eff@10 | C_W | Overall r |")
    emit("|---|---|---:|---:|---:|---:|")
    D4 = {"y": y, "pdb": g, "fam": fam,
          "seq40": can_i.loc[rows, "seq40_cluster"].astype(str).to_numpy(),
          "obs": can_i.loc[rows, "obs_key"].astype(str).to_numpy(),
          "state": can_i.loc[rows, "state_id"].astype(str).to_numpy()}
    res = []
    for name, blocks in CONTROLS.items():
        acc = []
        for sd in SEEDS:
            outer = s4.make_folds(D4["fam"], sd, 5)
            yhat = np.full(len(y), np.nan)
            for tr, te in outer:
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
            acc.append(mset(y[ok], yhat[ok], g[ok]))
        mm = {k: float(np.mean([a[k] for a in acc])) for k in acc[0]}
        emit(f"| `{name}` | {'+'.join(blocks)} | {mm['macro_spearman']:+.4f} | "
             f"{mm['eff10']:+.4f} | {mm['C_W']:+.4f} | {mm['overall_r']:+.4f} |")
        res.append({"control": name, "blocks": "+".join(blocks), **mm})
    emit("")
    emit("## 判定要点（§46）")
    emit("")
    emit("- **NC1**：若 `B1_wt_only`（仅 `r_wt`）已与 `B2_dr` 相当 → "
         "**扰动本身不必要**，增益来自结构身份")
    emit("- **NC2**：若 `mut_identity` 与 B2 相当 → 增益只是**氨基酸替换编码**")
    emit("- **NC3**：若 `dr_norm_only` 已接近 B2 → 增益只是**幅度**信息，"
         "方向/内容无贡献")
    emit("- **NC4**：若 `random_partner` 与 B3 相当 → **不得**声称 "
         "partner-specific 信号")
    emit("")
    pd.DataFrame(res).to_csv(PH7 / "negative_controls.csv", index=False,
                             encoding="utf-8")
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (PH7 / "negative_controls.md").write_text("\n".join(L) + "\n",
                                              encoding="utf-8")
    print(f"\n写出 {PH7 / 'negative_controls.csv'} 与 negative_controls.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
