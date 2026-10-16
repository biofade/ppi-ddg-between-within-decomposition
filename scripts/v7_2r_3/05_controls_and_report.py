# -*- coding: utf-8 -*-
"""scripts/v7_2r_3/05_controls_and_report.py -- V7.2R-3 PHASE 7：负对照（§19/§20）。

在原 U_ABC 上补充以下对照（与主结果同一 learner、同一冻结嵌套协议）：

| 对照 | 构成 | 判定 |
|---|---|---|
| **NC0** | `base` only | 基线 |
| **NC1 WT-only** | base + **强表示的 `r_wt`**（**无差分**） | 是否真的需要扰动？ |
| **NC2 identity** | base + mutation identity | 是否只是替换身份编码？ |
| **NC3 norm-only** | base + **`‖Δr‖`** | 是否只是幅度信息？ |
| **NC4 patch-shuffle** | base + **同 target 内逐行打乱的 Δr** | Δr 是否有排他性信息？（§19 变体） |

NC4 的动机：若把 `Δr` 在**同一 target 内**随机置换后性能不变，
说明 Δr 的**逐突变特异性**没有贡献（只是容量）。
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
V723 = C.ROOT / "experiments" / "v7_2r_3"
PH7 = V723 / "phase7_controls"
PH7.mkdir(parents=True, exist_ok=True)
IMP_B = C.ROOT / "features" / "v7_2r_2" / "trackB_imported"
IMP_P = C.ROOT / "features" / "v7_2r_3" / "ppi_imported"
TREES = 300
PCA_NEW = 32
MIN_PS_N = 5
SEEDS = [42, 43, 44]
L: list[str] = []

CONTROLS = {
    "NC0_base_only":      ["base"],
    "NC1_wt_only":        ["base", "P_wt"],
    "NC2_identity":       ["base", "F0"],
    "NC3_dr_norm":        ["base", "DRnorm"],
    "NC4_dr_shuffled":    ["base", "B_dr_shuf"],
    "ref_B2_dr":          ["base", "B_dr"],
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
    return {"macro_spearman": FM.macro_spearman(y, p, g, min_n=MIN_PS_N),
            "eff10": FM.eff10(y, p, g, n_rand=20),
            "C_W": float((cnt / cnt.sum() * np.array(
                [np.mean((y[g == c] - By[i]) * (p[g == c] - Bp[i]))
                 for i, c in enumerate(uq)])).sum()) / den,
            "overall_r": float(C.pearson(y, p)),
            "rmse": float(np.sqrt(np.mean((y - p) ** 2)))}


def main() -> int:
    t0 = time.time()
    emit("# V7.2R-3 PHASE 7 — 负对照（§19/§20）")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    pred_p = V723 / "phase5_eval" / "ppi_predictions.npz"
    if not pred_p.exists():
        emit("## 状态：**WAITING_FOR_V7_2R_3_EVAL**")
        emit("请先运行 `03_matched_eval.py`。")
        (PH7 / "negative_controls.md").write_text("\n".join(L) + "\n",
                                                  encoding="utf-8")
        return 2
    z = np.load(pred_p, allow_pickle=True)
    rows = z["row_id"].astype(np.int64)
    y = z["y"].astype(float)
    g = np.array([str(x) for x in z["pdb"]])
    fam = np.array([str(x) for x in z["family"]])

    d = C.build_pair_features(verbose=False)
    can = pd.read_csv(C.ROOT / "experiments" / "v6_2" / "phase2_split4" /
                      "matched_analysis_set.csv", low_memory=False)
    can_i = can.set_index("row")
    blk = d["X_blocks"]
    base = np.hstack([blk["phys_diff"], blk["Z_diff"]])[rows].astype(np.float32)

    # 强表示 r_wt / Δr（V7.2R-2 链隔离）与 PPI r_wt
    def load(dirp, keys):
        out = {}
        for p in sorted(dirp.glob("*.npz")):
            try:
                q = np.load(p, allow_pickle=True)
            except Exception:  # noqa: BLE001
                continue
            if not all(k in q.files for k in keys):
                continue
            for i, r in enumerate(q["row_id"].astype(np.int64)):
                out[int(r)] = {k: q[k][i].astype(np.float32) for k in keys}
        return out
    B = load(IMP_B, ["r_wt", "delta_r"])
    P = load(IMP_P, ["r_wt"])
    Bwt = np.stack([B[int(r)]["r_wt"] for r in rows]).astype(np.float32)
    Bdr = np.stack([B[int(r)]["delta_r"] for r in rows]).astype(np.float32)
    Pwt = np.stack([P[int(r)]["r_wt"] for r in rows]).astype(np.float32)
    DRnorm = np.linalg.norm(Bdr, axis=1, keepdims=True).astype(np.float32)
    # NC4：**同 target 内**逐行打乱 Δr（确定性）
    BdrS = np.zeros_like(Bdr)
    rs = np.random.RandomState(42)
    for c in np.unique(g):
        m = np.where(g == c)[0]
        BdrS[m] = Bdr[m][rs.permutation(len(m))]
    # NC2 identity
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
            for mm in re.finditer(r"(?<![A-Z])([A-Z])(\d+)([A-Z])(?![A-Z])", body):
                ms.append((mm.group(1), mm.group(3)))
        if not ms:
            ms = [("X", "X")]
        vs = [v7f.mut_identity(a2, b2) for a2, b2 in ms]
        ks = list(vs[0].keys())
        f0.append([float(np.mean([v[k] for v in vs])) for k in ks])
    F0 = np.array(f0, np.float32)
    RAW = {"base": base, "P_wt": Pwt, "F0": F0, "DRnorm": DRnorm,
           "B_dr": Bdr, "B_dr_shuf": BdrS}
    emit(f"- U_ABC: rows **{len(y)}** / targets **{len(set(g))}** / "
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
                    X = RAW[b]
                    if b == "base":
                        mtr.append(X[tr]); mte.append(X[te]); continue
                    from sklearn.decomposition import PCA
                    mu = np.nanmean(X[tr], axis=0)
                    sv = np.nanstd(X[tr], axis=0); sv[sv == 0] = 1
                    Ztr = (X[tr] - mu) / sv; Zte = (X[te] - mu) / sv
                    kk = min(PCA_NEW, Ztr.shape[1], max(2, len(tr) - 1))
                    pc = PCA(n_components=kk, random_state=sd).fit(Ztr)
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
        mm = {k: float(np.mean([a2[k] for a2 in acc])) for k in acc[0]}
        emit(f"| `{name}` | {'+'.join(blocks)} | {mm['macro_spearman']:+.4f} | "
             f"{mm['eff10']:+.4f} | {mm['C_W']:+.4f} | {mm['overall_r']:+.4f} |")
        res.append({"control": name, "blocks": "+".join(blocks), **mm})
    emit("")
    emit("## 判定要点")
    emit("")
    emit("- **NC1**：若 `NC1_wt_only` ≥ `ref_B2_dr` → **扰动本身不必要**，"
         "信号来自结构身份")
    emit("- **NC2**：若 `NC2_identity` 接近 B2 → 只是**替换身份编码**")
    emit("- **NC3**：若 `NC3_dr_norm` 接近 B2 → 只是**幅度**信息，方向无贡献")
    emit("- **NC4**：若 `NC4_dr_shuffled` 接近 B2 → Δr 的**逐突变特异性**无贡献"
         "（只是容量）")
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
