# -*- coding: utf-8 -*-
"""scripts/v7_2r_2/04_bootstrap.py -- V7.2R-2：B2−A2 / B2−B0 / B1−B0 / B3−B2.

**5000** paired target-cluster bootstrap（cluster = canonical target）。
若 Track B 预测尚未生成 → 输出 `WAITING_FOR_KAGGLE_ZIP` 并停止。
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v6_1"))
import v5_common as C  # noqa: E402
import fastmetrics as FM  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

V722 = C.ROOT / "experiments" / "v7_2r_2"
PH5 = V722 / "phase5_eval"
PH6 = V722 / "phase6_bootstrap"
PH6.mkdir(parents=True, exist_ok=True)
NPZ = PH5 / "trackB_predictions.npz"
N_BOOT = 5000
MIN_PS_N = 5
SEEDS = [42, 43, 44]
L: list[str] = []

CONTRASTS = [("B2_strong_dr", "A2_smallPLM_dh", "B2−A2  强扰动 vs 小PLM扰动（核心）"),
             ("B2_strong_dr", "B0_Base", "B2−B0  强扰动 vs 基线"),
             ("B1_strong_repr", "B0_Base", "B1−B0  强表示(无差分) vs 基线"),
             ("B3_strong_dr_partner", "B2_strong_dr", "B3−B2  partner 增量")]


def emit(m=""):
    print(m); L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def met(y, p, g):
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
    }


def main() -> int:
    t0 = time.time()
    emit("# V7.2R-2 PHASE 6 — 5000 paired target-cluster bootstrap")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    if not NPZ.exists():
        emit("")
        emit("## 状态：**WAITING_FOR_KAGGLE_ZIP**")
        emit("")
        emit(f"未找到 `{NPZ}`。请先完成 Kaggle 侧并运行：")
        emit("```")
        emit("python scripts/v7_2r_2/02_import_kaggle_features.py")
        emit("python scripts/v7_2r_2/03_matched_trackB_eval.py")
        emit("```")
        (PH6 / "bootstrap_report.md").write_text("\n".join(L) + "\n",
                                                 encoding="utf-8")
        return 2
    emit(f"- bootstrap **{N_BOOT}/{N_BOOT}**（cluster = canonical target）")
    emit("")
    z = np.load(NPZ, allow_pickle=True)
    y = z["y"].astype(float); g = np.array([str(x) for x in z["pdb"]])

    def avg(base):
        ks = [f"{base}_s{s}" for s in SEEDS if f"{base}_s{s}" in z.files]
        return np.mean([z[k].astype(float) for k in ks], axis=0) if ks else None

    emit(f"- rows **{len(y)}** / targets **{len(set(g))}**（U_AB）")
    n5 = sum(1 for c in np.unique(g) if (g == c).sum() >= MIN_PS_N)
    emit(f"- n≥{MIN_PS_N} target: **{n5}**")
    emit("")
    uq = np.unique(g)
    cidx = [np.where(g == c)[0] for c in uq]
    rows_out = []
    for a, b, lab in CONTRASTS:
        PA, PB = avg(a), avg(b)
        if PA is None or PB is None:
            emit(f"### {lab} → **NOT_RUN**（缺预测）"); emit(""); continue
        m = np.isfinite(PA) & np.isfinite(PB)
        yy, pa, pb, gg = y[m], PA[m], PB[m], g[m]
        rs = np.random.RandomState(42)
        pt = {k: met(yy, pa, gg)[k] - met(yy, pb, gg)[k] for k in
              ("macro_spearman", "eff10", "C_W", "overall_r", "rmse")}
        acc = {k: [] for k in pt}
        for _ in range(N_BOOT):
            pick = rs.randint(0, len(uq), len(uq))
            r_ = np.concatenate([cidx[j] for j in pick])
            rep = np.concatenate([np.full(len(cidx[j]), i)
                                  for i, j in enumerate(pick)])
            try:
                ma = met(y[r_], pa[r_], rep); mb = met(y[r_], pb[r_], rep)
            except Exception:  # noqa: BLE001
                continue
            for k in pt:
                v = ma[k] - mb[k]
                if v == v:
                    acc[k].append(v)
        emit(f"### {lab}")
        emit("")
        emit("| 指标 | point | 95% CI | 有效 | 排除 0? |")
        emit("|---|---:|---|---:|---|")
        rr = []
        for k in pt:
            arr = np.asarray(acc[k], float)
            if len(arr) < 20:
                lo = hi = np.nan; ex = False
            else:
                lo, hi = np.percentile(arr, [2.5, 97.5])
                ex = bool((lo > 0) or (hi < 0))
            emit(f"| Δ{k} | {pt[k]:+.4f} | [{lo:+.4f}, {hi:+.4f}] | "
                 f"{len(arr)} | {'✅' if ex else '❌ 含 0'} |")
            rr.append({"contrast": lab, "metric": f"d_{k}", "point": pt[k],
                       "ci_low": float(lo), "ci_high": float(hi),
                       "n_valid": len(arr), "ci_excludes_zero": ex})
        dmap = {r["metric"]: r for r in rr}
        sup = (dmap["d_macro_spearman"]["point"] > 0
               and dmap["d_macro_spearman"]["ci_excludes_zero"]
               and (dmap["d_eff10"]["point"] > 0 or dmap["d_C_W"]["point"] > 0)
               and dmap["d_rmse"]["ci_high"] < 0.5)
        # §47 TOPK_ONLY_SIGNAL
        topk_only = (not dmap["d_macro_spearman"]["ci_excludes_zero"]
                     and dmap["d_eff10"]["point"] > 0
                     and dmap["d_eff10"]["ci_excludes_zero"])
        emit("")
        emit(f"- **§23 成功判据：{'SUPPORTED ✅' if sup else 'NOT SUPPORTED ❌'}**")
        if topk_only:
            emit(f"- **§47 `TOPK_ONLY_SIGNAL`**（Macro ρ CI 含 0 但 Eff@10 显著 > 0）")
        emit("")
        rows_out += rr
        pd.DataFrame(rows_out).to_csv(PH6 / "contrasts.csv", index=False,
                                      encoding="utf-8")
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (PH6 / "bootstrap_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {PH6 / 'contrasts.csv'} 与 bootstrap_report.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
