# -*- coding: utf-8 -*-
"""scripts/v7_2r_3/04_bootstrap.py -- V7.2R-3：B3a/B3b/B4/B_shuf − B2 的 bootstrap.

**核心对比（§17）**：`B3/B4 − B2`（在控制强通用表示之后，真伙伴 PPI 上下文
是否增加靶点内排序信号）。
`B4 − B3a`（§18）：突变诱导的 PPI 变化是否超越静态复合物身份。
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

V723 = C.ROOT / "experiments" / "v7_2r_3"
PH5, PH6 = V723 / "phase5_eval", V723 / "phase6_bootstrap"
PH6.mkdir(parents=True, exist_ok=True)
NPZ = PH5 / "ppi_predictions.npz"
N_BOOT = 1500
MIN_PS_N = 5
SEEDS = [42, 43, 44]
CONTRASTS = [("B3b_dr_x_PPI", "B2_chain_isolated", "B3b-B2  Δr×PPI 交互 vs 链隔离（核心）"),
             ("B3a_static_PPI", "B2_chain_isolated", "B3a-B2  静态 PPI vs 链隔离"),
             ("B4_PPI_perturb", "B2_chain_isolated", "B4-B2   PPI 扰动 vs 链隔离"),
             ("B4_PPI_perturb", "B3a_static_PPI", "B4-B3a  PPI 扰动 vs 静态 PPI"),
             ("B_shuf_PPI", "B3a_static_PPI", "Shuf-真实  打乱 vs 真实 PPI")]
L: list[str] = []


def emit(m=""):
    print(m); L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def allm(y, p, g):
    uq = np.unique(g)
    By = np.array([y[g == c].mean() for c in uq])
    Bp = np.array([p[g == c].mean() for c in uq])
    cnt = np.array([(g == c).sum() for c in uq], float)
    den = float(y.std()) * float(p.std())
    if den <= 0:
        den = np.nan
    cw = float((cnt / cnt.sum() * np.array(
        [np.mean((y[g == c] - By[i]) * (p[g == c] - Bp[i]))
         for i, c in enumerate(uq)])).sum()) / den
    ken = np.nan
    pw = np.nan
    return (FM.macro_spearman(y, p, g, min_n=MIN_PS_N),
            FM.eff10(y, p, g, n_rand=20), cw, float(C.pearson(y, p)),
            float(np.sqrt(np.mean((y - p) ** 2))), ken, pw)


NAMES = ["macro_spearman", "eff10", "C_W", "overall_r", "rmse"]


def main() -> int:
    t0 = time.time()
    emit("# V7.2R-3 PHASE 6 — paired target-cluster bootstrap")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    if not NPZ.exists():
        emit("")
        emit("## 状态：**WAITING_FOR_V7_2R_3_ZIP**")
        emit("")
        emit("请先运行 `02_import_ppi_features.py` 与 `03_matched_eval.py`。")
        (PH6 / "bootstrap_report.md").write_text("\n".join(L) + "\n",
                                                 encoding="utf-8")
        return 2
    emit(f"- bootstrap **{N_BOOT}/{N_BOOT}**（cluster = canonical target）")
    emit("")
    z = np.load(NPZ, allow_pickle=True)
    y = z["y"].astype(float); g = np.array([str(x) for x in z["pdb"]])

    def avg(b):
        ks = [f"{b}_s{s}" for s in SEEDS if f"{b}_s{s}" in z.files]
        return np.mean([z[k].astype(float) for k in ks], axis=0) if ks else None
    uq = np.unique(g)
    cidx = [np.where(g == c)[0] for c in uq]
    emit(f"- rows **{len(y)}** / targets **{len(set(g))}**")
    emit(f"- n>={MIN_PS_N} target: "
         f"**{sum(1 for c in uq if (g == c).sum() >= MIN_PS_N)}**")
    emit("")
    out = []
    for a, b, lab in CONTRASTS:
        PA, PB = avg(a), avg(b)
        if PA is None or PB is None:
            emit(f"### {lab} → **NOT_RUN**（缺预测）"); emit(""); continue
        m = np.isfinite(PA) & np.isfinite(PB)
        yy, pa, pb, gg = y[m], PA[m], PB[m], g[m]
        A0, B0 = allm(yy, pa, gg), allm(yy, pb, gg)
        pt = {n: A0[i] - B0[i] for i, n in enumerate(NAMES)}
        acc = {n: [] for n in NAMES}
        rs = np.random.RandomState(42)
        for it in range(N_BOOT):
            pick = rs.randint(0, len(uq), len(uq))
            r_ = np.concatenate([cidx[j] for j in pick])
            rep = np.concatenate([np.full(len(cidx[j]), i)
                                  for i, j in enumerate(pick)])
            try:
                A = allm(y[r_], pa[r_], rep); Bx = allm(y[r_], pb[r_], rep)
            except Exception:  # noqa: BLE001
                continue
            for i, n in enumerate(NAMES):
                v = A[i] - Bx[i]
                if v == v:
                    acc[n].append(v)
        emit("")
        emit(f"### {lab}")
        emit("")
        emit("| 指标 | point | 95% CI | 有效 | 排除 0? |")
        emit("|---|---:|---|---:|---|")
        rr = []
        for n in NAMES:
            arr = np.asarray(acc[n], float)
            if len(arr) < 20:
                lo = hi = np.nan; ex = False
            else:
                lo, hi = np.percentile(arr, [2.5, 97.5])
                ex = bool((lo > 0) or (hi < 0))
            emit(f"| Δ{n} | {pt[n]:+.4f} | [{lo:+.4f}, {hi:+.4f}] | "
                 f"{len(arr)} | {'✅ 排除 0' if ex else '❌ 含 0'} |")
            rr.append({"contrast": lab, "metric": f"d_{n}", "point": pt[n],
                       "ci_low": float(lo), "ci_high": float(hi),
                       "n_valid": len(arr), "ci_excludes_zero": ex})
        dm = {r["metric"]: r for r in rr}
        sup = (dm["d_macro_spearman"]["point"] > 0
               and dm["d_macro_spearman"]["ci_excludes_zero"]
               and (dm["d_eff10"]["point"] > 0 or dm["d_C_W"]["point"] > 0)
               and dm["d_rmse"]["ci_high"] < 0.5)
        topk = (not dm["d_macro_spearman"]["ci_excludes_zero"]
                and dm["d_eff10"]["point"] > 0
                and dm["d_eff10"]["ci_excludes_zero"])
        emit("")
        emit(f"- **§23 成功判据：{'SUPPORTED ✅' if sup else 'NOT SUPPORTED ❌'}**")
        if topk:
            emit("- **§47 `TOPK_ONLY_SIGNAL`**")
        emit("")
        out += rr
        pd.DataFrame(out).to_csv(PH6 / "contrasts.csv", index=False,
                                 encoding="utf-8")
        emit(f"  ({(time.time() - t0) / 60:.1f} min)")
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (PH6 / "bootstrap_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {PH6 / 'contrasts.csv'} 与 bootstrap_report.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
