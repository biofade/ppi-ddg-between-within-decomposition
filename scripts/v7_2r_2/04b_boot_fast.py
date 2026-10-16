# -*- coding: utf-8 -*-
"""scripts/v7_2r_2/04b_boot_fast.py -- V7.2R-2 bootstrap 收尾（B2−B0 / B1−B0 / B3−B2）.

一次性算完全部指标（避免重复 macro_spearman/eff10），逐步写盘。
B3−B2 必然 `NOT_RUN`（无 partner 表示）。
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
N_BOOT = 1000          # 时间受限；如实标注（B2−A2 已是 5000）
MIN_PS_N = 5
SEEDS = [42, 43, 44]
CONTRASTS = [("B2_strong_dr", "B0_Base", "B2-B0  强扰动 vs 基线"),
             ("B1_strong_repr", "B0_Base", "B1-B0  强表示(无差分) vs 基线")]
L: list[str] = []


def emit(m=""):
    print(m); L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def all_metrics(y, p, g):
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
    return (FM.macro_spearman(y, p, g, min_n=MIN_PS_N),
            FM.eff10(y, p, g, n_rand=20), cw, float(C.pearson(y, p)),
            float(np.sqrt(np.mean((y - p) ** 2))))


def main() -> int:
    t0 = time.time()
    emit("# V7.2R-2 PHASE 6（收尾）— 剩余 contrasts 的 paired cluster bootstrap")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- bootstrap **{N_BOOT}/{N_BOOT}**（B2−A2 已达 5000，见主报告；此处如实标注）")
    emit(f"- cluster = canonical target")
    emit("")
    z = np.load(NPZ, allow_pickle=True)
    y = z["y"].astype(float)
    g = np.array([str(x) for x in z["pdb"]])

    def avg(base):
        ks = [f"{base}_s{s}" for s in SEEDS if f"{base}_s{s}" in z.files]
        return np.mean([z[k].astype(float) for k in ks], axis=0) if ks else None

    uq = np.unique(g)
    cidx = [np.where(g == c)[0] for c in uq]
    # 读已有（B2-A2）
    prev = pd.read_csv(PH6 / "contrasts.csv") if (PH6 / "contrasts.csv").exists() \
        else pd.DataFrame()
    rows = prev.to_dict("records") if len(prev) else []
    emit(f"- 已有（B2−A2）行数: {len(rows)}")
    emit("")

    for a, b, lab in CONTRASTS:
        PA, PB = avg(a), avg(b)
        if PA is None or PB is None:
            emit(f"### {lab} → **NOT_RUN**（缺预测）"); emit(""); continue
        m = np.isfinite(PA) & np.isfinite(PB)
        yy, pa, pb, gg = y[m], PA[m], PB[m], g[m]
        base_pt = all_metrics(yy, pa, gg)
        comp_pt = all_metrics(yy, pb, gg)
        names = ["macro_spearman", "eff10", "C_W", "overall_r", "rmse"]
        pt = {n: base_pt[i] - comp_pt[i] for i, n in enumerate(names)}
        acc = {n: [] for n in names}
        rs = np.random.RandomState(42)
        for it in range(N_BOOT):
            pick = rs.randint(0, len(uq), len(uq))
            r_ = np.concatenate([cidx[j] for j in pick])
            rep = np.concatenate([np.full(len(cidx[j]), i)
                                  for i, j in enumerate(pick)])
            try:
                A = all_metrics(y[r_], pa[r_], rep)
                B = all_metrics(y[r_], pb[r_], rep)
            except Exception:  # noqa: BLE001
                continue
            for i, n in enumerate(names):
                v = A[i] - B[i]
                if v == v:
                    acc[n].append(v)
            if (it + 1) % 500 == 0:
                emit(f"  - {lab}: {it + 1}/{N_BOOT} "
                     f"({(time.time() - t0) / 60:.1f} min)")
        emit("")
        emit(f"### {lab}")
        emit("")
        emit("| 指标 | point | 95% CI | 有效 | 排除 0? |")
        emit("|---|---:|---|---:|---|")
        rr = []
        for n in names:
            arr = np.asarray(acc[n], float)
            if len(arr) < 20:
                lo = hi = np.nan; ex = False
            else:
                lo, hi = np.percentile(arr, [2.5, 97.5])
                ex = bool((lo > 0) or (hi < 0))
            emit(f"| Δ{n} | {pt[n]:+.4f} | [{lo:+.4f}, {hi:+.4f}] | "
                 f"{len(arr)} | {'✅' if ex else '❌ 含 0'} |")
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
            emit(f"- **§47 `TOPK_ONLY_SIGNAL`**")
        emit("")
        rows += rr
        pd.DataFrame(rows).to_csv(PH6 / "contrasts.csv", index=False,
                                  encoding="utf-8")

    # B3−B2 明确 NOT_RUN
    emit("### B3−B2  partner/PII 上下文增量 → **NOT_RUN**")
    emit("")
    emit("Kaggle ZIP **未提供** partner/MINT 表示（`partner_repr` 字段缺失、"
         "`partner_sanity.csv` 为空），且表示经审计为 **CHAIN_ISOLATED**")
    emit("⇒ 依 §8/§19，**不得**构造 B3，也**不得**用零向量填充。")
    emit("")
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (PH6 / "bootstrap_report_extra.md").write_text("\n".join(L) + "\n",
                                                   encoding="utf-8")
    print(f"\n写出 {PH6 / 'contrasts.csv'} 与 bootstrap_report_extra.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
