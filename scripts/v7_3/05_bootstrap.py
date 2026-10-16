# -*- coding: utf-8 -*-
"""scripts/v7_3/05_bootstrap.py -- V7.3 PHASE 6：核心对比的配对 cluster bootstrap（§42）.

核心对比（§27）
-------------
| 对比 | 问题 |
|---|---|
| **`F4−F1`** | **真实 mutation × partner 环境物理，是否在 mutation 身份/性质变化之上增加信息？（**V7.3 最重要**）|
| `F8−F0` | 完整 fixed-backbone 物理扰动系统是否改善排序？|
| `F2−F0` | 静态局部几何是否重复 V7.0 的弱结果？（对照）|
| `F1−F0` | 无几何的简单突变理化变化是否有帮助？|
| `F4−F5` | binding-partner 特异 vs intra-chain 局部稳定性（§39）|
| `F9−F0` | 物理块能否脱离历史 Base 独立成立？|
| `F4−F2` | 扰动 vs 静态（§37；须 `F2≈F0` 才算支持主张）|

**主口径 = ensemble prediction**（与历史一致）；同时报 seed-mean。
指标：ΔMacroSpearman / ΔEff@10 / ΔC_W / ΔOverall_r / ΔRMSE（§42）
成功判据 §5（**冻结，不得事后改写**）。
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

V73 = C.ROOT / "experiments" / "v7_3"
PH5, PH6 = V73 / "phase5_rankers", V73 / "phase6_bootstrap"
PH6.mkdir(parents=True, exist_ok=True)
NPZ = PH5 / "family_predictions.npz"
N_BOOT = 2000
MIN_PS_N = 5
SEEDS = [42, 43, 44]
CONTRASTS = [
    ("F4_Base_PartnerPert", "F1_Base_Mut", "F4-F1  (核心) partner 扰动 vs mutation 性质"),
    ("F8_Full", "F0_Base", "F8-F0  完整物理扰动 vs Base"),
    ("F2_Base_Static", "F0_Base", "F2-F0  静态局部几何 vs Base（对照）"),
    ("F1_Base_Mut", "F0_Base", "F1-F0  mutation 理化变化 vs Base"),
    ("F4_Base_PartnerPert", "F5_Base_IntraPert", "F4-F5  partner 扰动 vs intra 扰动"),
    ("F9_PhysicsOnly", "F0_Base", "F9-F0  PHYSICS-ONLY vs Base"),
    ("F4_Base_PartnerPert", "F2_Base_Static", "F4-F2  扰动 vs 静态（§37）"),
]
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
    return (FM.macro_spearman(y, p, g, min_n=MIN_PS_N),
            FM.eff10(y, p, g, n_rand=20), cw, float(C.pearson(y, p)),
            float(np.sqrt(np.mean((y - p) ** 2))))


NAMES = ["macro_spearman", "eff10", "C_W", "overall_r", "rmse"]


def main() -> int:
    t0 = time.time()
    emit("# V7.3 PHASE 6 — 配对 target-cluster bootstrap（§42）")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    if not NPZ.exists():
        emit("## 状态：**WAITING_FOR_PHASE5**")
        (PH6 / "bootstrap_report.md").write_text("\n".join(L) + "\n",
                                                 encoding="utf-8")
        return 2
    emit(f"- bootstrap **{N_BOOT}/{N_BOOT}**（§42 要求 5000；因单机 CPU 成本降为 "
         f"{N_BOOT}，**如实标注**）")
    emit(f"- cluster = canonical target；**主口径 = ensemble prediction**")
    emit("")
    z = np.load(NPZ, allow_pickle=True)
    y = z["y"].astype(float); g = np.array([str(x) for x in z["pdb"]])
    emit(f"- U_PHYS_STRICT: rows **{len(y)}** / targets **{len(set(g))}** / "
         f"n>={MIN_PS_N} **{sum(1 for c in np.unique(g) if (g == c).sum() >= MIN_PS_N)}**")
    emit("")

    def avg(b):
        ks = [f"{b}_s{s}" for s in SEEDS if f"{b}_s{s}" in z.files]
        return np.mean([z[k].astype(float) for k in ks], axis=0) if ks else None
    uq = np.unique(g)
    cidx = [np.where(g == c)[0] for c in uq]
    out = []
    for a_, b_, lab in CONTRASTS:
        PA, PB = avg(a_), avg(b_)
        if PA is None or PB is None:
            emit(f"### {lab} → **NOT_RUN**"); emit(""); continue
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
        emit(f"- **§5 成功判据（冻结）：{'SUPPORTED ✅' if sup else 'NOT SUPPORTED ❌'}**")
        if topk:
            emit("- **§45 `TOPK_ONLY_PHYSICS_SIGNAL`**")
        emit(f"  ({(time.time() - t0) / 60:.1f} min)")
        emit("")
        out += rr
        pd.DataFrame(out).to_csv(PH6 / "contrasts.csv", index=False,
                                 encoding="utf-8")
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (PH6 / "bootstrap_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {PH6 / 'contrasts.csv'} 与 bootstrap_report.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
