# -*- coding: utf-8 -*-
"""scripts/v7_2r/02_bootstrap.py -- V7.2R PHASE 6: 关键 contrast 的配对 cluster bootstrap.

cluster = **canonical target**；5000/5000（§22）
对比（§17）：
    A2-B0  (R3_Base_dh   - R0_Base)
    A1-B0  (R2_Base_Residue - R0_Base)
    R1-B0  (R1_Base_F0   - R0_Base)   身份对照
指标：ΔMacroSpearman / ΔEff@10 / ΔC_W / ΔOverall_r / ΔRMSE

成功判据（§23，**不得事后改写**）：
    ΔMacroSpearman > 0 且 95% CI 排除 0
    且（ΔEff@10 > 0 或 ΔC_W > 0）
    且无灾难性 RMSE 退化
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

V72 = C.ROOT / "experiments" / "v7_2"
PH5 = V72 / "phase5_nested_cv"
PH6 = V72 / "phase7_bootstrap"
PH6.mkdir(parents=True, exist_ok=True)
NPZ = PH5 / "ranker_predictions.npz"
N_BOOT = 1500
MIN_PS_N = 5
SEEDS = [42, 43, 44]
CONTRASTS = [("R3_Base_dh", "R0_Base", "A2-B0  残基扰动增益"),
             ("R2_Base_Residue", "R0_Base", "A1-B0  残基表示增益"),
             ("R1_Base_F0", "R0_Base", "R1-B0  mutation identity 对照")]
L: list[str] = []


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
    # 基线/候选必须用**同一 target 集合**：Macro ρ 要求在双方都有 n>=MIN 的 target 上
    return {
        "macro_spearman": FM.macro_spearman(y, p, g, min_n=MIN_PS_N),
        "eff10": FM.eff10(y, p, g, n_rand=10),
        "C_W": float((cnt / cnt.sum() * np.array(
            [np.mean((y[g == c] - By[i]) * (p[g == c] - Bp[i]))
             for i, c in enumerate(uq)])).sum()) / den,
        "overall_r": float(C.pearson(y, p)),
        "rmse": float(np.sqrt(np.mean((y - p) ** 2))),
    }


def boot(y, pa, pb, g, n=N_BOOT, seed=42):
    uq = np.unique(g)
    cidx = [np.where(g == c)[0] for c in uq]
    rs = np.random.RandomState(seed)
    pt = {k: met(y, pa, g)[k] - met(y, pb, g)[k] for k in
          ("macro_spearman", "eff10", "C_W", "overall_r", "rmse")}
    acc = {k: [] for k in pt}
    for _ in range(n):
        pick = rs.randint(0, len(uq), len(uq))
        r_ = np.concatenate([cidx[j] for j in pick])
        rep = np.concatenate([np.full(len(cidx[j]), i)
                              for i, j in enumerate(pick)])
        try:
            ma = met(y[r_], pa[r_], rep)
            mb = met(y[r_], pb[r_], rep)
        except Exception:  # noqa: BLE001
            continue
        for k in pt:
            v = ma[k] - mb[k]
            if v == v:
                acc[k].append(v)
    out = []
    for k in pt:
        a = np.asarray(acc[k], float)
        if len(a) < 20:
            out.append({"metric": f"d_{k}", "point": pt[k], "ci_low": np.nan,
                        "ci_high": np.nan, "n_valid": len(a),
                        "ci_excludes_zero": False})
            continue
        lo, hi = np.percentile(a, [2.5, 97.5])
        out.append({"metric": f"d_{k}", "point": float(pt[k]),
                    "ci_low": float(lo), "ci_high": float(hi),
                    "n_valid": len(a),
                    "ci_excludes_zero": bool((lo > 0) or (hi < 0))})
    return out


def main() -> int:
    t0 = time.time()
    emit("# V7.2R PHASE 6 — 配对 cluster bootstrap（cluster = canonical target）")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- bootstrap **{N_BOOT}/{N_BOOT}**（§22）")
    emit(f"- 成功判据（§23，**冻结，不得事后改写**）：")
    emit(f"  `ΔMacroSpearman > 0` 且 95% CI 排除 0，")
    emit(f"  且（`ΔEff@10 > 0` 或 `ΔC_W > 0`），且无灾难性 RMSE 退化")
    emit("")
    z = np.load(NPZ, allow_pickle=True)
    y = z["y"].astype(float)
    g = np.array([str(x) for x in z["pdb"]])
    # 用 3 seed 的平均预测作为主结果
    def avg(base):
        ks = [f"{base}_s{s}" for s in SEEDS if f"{base}_s{s}" in z.files]
        if not ks:
            return None
        return np.mean([z[k].astype(float) for k in ks], axis=0)
    emit(f"- rows **{len(y)}** / targets **{len(set(g))}**")
    n5 = sum(1 for c in np.unique(g) if (g == c).sum() >= MIN_PS_N)
    emit(f"- n≥{MIN_PS_N} 的 target（Macro ρ 资格）: **{n5}**")
    emit("")
    if n5 < 5:
        emit(f"!! 合格 target 仅 {n5} 个，bootstrap 无统计意义 → **NOT_RUN**")
        (PH6 / "bootstrap_report.md").write_text("\n".join(L) + "\n",
                                                 encoding="utf-8")
        return 1
    rows = []
    for a, b, lab in CONTRASTS:
        pa, pb = avg(a), avg(b)
        if pa is None or pb is None:
            emit(f"### {lab} → **NOT_RUN**（缺预测）")
            emit("")
            continue
        m = np.isfinite(y) & np.isfinite(pa) & np.isfinite(pb)
        emit(f"### {lab}")
        emit("")
        emit("| 指标 | point | 95% CI | 有效 | 排除 0? |")
        emit("|---|---:|---|---:|---|")
        res = boot(y[m], pa[m], pb[m], g[m])
        for r in res:
            emit(f"| Δ{r['metric'][2:]} | {r['point']:+.4f} | "
                 f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}] | "
                 f"{r['n_valid']} | {'✅' if r['ci_excludes_zero'] else '❌ 含 0'} |")
            rows.append({"contrast": lab, **r})
        emit("")
        pd.DataFrame(rows).to_csv(PH6 / "contrasts.csv", index=False,
                                  encoding="utf-8")
        # 成功判据
        d = {r["metric"]: r for r in res}
        ok = (d["d_macro_spearman"]["point"] > 0
              and d["d_macro_spearman"]["ci_excludes_zero"]
              and (d["d_eff10"]["point"] > 0 or d["d_C_W"]["point"] > 0)
              and d["d_rmse"]["ci_high"] < 0.5)
        emit(f"- **§23 成功判据：{'SUPPORTED ✅' if ok else 'NOT SUPPORTED ❌'}**")
        emit("")
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (PH6 / "bootstrap_report.md").write_text("\n".join(L) + "\n",
                                             encoding="utf-8")
    print(f"\n写出 {PH6/'contrasts.csv'} 与 bootstrap_report.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
