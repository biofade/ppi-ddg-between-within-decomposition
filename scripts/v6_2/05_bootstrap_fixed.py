# -*- coding: utf-8 -*-
"""scripts/v6_2/05_bootstrap_fixed.py -- V6.2 PHASE 5: Bootstrap 修正.

V6.1 的 bug（本 PHASE 修复）
---------------------------
V6.1 的 `01_cluster_bootstrap_ci.py` 设了 **`if unique complex < 15: skip`**，
导致 **5000 次里只有 1151 次有效（23%）**，且被跳过的恰是
**抽样重复度高**（更极端）的迭代 → **CI 系统性偏窄**。

修正（V6.2 §5.1）
-----------------
* 5000 次有放回抽样
* **不设 unique 阈值**
* **不跳过任何迭代**（仅跳过 metric 返回 NaN 的退化迭代，并如实计数）

Seeds 处理（V6.2 §5.2）
-----------------------
* **Primary**：先对 3 个 seed 的预测**取平均**，再算 metric
* **Secondary**：逐 seed 的 mean ± SD
* bootstrap 抽样单位 = **complex**

输出
----
* `phase5_bootstrap/bootstrap_comparison.csv`（V6.1 vs V6.2 并排）
* `phase5_bootstrap/primary_statistics_ci.csv`
* `phase5_bootstrap/bootstrap_report.md`
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v6_1"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v6_2"))
import v5_common as C  # noqa: E402
import fastmetrics as FM  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

STD = C.ROOT / "experiments" / "v6_sota_audit" / "standardized_preds.csv"
PH2 = C.ROOT / "experiments" / "v6_2" / "phase2_split4"
PH5 = C.ROOT / "experiments" / "v6_2" / "phase5_bootstrap"
PH5.mkdir(parents=True, exist_ok=True)
N_BOOT = 5000
SEED = 42
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def boot_ci(fn, y, p, g, n_boot=N_BOOT, seed=SEED):
    """**修正版** cluster bootstrap：无 unique 阈值、不跳过迭代。"""
    uniq = np.unique(g)
    idx = {u: np.where(g == u)[0] for u in uniq}
    pt = fn(y, p, g)
    rs = np.random.RandomState(seed)
    vals, nan_cnt = [], 0
    for b in range(n_boot):
        pick = rs.randint(0, len(uniq), len(uniq))
        rows = np.concatenate([idx[uniq[j]] for j in pick])
        rep = np.concatenate([np.full(len(idx[uniq[j]]), i)
                              for i, j in enumerate(pick)])
        v = fn(y[rows], p[rows], rep)
        if v == v:
            vals.append(v)
        else:
            nan_cnt += 1
    a = np.asarray(vals, float)
    if len(a) < 20:
        return {"point": pt, "ci_low": np.nan, "ci_high": np.nan, "sd": np.nan,
                "n_valid": len(a), "n_nan": nan_cnt, "n_requested": n_boot}
    lo, hi = np.percentile(a, [2.5, 97.5])
    return {"point": float(pt), "ci_low": float(lo), "ci_high": float(hi),
            "sd": float(a.std(ddof=1)), "n_valid": len(a), "n_nan": nan_cnt,
            "n_requested": n_boot}


def macro_spearman(y, p, g):
    return FM.macro_spearman(y, p, g)


def main() -> int:
    t0 = time.time()
    emit("# V6.2 PHASE 5 — Bootstrap 修正")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- **修正**：5000 次抽样，**无 unique 阈值、不跳过迭代**")
    emit(f"- **Primary**：3 seed 预测先平均（ensemble），再算 metric")
    emit(f"- Secondary：逐 seed mean ± SD")
    emit(f"- 抽样单位: **complex**")
    emit("")

    cmp_rows, prim_rows = [], []

    # ================= CATH 单模型 CI =================
    emit("## 5.1 CATH 单模型 CI（修正版）")
    emit("")
    if STD.exists():
        d = pd.read_csv(STD, low_memory=False)
        ca = d[d["dataset_name"] == "CATH"]
        emit("| 模型 | n | n_complex | macro Spearman | 95% CI | SD | "
             "有效/请求 |")
        emit("|---|---:|---:|---:|---|---:|---|")
        for mo in sorted(ca["model"].unique()):
            g0 = ca[ca["model"] == mo]
            y = g0["true_ddG"].to_numpy(float)
            p = g0["pred_ddG"].to_numpy(float)
            g = g0["complex_id"].astype(str).to_numpy()
            r = boot_ci(macro_spearman, y, p, g)
            emit(f"| {mo} | {len(y)} | {g0['complex_id'].nunique()} | "
                 f"{r['point']:+.4f} | [{r['ci_low']:+.4f}, {r['ci_high']:+.4f}] | "
                 f"{r['sd']:.4f} | {r['n_valid']}/{r['n_requested']} |")
            cmp_rows.append({"source": "CATH", "model": mo,
                             "metric": "macro_pdb_spearman",
                             "point": r["point"], "ci_low": r["ci_low"],
                             "ci_high": r["ci_high"], "sd": r["sd"],
                             "n_valid": r["n_valid"], "n_requested": r["n_requested"]})
        emit("")

    # ================= 5.4 与 V6.1 并排 =================
    emit("## 5.4 与 V6.1 并排对比（关键：CI 是否变宽、结论是否改变）")
    emit("")
    v61 = C.ROOT / "experiments" / "v6_1" / "phase1_ci" / "bootstrap_ci_models.csv"
    if v61.exists():
        v = pd.read_csv(v61)
        v = v[v["metric"] == "macro_pdb_spearman"]
        emit("| 模型 | V6.1 point | V6.1 95% CI | V6.1 有效 | V6.2 point | "
             "V6.2 95% CI | V6.2 有效 | CI 宽度变化 | 结论改变? |")
        emit("|---|---:|---|---:|---:|---|---:|---|---|")
        for mo in sorted(set(v["model"]) & set(r["model"] for r in cmp_rows)):
            a = v[v["model"] == mo].iloc[0]
            b = [r for r in cmp_rows
                 if r["model"] == mo and r["source"] == "CATH"][0]
            wa = a["ci_high"] - a["ci_low"]
            wb = b["ci_high"] - b["ci_low"]
            sign_a = "含0" if a["ci_low"] <= 0 <= a["ci_high"] else "排除0"
            sign_b = "含0" if b["ci_low"] <= 0 <= b["ci_high"] else "排除0"
            emit(f"| {mo} | {a['point']:+.4f} | "
                 f"[{a['ci_low']:+.4f}, {a['ci_high']:+.4f}] | {a['n_valid']} | "
                 f"{b['point']:+.4f} | [{b['ci_low']:+.4f}, {b['ci_high']:+.4f}] | "
                 f"{b['n_valid']} | {wa:.4f} → {wb:.4f} "
                 f"({'**变宽**' if wb > wa else '变窄'}) | "
                 f"{'**是**' if sign_a != sign_b else '否'} |")
            cmp_rows.append({"source": "CATH_v61_vs_v62", "model": mo,
                             "metric": "macro_pdb_spearman",
                             "v61_point": a["point"], "v61_width": wa,
                             "v62_point": b["point"], "v62_width": wb,
                             "v61_n_valid": a["n_valid"], "v62_n_valid": b["n_valid"],
                             "verdict_changed": bool(sign_a != sign_b)})
        emit("")
        emit("> **V6.1 的 bug**：`if unique complex < 15: skip` → "
             "5000 次里仅 1151 次有效，被跳过的恰是抽样重复度更高的迭代，")
        emit("> 因此 **V6.1 的 CI 系统性偏窄**。上表给出宽度与「是否含 0」的实际变化。")
        emit("")

    # ================= split4 的 CI =================
    npz = PH2 / "split4_predictions.npz"
    if npz.exists():
        z = np.load(npz, allow_pickle=True)
        levels = [k[:-2] for k in sorted(z.files) if k.endswith("_y")]
        emit("## 5.3 split4 每档 × 每模型 CI（primary = seed ensemble）")
        emit("")
        emit("| 档 | 模型 | 指标 | point | 95% CI | SD | 有效/请求 |")
        emit("|---|---|---|---:|---|---:|---|")
        for lv in levels:
            if f"{lv}_y" not in z.files:
                continue
            y = z[f"{lv}_y"]; g = np.array([str(x) for x in z[f"{lv}_pdb"]])
            for mk in ("M0", "M1"):
                ks = [f"{lv}_{mk}_s{s}" for s in (42, 43, 44)
                      if f"{lv}_{mk}_s{s}" in z.files]
                if not ks:
                    continue
                m = np.isfinite(y)
                for k in ks:
                    m &= np.isfinite(z[k])
                p_ens = np.mean([z[k][m] for k in ks], axis=0)
                yy, gg = y[m], g[m]
                for name, fn in (("macro_pdb_spearman", macro_spearman),
                                 ("C_W", None)):
                    if fn is None:
                        continue
                    r = boot_ci(fn, yy, p_ens, gg)
                    emit(f"| {lv} | {mk} | {name} | {r['point']:+.4f} | "
                         f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}] | "
                         f"{r['sd']:.4f} | {r['n_valid']}/{r['n_requested']} |")
                    prim_rows.append({"level": lv, "model": mk, "metric": name,
                                      "point": r["point"], "ci_low": r["ci_low"],
                                      "ci_high": r["ci_high"], "sd": r["sd"],
                                      "n_valid": r["n_valid"],
                                      "n_requested": r["n_requested"],
                                      "seeds": "ensemble"})
                # secondary: 逐 seed
                ss = [macro_spearman(y[m], z[k][m], g[m]) for k in ks]
                emit(f"| {lv} | {mk} | （逐 seed） | "
                     f"{np.mean(ss):+.4f} ± {np.std(ss, ddof=1):.4f} | — | — | "
                     f"{len(ss)} seeds |")
        emit("")
        # primary statistics CI 在 PHASE 3 已算；此处汇总对照
        p3 = C.ROOT / "experiments" / "v6_2" / "phase3_decomposition" / "primary_statistics.csv"
        if p3.exists():
            pp = pd.read_csv(p3)
            emit("## 5.3b Primary Statistics 的 CI（来自 PHASE 3，同为无阈值 bootstrap）")
            emit("")
            emit("| 统计量 | 档 | 模型 | point | 95% CI | 有效 | 排除 0? |")
            emit("|---|---|---|---:|---|---:|---|")
            for _, r in pp.iterrows():
                emit(f"| {r['statistic']} | {r['level']} | {r['model']} | "
                     f"{r['point']:+.4f} | [{r['ci_low']:+.4f}, {r['ci_high']:+.4f}] | "
                     f"{r['n_valid']} | "
                     f"{'✅' if r['ci_excludes_zero'] else '❌ 含 0'} |")
            emit("")

    pd.DataFrame(cmp_rows).to_csv(PH5 / "bootstrap_comparison.csv", index=False,
                                  encoding="utf-8")
    pd.DataFrame(prim_rows).to_csv(PH5 / "primary_statistics_ci.csv", index=False,
                                   encoding="utf-8")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    (PH5 / "bootstrap_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {PH5 / 'bootstrap_comparison.csv'}")
    print(f"写出 {PH5 / 'primary_statistics_ci.csv'}")
    print(f"写出 {PH5 / 'bootstrap_report.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
