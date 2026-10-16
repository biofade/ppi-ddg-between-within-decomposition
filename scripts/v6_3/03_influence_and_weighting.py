# -*- coding: utf-8 -*-
"""scripts/v6_3/03_influence_and_weighting.py -- V6.3 PHASE 4 + PHASE 5.

PHASE 4：family / source 集中度、leave-one-family-out、leave-top-targets-out
PHASE 5：weighting 敏感性（row / equal-target / sqrt(n) / capped30 / log(1+n)）

关键判据
--------
* 若单个 family 移除使现象消失 → `FAMILY_DRIVEN`
* 若 1–5 个超级大 target 驱动 → 记录
* `source` 在本项目为**常量**（全部 SKEMPI v2.0）→ 标 `UNAVAILABLE`
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

V62P2 = C.ROOT / "experiments" / "v6_2" / "phase2_split4"
PH1 = C.ROOT / "experiments" / "v6_3" / "phase1_target_table"
PH4 = C.ROOT / "experiments" / "v6_3" / "phase4_family_source"
PH5 = C.ROOT / "experiments" / "v6_3" / "phase5_influence"
for d in (PH4, PH5):
    d.mkdir(parents=True, exist_ok=True)
L: list[str] = []
W: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def emit_w(m: str = "") -> None:
    print(m)
    W.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def between_stats(y, p, g, w=None):
    """加权 between 分解。w=None → row-weighted（n_c/N）。"""
    uq = np.unique(g)
    By = np.array([y[g == c].mean() for c in uq])
    Bp = np.array([p[g == c].mean() for c in uq])
    cnt = np.array([(g == c).sum() for c in uq], float)
    if w is None:
        ww = cnt / cnt.sum()
    elif w == "equal":
        ww = np.ones(len(uq)) / len(uq)
    elif w == "sqrt":
        ww = np.sqrt(cnt); ww /= ww.sum()
    elif w == "cap30":
        cc = np.minimum(cnt, 30.0); ww = cc / cc.sum()
    elif w == "log1p":
        cc = np.log1p(cnt); ww = cc / cc.sum()
    else:
        ww = np.ones(len(uq)) / len(uq)
    yb, pb = y.mean(), p.mean()
    C_B = float((ww * (By - yb) * (Bp - pb)).sum()) / (y.std() * p.std())
    # within：用同样权重对每 target 的 within 协方差加权
    vy = np.array([np.mean((y[g == c] - By[i]) * (p[g == c] - Bp[i]))
                   for i, c in enumerate(uq)])
    C_W = float((ww * vy).sum()) / (y.std() * p.std())
    By2 = np.array([y[g == c].mean() for c in uq])
    Bp2 = np.array([p[g == c].mean() for c in uq])
    mr = C.pearson(By2, Bp2) if len(By2) >= 3 else float("nan")
    wmr = (float(np.corrcoef(By2, Bp2)[0, 1]) if len(By2) >= 3 else float("nan"))
    return dict(C_B=C_B, C_W=C_W, macro_mean_pearson=mr, n_target=len(uq))


def main() -> int:
    t0 = time.time()
    emit("# V6.3 PHASE 4 + 5 — Family/Source 集中度 + Influence + Weighting")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("")
    T = pd.read_csv(PH1 / "target_mechanism_table.csv", low_memory=False)
    z6 = np.load(V62P2 / "split4_predictions.npz", allow_pickle=True)
    lv = "L4_family"
    y = z6[f"{lv}_y"].astype(float)
    g = np.array([str(x) for x in z6[f"{lv}_pdb"]])
    ks = [f"{lv}_M1_s{s}" for s in (42, 43, 44) if f"{lv}_M1_s{s}" in z6.files]
    p = np.mean([z6[k] for k in ks], axis=0).astype(float)
    m = np.isfinite(y) & np.isfinite(p)
    y, p, g = y[m], p[m], g[m]
    emit(f"- L4 M1 ensemble: {len(y)} 行 / {len(set(g))} targets")
    emit("")

    # ================= PHASE 4.1 family 集中度 =================
    emit("## 4.1 Family 集中度")
    emit("")
    T["grp"] = pd.cut(T.n_mutations, bins=[0, 9, 29, 10 ** 9],
                      labels=["2-9", "10-29", ">=30"])
    rich = T[T.grp == ">=30"]
    poor = T[T.grp != ">=30"]
    emit("| 集合 | n_target | n_family | rows | top1 family 占比(targets) | "
         "top1 family 占比(rows) |")
    emit("|---|---:|---:|---:|---:|---:|")
    for lab, s in ((">=30", rich), ("<30", poor), ("ALL", T)):
        if not len(s):
            continue
        fc = s.family.value_counts()
        rowsum = s.groupby("family").n_rows.sum()
        emit(f"| {lab} | {len(s)} | {s.family.nunique()} | {int(s.n_rows.sum())} | "
             f"{fc.iloc[0] / len(s):.3f} | {rowsum.iloc[0] / s.n_rows.sum():.3f} |")
    emit("")
    emit(f"- `>=30` 组 top-5 family：")
    for f, n in rich.family.value_counts().head(5).items():
        rr = int(rich.loc[rich.family == f, "n_rows"].sum())
        emit(f"  - `{f}`：{n} targets / {rr} rows")
    emit("")

    # ================= 4.2 source 集中度 =================
    emit("## 4.2 Source 集中度")
    emit("")
    emit(f"- `source_top` 唯一值数: **{T.source_top.nunique()}** → "
         f"值 = `{T.source_top.iloc[0]}`")
    emit(f"- **判定：`UNAVAILABLE`** —— 全部 target 来自同一 source，"
         f"无法做 source 分层或 leave-one-source-out")
    emit(f"- `assay_top` 唯一值数: **{T.assay_top.nunique()}** → 可用作 assay 分层")
    emit("")
    emit("| assay | n_target | n_rows | mean abs_err | mean signed_err |")
    emit("|---|---:|---:|---:|---:|")
    for a, s in T.groupby("assay_top"):
        if len(s) < 3:
            continue
        emit(f"| `{a}` | {len(s)} | {int(s.n_rows.sum())} | "
             f"{s.abs_mean_error.mean():.4f} | {s.signed_mean_error.mean():+.4f} |")
    emit("")

    # ================= 4.3 leave-one-family-out =================
    emit("## 4.3 Leave-one-largest-family-out（在 >=30 组内）")
    emit("")
    base = between_stats(y, p, g)
    emit("| 移除的 family | 剩余 target | 移除 rows | macro mean-Pearson | row-w C_B |")
    emit("|---|---:|---:|---:|---:|")
    infl = []
    top_fams = rich.family.value_counts().head(8).index.tolist()
    for f in top_fams:
        drop_t = set(rich.loc[rich.family == f, "pdb_id"])
        mm = ~np.isin(g, list(drop_t))
        if mm.sum() < 100:
            continue
        st = between_stats(y[mm], p[mm], g[mm])
        emit(f"| `{f}` | {len(set(g[mm]))} | {int((~mm).sum())} | "
             f"{st['macro_mean_pearson']:+.4f} | {st['C_B']:+.4f} |")
        infl.append({"analysis": "leave_family_out", "family": f,
                     "n_target_removed": len(drop_t),
                     "macro_mean_pearson": st["macro_mean_pearson"],
                     "C_B": st["C_B"]})
    emit("")
    emit(f"- **基线**：macro mean-Pearson = **{base['macro_mean_pearson']:+.4f}**，"
         f"row-w C_B = **{base['C_B']:+.4f}**")
    if infl:
        mr = [x["macro_mean_pearson"] for x in infl]
        emit(f"- 移除任一 top family 后 macro mean-Pearson 范围: "
             f"[{min(mr):+.4f}, {max(mr):+.4f}]")
        driven = (max(mr) > 0.3 and base["macro_mean_pearson"] < 0.1)
        verdict4 = "FAMILY_DRIVEN" if driven else "NOT_SINGLE_FAMILY_DRIVEN"
        emit(f"- **单 family 移除是否使现象消失**: "
             f"{'是 → FAMILY_DRIVEN' if driven else '否'}")
    else:
        verdict4 = "INCONCLUSIVE"
    emit("")

    # ================= 4.4 leave-top-targets-out =================
    emit("## 4.4 Leave-top-targets-out（按 n_mutations 降序）")
    emit("")
    emit("| 移除 | 剩余 targets | macro mean-Pearson | row-w C_B | overall r |")
    emit("|---|---:|---:|---:|---:|")
    order = T.sort_values("n_mutations", ascending=False).pdb_id.tolist()
    for k in (1, 3, 5, 10, 20):
        drop_t = set(order[:k])
        mm = ~np.isin(g, list(drop_t))
        st = between_stats(y[mm], p[mm], g[mm])
        r_all = C.pearson(y[mm], p[mm])
        emit(f"| top {k} | {len(set(g[mm]))} | {st['macro_mean_pearson']:+.4f} | "
             f"{st['C_B']:+.4f} | {r_all:+.4f} |")
        infl.append({"analysis": "leave_top_targets", "k": k,
                     "macro_mean_pearson": st["macro_mean_pearson"],
                     "C_B": st["C_B"], "overall_r": r_all})
    emit("")
    pd.DataFrame(infl).to_csv(PH4 / "influence_analysis.csv", index=False,
                              encoding="utf-8")
    pd.DataFrame([{"family": f, "n_target": int(n),
                   "n_rows": int(rich.loc[rich.family == f, "n_rows"].sum())}
                  for f, n in rich.family.value_counts().items()]
                 ).to_csv(PH4 / "family_concentration.csv", index=False,
                          encoding="utf-8")
    pd.DataFrame([{"source": T.source_top.iloc[0],
                   "n_target": len(T), "status": "UNAVAILABLE_SINGLE_SOURCE"}]
                 ).to_csv(PH4 / "source_concentration.csv", index=False,
                          encoding="utf-8")
    emit("### 4.3 判定")
    emit("")
    emit(f"> **`{verdict4}`**")
    emit("")

    # ================= PHASE 5 weighting =================
    emit("## PHASE 5 — Weighting 敏感性")
    emit("")
    emit("> 目的**不是**挑最好看的 weighting，而是回答：")
    emit("> **performance conclusion 对 target-size weighting 有多敏感？**")
    emit("")
    emit("| weighting | L4 C_B | L4 C_W | L4 macro mean-Pearson | n_target |")
    emit("|---|---:|---:|---:|---:|")
    wrows = []
    for wname, warg in (("row (n_c/N)", None), ("equal-target", "equal"),
                        ("sqrt(n)", "sqrt"), ("capped 30", "cap30"),
                        ("log(1+n)", "log1p")):
        st = between_stats(y, p, g, w=warg)
        emit(f"| {wname} | {st['C_B']:+.4f} | {st['C_W']:+.4f} | "
             f"{st['macro_mean_pearson']:+.4f} | {st['n_target']} |")
        wrows.append({"weighting": wname, "C_B": st["C_B"], "C_W": st["C_W"],
                      "macro_mean_pearson": st["macro_mean_pearson"]})
    emit("")
    # L1 对照
    y1 = z6["L1_random_y"].astype(float)
    g1 = np.array([str(x) for x in z6["L1_random_pdb"]])
    k1 = [f"L1_random_M1_s{s}" for s in (42, 43, 44)
          if f"L1_random_M1_s{s}" in z6.files]
    p1 = np.mean([z6[k] for k in k1], axis=0).astype(float)
    m1 = np.isfinite(y1) & np.isfinite(p1)
    emit("| weighting | L1 C_B | L1 macro mean-Pearson |")
    emit("|---|---:|---:|")
    w1rows = []
    for wname, warg in (("row (n_c/N)", None), ("equal-target", "equal"),
                        ("sqrt(n)", "sqrt"), ("capped 30", "cap30"),
                        ("log(1+n)", "log1p")):
        st = between_stats(y1[m1], p1[m1], g1[m1], w=warg)
        emit(f"| {wname} | {st['C_B']:+.4f} | {st['macro_mean_pearson']:+.4f} |")
        w1rows.append({"weighting": wname, "L1_C_B": st["C_B"],
                       "L1_macro_mean_pearson": st["macro_mean_pearson"]})
    emit("")
    W5 = pd.DataFrame(wrows).merge(pd.DataFrame(w1rows), on="weighting")
    W5.to_csv(PH5 / "weighting_sensitivity.csv", index=False, encoding="utf-8")
    spread = W5.C_B.max() - W5.C_B.min()
    emit(f"- **L4 `C_B` 跨 weighting 的极差 = {spread:.4f}**"
         f"（min {W5.C_B.min():+.4f} / max {W5.C_B.max():+.4f}）")
    emit(f"- 符号是否随 weighting 变化: "
         f"**{'是（结论不稳定）' if (W5.C_B.min() < 0 < W5.C_B.max()) else '否'}**")
    emit("")
    emit_w("# V6.3 PHASE 5 — Weighting 敏感性")
    emit_w("")
    emit_w("| weighting | L4 C_B | L4 C_W | L4 macro mean-Pearson | L1 C_B |")
    emit_w("|---|---:|---:|---:|---:|")
    for _, r in W5.iterrows():
        emit_w(f"| {r['weighting']} | {r['C_B']:+.4f} | {r['C_W']:+.4f} | "
               f"{r['macro_mean_pearson']:+.4f} | {r['L1_C_B']:+.4f} |")
    emit_w("")
    emit_w(f"- L4 `C_B` 极差 **{spread:.4f}**；符号不随 weighting 翻转时结论稳健。")
    emit_w("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    (PH4 / "family_source_report.md").write_text("\n".join(L) + "\n",
                                                encoding="utf-8")
    (PH5 / "weighting_report.md").write_text("\n".join(W) + "\n",
                                             encoding="utf-8")
    print(f"\n写出 {PH4}/ 与 {PH5}/")
    print(f"\n>>> PHASE4 判定: {verdict4}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
