# -*- coding: utf-8 -*-
"""scripts/v7_3/06_negative_controls.py -- V7.3 PHASE 7：负对照（§35–§39）.

| 对照 | 设计 | 判据 |
|---|---|---|
| **`F4_SHUFFLED_ENV`**（§35） | **只在 train 内**打乱「突变位点 ↔ 局部环境」的配对（保留维度与边缘分布、保留 mutation deltas） | 若真实 F4 不优于它 → **不得**声称局部结构特异性 |
| **`F4_SCRAMBLED_CHEM`**（§36） | 保留接触数/加权总量，但**置换伙伴的化学类别**（在 partner 环境列之间重排 chem 类） | 检验增益来自**真实伙伴化学**还是仅接触几何 |
| **`F2`（静态）**（§37） | 已在 PHASE 6 完成 | `F4>F2` 且 `F2≈F0` → 支持「扰动而非静态几何」 |
| **`F1`（mutation only）**（§38） | 已在 PHASE 6 完成 | `F4>F1` → 支持几何条件化的扰动 |
| **`F4 vs F5`**（§39） | 已在 PHASE 6 完成 | partner vs intra |

实现要点
--------
* 打乱**只在 outer-train 行内**进行（绝不触碰 outer-test 的配对），
  这样 `F4_SHUFFLED_ENV` 与 `F4` 唯一差别就是**配对正确性**。
* 化学类 scramble：把 partner 环境里的 `wc_*/c_*`（按化学类别）在**该 target 内**
  的位点间重排，保留总量的量级与维度。
"""

from __future__ import annotations

import argparse
import gc
import importlib.util
import json
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

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

FEAT = C.ROOT / "features" / "v7_3"
V73 = C.ROOT / "experiments" / "v7_3"
PH7 = V73 / "phase7_controls"
PH7.mkdir(parents=True, exist_ok=True)
TREES = 300
MIN_PS_N = 5
SEEDS = [42, 43, 44]
N_BOOT = 1500
L: list[str] = []


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
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-minutes", type=float, default=60.0)
    a = ap.parse_args()
    t0 = time.time()
    emit("# V7.3 PHASE 7 — 负对照（§35–§39）")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("")

    PH = pd.read_parquet(FEAT / "physics_features.parquet")
    FD = pd.read_csv(V73 / "phase3_features" / "feature_dict.csv")
    fam_spec = json.loads((V73 / "phase5_rankers" / "family_spec.json")
                          .read_text(encoding="utf-8"))
    lab = pd.read_csv(C.ROOT / "experiments" / "v7_1_1" / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    _m = lab.set_index("pair_row")[["target_id_canonical",
                                    "family_id_canonical"]]
    PH = PH.merge(_m, left_on="row_id", right_index=True, how="left")
    PH["family_id_canonical"] = PH["family_id_canonical"].fillna("NA")
    can = pd.read_csv(C.ROOT / "experiments" / "v6_2" / "phase2_split4" /
                      "matched_analysis_set.csv", low_memory=False).set_index("row")
    d = C.build_pair_features(verbose=False)
    rows = PH.row_id.to_numpy(int)
    y = d["y"][rows].astype(float)
    g = PH.target_id_canonical.astype(str).to_numpy()
    fam = PH.family_id_canonical.astype(str).to_numpy()
    emit(f"- U_PHYS_STRICT: rows **{len(y)}** / targets **{len(set(g))}** / "
         f"families **{len(set(fam))}**")
    emit("")

    F4 = fam_spec["F4_Base_PartnerPert"]
    # partner 环境/扰动列（用于打乱与 scramble）
    env_cols = [c for c in F4 if (c.startswith("pt_") or
                                  "partner" in c or "_x_" in c)]
    chem_cols = [c for c in env_cols if ("wc_" in c or c.startswith("pt_c_")
                                         or "_x_" in c)]
    emit(f"- `F4` 列数 **{len(F4)}**；其中 partner 环境/扰动列 **{len(env_cols)}**；"
         f"化学类列 **{len(chem_cols)}**")
    emit("")

    D4 = {"y": y, "pdb": g, "fam": fam,
          "seq40": can.loc[rows, "seq40_cluster"].astype(str).to_numpy(),
          "obs": can.loc[rows, "obs_key"].astype(str).to_numpy(),
          "state": can.loc[rows, "state_id"].astype(str).to_numpy()}
    X4 = PH[F4].apply(pd.to_numeric, errors="coerce").to_numpy(np.float32)
    env_idx = [F4.index(c) for c in env_cols]
    chem_idx = [F4.index(c) for c in chem_cols]

    def fit_predict(mode):
        """mode: 'real' | 'shuffle_env' | 'scramble_chem'"""
        preds = {}
        for sd in SEEDS:
            outer = s4.make_folds(D4["fam"], sd, 5)
            yhat = np.full(len(y), np.nan)
            for tr, te in outer:
                Xtr = np.array(X4[tr], np.float32, copy=True)
                Xte = np.array(X4[te], np.float32, copy=True)
                if mode != "real":
                    rs = np.random.RandomState(sd)
                    if mode == "shuffle_env":
                        # 只在 train 内打乱环境列的行顺序（保留 mutation deltas）
                        perm = rs.permutation(len(tr))
                        Xtr[:, env_idx] = Xtr[np.ix_(perm, env_idx)]
                    elif mode == "scramble_chem":
                        # 在同 target 内重排化学类列（保留每行总量量级）
                        gt = g[tr]
                        for c in np.unique(gt):
                            m = np.where(gt == c)[0]
                            if len(m) < 2:
                                continue
                            p2 = rs.permutation(len(m))
                            Xtr[np.ix_(m, chem_idx)] = Xtr[np.ix_(m[p2],
                                                                 chem_idx)]
                med = np.nanmedian(Xtr, axis=0)
                med = np.where(np.isfinite(med), med, 0.0)
                for A in (Xtr, Xte):
                    bad = ~np.isfinite(A)
                    if bad.any():
                        A[bad] = np.take(med, np.where(bad)[1])
                lgb = C.load_lgb()
                mdl = lgb.LGBMRegressor(**{**C.LGB_PARAMS,
                                           "n_estimators": TREES,
                                           "random_state": sd})
                mdl.fit(np.vstack([Xtr, -Xtr]), np.r_[y[tr], -y[tr]])
                yhat[te] = mdl.predict(Xte)
                del mdl; gc.collect()
            preds[sd] = yhat
        return preds

    results = {}
    for mode in ("real", "shuffle_env", "scramble_chem"):
        if (time.time() - t0) / 60 > a.max_minutes:
            emit(f"- ⏱ 软时限，停止于 {mode}"); break
        pr = fit_predict(mode)
        P = np.mean([pr[s] for s in SEEDS], axis=0)
        ok = np.isfinite(P)
        mm = mset(y[ok], P[ok], g[ok])
        results[mode] = {"pred": P, **mm}
        emit(f"- `{mode}`: Macro ρ = **{mm['macro_spearman']:+.4f}** | "
             f"Eff@10 = {mm['eff10']:+.4f} | C_W = {mm['C_W']:+.4f}  "
             f"({(time.time() - t0) / 60:.1f} min)")

    emit("")
    emit("## 7.1 负对照结果（ensemble 口径）")
    emit("")
    emit("| 变体 | Macro ρ | Eff@10 | C_W | Overall r |")
    emit("|---|---:|---:|---:|---:|")
    NAME = {"real": "**真实 F4**", "shuffle_env": "`F4_SHUFFLED_ENV`",
            "scramble_chem": "`F4_SCRAMBLED_CHEM`"}
    for k in ("real", "shuffle_env", "scramble_chem"):
        if k not in results:
            continue
        r = results[k]
        emit(f"| {NAME[k]} | {r['macro_spearman']:+.4f} | {r['eff10']:+.4f} | "
             f"{r['C_W']:+.4f} | {r['overall_r']:+.4f} |")
    emit("")
    pd.DataFrame([{"variant": k, **{kk: vv for kk, vv in v.items()
                                    if kk != "pred"}}
                  for k, v in results.items()]).to_csv(
        PH7 / "negative_controls.csv", index=False, encoding="utf-8")

    # ---------- bootstrap：real − shuffled / real − scrambled ----------
    emit("## 7.2 配对 bootstrap（real − 对照）")
    emit("")
    uq = np.unique(g)
    cidx = [np.where(g == c)[0] for c in uq]
    out = []
    for other, lab in (("shuffle_env", "真实 F4 − F4_SHUFFLED_ENV"),
                       ("scramble_chem", "真实 F4 − F4_SCRAMBLED_CHEM")):
        if "real" not in results or other not in results:
            continue
        PA, PB = results["real"]["pred"], results[other]["pred"]
        m = np.isfinite(PA) & np.isfinite(PB)
        yy, pa, pb, gg = y[m], PA[m], PB[m], g[m]
        A0, B0 = mset(yy, pa, gg), mset(yy, pb, gg)
        names = ["macro_spearman", "eff10", "C_W", "overall_r", "rmse"]
        pt = {n: A0[n] - B0[n] for n in names}
        acc = {n: [] for n in names}
        rs = np.random.RandomState(42)
        for _ in range(N_BOOT):
            pick = rs.randint(0, len(uq), len(uq))
            r_ = np.concatenate([cidx[j] for j in pick])
            rep = np.concatenate([np.full(len(cidx[j]), i)
                                  for i, j in enumerate(pick)])
            try:
                A = mset(y[r_], pa[r_], rep); Bx = mset(y[r_], pb[r_], rep)
            except Exception:  # noqa: BLE001
                continue
            for n in names:
                v = A[n] - Bx[n]
                if v == v:
                    acc[n].append(v)
        emit(f"### {lab}")
        emit("")
        emit("| 指标 | point | 95% CI | 排除 0? |")
        emit("|---|---:|---|---|")
        for n in names:
            arr = np.asarray(acc[n], float)
            if len(arr) < 20:
                lo = hi = np.nan; ex = False
            else:
                lo, hi = np.percentile(arr, [2.5, 97.5])
                ex = bool((lo > 0) or (hi < 0))
            emit(f"| Δ{n} | {pt[n]:+.4f} | [{lo:+.4f}, {hi:+.4f}] | "
                 f"{'✅ 排除 0' if ex else '❌ 含 0'} |")
            out.append({"contrast": lab, "metric": f"d_{n}", "point": pt[n],
                        "ci_low": float(lo), "ci_high": float(hi),
                        "ci_excludes_zero": ex})
        emit("")
        pd.DataFrame(out).to_csv(PH7 / "control_contrasts.csv", index=False,
                                 encoding="utf-8")

    # ---------- 判定 ----------
    emit("## 7.3 负对照判定（§35/§36）")
    emit("")
    if "real" in results and "shuffle_env" in results:
        dr = results["real"]["macro_spearman"] - results["shuffle_env"]["macro_spearman"]
        emit(f"- 真实 F4 − `F4_SHUFFLED_ENV` 的 Macro ρ 差 = **{dr:+.4f}**")
        sh = [o for o in out if o["contrast"].startswith("真实 F4 − F4_SHUFFLED")
              and o["metric"] == "d_macro_spearman"]
        if sh:
            r = sh[0]
            if r["ci_excludes_zero"] and r["point"] > 0:
                emit(f"  ⇒ **环境特异性 SUPPORTED ✅**"
                     f"（Δ={r['point']:+.4f} [{r['ci_low']:+.4f}, {r['ci_high']:+.4f}]）")
                emit(f"  ⇒ 增益**依赖正确的局部相互作用环境**（§35 通过）")
            else:
                emit(f"  ⇒ **环境特异性 NOT SUPPORTED ❌**"
                     f"（Δ={r['point']:+.4f} [{r['ci_low']:+.4f}, {r['ci_high']:+.4f}] 含 0）")
                emit(f"  ⇒ **不得**声称局部结构特异性")
    if "real" in results and "scramble_chem" in results:
        sc = [o for o in out if o["contrast"].startswith("真实 F4 − F4_SCRAMBLED")
              and o["metric"] == "d_macro_spearman"]
        if sc:
            r = sc[0]
            emit(f"- 真实 F4 − `F4_SCRAMBLED_CHEM`：Δ={r['point']:+.4f} "
                 f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}] → "
                 f"{'✅ 伙伴化学重要' if r['ci_excludes_zero'] and r['point']>0 else '❌ 不显著（可能是接触几何而非化学）'}")
    emit("")
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (PH7 / "negative_controls.md").write_text("\n".join(L) + "\n",
                                              encoding="utf-8")
    print(f"\n写出 negative_controls.csv / control_contrasts.csv / "
          f"negative_controls.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
