# -*- coding: utf-8 -*-
"""scripts/v6_3/02_univariate_and_multivariable.py -- V6.3 PHASE 2 + PHASE 3.

分析单位 = **target**（N ≈ 203），**禁止**把 7523 rows 当独立样本。

PHASE 2：单变量关联（Spearman + BH-FDR；分组变量用 Kruskal-Wallis / Mann-Whitney）
PHASE 3：多变量 + 混杂控制（family-level CV；partial association）

**关键判定**：`log_n_mutations` 在控制 library composition 后是否仍有独立关联 →
`SIZE_IS_PROXY` / `SIZE_HAS_INDEPENDENT_ASSOCIATION` / `INCONCLUSIVE`

**只能说 association，不能说 causation。**
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

PH1 = C.ROOT / "experiments" / "v6_3" / "phase1_target_table"
PH2 = C.ROOT / "experiments" / "v6_3" / "phase2_univariate"
PH3 = C.ROOT / "experiments" / "v6_3" / "phase3_multivariable"
for d in (PH2, PH3):
    d.mkdir(parents=True, exist_ok=True)
TBL = PH1 / "target_mechanism_table.csv"
O_UNI = PH2 / "univariate_associations.csv"
O_TOP = PH2 / "top_associations.md"
O_MOD = PH3 / "model_comparison.csv"
O_COEF = PH3 / "coefficients.csv"
O_CV = PH3 / "cv_results.csv"
O_PAR = PH3 / "partial_associations.csv"
O_MD = PH3 / "multivariable_report.md"
L: list[str] = []
U: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def emit_u(m: str = "") -> None:
    print(m)
    U.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


# 部署可用 vs label-derived（与 PHASE 1 字典一致）
LABEL_DERIVED = {"true_mean_ddg", "signed_mean_error", "abs_mean_error",
                 "squared_mean_error", "ddg_std", "ddg_range", "ddg_IQR",
                 "beneficial_fraction", "deleterious_fraction",
                 "extreme_fraction", "within_spearman", "within_kendall"}
ID_COLS = {"pdb_id", "family", "source_top", "assay_top"}

OUTCOMES = [("abs_mean_error", "primary"),
            ("signed_mean_error", "secondary"),
            ("within_spearman", "secondary")]


def main() -> int:
    t0 = time.time()
    emit("# V6.3 PHASE 2 + 3 — 单变量与多变量机制分析")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("- **分析单位 = target**（N ≈ 203）；禁止 row-level 当独立样本")
    emit("")
    T = pd.read_csv(TBL, low_memory=False)
    emit(f"- 载入: **{T.shape[0]}** targets × {T.shape[1]} cols")
    emit("")
    # 数值特征
    num_cols = [c for c in T.columns
                if c not in ID_COLS
                and pd.api.types.is_numeric_dtype(pd.to_numeric(T[c],
                                                               errors="coerce"))]
    num_cols = [c for c in num_cols if T[c].notna().sum() >= 50]
    emit(f"- 可分析数值特征: **{len(num_cols)}**"
         f"（非空 ≥50 个 target）")
    emit("")

    # ================= PHASE 2 =================
    from scipy.stats import spearmanr
    emit("## PHASE 2 — 单变量关联（Spearman + BH-FDR）")
    emit("")
    rows = []
    for out, kind in OUTCOMES:
        yo = pd.to_numeric(T[out], errors="coerce")
        for c in num_cols:
            if c == out:
                continue
            x = pd.to_numeric(T[c], errors="coerce")
            m = x.notna() & yo.notna()
            if m.sum() < 30 or x[m].nunique() < 3:
                continue
            rho, p = spearmanr(x[m], yo[m])
            if rho != rho:
                continue
            rows.append({"outcome": out, "outcome_kind": kind, "feature": c,
                         "feature_kind": ("LABEL_DERIVED"
                                          if c in LABEL_DERIVED else "DEPLOYABLE"),
                         "spearman_rho": float(rho), "p": float(p),
                         "n": int(m.sum())})
    U1 = pd.DataFrame(rows)
    # BH-FDR per outcome
    parts = []
    for out, _ in OUTCOMES:
        s = U1[U1.outcome == out].copy()
        if not len(s):
            continue
        s = s.sort_values("p")
        q = C.bh_fdr(list(s["p"]))
        s["bh_q"] = q
        s["abs_rho"] = s.spearman_rho.abs()
        parts.append(s)
    UA = pd.concat(parts, ignore_index=True) if parts else U1
    UA.to_csv(O_UNI, index=False, encoding="utf-8")

    # 输出重点
    emit("### 2.3 必须报告的变量（针对 primary outcome `abs_mean_error`）")
    emit("")
    emit("| 特征 | 类型 | Spearman ρ | p | BH q | n |")
    emit("|---|---|---:|---:|---:|---:|")
    FOCUS = ["log_n_mutations", "n_mutations", "n_unique_mutation_sites",
             "site_coverage_all", "site_coverage_interface", "site_entropy",
             "mutation_pair_entropy", "alanine_fraction", "single_mut_fraction",
             "ddg_std", "ddg_range", "ddg_IQR",
             "struct_n_iface_res_lig", "struct_n_iface_res_rec",
             "struct_n_contact_pairs", "struct_dsasa_total",
             "struct_interface_residue_ratio", "mutations_per_site",
             "mutation_count_gini", "conservative_fraction"]
    prim = UA[(UA.outcome == "abs_mean_error")].set_index("feature")
    focus_rows = []
    for f in FOCUS:
        if f not in prim.index:
            emit(f"| `{f}` | — | — | — | — | UNAVAILABLE |")
            continue
        r = prim.loc[f]
        emit(f"| `{f}` | {r.feature_kind} | {r.spearman_rho:+.4f} | {r.p:.3g} | "
             f"{r.bh_q:.3g} | {int(r.n)} |")
        focus_rows.append({"feature": f, "kind": r.feature_kind,
                           "rho": r.spearman_rho, "p": r.p, "q": r.bh_q,
                           "n": int(r.n)})
    emit("")

    # top deployable / label-derived
    emit("### Top associations（按 |ρ|）")
    emit("")
    for kind in ("DEPLOYABLE", "LABEL_DERIVED"):
        s = UA[(UA.outcome == "abs_mean_error") & (UA.feature_kind == kind)]
        s = s.reindex(s.abs_rho.sort_values(ascending=False).index).head(10)
        emit(f"**{kind}**（outcome = `abs_mean_error`）")
        emit("")
        emit("| 特征 | ρ | p | BH q |")
        emit("|---|---:|---:|---:|")
        for _, r in s.iterrows():
            emit(f"| `{r.feature}` | {r.spearman_rho:+.4f} | {r.p:.3g} | "
                 f"{r.bh_q:.3g} |")
        emit("")

    # 分组变量
    emit("### 2.2 分组变量")
    emit("")
    emit("| 分组变量 | 检验 | 统计量 | p | 组数 | effect size |")
    emit("|---|---|---:|---:|---:|---:|")
    from scipy.stats import kruskal, mannwhitneyu
    grp_rows = []
    yo = pd.to_numeric(T.abs_mean_error, errors="coerce")
    for gc, lab in (("assay_top", "assay 类型"),):
        if gc not in T.columns:
            continue
        groups = [pd.to_numeric(T.loc[T[gc] == v, "abs_mean_error"],
                                errors="coerce").dropna().to_numpy()
                  for v in T[gc].dropna().unique()]
        groups = [g for g in groups if len(g) >= 3]
        if len(groups) < 2:
            emit(f"| {lab} | Kruskal-Wallis | — | — | {len(groups)} | "
                 f"组不足 |")
            continue
        h, p = kruskal(*groups)
        # effect size: epsilon-squared
        n = sum(len(g) for g in groups)
        eps2 = (h - len(groups) + 1) / (n - len(groups)) if n > len(groups) else np.nan
        emit(f"| {lab} | Kruskal-Wallis | {h:.3f} | {p:.4g} | {len(groups)} | "
             f"ε²={eps2:.4f} |")
        grp_rows.append({"group_var": gc, "test": "kruskal", "stat": h,
                         "p": p, "n_groups": len(groups), "effect": eps2})
    # 二分组：mutation-rich vs not（预注册 ≥30）
    rich = T.n_mutations >= 30
    a = pd.to_numeric(T.loc[rich, "abs_mean_error"], errors="coerce").dropna()
    b = pd.to_numeric(T.loc[~rich, "abs_mean_error"], errors="coerce").dropna()
    if len(a) >= 3 and len(b) >= 3:
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        # rank-biserial
        rb = 1 - (2 * u) / (len(a) * len(b))
        emit(f"| `n_mutations ≥30` vs `<30` | Mann-Whitney | {u:.1f} | {p:.4g} | "
             f"{len(a)} vs {len(b)} | rank-biserial={abs(rb):.4f} |")
        emit("")
        emit(f"- ≥30 组 |mean error| 中位 **{a.median():.4f}**（n={len(a)}）"
             f"vs <30 组 **{b.median():.4f}**（n={len(b)}）")
        grp_rows.append({"group_var": "n_mutations>=30", "test": "mannwhitney",
                         "stat": u, "p": p, "n_groups": 2, "effect": abs(rb)})
    emit("")
    pd.DataFrame(grp_rows).to_csv(PH2 / "group_tests.csv", index=False,
                                  encoding="utf-8")

    # top_associations.md
    emit_u("# V6.3 PHASE 2 — Top associations")
    emit_u("")
    emit_u("Primary outcome：`abs_mean_error`（target-level，N=203）")
    emit_u("")
    for kind in ("DEPLOYABLE", "LABEL_DERIVED"):
        s = UA[(UA.outcome == "abs_mean_error") & (UA.feature_kind == kind)]
        s = s.reindex(s.abs_rho.sort_values(ascending=False).index).head(12)
        emit_u(f"## {kind}")
        emit_u("")
        emit_u("| 特征 | ρ | p | BH q |")
        emit_u("|---|---:|---:|---:|")
        for _, r in s.iterrows():
            emit_u(f"| `{r.feature}` | {r.spearman_rho:+.4f} | {r.p:.3g} | "
                   f"{r.bh_q:.3g} |")
        emit_u("")
    emit_u("## 关键读法")
    emit_u("")
    emit_u("- `SIZE_HAS_INDEPENDENT_ASSOCIATION` 的前提是 `log_n_mutations` "
            "在**多变量控制后**仍保持强关联（见 PHASE 3）。")
    emit_u("- 单变量显著**不等于**因果；本项目**没有** causal design。")
    O_TOP.write_text("\n".join(U) + "\n", encoding="utf-8")

    # ================= PHASE 3 =================
    emit("## PHASE 3 — 多变量 / 混杂控制")
    emit("")

    def zscore(s):
        s = pd.to_numeric(s, errors="coerce").astype(float)
        sd = s.std(ddof=0)
        return (s - s.mean()) / sd if sd > 0 else s * 0.0

    # 候选特征集（只保留非高度共线）
    SET_B = ["log_n_mutations", "n_unique_mutation_sites", "site_entropy",
             "mutation_pair_entropy", "alanine_fraction", "single_mut_fraction"]
    SET_C = SET_B + ["struct_n_iface_res_lig", "struct_n_iface_res_rec",
                     "struct_n_contact_pairs", "struct_dsasa_total",
                     "struct_interface_residue_ratio", "protein_length"]
    SET_D = ["ddg_std", "ddg_range", "beneficial_fraction"]
    MODELS = {
        "A_size_only": ["log_n_mutations"],
        "B_library_composition": SET_B,
        "C_context": SET_C,
        "D_explanatory_only": SET_B + SET_D,
    }

    y = pd.to_numeric(T.abs_mean_error, errors="coerce").to_numpy(float)
    fam = T["family"].astype(str).to_numpy()
    cv_rows, coef_rows, mod_rows = [], [], []

    from sklearn.linear_model import RidgeCV, LinearRegression
    from sklearn.model_selection import GroupKFold
    from sklearn.metrics import mean_absolute_error, r2_score
    from scipy.stats import spearmanr as spr

    for name, feats in MODELS.items():
        feats = [f for f in feats if f in T.columns
                 and T[f].notna().sum() >= 100]
        if not feats:
            continue
        X = np.column_stack([zscore(T[f]).to_numpy(float) for f in feats])
        m = np.isfinite(y) & np.isfinite(X).all(axis=1)
        Xm, ym, fm = X[m], y[m], fam[m]
        # 全样本拟合（报告系数）
        lr = LinearRegression().fit(Xm, ym)
        for f, c in zip(feats, lr.coef_):
            coef_rows.append({"model": name, "feature": f,
                              "std_coef": float(c),
                              "kind": ("LABEL_DERIVED" if f in LABEL_DERIVED
                                       else "DEPLOYABLE"),
                              "n": int(m.sum())})
        # family-level CV
        n_g = len(set(fm))
        k = min(5, max(2, n_g))
        gkf = GroupKFold(n_splits=k)
        pred = np.full(len(ym), np.nan)
        for tr, te in gkf.split(Xm, ym, groups=fm):
            mu = LinearRegression().fit(Xm[tr], ym[tr])
            pred[te] = mu.predict(Xm[te])
        ok = np.isfinite(pred)
        r2 = r2_score(ym[ok], pred[ok])
        mae = mean_absolute_error(ym[ok], pred[ok])
        rho = spr(ym[ok], pred[ok]).statistic
        cv_rows.append({"model": name, "n_features": len(feats),
                        "n": int(ok.sum()), "cv_r2": float(r2),
                        "cv_mae": float(mae), "cv_spearman": float(rho),
                        "features": "|".join(feats)})
        emit(f"- **{name}**（{len(feats)} 特征）：family-CV R² = **{r2:+.4f}**，"
             f"MAE = {mae:.4f}，Spearman(pred,obs) = {rho:+.4f}")
    emit("")
    emit("## 3.7 Cross-validation 汇总（GroupKFold by family）")
    emit("")
    emit("| 模型 | 特征数 | CV R² | CV MAE | CV Spearman |")
    emit("|---|---:|---:|---:|---:|")
    for r in cv_rows:
        emit(f"| {r['model']} | {r['n_features']} | {r['cv_r2']:+.4f} | "
             f"{r['cv_mae']:.4f} | {r['cv_spearman']:+.4f} |")
    emit("")
    pd.DataFrame(cv_rows).to_csv(O_CV, index=False, encoding="utf-8")
    pd.DataFrame(coef_rows).to_csv(O_COEF, index=False, encoding="utf-8")
    pd.DataFrame(cv_rows).to_csv(O_MOD, index=False, encoding="utf-8")

    # ---------- 3.8 Partial association ----------
    emit("## 3.8 Partial association（残差化）")
    emit("")
    emit("**方法**：对 `log_n_mutations` 与 `abs_mean_error` 分别回归掉控制变量，")
    emit("再对**残差**求 Spearman。同时报告普通 Spearman。")
    emit("")
    CTRL_SETS = {
        "none": [],
        "n_unique_sites": ["n_unique_mutation_sites"],
        "diversity": ["n_unique_mutation_sites", "site_entropy",
                      "mutation_pair_entropy"],
        "diversity+alanine": ["n_unique_mutation_sites", "site_entropy",
                              "mutation_pair_entropy", "alanine_fraction"],
        "full_B": SET_B[1:],
        "full_C": SET_C[1:],
    }
    par_rows = []
    emit("| 控制集 | 普通 ρ | partial ρ | p | n |")
    emit("|---|---:|---:|---:|---:|")
    x0 = pd.to_numeric(T["log_n_mutations"], errors="coerce")
    for cname, ctrls in CTRL_SETS.items():
        ctrls = [c for c in ctrls if c in T.columns
                 and T[c].notna().sum() >= 100]
        if not ctrls:
            rho0 = spr(x0, pd.Series(y)).statistic
            par_rows.append({"control_set": cname, "plain_rho": float(rho0),
                             "partial_rho": float(rho0), "p": np.nan,
                             "n": int(np.isfinite(y).sum())})
            emit(f"| {cname} | {rho0:+.4f} | {rho0:+.4f}（= 普通）| — | "
                 f"{int(np.isfinite(y).sum())} |")
            continue
        Xc = np.column_stack([zscore(T[c]).to_numpy(float) for c in ctrls])
        mm = np.isfinite(y) & np.isfinite(Xc).all(axis=1) & x0.notna().to_numpy()
        Xc, ym, xm = Xc[mm], y[mm], x0.to_numpy(float)[mm]
        # 残差化
        rx = xm - LinearRegression().fit(Xc, xm).predict(Xc)
        ry = ym - LinearRegression().fit(Xc, ym).predict(Xc)
        pr, pp = spr(rx, ry)
        pl, _ = spr(xm, ym)
        emit(f"| {cname} | {pl:+.4f} | **{pr:+.4f}** | {pp:.3g} | {len(ym)} |")
        par_rows.append({"control_set": cname, "plain_rho": float(pl),
                         "partial_rho": float(pr), "p": float(pp),
                         "n": int(len(ym))})
    emit("")
    pd.DataFrame(par_rows).to_csv(O_PAR, index=False, encoding="utf-8")

    # ---------- 判定 ----------
    emit("## Primary mechanism decision")
    emit("")
    pr_full = [r for r in par_rows if r["control_set"] in ("diversity", "full_B",
                                                           "full_C")]
    plain = [r["plain_rho"] for r in par_rows
             if r["control_set"] == "none"][0]
    verdict = "INCONCLUSIVE"
    if pr_full:
        worst = max(abs(r["partial_rho"]) for r in pr_full)
        best_p = min(r["p"] for r in pr_full)
        emit(f"- 普通 Spearman(log_n_mut, |error|) = **{plain:+.4f}**")
        emit(f"- 控制 diversity 后 |partial ρ| 最大 = **{worst:.4f}**"
             f"（最小 p = {best_p:.3g}）")
        emit(f"- 衰减比例 = {1 - worst / abs(plain):.2%}" if plain else "")
        emit("")
        if worst < 0.15 or best_p > 0.05:
            verdict = "SIZE_IS_PROXY"
            emit("> ### **`SIZE_IS_PROXY`**")
            emit(">")
            emit("> `log_n_mutations` 的关联在控制 mutation-library diversity 后"
                 "**基本消失** → size 本身不是独立解释变量。")
        elif worst > 0.3 and best_p < 0.01:
            verdict = "SIZE_HAS_INDEPENDENT_ASSOCIATION"
            emit("> ### **`SIZE_HAS_INDEPENDENT_ASSOCIATION`**")
            emit(">")
            emit("> 控制 diversity 后仍有强独立关联。**注意：这只是 association，"
                 "不是 causation。**")
        else:
            verdict = "INCONCLUSIVE"
            emit("> ### **`INCONCLUSIVE`**")
            emit(">")
            emit("> 控制后关联衰减到中间水平，无明确判定。")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    O_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {O_UNI}")
    print(f"写出 {O_TOP}")
    print(f"写出 {O_MOD} / {O_CV} / {O_COEF} / {O_PAR}")
    print(f"写出 {O_MD}")
    print(f"\n>>> {verdict}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
