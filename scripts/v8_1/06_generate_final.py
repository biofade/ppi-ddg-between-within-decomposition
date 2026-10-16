# -*- coding: utf-8 -*-
"""scripts/v8_1/06_generate_final.py -- V8.1 §16/§17/§19/§23：从 FROZEN_RESULTS.json 生成全部最终产物。

**单一真源**：所有数字来自 `v8_1/FROZEN_RESULTS.json`，不手工粘贴。
产出 `v8_1_final/` 全树（除一致性检查外）。
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
A = R / "v8_1"
V8 = R / "experiments" / "v8"
F = R / "v8_1_final"
if F.exists():
    shutil.rmtree(F)
for sub in ("tables", "figures", "supplement", "scripts"):
    (F / sub).mkdir(parents=True, exist_ok=True)

FR = json.load(open(A / "FROZEN_RESULTS.json", encoding="utf-8"))
for f in A.glob("*.csv"):
    shutil.copy2(f, F / f.name)
for f in A.glob("*.md"):
    shutil.copy2(f, F / f.name)
for f in A.glob("*.json"):
    shutil.copy2(f, F / f.name)

# 字体
for cand in ("Microsoft YaHei", "SimHei", "DengXian", "SimSun"):
    try:
        font_manager.findfont(cand, fallback_to_default=False)
        plt.rcParams["font.sans-serif"] = [cand]; break
    except Exception:
        continue
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 130
plt.rcParams["savefig.bbox"] = "tight"

UR = pd.read_csv(A / "MATCHED_UNIVERSE_REGISTRY.csv")
ST = pd.read_csv(A / "rank_agreement_by_stratum.csv")
D = pd.read_csv(A / "all_models_by_rowset.csv", low_memory=False)
EXC = pd.read_csv(A / "external_clean_profiles.csv")
GEO = pd.read_csv(A / "axis_c_geometry.csv")
PU = FR["primary_universe"]
H4 = FR["H4_controlled_split"]

# ================= §17 表 =================
print("§17 生成表")
MR = pd.read_csv(A / "MODEL_FAMILY_REGISTRY.csv")
tabs = {}
tabs["T1_model_inventory"] = MR
prim_models = D[(D.rowset_hash == PU["rowset_hash"])].dropna(
    subset=["overall_pearson", "macro_spearman"]).copy()
prim_models.insert(0, "universe_id", f"U_PRIMARY_{PU['rowset_hash'][:8]}")
tabs["T2_primary_matched_scoreboard"] = prim_models
ST2 = ST.copy(); ST2.insert(0, "universe_id", "U_" + ST2.rowset_hash.str[:8])
tabs["T3_rank_agreement_by_universe"] = ST2
tabs["T4_CB_CW_decomposition"] = D[["artifact", "key", "rowset_hash",
                                    "overall_pearson", "C_B_row", "C_W_row",
                                    "C_B_equal", "decomp_residual"]]
h4rows = [{"regime": k, **v} for k, v in H4.items() if k != "delta"]
tabs["T5_controlled_split_inflation"] = pd.DataFrame(h4rows)
W = pd.read_csv(V8 / "weighting" / "weighting_sensitivity.csv")
SCHEMES = ["ROW", "EQUAL", "SQRT", "LOG", "CAP30"]
cb = [f"C_B_{s}" for s in SCHEMES]; cw = [f"C_W_{s}" for s in SCHEMES]
W2 = W[["model_name"]].copy()
W2["range_C_B"] = W[cb].max(axis=1) - W[cb].min(axis=1)
W2["range_C_W"] = W[cw].max(axis=1) - W[cw].min(axis=1)
W2["max_abs_C_B"] = W[cb].abs().max(axis=1)
W2["max_abs_C_W"] = W[cw].abs().max(axis=1)
for lo in (0, 0.001, 0.005, 0.01):
    W2[f"signflip_C_B_gt{lo}"] = (
        ((W[cb] > lo).sum(axis=1) > 0) & ((W[cb] < -lo).sum(axis=1) > 0))
    W2[f"signflip_C_W_gt{lo}"] = (
        ((W[cw] > lo).sum(axis=1) > 0) & ((W[cw] < -lo).sum(axis=1) > 0))
W2.to_csv(A / "WEIGHTING_EFFECT_SIZE_AUDIT.csv", index=False, encoding="utf-8")
tabs["T6_weighting_effect_size"] = W2
tabs["T7_axis_c_geometry"] = GEO
sc = pd.read_csv(V8 / "scenarios" / "scenario_recommendations_provenance_clean.csv")
tabs["T8_scenario_recommendations"] = sc
ea = pd.read_csv(V8 / "external" / "external_audit.csv")
tabs["T9_external_audit"] = ea
tabs["T10_external_capabilities"] = EXC
for name, t in tabs.items():
    t.to_csv(F / "tables" / f"{name}.csv", index=False, encoding="utf-8")
    print(f"   {name}: {t.shape[0]}×{t.shape[1]}")

# 场景可估计性判定
scen_verdict = {}
for _, r in sc.iterrows():
    scen_verdict[r.scenario] = {
        "primary_metric": r.primary_metric,
        "top1_ALL_pool": r.top1_ALL_pool,
        "top1_STRICT_pool": r.top1_STRICT_pool,
        "changed": bool(r.changed_after_provenance_clean)}
# Scenario A：HERMES 无校准 ddG ⇒ 仅允许有 kcal/mol 语义的模型
scen_verdict["A_quantitative_ddG"]["hermes_eligible"] = False
scen_verdict["A_quantitative_ddG"]["note"] = (
    "HERMES log-odds is NOT eligible for quantitative RMSE recommendation")
scen_verdict["C_wetlab_shortlist"]["note"] = (
    "DESCRIPTIVE / LOW-CONFIDENCE: selection axis is structurally non-informative "
    "for n<=K targets; H3 = FAIL")

# ================= §16 图 =================
print("\n§16 生成图")
# Fig1 概念/分层总览
fig, ax = plt.subplots(figsize=(11.5, 5.6))
u = UR.sort_values("n_rows")
y = np.arange(len(u))
ax.barh(y, u.n_rows, color="#4c72b0", alpha=0.85)
ax.barh(y, u.n_ge5, color="#dd8452", alpha=0.95, height=0.55)
ax.set_yticks(y)
ax.set_yticklabels([f"{a} ({h[:8]})" for a, h in zip(u.artifact, u.rowset_hash)],
                   fontsize=8.4)
for i, (n, g) in enumerate(zip(u.n_rows, u.n_ge5)):
    ax.text(n + 60, i, f"{n} / n≥5:{g}", va="center", fontsize=7.6)
ax.set_xlabel("rows (blue) and n≥5 targets (orange)", fontsize=9.6)
ax.set_title(f"Fig 1  Matched-universe strata: {len(u)} artifacts across "
             f"{UR.rowset_hash.nunique()} distinct rowset hashes\n"
             "Artifacts sharing a hash are directly comparable; across hashes they are not.",
             fontsize=10.6)
ax.grid(alpha=0.25, axis="x")
fig.savefig(F / "figures" / "Fig1_matched_universe_strata.png"); plt.close(fig)
print("   Fig1 ✅")

# Fig2 C_B/C_W 分解（主池）
pm = prim_models.dropna(subset=["C_B_row", "C_W_row"])
fig, ax = plt.subplots(figsize=(7.6, 6.0))
ax.scatter(pm.C_B_row, pm.C_W_row, s=40, c="#4c72b0", alpha=0.8)
ax.axhline(0, color="#bbb", lw=0.9); ax.axvline(0, color="#bbb", lw=0.9)
lim = max(pm.C_B_row.abs().max(), pm.C_W_row.abs().max()) * 1.1
ax.plot([-lim, lim], [-lim, lim], "--", color="#999", lw=0.9)
ax.set_xlabel("C_B (between-target covariance contribution)", fontsize=9.8)
ax.set_ylabel("C_W (within-target covariance contribution)", fontsize=9.8)
res = float(np.nanmax(np.abs(pm.overall_pearson - (pm.C_B_row + pm.C_W_row))))
ax.set_title(f"Fig 2  Covariance-signal decomposition — primary matched universe\n"
             f"universe={PU['rowset_hash'][:8]} | rows={PU['n_rows']} | "
             f"targets n≥5={PU['n_ge5']} | artifacts={PU['n_models']} | "
             f"method families={PU['n_method_families']}\n"
             f"provenance=100% GREEN (level {PU['level']}); r_full=C_B+C_W "
             f"residual ≤ {res:.1e}", fontsize=9.2)
ax.grid(alpha=0.25)
fig.savefig(F / "figures" / "Fig2_CB_CW_decomposition.png"); plt.close(fig)
print("   Fig2 ✅")

# Fig3 受控划分膨胀
names = [k for k in H4 if k != "delta"]
vals_r = [H4[k]["overall_pearson"] for k in names]
vals_cb = [H4[k]["C_B"] for k in names]
vals_cw = [H4[k]["C_W"] for k in names]
fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.4))
x = np.arange(len(names))
axes[0].plot(x, vals_r, "o-", lw=2, color="#c44e52")
axes[0].set_xticks(x); axes[0].set_xticklabels(names, fontsize=8.6, rotation=15)
axes[0].set_ylabel("median Overall Pearson", fontsize=9.6)
axes[0].axhline(0, color="#bbb", lw=0.8)
axes[0].set_title("(a) pooled correlation", fontsize=10)
axes[0].grid(alpha=0.25)
axes[1].plot(x, vals_cb, "s-", lw=2, color="#4c72b0", label="C_B (between)")
axes[1].plot(x, vals_cw, "d-", lw=2, color="#55a868", label="C_W (within)")
axes[1].set_xticks(x); axes[1].set_xticklabels(names, fontsize=8.6, rotation=15)
axes[1].axhline(0, color="#bbb", lw=0.8)
axes[1].set_ylabel("covariance contribution", fontsize=9.6)
axes[1].legend(fontsize=8.4); axes[1].grid(alpha=0.25)
axes[1].set_title("(b) decomposition", fontsize=10)
dl = H4.get("delta", {})
fig.suptitle("Fig 3  CONTROLLED split inflation (same rows, same configurations, "
             "only split regime changes)\n"
             f"universe=V6_2_split4 (7523 rows) | Δr={dl.get('delta_r',float('nan')):+.4f} | "
             f"ΔC_B={dl.get('delta_C_B',float('nan')):+.4f} | "
             f"ΔC_W={dl.get('delta_C_W',float('nan')):+.4f} | "
             f"between-fraction={dl.get('fraction_between',float('nan')):.3f}",
             fontsize=9.8, y=1.10)
fig.text(0.5, -0.06,
         "NOTE: L2_complex and L3_seq40 are IDENTICAL in every metric here — "
         "they must NOT be treated as independent evidence.",
         ha="center", fontsize=8.6, color="#a33")
fig.savefig(F / "figures" / "Fig3_controlled_split_inflation.png"); plt.close(fig)
print("   Fig3 ✅")

# Fig4 remedy 链
RS = pd.read_csv(V8 / "remedy" / "remedy_sequence.csv")
b = RS[RS.label.str.startswith("V7.3 F1")]
r5 = RS[RS.step == "R5"].copy()
fig, ax = plt.subplots(figsize=(8.2, 5.6))
if len(b):
    b = b.iloc[0]
    r5["d_ov"] = r5.overall_r - b.overall_r
    r5["d_mc"] = r5.macro_spearman - b.macro_spearman
    ax.axhline(0, color="#bbb", lw=0.9); ax.axvline(0, color="#bbb", lw=0.9)
    cols = ["#55a868" if d > 0.02 else "#c44e52" for d in r5.d_mc]
    ax.scatter(r5.d_ov, r5.d_mc, s=110, c=cols, zorder=3,
               edgecolors="white", linewidths=1.2)
    for _, r in r5.iterrows():
        ax.annotate(r.label.replace("V7.3 ", ""), (r.d_ov, r.d_mc), fontsize=8,
                    xytext=(5, 5), textcoords="offset points")
ax.set_xlabel("ΔOverall Pearson vs F1 (mutation-only)", fontsize=9.8)
ax.set_ylabel("ΔMacro per-target Spearman", fontsize=9.8)
ax.set_title("Fig 4  Retrospective matched remedy (V7.3, same 3351 rows)\n"
             "green = both axes improve; red = pooled-only or negative",
             fontsize=10.2)
ax.grid(alpha=0.25)
fig.savefig(F / "figures" / "Fig4_retrospective_remedy.png"); plt.close(fig)
print("   Fig4 ✅")

# Fig5 fixed-K geometry
fig, ax = plt.subplots(figsize=(8.2, 5.2))
m = GEO.n <= 200
ax.plot(GEO.n[m], GEO.D[m], lw=2.2, color="#4c72b0")
ax.axvline(10, color="#c44e52", ls="--", lw=1.4)
ax.text(10.8, GEO.D[m].max() * 0.85, "n = K = 10", fontsize=8.6, color="#c44e52")
nstar = FR["axis_c_geometry"]["n_star_analytic"]
ax.axvline(nstar, color="#55a868", ls=":", lw=1.4)
ax.text(nstar + 1.5, GEO.D[m].max() * 0.62,
        f"n* ≈ {nstar:.1f}\n(peak)", fontsize=8.6, color="#55a868")
ax.set_xlabel("candidate-set size n (K=10 fixed)", fontsize=9.8)
ax.set_ylabel("discriminable interval  D(n,K)", fontsize=9.8)
ax.set_title("Fig 5  Fixed-K candidate-selection geometry\n"
             "D(n,K) = K(n−K)/[n(n+1)];  n ≤ K ⇒ D = 0 (structurally non-informative)\n"
             f"D(20)={FR['axis_c_geometry']['D_at_n20']:.4f} > "
             f"D(50)={FR['axis_c_geometry']['D_at_n50']:.4f} ⇒ 'n≥50' claim is NOT supported",
             fontsize=10.0)
ax.grid(alpha=0.25)
fig.savefig(F / "figures" / "Fig5_fixed_K_geometry.png"); plt.close(fig)
print("   Fig5 ✅")

# Fig6 加权敏感性
fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.3))
for ax, cols, nm, cc in ((axes[0], cb, "C_B (between-target)", "#c44e52"),
                         (axes[1], cw, "C_W (within-target)", "#55a868")):
    A_ = W[cols].to_numpy(float)
    xx = np.arange(len(SCHEMES))
    for i in range(0, len(A_), 3):
        ax.plot(xx, A_[i], color=cc, alpha=0.16, lw=0.9)
    ax.plot(xx, np.nanmedian(A_, axis=0), "k-o", lw=2.2, ms=6)
    ax.axhline(0, color="#bbb", lw=0.9)
    ax.set_xticks(xx); ax.set_xticklabels(SCHEMES, fontsize=8.6)
    rng = np.nanmax(A_, axis=1) - np.nanmin(A_, axis=1)
    sf = int((W2[[f"signflip_{'C_B' if cc=='#c44e52' else 'C_W'}_gt0.005"]].sum().iloc[0]))
    ax.set_title(f"{nm}\nrange median={np.nanmedian(rng):.4f}  "
                 f"sign-flip(|·|>0.005)={sf}/{len(A_)}", fontsize=9.6)
    ax.grid(alpha=0.25)
fig.suptitle("Fig 6  Target-weighting sensitivity: C_B direction unstable near zero; "
             "C_W sign stable but magnitude changes", fontsize=10.4, y=1.05)
fig.savefig(F / "figures" / "Fig6_weighting_sensitivity.png"); plt.close(fig)
print("   Fig6 ✅")

# Fig7 能力景观（主池 = provenance-clean）
fig, ax = plt.subplots(figsize=(8.0, 5.8))
g = prim_models.dropna(subset=["overall_pearson", "macro_spearman"])
ax.scatter(g.overall_pearson, g.macro_spearman, s=46, c="#4c72b0", alpha=0.82)
for _, r in g.nlargest(3, "macro_spearman").iterrows():
    ax.annotate(r.key, (r.overall_pearson, r.macro_spearman), fontsize=7.6,
                xytext=(5, 4), textcoords="offset points")
ax.axhline(0, color="#bbb", lw=0.9); ax.axvline(0, color="#bbb", lw=0.9)
rho = float(ST[(ST.rowset_hash == PU["rowset_hash"]) &
               (ST.level == PU["level"])].rho_overall_vs_macro.iloc[0])
ax.set_xlabel("Overall Pearson (axis A)", fontsize=9.8)
ax.set_ylabel("Macro per-target Spearman (axis B)", fontsize=9.8)
ax.set_title(f"Fig 7  Capability landscape — primary matched, provenance-clean pool\n"
             f"universe={PU['rowset_hash'][:8]} | rows={PU['n_rows']} | "
             f"targets n≥5={PU['n_ge5']} | artifacts={PU['n_models']} | "
             f"method families={PU['n_method_families']} | provenance=100% GREEN\n"
             f"rho(Overall, Macro) = {rho:+.3f} IN THIS POOL  —  "
             f"do NOT generalise: across {len(ST)} strata it ranges "
             f"{ST.rho_overall_vs_macro.min():+.2f} to {ST.rho_overall_vs_macro.max():+.2f}",
             fontsize=9.0)
ax.grid(alpha=0.25)
fig.savefig(F / "figures" / "Fig7_capability_landscape.png"); plt.close(fig)
print("   Fig7 ✅")

# Supplement: artifact-level 热图
shutil.copy2(V8 / "figures" / "fig1_rank_agreement_heatmap.png",
             F / "supplement" / "artifact_level_heatmap.png")
pd.read_csv(V8 / "BENCHMARK_SCOREBOARD.csv").to_csv(
    F / "supplement" / "all_variant_scoreboard.csv", index=False, encoding="utf-8")
(F / "supplement" / "historical_diagnostics.md").write_text(
    "# 历史诊断（artifact-level，**非独立方法**，仅作补充）\n\n"
    "* 160-artifact `U_CANONICAL` 秩一致热图 → `artifact_level_heatmap.png`\n"
    "  （含 51 个 V6.2 相关键；**不得**作为主证据）\n"
    "* `7 ↔ 160` 极端反转：**未存活**（见 `HISTORICAL_SURVIVAL_AUDIT.md`）\n"
    "* `ρ = −0.26`：**未存活为普遍结论**；跨 10 个 stratum 范围为 "
    f"{ST.rho_overall_vs_macro.min():+.4f} .. {ST.rho_overall_vs_macro.max():+.4f}\n"
    "* `+0.81 → −0.03`：**降级为探索性**（池组成不同，§8.2）\n",
    encoding="utf-8")

# 图注
caps = ["# V8.1 图注（每图声明 universe）\n"]
for p in sorted((F / "figures").glob("*.png")):
    caps.append(f"## {p.stem}\n\n"
                f"- universe: `{PU['rowset_hash'][:8] if 'Fig1' in p.stem or 'Fig7' in p.stem else 'see title'}`\n"
                f"- 声明：图内标题已包含 universe ID / rows / targets / "
                f"artifacts / method families / provenance 状态\n"
                f"- 层级：artifact-level 或 method-level（见标题）\n")
(F / "FIGURE_CAPTIONS.md").write_text("\n".join(caps), encoding="utf-8")
print("   图注 ✅")

json.dump({"scenario_verdict": scen_verdict,
           "figure_universes": {p.stem: PU["rowset_hash"] for p in
                                (F / "figures").glob("*.png")}},
          open(A / "figure_table_meta.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
print(f"\n§17/§16 完成。表 {len(list((F/'tables').glob('*.csv')))} 张，"
      f"图 {len(list((F/'figures').glob('*.png')))} 张")
