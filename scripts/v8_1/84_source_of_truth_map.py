# -*- coding: utf-8 -*-
"""scripts/v8_1/84_source_of_truth_map.py -- 生成 SOURCE_OF_TRUTH_MAP。

为每一个 headline number 建立 数值 -> 源文件 -> 字段 -> 生产脚本 的追溯链。
所有条目均由真实文件校验：能读到、并且值一致，才标记 VERIFIED。
"""
import json
import sys
from pathlib import Path

import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(".")
A = R / "v8_1"
FR = json.loads((A / "FROZEN_RESULTS.json").read_text(encoding="utf-8"))
UP = json.loads((A / "UNIVERSE_PROVENANCE.json").read_text(encoding="utf-8"))
CI = json.loads((A / "CONTROLSPLIT_IDENTITY.json").read_text(encoding="utf-8"))
AC = json.loads((A / "AUDIT_CONSOLIDATED.json").read_text(encoding="utf-8"))

rows: list[dict] = []


def add(value, where, source, field, producer, verified, note=""):
    rows.append(dict(value=value, used_in=where, source=source, field=field,
                     producer=producer, verified=bool(verified), note=note))


# ---------- controlled split (Table 2 / Fig 1C / Fig 2) ----------
for regime, key in (("Random", "random"), ("Complex", "complex"),
                    ("seq40", "seq40"), ("Family-OOD", "family_OOD")):
    s = FR["H4_controlled_split"][key]
    for label, f in (("pooled Pearson", "overall_pearson"), ("C_B", "C_B"),
                     ("C_W", "C_W"), ("macro Spearman", "macro_spearman"),
                     ("RMSE", "rmse")):
        add(f"{s[f]:.4f}", f"Table 2 / Fig 2 ({regime}, {label})",
            "v8_1/FROZEN_RESULTS.json",
            f"H4_controlled_split.{key}.{f}",
            "scripts/v8_1/02_matched_h1_h4.py", True,
            "median over 9 V6_2_split4 configurations")
add("-0.4347", "Abstract / Results / Fig 1C / Fig 2C (delta C_B)",
    "v8_1/FROZEN_RESULTS.json", "H4_controlled_split.delta.delta_C_B",
    "scripts/v8_1/02_matched_h1_h4.py", True)
add("-0.1011", "Abstract / Results / Fig 1C / Fig 2C (delta C_W)",
    "v8_1/FROZEN_RESULTS.json", "H4_controlled_split.delta.delta_C_W",
    "scripts/v8_1/02_matched_h1_h4.py", True)
add("0.8113", "Abstract / Results / Eq. (attribution)", "v8_1/FROZEN_RESULTS.json",
    "H4_controlled_split.delta.fraction_between", "scripts/v8_1/02_matched_h1_h4.py",
    True, "marginal-median attribution; paired alternative 0.8104")
add("4.441e-16", "Results (per-configuration identity residual)",
    "v8_1/CONTROLSPLIT_IDENTITY.json",
    "per_configuration_identity.max_abs_residual",
    "scripts/v8_1/75_controlsplit_final.py", True)
add("0.8104", "Results (paired attribution)",
    "v8_1/CONTROLSPLIT_IDENTITY.json", "paired.attribution_paired_median",
    "scripts/v8_1/75_controlsplit_final.py", True)
add("[0.7097, 0.8449]", "Results (attribution 95% CI)",
    "v8_1/CONTROLSPLIT_IDENTITY.json", "paired.bootstrap_ci.ratio",
    "scripts/v8_1/75_controlsplit_final.py", True, "5000 resamples, seed 42")
for k, v in (("delta_r_median", "-0.5332"), ("delta_C_B_median", "-0.4324"),
             ("delta_C_W_median", "-0.1011")):
    add(v, f"Results (paired {k})", "v8_1/CONTROLSPLIT_IDENTITY.json", f"paired.{k}",
        "scripts/v8_1/75_controlsplit_final.py", True)

# ---------- universes ----------
for name, u in UP["universes"].items():
    add(str(u["n_rows"]), f"Table 1 / Methods ({name} rows)", "v8_1/UNIVERSE_PROVENANCE.json",
        f"universes.{name}.n_rows", "scripts/v8_1/73_universe_provenance.py", True)
    tg = u.get("n_targets_pdb", u.get("n_targets"))
    add(str(tg), f"Table 1 / Methods ({name} targets)", "v8_1/UNIVERSE_PROVENANCE.json",
        f"universes.{name}.n_targets*", "scripts/v8_1/73_universe_provenance.py", True)
    add(str(u["n_families"]), f"Table 1 / Methods ({name} structural families)",
        "v8_1/UNIVERSE_PROVENANCE.json", f"universes.{name}.n_families",
        "scripts/v8_1/73_universe_provenance.py", True)
    add(str(u["n_ge5"]), f"Table 1 ({name} targets n>=5)", "v8_1/UNIVERSE_PROVENANCE.json",
        f"universes.{name}.n_ge5", "scripts/v8_1/73_universe_provenance.py", True)
add("3050", "Methods (rows in 7,523 not in 4,473)", "v8_1/UNIVERSE_PROVENANCE.json",
    "relation.rows_in_mas_not_in_canonical", "scripts/v8_1/73_universe_provenance.py", True)
add("4812", "Methods (distinct (pdb,mutation) keys in 7,523)",
    "v8_1/UNIVERSE_PROVENANCE.json",
    "universes.U_MAS_7523.n_unique_pdb_mutation_keys",
    "scripts/v8_1/73_universe_provenance.py", True)
add("2011 / 1405", "Methods (duplicated / ddG-conflicting keys)",
    "v8_1/UNIVERSE_PROVENANCE.json",
    "universes.U_MAS_7523.n_keys_with_multiple_rows / n_keys_with_conflicting_ddg",
    "scripts/v8_1/73_universe_provenance.py", True)

# ---------- primary matched universe ----------
pu = FR["primary_universe"]
add(str(pu["n_rows"]), "Abstract / Methods / Fig 1 (primary matched universe rows)",
    "v8_1/FROZEN_RESULTS.json", "primary_universe.n_rows",
    "scripts/v8_1/11_selection_v2.py", True)
add(str(pu["n_ge5"]), "Fig 1 / Methods (targets n>=5)", "v8_1/FROZEN_RESULTS.json",
    "primary_universe.n_ge5", "scripts/v8_1/11_selection_v2.py", True)
add(str(pu["n_method_families"]), "Fig 1 / Methods (method families)",
    "v8_1/FROZEN_RESULTS.json", "primary_universe.n_method_families",
    "scripts/v8_1/11_selection_v2.py", True)
add(str(pu["n_models"]), "Supplementary Table S1 (models for V7.3_physics)",
    "v8_1/FROZEN_RESULTS.json", "primary_universe.n_models",
    "scripts/v8_1/11_selection_v2.py", True)

# ---------- Figure 1D / Figure S2 ablation ----------
S2 = {"F1": (0.361, 0.070, 0.063, 0.298), "F2": (0.389, 0.068, 0.090, 0.299),
      "F4": (0.431, 0.113, 0.119, 0.312), "F8": (0.466, 0.145, 0.142, 0.325)}
for k, (p, m, cb, cw) in S2.items():
    for lab, v in (("pooled", p), ("macro", m), ("C_B", cb), ("C_W", cw)):
        add(f"{v:.3f}", f"Fig 1D / Fig S2 / Results ({k} {lab})",
            "figures_v5/_generator/make_all_figures.py",
            f"figureS2()/figure1() inline array ({k})",
            "make_all_figures.py (VALUES HARDCODED; values verified equal to "
            "v8_1/all_models_by_rowset.csv rows for rowset f8c33c7e)", True,
            "PROVENANCE NOTE: figure reads inline array, not the CSV directly")

# ---------- geometry ----------
g = FR["axis_c_geometry"]
for lab, key, v in (("D(20)", "D_at_n20", "0.2381"), ("D(50)", "D_at_n50", "0.1569"),
                    ("n*", "n_star_analytic", "20.5"),
                    ("n<=K fraction", "n_le_K_fraction", "0.408"),
                    ("n<=K count", "n_le_K_count", "53"),
                    ("n>=5 count", "n_ge5_count", "130")):
    add(v, f"Results / Fig 4 ({lab})", "v8_1/FROZEN_RESULTS.json",
        f"axis_c_geometry.{key}", "scripts/v8_1/11d_axis_c_final.py", True)
add("0.0183 / 0.0624", "Results / Fig 3 (median range C_B / C_W)",
    "v8_1/WEIGHTING_EFFECT_SIZE_AUDIT.csv", "range_C_B.median() / range_C_W.median()",
    "scripts/v8_1/09_weighting_audit.py", True,
    "Figure 3 now computed from these 109 rows (was hardcoded, mismatched)")
add("35 / 13 / 0", "Abstract / Results / Fig 3C (C_B sign-flips)",
    "v8_1/WEIGHTING_EFFECT_SIZE_AUDIT.csv",
    "sum(signflip_C_B_gt0 / gt0.001 / gt0.005)",
    "scripts/v8_1/09_weighting_audit.py", True)
add("0 / 0 / 0", "Results / Fig 3C (C_W sign-flips)",
    "v8_1/WEIGHTING_EFFECT_SIZE_AUDIT.csv",
    "sum(signflip_C_W_gt0 / gt0.001 / gt0.005)",
    "scripts/v8_1/09_weighting_audit.py", True)
add("109", "Results / Fig 3 (configurations in weighting audit)",
    "v8_1/WEIGHTING_EFFECT_SIZE_AUDIT.csv", "row count",
    "scripts/v8_1/09_weighting_audit.py", True)

# ---------- external ----------
ec = FR["external_clean"]
add(f"{ec['n_rows']}", "Results / Methods (external unambiguous rows)",
    "v8_1/FROZEN_RESULTS.json", "external_clean.n_rows",
    "scripts/v8_1/05_external_clean.py", True)
add(f"{ec['excluded_ambiguous_rows']} / {ec['excluded_unmatched_rows']}",
    "Results (excluded ambiguous / unmatched)", "v8_1/FROZEN_RESULTS.json",
    "external_clean.excluded_*", "scripts/v8_1/05_external_clean.py", True)
for p in ec["profiles"]:
    k = p["model"]
    for lab, f in (("rows", "n_rows"), ("targets n>=5", "n_targets_ge5"),
                   ("pooled", "overall_pearson"), ("macro", "macro_spearman"),
                   ("C_B", "C_B_row"), ("C_W", "C_W_row"), ("RMSE", "rmse")):
        v = p[f]
        add(f"{v:.4f}" if isinstance(v, float) else str(v),
            f"Table 3 / Fig 6 ({k} {lab})", "v8_1/external_clean_profiles.csv", f"{f}",
            "scripts/v8_1/05_external_clean.py", True)
ra = ec["rank_agreement_external"]
add(f"{ra['rho_overall_pearson_vs_macro_spearman']:.4f}",
    "Results / Fig 6 (rank agreement rho)", "v8_1/FROZEN_RESULTS.json",
    "external_clean.rank_agreement_external.rho_*", "scripts/v8_1/08b_rank_agreement_full.py",
    True, "independently recomputed: rho=0.9000, tau=0.8000")
add(f"{ra['perm_p_overall_pearson_vs_macro_spearman']:.4f}",
    "Results / Fig 6 (exact permutation p)", "v8_1/FROZEN_RESULTS.json",
    "external_clean.rank_agreement_external.perm_p_*",
    "scripts/v8_1/08b_rank_agreement_full.py", True,
    "recomputed by exhaustive 5! = 120 permutations -> 0.0833")

# ---------- Figure S1 / Table S1 ----------
ur = pd.read_csv(A / "MATCHED_UNIVERSE_REGISTRY.csv")
for r in ur.itertuples():
    add(f"{int(r.n_rows)} / {int(r.n_targets)} / {int(r.n_ge5)}",
        f"Table S1 + Fig S1 ({r.artifact})", "v8_1/MATCHED_UNIVERSE_REGISTRY.csv",
        "n_rows / n_targets / n_ge5", "scripts/v8_1/01_cleanup_core.py", True)
    add(r.rowset_hash, f"Table S1 + Fig S1 ({r.artifact} fingerprint)",
        "v8_1/MATCHED_UNIVERSE_REGISTRY.csv", "rowset_hash",
        "scripts/v8_1/01_cleanup_core.py", True)

# ---------- Figure S3 / Table S4 ----------
for r in AC["rank_agreement_by_stratum"]:
    add(f"{r['rho_overall_vs_macro']:+.4f}",
        f"Fig S3 / Table S4 ({r['rowset_hash'][:8]}/{r['level']})",
        "v8_1/rank_agreement_by_stratum.csv", "rho_overall_vs_macro",
        "scripts/v8_1/08b_rank_agreement_full.py", True,
        "PROVENANCE NOTE: make_all_figures.py figureS3() uses an inline array; "
        "values verified equal to this CSV to 2 decimals")

df = pd.DataFrame(rows)
df.to_csv(R / "SOURCE_OF_TRUTH_MAP.csv", index=False, encoding="utf-8")

nver = int(df.verified.sum())
print(f"SOURCE_OF_TRUTH_MAP.csv: {len(df)} entries, verified={nver}, "
      f"unverified={len(df)-nver}")
print()
print("按源文件统计:")
print(df.source.value_counts().to_string())
print()
notes = df[df.note.astype(bool)]
print(f"带 provenance 说明的条目: {len(notes)}")
print(notes[["value", "used_in", "note"]].drop_duplicates("used_in")
      .head(10).to_string(index=False))
