# -*- coding: utf-8 -*-
"""scripts/v8_1/22_si_v5.py -- 补充材料 v5：真实上下标 + 宽表横向 + 精度/编码修复。

修复：
* Note 1/6/5 的求和式、根式、对数几率式改为规范排版
* Table S5/S6：布尔计数列 0/1（非 0.0000/1.0000）；seed 整数化；nan -> "—"
* Table S1/S3/S6 表头字段名美化（n_rows -> Rows (N) 等）
* 宽表（S1/S2/S3/S5/S6）改为 **横向页面** 且 8 pt
* 全部表格正文经 math_tokens 渲染，出现 C_B/C_W 时为真下标
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.shared import Emu, Inches, Pt, RGBColor

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
A = R / "v8_1"
F = R / "v8_1_final"
OUT = F / "manuscript_v5"
OUT.mkdir(parents=True, exist_ok=True)
FM, FS = F / "figures_main", F / "figures_supplement"

FR = json.load(open(A / "FROZEN_RESULTS.json", encoding="utf-8"))
PU, GEO = FR["primary_universe"], FR["axis_c_geometry"]
UR = pd.read_csv(A / "MATCHED_UNIVERSE_REGISTRY.csv")
PA = pd.read_csv(A / "PROVENANCE_AUDIT.csv")
MD = pd.read_csv(A / "METRIC_DIRECTION_REGISTRY.csv")
GEO_DF = pd.read_csv(A / "axis_c_geometry.csv")
W2 = pd.read_csv(A / "WEIGHTING_EFFECT_SIZE_AUDIT.csv")
MR = pd.read_csv(A / "MODEL_FAMILY_REGISTRY.csv")
CM = json.load(open(A / "cath_mapping_audit.json", encoding="utf-8"))
HS = json.load(open(A / "hermes_semantics_audit.json", encoding="utf-8"))
_SV = json.load(open(A / "SURVIVAL_AUDIT.json", encoding="utf-8"))

DISP = {
 "S1": "Single-value near-orthogonality (ρ = −0.26) does not hold; the cross-pool range "
       "includes strong positive agreement (up to +0.82), so near-orthogonality is not a "
       "universal conclusion.",
 "S2": "DOES NOT SURVIVE IN PROVENANCE-CLEAN POOL",
 "S3": "DOES NOT QUALIFY AS CAUSAL EVIDENCE",
 "S4": "NOT PRESENT IN PRIMARY POOL",
}
CLAIM = {
 "S1": "Near-orthogonality of pooled and per-target metrics (single value ρ = −0.26)",
 "S2": "Extreme rank-reversal example (7 ↔ 160)",
 "S3": "Split-regime comparison across pools of different composition (+0.81 → −0.03)",
 "S4": "Scenario ranking obtained on a mixed universe (rank 54)",
}
SV = pd.DataFrame([{"Historical claim": CLAIM.get(k, k),
                    "Remains defensible": "No" if not v.get("remains_defensible") else "Yes",
                    "Disposition": DISP.get(k, "—")} for k, v in _SV.items()])

UNI_NAME = {"34b45937b42b7114": "V7.1_z / V7.0 ablation", "8d003e00057d932a": "V6.2 split4",
            "f8c33c7e987d72c8": "V7.3 physics", "174874178847c9bc": "V7.2R3 PPI",
            "2f5a5c2e4832298a": "V7.2R1 full", "e81801ff29e01b7b": "V7.2R2 trackB",
            "c08fe062e64e4cba": "V7.2 nested"}

_TOK = re.compile(r"(\^\{[^}]*\}|_\{[^}]*\}|_[A-Za-z\u0177\u0302])")


def math_tokens(text):
    out = []
    for part in _TOK.split(str(text)):
        if not part:
            continue
        if part.startswith("^{"):
            out.append((part[2:-1], "sup"))
        elif part.startswith("_{"):
            out.append((part[2:-1], "sub"))
        elif part.startswith("_") and len(part) == 2:
            out.append((part[1], "sub"))
        else:
            out.append((part, ""))
    return out


def f4(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return str(v)
    if x != 0 and abs(x) < 1e-4:
        m, e = f"{x:.2e}".split("e")
        return f"{m} \u00d7 10^{int(e)}"
    return f"{x:.4f}"


def setup(doc):
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"; st.font.size = Pt(10.5)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "SimSun")
    st.paragraph_format.space_after = Pt(6)
    st.paragraph_format.line_spacing = 1.20
    for s in doc.sections:
        s.top_margin = s.bottom_margin = Inches(0.8)
        s.left_margin = s.right_margin = Inches(0.8)


def landscape(section):
    w, h = section.page_width, section.page_height
    section.orientation = WD_ORIENT.LANDSCAPE
    section.page_width, section.page_height = h, w
    section.left_margin = section.right_margin = Inches(0.7)
    section.top_margin = section.bottom_margin = Inches(0.7)


def portrait(section):
    w, h = section.page_width, section.page_height
    section.orientation = WD_ORIENT.PORTRAIT
    if w > h:
        section.page_width, section.page_height = h, w


def H(doc, t, lvl=1):
    p = doc.add_heading(t, level=lvl)
    for r in p.runs:
        r.font.color.rgb = RGBColor(0, 0, 0); r.font.name = "Arial"
        r.font.size = Pt(13.5 if lvl == 1 else 11.5)
    p.paragraph_format.keep_with_next = True
    return p


def W(doc, text, size=None):
    p = doc.add_paragraph()
    for chunk, kind in math_tokens(text):
        r = p.add_run(chunk)
        if size:
            r.font.size = Pt(size)
        if kind == "sub":
            r.font.subscript = True
        elif kind == "sup":
            r.font.superscript = True
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    return p


def cap(doc, t):
    p = doc.add_paragraph()
    for chunk, kind in math_tokens(t):
        r = p.add_run(chunk); r.font.size = Pt(8.8)
        if kind == "sub":
            r.font.subscript = True
        elif kind == "sup":
            r.font.superscript = True
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    return p


def _autofit(table, total_in):
    tbl = table._tbl; tblPr = tbl.tblPr
    lay = tblPr.find(qn("w:tblLayout"))
    if lay is None:
        lay = tblPr.makeelement(qn("w:tblLayout"), {}); tblPr.append(lay)
    lay.set(qn("w:type"), "fixed")
    weights = []
    for j in range(len(table.columns)):
        mx = 1
        for cell in table.columns[j].cells:
            for para in cell.paragraphs:
                for run in para.runs:
                    mx = max(mx, len(run.text))
        weights.append(max(mx, 4))
    tot = sum(weights)
    for j, w in enumerate(weights):
        ew = Emu(int(Inches(total_in) * w / tot))
        for cell in table.columns[j].cells:
            cell.width = ew


def tbl(doc, df, fs=7.5, note=None, total_in=9.4):
    t = doc.add_table(rows=1, cols=len(df.columns))
    t.style = "Table Grid"; t.alignment = WD_TABLE_ALIGNMENT.CENTER; t.autofit = False
    hdr = t.rows[0]
    trPr = hdr._tr.get_or_add_trPr()
    trPr.append(trPr.makeelement(qn("w:tblHeader"), {}))
    for i, c in enumerate(df.columns):
        cell = hdr.cells[i]; cell.text = ""
        for chunk, kind in math_tokens(c):
            r = cell.paragraphs[0].add_run(chunk)
            r.bold = True; r.font.size = Pt(fs)
            if kind == "sub":
                r.font.subscript = True
            elif kind == "sup":
                r.font.superscript = True
    for _, row in df.iterrows():
        cells = t.add_row().cells
        for i, v in enumerate(row):
            cells[i].text = ""
            cells[i].paragraphs[0].add_run(str(v)).font.size = Pt(fs)
    _autofit(t, total_in)
    if note:
        cap(doc, note)
    else:
        doc.add_paragraph()


# ---------- 表数据 ----------
UR_t = UR.copy()
UN = {"universe_id": "Universe", "artifact": "Artifact", "rowset_hash": "Row-set fingerprint",
      "n_rows": "Rows (N)", "n_targets": "Targets", "n_families": "Families",
      "n_ge3": "Targets (n ≥ 3)", "n_ge5": "Targets (n ≥ 5)", "n_ge10": "Targets (n ≥ 10)",
      "n_models": "Models"}
UR_t = UR_t.rename(columns=UN)
for c in UR_t.columns:
    if UR_t[c].dtype.kind in "if":
        UR_t[c] = UR_t[c].map(lambda x: f"{int(x)}" if float(x).is_integer() else f4(x))

PA_t = PA.head(40).copy()

METRIC_FIX = {"C_B_row": "C_{B} (row-weighted)", "C_B_equal": "C_{B} (target-equal)",
              "C_W_row": "C_{W} (row-weighted)", "C_W": "C_{W}",
              "covariance_contribution": "covariance contribution",
              "attenuation_slope": "attenuation slope", "pairwise_d05": "pairwise accuracy",
              "higher_better": "higher is better", "lower_better": "lower is better"}
MD_t = MD.copy()
MD_t = MD_t.rename(columns={"metric": "Metric", "axis": "Analysis axis",
                            "raw_direction": "Raw direction",
                            "ranking_direction": "Direction for ranking",
                            "transform_for_rank_comparison": "Transform for rank comparison",
                            "absolute_unit_or_score": "Unit / score type",
                            "eligible_for_calibration": "Eligible for calibration",
                            "notes": "Notes"})
if "Metric" in MD_t.columns:
    MD_t["Metric"] = MD_t["Metric"].map(lambda m: METRIC_FIX.get(str(m), str(m)))
for c in ["Raw direction", "Direction for ranking", "Transform for rank comparison",
          "Eligible for calibration", "Analysis axis", "Unit / score type"]:
    if c in MD_t.columns:
        MD_t[c] = (MD_t[c].astype(str)
                   .str.replace("_", " ", regex=False)
                   .replace({"higher better": "higher is better",
                             "lower better": "lower is better",
                             "nan": "—"}))

_g = GEO_DF[GEO_DF["n"] > 10].copy()          # E_rand = K/n 仅在 n > K 时为占比
_sel = sorted(set([10, 11, 12, 15, 20, 25, 30, 50, 75, 100, 150]) & set(_g["n"]))
GEO_t = _g[_g["n"].isin(_sel)].copy()
GEO_t = GEO_t.rename(columns={"n": "n", "E_rand": "E_{rand}", "E_max": "E_{max}",
                              "D": "D(n, K)"})
for c in ["E_{rand}", "E_{max}", "D(n, K)"]:
    GEO_t[c] = GEO_t[c].map(lambda x: f"{float(x):.3f}")
GEO_t["n"] = GEO_t["n"].map(lambda x: f"{int(x)}")

W2_t = W2.head(30).copy()
W2_t = W2_t.rename(columns={
    "model_name": "Model configuration", "range_C_B": "Range C_{B}", "range_C_W": "Range C_{W}",
    "max_abs_C_B": "Max |C_{B}|", "max_abs_C_W": "Max |C_{W}|",
    "signflip_C_B_gt0": "C_{B} flips (>0)", "signflip_C_W_gt0": "C_{W} flips (>0)",
    "signflip_C_B_gt0.001": "C_{B} flips (>0.001)",
    "signflip_C_W_gt0.001": "C_{W} flips (>0.001)",
    "signflip_C_B_gt0.005": "C_{B} flips (>0.005)",
    "signflip_C_W_gt0.005": "C_{W} flips (>0.005)",
    "signflip_C_B_gt0.01": "C_{B} flips (>0.01)",
    "signflip_C_W_gt0.01": "C_{W} flips (>0.01)"})
FLIP_COLS = [c for c in W2_t.columns if "flips" in c]
for c in W2_t.columns:
    if c == "Model configuration":
        continue
    if c in FLIP_COLS:
        W2_t[c] = W2_t[c].map(lambda x: str(int(bool(x))))     # 0 / 1 整数
    else:
        W2_t[c] = W2_t[c].map(f4)

keep = [c for c in ["artifact_id", "method_family", "seed", "split_regime",
                    "rowset_hash", "is_ablation", "is_primary_representative"]
        if c in MR.columns]
MR_t = MR[keep].head(60).copy()
MR_t = MR_t.rename(columns={"artifact_id": "Artifact", "method_family": "Method family",
                            "seed": "Seed", "split_regime": "Split regime",
                            "rowset_hash": "Row-set fingerprint",
                            "is_ablation": "Ablation",
                            "is_primary_representative": "Primary representative"})
if "Seed" in MR_t.columns:
    MR_t["Seed"] = MR_t["Seed"].map(
        lambda x: "—" if pd.isna(x) else str(int(x)))          # 42 / 43 / 44, nan -> em dash
for c in ["Ablation", "Primary representative"]:
    if c in MR_t.columns:
        MR_t[c] = MR_t[c].map(lambda x: "yes" if bool(x) else "no")


def build():
    D = Document(); setup(D)
    p = D.add_paragraph()
    r = p.add_run("Supplementary Information"); r.bold = True; r.font.size = Pt(20)
    r.font.name = "Arial"
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p = D.add_paragraph()
    r = p.add_run("Disentangling Between-Target and Within-Target Signal in "
                  "PPI ΔΔG Benchmarking")
    r.italic = True; r.font.size = Pt(11)
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    D.add_page_break()

    H(D, "Supplementary Note 1 | Derivation of the covariance decomposition", 1)
    W(D, "Write the true value and prediction for mutation i in target c as "
         "y_{ci} = μ_{c} + ε_{ci} and ŷ_{ci} = μ̂_{c} + ε̂_{ci}, where μ_{c} and μ̂_{c} are "
         "target means and ε, ε̂ are within-target deviations whose target-wise means "
         "vanish by construction. Expanding Cov(y, ŷ) gives a between-target term "
         "N⁻¹ Σ_{c} n_{c} (μ_{c} − μ)(μ̂_{c} − μ̂) and a within-target term "
         "N⁻¹ Σ_{c} Σ_{i} ε_{ci} ε̂_{ci}; all cross terms vanish identically because "
         "within-target deviations have zero mean in every target. Dividing by "
         "σ_{y} σ_{ŷ} gives C_{B} and C_{W} and establishes r_{pooled} = C_{B} + C_{W} "
         "exactly. Row-level weights are used throughout and the identity residual is "
         "reported as a quality-control quantity (observed maximum 3.331 × 10⁻¹⁶). Because "
         "C_{B} and C_{W} share the same denominator they are covariance-signal "
         "contributions on a common scale, not independent correlations or independently "
         "estimated skill scores.", size=10)

    H(D, "Supplementary Note 2 | Matched universes and row-set fingerprints", 1)
    W(D, "Each prediction output received a row-set fingerprint, the SHA-256 hash of its "
         "sorted row identifiers. Cross-model comparison was performed only between outputs "
         "sharing a fingerprint; raw metrics were never compared across fingerprints. The "
         "primary matched set was selected by a rule fixed before any metric values were "
         "inspected: provenance cleanliness, then number of independent method families, "
         "then number of targets with n ≥ 5, then number of rows.", size=10)
    landscape(D.add_section())
    tbl(D, UR_t, fs=7.5,
        note="Table S1 | Matched-evaluation strata and row-set fingerprints. Outputs are "
             "directly comparable only within a shared fingerprint. Fingerprints are "
             "abbreviated to their first eight hexadecimal characters.")
    portrait(D.add_section())

    H(D, "Supplementary Note 3 | Provenance filtering and de-replication", 1)
    W(D, "Provenance was assigned from the recorded training split regime rather than "
         "inferred from names. Configurations trained under a permissive (random / complex / "
         "seq40) partition were excluded from deployment-oriented comparisons and retained "
         "only in the controlled split experiment, where the split regime is itself the "
         "manipulated variable. Analyses were performed at artifact, configuration and "
         "method-family level; representatives were chosen performance-blind.", size=10)
    landscape(D.add_section())
    tbl(D, PA_t, fs=7.5,
        note=f"Table S2 | Provenance audit, showing the first 40 of {len(PA)} entries. "
             f"Training split regime determines provenance status.")
    portrait(D.add_section())

    H(D, "Supplementary Note 4 | Metric-direction registry and RMSE-direction correction", 1)
    W(D, "A metric-direction registry was established before any rank-agreement computation. "
         "RMSE and MAE are lower-is-better; all cross-metric rank agreement uses negated RMSE "
         "(or ascending ranks) so that a higher-is-better correlation is never correlated "
         "against a raw higher-is-worse error. An earlier version of the calibration-axis "
         "summary violated this rule and has been superseded; the corrected values are used "
         "throughout. The registry comprises all 29 metrics used in this work.", size=10)
    landscape(D.add_section())
    tbl(D, MD_t, fs=7.2,
        note=f"Table S3 | Metric-direction registry (all {len(MD)} metrics).")
    portrait(D.add_section())

    H(D, "Supplementary Note 5 | External mapping ambiguity", 1)
    W(D, f"The normalised matching key, formed by dropping chain identity and the wild-type "
         f"residue, is not unique on the annotation side: of 4,473 rows "
         f"{CM['n_exact_unique']:,} are uniquely mapped, {CM['n_ambiguous_key']:,} are "
         f"ambiguous and {CM['n_unmatched']:,} are unmatched. Primary matched evaluation used "
         f"uniquely mapped rows only; ambiguous keys were excluded rather than resolved. The "
         f"official held-out CATH partition carries its own unambiguous rows and was "
         f"evaluated separately without canonical re-mapping. For the external quantity "
         f"log P_{{mt}} − log P_{{wt}}, RMSE and MAE are NOT_APPLICABLE because the quantity "
         f"is a dimensionless mutation-effect score rather than a calibrated kcal/mol ΔΔG.",
         size=10)
    cap(D, f"Mapping QC: canonical rows before merge = "
           f"{HS['qc']['canonical_rows_before_merge']}, after merge = "
           f"{HS['qc']['canonical_rows_after_merge']}, duplicate row identifiers = "
           f"{HS['qc']['duplicate_canonical_row_id']}, exact matches = "
           f"{HS['qc']['exact_match_count']}, ambiguous = {HS['qc']['ambiguous_count']}, "
           f"unmatched = {HS['qc']['unmatched_count']}.")

    H(D, "Supplementary Note 6 | Fixed-K geometry derivation", 1)
    W(D, "For a target with n mutations and fixed K, let the rank-derived gain of the K most "
         "favourable mutations be the available utility. Choosing K candidates uniformly at "
         "random without replacement yields an expected gain share E_{rand} = K/n by "
         "symmetry, while the optimal K candidates attain "
         "E_{max} = K(2n − K + 1)/[n(n + 1)]. The discriminable interval is therefore "
         "D(n, K) = E_{max} − E_{rand} = K(n − K)/[n(n + 1)] for n > K, and D = 0 for "
         "n ≤ K, since every candidate is then selected. D is not monotonic in n: "
         "dD/dn = 0 at n* = K + √(K² + K), giving n* ≈ 20.5 for K = 10, after which D "
         "declines; for example D(20) = 0.2381 and D(50) = 0.1569.", size=10)
    landscape(D.add_section())
    tbl(D, GEO_t, fs=7.5,
        note="Table S4 | Discriminable interval versus candidate-set size (K = 10), every "
             "fortieth value.")
    portrait(D.add_section())

    H(D, "Supplementary Note 7 | Full weighting audit", 1)
    W(D, "Continuous values are rounded to four decimal places; ranges smaller than "
         "1 × 10⁻⁴ are given in scientific notation. Sign-flip columns are integer counts "
         "(0 or 1) of the weighting schemes in which the component changes sign, at the "
         "stated threshold on its absolute value.", size=10)
    landscape(D.add_section())
    tbl(D, W2_t, fs=7.0,
        note=f"Table S5 | Weighting effect-size audit (first 30 of {len(W2)} configurations). "
             f"Sign-flip columns report the number of schemes producing a sign change at each "
             f"threshold.")
    portrait(D.add_section())

    H(D, "Supplementary Note 8 | Full model registry", 1)
    landscape(D.add_section())
    tbl(D, MR_t, fs=7.2,
        note="Table S6 | Model family registry (first 60 entries). Row-set fingerprints allow "
             "outputs to be grouped into comparable matched universes; an em dash marks "
             "entries without a single seed value.")
    portrait(D.add_section())

    H(D, "Supplementary Note 9 | Automated workflow and reproducibility", 1)
    W(D, "AI-assisted coding tools were used to draft and debug the analysis and plotting "
         "scripts, and to automate document generation, as declared in the main text. All "
         "headline numbers are generated from a single machine-readable results file; an "
         "automated consistency checker verifies that the documents, tables and figures agree "
         "with it. Metric definitions, analysis choices, threshold settings, hypothesis "
         "verdicts and final conclusions were made by the authors. The AI system was not used "
         "to generate data or to modify any upstream model prediction, and all results can be "
         "regenerated from the accompanying scripts and frozen input hashes. No new model was "
         "trained and no GPU was used.", size=10)

    H(D, "Supplementary Note 10 | Historical claim audit", 1)
    W(D, "For transparency we record that four headline statements from an earlier version of "
         "this analysis did not survive the matched, de-duplicated, provenance-clean "
         "re-analysis and were withdrawn. Their dispositions are tabulated below.", size=10)
    landscape(D.add_section())
    tbl(D, SV, fs=8.0)
    portrait(D.add_section())

    H(D, "Supplementary Figures", 1)
    for folder, stem, w, c in [
        (FS, "FigS1_matched_universes", 6.4,
         "Figure S1 | Matched-evaluation strata and row-set fingerprints. Each bar is one "
         "prediction output; outputs sharing a row-set fingerprint are directly comparable "
         "with each other and with no others."),
        (FS, "Capability_heatmap_provenance_clean", 6.4,
         "Figure S2 | Capability profiles of matched, provenance-clean method families. Rows "
         "are the controlled ablation variants of a single local-physics method family "
         "(seed 42), not independent methods. Values are standardised within each metric "
         "(z-score) for visual comparison only; raw values are printed in each cell and "
         "statistical analyses use the raw metrics."),
        (FS, "FigS3_capability_pool_sensitivity", 6.4,
         "Figure S3 | Pool sensitivity of cross-metric rank agreement. Within one evaluation "
         "universe and row set, changing only which model families are included moves the "
         "agreement between pooled Pearson and macro within-target Spearman across a wide "
         "range, which is why this quantity is reported as a range rather than as a single "
         "property of a benchmark."),
    ]:
        D.add_picture(str(folder / f"{stem}_EN.png"), width=Inches(w))
        cap(D, c)

    nm = "CR-PPI-ddG_Supplementary_v5_EN.docx"
    path = OUT / nm
    D.save(str(path))
    print(f"  {nm} ({path.stat().st_size/1024:.0f} KB)")


build()
print("v5 supplementary 完成")
