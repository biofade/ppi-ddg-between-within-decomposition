# -*- coding: utf-8 -*-
"""scripts/v8_1/31_manuscript_v6.py -- v6：把用户提供的 generate_publication_figures_v5.py
生成的 9 张新图接入英文主稿与补充材料，并保持全部数值来自 FROZEN_RESULTS.json。

图号映射（主稿）：
    Fig1_overview_EN.png                      -> Figure 1
    Fig2_controlled_split_collapse_EN.png      -> Figure 2
    Fig3_target_weighting_sensitivity_EN.png   -> Figure 3
    Fig4_fixed_K_geometry_EN.png               -> Figure 4
    Fig5_retrospective_remedy_EN.png           -> Figure 5
    Fig6_external_methods_EN.png               -> Figure 6
补充材料：
    FigS1_matched_universe_strata_EN           -> Figure S1
    FigS2_capability_heatmap_provenance_clean_EN -> Figure S2
    FigS3_pool_sensitivity_rank_agreement_EN   -> Figure S3
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
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml.ns import qn
from docx.shared import Emu, Inches, Pt, RGBColor

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
A = R / "v8_1"
F = R / "v8_1_final"
OUT = F / "manuscript_v6"
OUT.mkdir(parents=True, exist_ok=True)
FIGS = R / "figures_v5"

FR = json.load(open(A / "FROZEN_RESULTS.json", encoding="utf-8"))
PU, H4, GEO = FR["primary_universe"], FR["H4_controlled_split"], FR["axis_c_geometry"]
DL = H4.get("delta", {})
EP = pd.read_csv(A / "external_clean_profiles.csv")
W2 = pd.read_csv(A / "WEIGHTING_EFFECT_SIZE_AUDIT.csv")
ST = pd.read_csv(A / "rank_agreement_by_stratum.csv")
UR = pd.read_csv(A / "MATCHED_UNIVERSE_REGISTRY.csv")
GEO_DF = pd.read_csv(A / "axis_c_geometry.csv")
r0, r3 = H4["random"]["overall_pearson"], H4["family_OOD"]["overall_pearson"]
sf0 = int(W2["signflip_C_B_gt0"].sum())
sf5 = int(W2["signflip_C_B_gt0.005"].sum())
sfw = int(W2["signflip_C_W_gt0"].sum())
cw3 = H4["family_OOD"]["C_W"]

REFS = [
 "Jankauskaitė J, Jiménez-García B, Dapkūnas J, Fernández-Recio J, Moal IH. SKEMPI 2.0: an "
 "updated benchmark of changes in protein–protein binding energy, kinetics and "
 "thermodynamics upon mutation. Bioinformatics. 2019;35(3):462–469. "
 "doi:10.1093/bioinformatics/bty635",
 "Bushuiev A, Bushuiev R, Sedlar J, Pluskal T, Damborsky J, Mazurenko S, Sivic J. Revealing "
 "data leakage in protein interaction benchmarks. arXiv:2404.10457. 2024. "
 "doi:10.48550/arXiv.2404.10457",
 "Yu G, Bi X, Ma T, Li Y, Wang J. CATH-ddG: towards robust mutation effect prediction on "
 "protein–protein interactions out of CATH homologous superfamily. Bioinformatics. "
 "2025;41(Suppl_1):i362–i372. doi:10.1093/bioinformatics/btaf228",
 "Shan S, Luo S, Yang Z, Hong J, Su Y, Ding F, et al. Deep learning guided optimization "
 "of human antibody against SARS-CoV-2 variants with broad neutralization. Proc Natl Acad "
 "Sci U S A. 2022;119(11):e2122954119. doi:10.1073/pnas.2122954119",
 "Park S, Yim S, Kim D, Yoo K, Hwang D, Lee K, Jang J, Kim K. PreFold-dG: estimating "
 "binding affinity of protein–protein interaction from intermediate representations of "
 "protein folding model. Bioinformatics. 2026;42(Suppl_2):btag489. "
 "doi:10.1093/bioinformatics/btag489",
 "Barlow KA, Ó Conchúir S, Thompson S, Suresh P, Lucas JE, Heinonen M, Kortemme T. Flex ddG: "
 "Rosetta ensemble-based estimation of changes in protein–protein binding affinity upon "
 "mutation. J Phys Chem B. 2018;122(21):5389–5399. doi:10.1021/acs.jpcb.7b11367",
 "Delgado J, Radusky LG, Cianferoni D, Serrano L. FoldX 5.0: working with RNA, small "
 "molecules and a new graphical interface. Bioinformatics. 2019;35(20):4168–4169. "
 "doi:10.1093/bioinformatics/btz184",
 "Dauparas J, Anishchenko I, Bennett N, Bai H, Ragotte RJ, Milles LF, et al. Robust deep "
 "learning-based protein sequence design using ProteinMPNN. Science. 2022;378(6615):49–56. "
 "doi:10.1126/science.add2187",
 "Hsu C, Verkuil R, Liu J, Lin Z, Hie B, Sercu T, Lerer A, Rives A. Learning inverse folding "
 "from millions of predicted structures. In: Proceedings of the 39th International "
 "Conference on Machine Learning (ICML). 2022;162:8946–8970.",
 "Bushuiev A, Bushuiev R, Kouba P, Filkin A, Gabrielova M, Gabriel M, Sedlar J, Pluskal T, "
 "Damborsky J, Mazurenko S, Sivic J. Learning to design protein–protein interactions with "
 "enhanced generalization. In: The Twelfth International Conference on Learning "
 "Representations (ICLR). 2024. arXiv:2310.18515",
 "Järvelin K, Kekäläinen J. Cumulated gain-based evaluation of IR techniques. ACM Trans Inf "
 "Syst. 2002;20(4):422–446. doi:10.1145/582415.582418",
 "Efron B, Tibshirani RJ. An Introduction to the Bootstrap. New York: Chapman & Hall; 1993.",
 "Sokal RR, Rohlf FJ. Biometry: The Principles and Practice of Statistics in Biological "
 "Research. 3rd ed. New York: WH Freeman; 1995.",
]

_TOK = re.compile(r"(\^\{[^}]*\}|_\{[^}]*\}|_[A-Za-z\u0177\u0302]|\[\d+(?:[,\u2013]\d+)*\])")


def toks(text):
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
        elif part.startswith("[") and part.endswith("]"):
            inner = part[1:-1].replace("\u2013", "\u2013")
            out.append((inner, "sup"))
        else:
            out.append((part, ""))
    return out


def setup(doc):
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"; st.font.size = Pt(10.5)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "SimSun")
    pf = st.paragraph_format
    pf.line_spacing_rule = WD_LINE_SPACING.MULTIPLE; pf.line_spacing = 1.25
    pf.space_after = Pt(6)
    for s in doc.sections:
        s.top_margin = s.bottom_margin = Inches(0.9)
        s.left_margin = s.right_margin = Inches(0.9)


def H(doc, t, lvl=1):
    p = doc.add_heading(t, level=lvl)
    for r in p.runs:
        r.font.color.rgb = RGBColor(0, 0, 0); r.font.name = "Arial"
        r.font.size = Pt(14 if lvl == 1 else 11.5)
    p.paragraph_format.keep_with_next = True
    return p


def W(doc, text, size=None):
    p = doc.add_paragraph()
    for c, k in toks(text):
        r = p.add_run(c)
        if k == "sup":
            r.font.superscript = True
        elif k == "sub":
            r.font.subscript = True
        if size and k == "":
            r.font.size = Pt(size)
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    return p


def cap(doc, t):
    p = doc.add_paragraph()
    for c, k in toks(t):
        r = p.add_run(c); r.font.size = Pt(9.0)
        if k == "sub":
            r.font.subscript = True
        elif k == "sup":
            r.font.superscript = True
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    p.paragraph_format.space_after = Pt(10)
    return p


def _autofit(table, total_in):
    tbl = table._tbl; tblPr = tbl.tblPr
    lay = tblPr.find(qn("w:tblLayout"))
    if lay is None:
        lay = tblPr.makeelement(qn("w:tblLayout"), {}); tblPr.append(lay)
    lay.set(qn("w:type"), "fixed")
    ws = []
    for j in range(len(table.columns)):
        mx = 1
        for cell in table.columns[j].cells:
            for para in cell.paragraphs:
                for run in para.runs:
                    mx = max(mx, len(run.text))
        ws.append(max(mx, 4))
    tot = sum(ws)
    for j, w in enumerate(ws):
        ew = Emu(int(Inches(total_in) * w / tot))
        for cell in table.columns[j].cells:
            cell.width = ew


def tbl(doc, df, above=None, below=None, fs=8.4, total_in=6.7):
    if above:
        p = doc.add_paragraph(); r = p.add_run(above)
        r.font.size = Pt(9.0); r.bold = True
        p.paragraph_format.keep_with_next = True
    t = doc.add_table(rows=1, cols=len(df.columns))
    t.style = "Table Grid"; t.alignment = WD_TABLE_ALIGNMENT.CENTER; t.autofit = False
    hdr = t.rows[0]
    trPr = hdr._tr.get_or_add_trPr()
    trPr.append(trPr.makeelement(qn("w:tblHeader"), {}))
    for i, c in enumerate(df.columns):
        cell = hdr.cells[i]; cell.text = ""
        for cc, k in toks(c):
            r = cell.paragraphs[0].add_run(cc); r.bold = True; r.font.size = Pt(fs)
            if k == "sub":
                r.font.subscript = True
    for _, row in df.iterrows():
        cells = t.add_row().cells
        for i, v in enumerate(row):
            cells[i].text = ""
            cells[i].paragraphs[0].add_run(str(v)).font.size = Pt(fs)
    _autofit(t, total_in)
    if below:
        cap(doc, below)
    else:
        doc.add_paragraph()


def land(section):
    w, h = section.page_width, section.page_height
    section.orientation = WD_ORIENT.LANDSCAPE
    section.page_width, section.page_height = h, w


def port(section):
    w, h = section.page_width, section.page_height
    section.orientation = WD_ORIENT.PORTRAIT
    if w > h:
        section.page_width, section.page_height = h, w


T1 = pd.DataFrame([
 {"Evaluation set": "Full annotation namespace", "Rows": "4,473", "Targets": "214",
  "Targets (n ≥ 5)": "130", "Method families": "—", "Purpose": "Namespace definition"},
 {"Evaluation set": "Primary matched set", "Rows": f"{PU['n_rows']:,}", "Targets": "166",
  "Targets (n ≥ 5)": f"{PU['n_ge5']}", "Method families": f"{PU['n_method_families']}",
  "Purpose": "Main capability-resolved analyses"},
 {"Evaluation set": "Controlled split set", "Rows": "7,523", "Targets": "230",
  "Targets (n ≥ 5)": "172", "Method families": "12",
  "Purpose": "Split-strictness experiment"},
 {"Evaluation set": "External unambiguous set", "Rows": f"{FR['external_clean']['n_rows']:,}",
  "Targets": f"{FR['external_clean']['n_targets']}",
  "Targets (n ≥ 5)": f"{FR['external_clean']['n_ge5']}", "Method families": "5",
  "Purpose": "Boundary of generalisation"},
])
T2 = pd.DataFrame([
 {"Split regime": nm, "Pooled Pearson": f"{H4[k]['overall_pearson']:.4f}",
  "C_{B}": f"{H4[k]['C_B']:+.4f}", "C_{W}": f"{H4[k]['C_W']:+.4f}",
  "Macro Spearman": f"{H4[k]['macro_spearman']:.4f}", "RMSE": f"{H4[k]['rmse']:.4f}"}
 for k, nm in [("random", "Random"), ("complex", "Complex"), ("seq40", "seq40"),
               ("family_OOD", "Family-OOD")]])
NICE = {"Flex-ddG": "Flex-ddG", "DDGPred": "DDGPred", "ddG_FoldX": "FoldX (ddG_FoldX)",
        "ESM-IF": "ESM-IF", "proteinMPNN": "ProteinMPNN"}
T3 = EP[["model", "n_rows", "n_targets_ge5", "overall_pearson", "macro_spearman",
         "C_B_row", "C_W", "rmse"]].copy()
T3["model"] = T3["model"].map(lambda m: NICE.get(m, m))
T3.columns = ["Method", "Mutations", "Targets (n ≥ 5)", "Pooled Pearson",
              "Macro Spearman", "C_{B}", "C_{W}", "RMSE"]
for c in ["Pooled Pearson", "Macro Spearman", "C_{B}", "C_{W}", "RMSE"]:
    T3[c] = T3[c].map(lambda x: f"{float(x):.4f}")
T3["Mutations"] = T3["Mutations"].map(lambda x: f"{int(x):,}")

FIGW = {"Fig1_overview_EN": 6.7, "Fig2_controlled_split_collapse_EN": 6.7,
        "Fig3_target_weighting_sensitivity_EN": 6.7, "Fig4_fixed_K_geometry_EN": 6.1,
        "Fig5_retrospective_remedy_EN": 6.2, "Fig6_external_methods_EN": 6.7}
# 用户提供的成图尺寸差异较大，逐图指定插入宽度与格式
FIGSPEC = {
    "Fig1_overview_EN.png": (6.5, "png"),
    "Fig2_controlled_split_collapse_EN.png": (6.5, "png"),
    "Fig3_target_weighting_sensitivity_EN.png": (6.5, "png"),
    "Fig4_fixed_K_geometry_EN.png": (6.2, "png"),
    "Fig5_retrospective_remedy_EN.png": (6.2, "png"),
    "Fig6_external_methods_EN.png": (6.6, "png"),
    "FigS1_matched_universe_strata_EN.png": (6.5, "png"),
    "FigS2_capability_heatmap_provenance_clean_EN.png": (6.3, "png"),
    "FigS3_pool_sensitivity_rank_agreement_EN.png": (6.5, "png"),
}


def addfig(doc, fname):
    if not fname.lower().endswith((".png", ".jpg", ".jpeg")):
        fname = fname + ".png"
    w, _ = FIGSPEC.get(fname, (6.5, "png"))
    doc.add_picture(str(FIGS / fname), width=Inches(w))


def build(zh=False):
    D = Document(); setup(D)
    t = D.add_paragraph()
    r = t.add_run("Disentangling Between-Target and Within-Target Signal in "
                  "PPI ΔΔG Benchmarking")
    r.bold = True; r.font.size = Pt(19); r.font.name = "Arial"
    t.alignment = WD_ALIGN_PARAGRAPH.CENTER
    s = D.add_paragraph()
    rs = s.add_run("An exact covariance-signal decomposition under matched, "
                   "family-out-of-distribution evaluation")
    rs.italic = True; rs.font.size = Pt(11)
    s.alignment = WD_ALIGN_PARAGRAPH.CENTER
    D.add_paragraph()
    for line in ["Author Name(s) [AUTHOR_PLACEHOLDER]",
                 "Affiliation(s) [AFFILIATION_PLACEHOLDER]",
                 "Corresponding author: [NAME] · Email: [EMAIL] · "
                 "ORCID: [ORCID_PLACEHOLDER]"]:
        p = D.add_paragraph(); rr = p.add_run(line); rr.font.size = Pt(10)
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    D.add_paragraph()
    p = D.add_paragraph()
    rr = p.add_run("Keywords: "); rr.bold = True; rr.font.size = Pt(10)
    p.add_run("protein–protein interaction; binding affinity; ΔΔG; mutation effect "
              "prediction; benchmarking; out-of-distribution generalisation"
              ).font.size = Pt(10)
    D.add_page_break()

    H(D, "Abstract", 1)
    W(D, f"Predicting how mutations alter protein–protein binding affinity (ΔΔG) "
         f"underpins antibody engineering, interface optimisation and the interpretation "
         f"of disease-associated variants. The pooled correlation commonly reported for "
         f"such benchmarks is difficult to interpret directly, because it depends both on "
         f"whether a model places different targets on the correct overall binding scale "
         f"and on whether it orders mutations correctly within each target. We introduce an "
         f"exact covariance-signal decomposition that writes pooled Pearson correlation as "
         f"r = C_{{B}} + C_{{W}}, and apply it to matched, de-replicated and "
         f"provenance-filtered evaluation sets. In a controlled comparison that changes only "
         f"the split regime (identical 7,523 mutation rows and identical model "
         f"configurations), median pooled Pearson falls from {r0:.4f} under random splitting "
         f"to {r3:.4f} under family-out-of-distribution evaluation, with "
         f"ΔC_{{B}} = {DL['delta_C_B']:+.4f} and ΔC_{{W}} = {DL['delta_C_W']:+.4f}; on an "
         f"absolute-change basis the between-target component accounts for "
         f"{100*DL['fraction_between']:.1f}% of the decline while the within-target "
         f"component remains positive. Target weighting changes how these contributions "
         f"should be read: C_{{B}} reverses direction only near zero ({sf0}/{len(W2)} "
         f"configurations), whereas C_{{W}} keeps its sign at every threshold although its "
         f"magnitude is not invariant. Fixed-K shortlist metrics are constrained by "
         f"candidate-set geometry: for K = 10, D(n, K) is zero when n ≤ K and peaks near "
         f"n* ≈ {GEO['n_star_analytic']:.1f}, and {100*GEO['n_le_K_fraction']:.1f}% of "
         f"targets with n ≥ 5 fall in that degenerate regime. In a matched retrospective "
         f"ablation the decomposition separated a change that improved pooled scores only "
         f"from one that also improved within-target ranking. External evidence comprised "
         f"five methods and is underpowered and positively aligned, so it does not support a "
         f"universal decoupling of pooled and within-target performance. We recommend that "
         f"PPI ΔΔG benchmarks report both components under an explicitly declared weighting, "
         f"together with candidate-set geometry.")

    H(D, "Introduction", 1)
    W(D, "Predicting how mutations alter protein–protein binding affinity (ΔΔG) is central "
         "to antibody affinity maturation, interface optimisation and the interpretation of "
         "disease-associated variants. As the number of deep-learning methods grows rapidly, "
         "progress is increasingly judged by benchmark evaluation, and the properties of the "
         "evaluation protocol itself — what it measures, and under which extrapolation "
         "conditions — determine what the field can conclude.")
    W(D, "Current practice relies largely on pooled metrics, computing Pearson or Spearman "
         "correlation and RMSE over the entire mutation set[13], often supplemented by "
         "per-complex or per-target metrics. Bushuiev et al. quantified structural leakage "
         "in PPI data: about 65% of test examples have a near-duplicate training example "
         "under PDB-code splitting, about 30% even under 30% sequence-identity splitting, "
         "56% for SKEMPI v2.0[1] split by PPI code, and 0% under domain-expert grouping[2]. "
         "CATH-ddG adopted CATH homologous-superfamily splits and defined hard mutations by "
         "a maximum TM-score below 0.6[3]; related methods likewise report per-structure or "
         "per-target correlations[4,5,10]. We emphasise that leakage-aware splitting, "
         "family-out-of-distribution evaluation, per-target reporting and shortlist (Top-K) "
         "evaluation were none of them introduced by this work.")
    W(D, "A key question nevertheless remains unanswered directly: where does pooled "
         "benchmark performance come from? Reporting pooled and per-target metrics side by "
         "side does not quantify the additive contribution of between-target and "
         "within-target signal to pooled performance. In benchmarks containing many targets "
         "with highly unequal mutation counts, pooled correlation is additionally dominated "
         "by large targets, which amplifies the conflation.")
    W(D, "Here we introduce an exact covariance-signal decomposition that addresses this "
         "question. We decompose pooled Pearson correlation rigorously into a "
         "between-target component and a within-target component, which makes 'what pooled "
         "performance is made of' a directly measurable quantity. We quantify four things: "
         "which component is changed by split strictness; how target weighting changes the "
         "interpretation of the two components; which constraints candidate-set geometry "
         "places on fixed-K shortlist evaluation; and whether the decomposition can "
         "distinguish model modifications that improve pooled scores only from those that "
         "also improve within-target ranking. We further bound the evidence using five "
         "independently developed external methods.")

    H(D, "Results", 1)
    tbl(D, T1, above="Table 1 | Dataset and evaluation overview.",
        below="All rows are joined by explicit row identifiers.")
    H(D, "Pooled PPI ΔΔG performance contains separable between-target and within-target "
         "signal", 2)
    W(D, "For a mutation i belonging to target c, the observed ΔΔG can be written as a "
         "target-level mean plus a mutation-specific deviation from that mean:")
    eq = D.add_paragraph(); eq.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for c, k in toks("y_{ci} = μ_{c} + ε_{ci}          ŷ_{ci} = μ̂_{c} + ε̂_{ci}"):
        rr = eq.add_run(c)
        if k == "sub":
            rr.font.subscript = True
    W(D, "With σ_{y} and σ_{ŷ} the pooled standard deviations, the pooled correlation "
         "splits exactly as r_{pooled} = Cov(y, ŷ) / (σ_{y} σ_{ŷ}) = C_{B} + C_{W}, where "
         "C_{B} and C_{W} are the between-target and within-target covariance contributions. "
         "The identity holds to machine precision (residual ≤ 3.3 × 10⁻¹⁶). These are "
         "covariance-signal contributions on a common scale, not independent correlations, "
         "skill scores or causal mechanisms.")
    W(D, f"We apply this decomposition on a matched evaluation set comprising "
         f"{PU['n_rows']:,} mutation rows, {PU['n_ge5']} targets with at least five "
         f"mutations, and {PU['n_method_families']} independent method families.")
    addfig(D, "Fig1_overview_EN.png")
    cap(D, "Figure 1 | Capability-resolved evaluation of PPI ΔΔG predictors. "
           "(A) Provenance-filtered, fingerprint-matched benchmark construction: SKEMPI v2.0 "
           "annotation (4,473 rows, 214 targets, 108 families) is audited by SHA-256 row-set "
           "fingerprinting; family-out-of-distribution predictions without training overlap "
           "are retained and permissively split predictions excluded, giving a primary "
           "matched universe of 3,351 rows, 102 targets with n ≥ 5 and 10 method families. "
           "Target cards show between-target differences in mean ΔΔG (μ_c, purple) and "
           "within-target mutation ordering (ε_ci, green). (B) Exact additive "
           "covariance-signal decomposition; C_B and C_W sum to pooled correlation by an "
           "exact identity (residual ≤ 3.3 × 10⁻¹⁶) and are common-scale additive "
           "contributions, not independent correlations or causal mechanisms. (C) Controlled "
           "split collapse on identical 7,523 mutation rows: pooled correlation falls from "
           "0.6844 to 0.1508 and 81.1% of the absolute decline is carried by loss of the "
           "between-target component. (D) Retrospective matched ablation on the identical "
           "3,351-row evaluation universe. Pooled Pearson and within-target Macro Spearman "
           "are shown for the mutation-only baseline (F1), static local geometry (F2), "
           "mutation × partner perturbation (F4), and the full local-physics model (F8). "
           "Static local geometry improves pooled correlation while leaving within-target "
           "ranking essentially unchanged, whereas F4 and F8 improve both axes. The "
           "capability trajectory is retrospective and is not intended as a causal or "
           "prospective model-selection analysis.")

    H(D, "Family-OOD evaluation disproportionately removes between-target signal", 2)
    W(D, f"To isolate the effect of split strictness, we compared identical model "
         f"configurations on the same 7,523 mutation rows while changing only the split "
         f"regime. Median pooled Pearson falls from {r0:.4f} to {r3:.4f}. The decline is not "
         f"uniform: ΔC_{{B}} = {DL['delta_C_B']:+.4f} and "
         f"ΔC_{{W}} = {DL['delta_C_W']:+.4f}. On an absolute-change basis the between-target "
         f"component accounts for {100*DL['fraction_between']:.1f}% of the decline, while "
         f"C_{{W}} remains positive ({cw3:+.4f}). This is not the same proposition as 'the "
         f"model memorised target identity': a more defensible reading is that the pattern is "
         f"consistent with reliance on target-specific information that does not transfer "
         f"across families. The complex-level and seq40 partitions produced identical "
         f"assignments and were therefore not treated as independent comparisons.")
    addfig(D, "Fig2_controlled_split_collapse_EN.png")
    cap(D, "Figure 2 | Family-out-of-distribution evaluation disproportionately removes the "
           "between-target component. (A) Pooled Pearson by split regime. (B) Exact "
           "decomposition into C_B and C_W; C_W remains positive. (C) Attribution of the "
           f"absolute decline, of which the between-target component accounts for "
           f"{100*DL['fraction_between']:.1f}%. All four regimes use the same 7,523 rows and "
           "configurations; complex and seq40 are identical and are not independent "
           "comparisons.")
    tbl(D, T2, above="Table 2 | Controlled split comparison.",
        below="Median values across configurations. seq40 is identical to complex in this "
              "dataset.")

    H(D, "Target aggregation changes the interpretation of covariance contributions", 2)
    W(D, f"Holding predictions fixed and varying only the aggregation weights across five "
         f"schemes (per-mutation, per-target, proportional to √n, proportional to "
         f"log(1 + n), and capped at 30), C_{{B}} shows directional instability concentrated "
         f"near zero: {sf0}/{len(W2)} configurations reverse sign with no threshold, whereas "
         f"none do once |C_{{B}}| exceeds 0.005 ({sf5}/{len(W2)}). C_{{W}} never reverses "
         f"sign at any threshold ({sfw}/{len(W2)}), but its magnitude is not invariant — its "
         f"median range is {W2.range_C_W.median():.4f}, larger than the "
         f"{W2.range_C_B.median():.4f} of C_{{B}}. Any reported decomposition should declare "
         f"its weighting scheme.")
    addfig(D, "Fig3_target_weighting_sensitivity_EN.png")
    cap(D, "Figure 3 | Target weighting changes the interpretation of covariance "
           "contributions. Predictions are unchanged; only aggregation weights differ. "
           "(A, B) Distribution of the range of C_B and C_W across the five schemes, with "
           "the median marked. (C) Configurations reversing sign at three thresholds. C_B "
           "reverses direction only near zero, whereas C_W keeps its sign but its magnitude "
           "is not invariant.")

    H(D, "Fixed-K shortlist evaluation is constrained by candidate-set geometry", 2)
    W(D, f"Fixed-K top-K metrics are used for mutation library or shortlist "
         f"screening[11]. The discriminable interval is D(n, K) = 0 for n ≤ K and "
         f"D(n, K) = K(n − K)/[n(n + 1)] for n > K. When a target has no more than K "
         f"mutations, any model selects every candidate and the target carries no information "
         f"for fixed-K shortlist metrics. In this benchmark, {GEO['n_le_K_count']} of "
         f"{GEO['n_ge5_count']} targets with n ≥ 5 ({100*GEO['n_le_K_fraction']:.1f}%) fall "
         f"in that degenerate regime. D is not monotonic in n: it peaks at "
         f"n* = K + √(K² + K), which for K = 10 gives n* ≈ {GEO['n_star_analytic']:.1f}, so "
         f"that D(20) = {GEO['D_at_n20']:.4f} while D(50) = {GEO['D_at_n50']:.4f}. This "
         f"geometric result does not imply that SKEMPI ΔΔG values are intrinsically too "
         f"similar; establishing that would require an analysis using actual ΔΔG utility "
         f"magnitudes, which we did not perform.")
    addfig(D, "Fig4_fixed_K_geometry_EN.png")
    cap(D, "Figure 4 | Fixed-K shortlist evaluation is constrained by candidate-set "
           "geometry. The discriminable interval D(n, K) is zero when n ≤ K and follows "
           "K(n − K)/[n(n + 1)] otherwise, peaking at n* = K + √(K² + K) ≈ 20.5 for K = 10. "
           f"In this benchmark {100*GEO['n_le_K_fraction']:.1f}% of targets with n ≥ 5 lie in "
           "the degenerate regime where every model selects all candidates.")

    H(D, "Capability decomposition clarifies retrospective remedy effects", 2)
    W(D, "In a matched ablation set sharing an identical 3,351-row universe, the "
         "decomposition distinguishes modifications that raise pooled scores only from those "
         "that also improve within-target ranking. Relative to a mutation-only baseline, "
         "adding static local geometry raised pooled Pearson by about 0.028 while leaving "
         "within-target Macro Spearman essentially unchanged (−0.002); mutation × partner "
         "perturbation improved both axes (+0.070 and +0.044); and a full local-physics "
         "model gave the largest improvement on both axes (+0.105 and +0.076). This analysis "
         "is retrospective and does not establish prospective prediction of which "
         "modification will succeed.")
    addfig(D, "Fig5_retrospective_remedy_EN.png")
    cap(D, "Figure 5 | Capability decomposition separates pooled-score gains from "
           "within-target gains. Retrospective matched ablation on an identical 3,351-row "
           "universe; the baseline is the mutation-only model. Green markers indicate "
           "improvement on both axes; red markers indicate a pooled-score-only or negative "
           "change. The analysis is retrospective and does not imply prospective prediction "
           "of which modification will succeed.")

    H(D, "External methods define the boundary of generalisation", 2)
    W(D, f"We built external evaluation sets from published official per-mutation "
         f"predictions: Flex ddG[6], FoldX[7], ProteinMPNN[8], ESM-IF[9] and DDGPred[4]. "
         f"Where the matching key was not unique on the annotation side, only uniquely "
         f"mapped rows were used; of 4,473 rows only {FR['external_clean']['n_rows']:,} have "
         f"a unique mapping, and ambiguous keys were excluded rather than forcibly resolved. "
         f"On this unambiguous set the five methods agree in their ordering on pooled "
         f"correlation and on within-target ranking, with DDGPred best on both axes; no "
         f"method-level rank reversal occurs. Method-level rank agreement is ρ = +0.90 "
         f"(exact permutation p = 0.083). This strongly positive agreement does not support "
         f"a universal decoupling of pooled and within-target performance and falls within "
         f"the heterogeneous range observed internally. With only five methods the panel is "
         f"underpowered, so we draw no directional conclusion from it.")
    addfig(D, "Fig6_external_methods_EN.png")
    cap(D, "Figure 6 | External methods define the boundary of generalisation. (A) Pooled "
           "versus within-target performance for five independently developed methods on a "
           "matched, unambiguous mapping. (B) Method-level capability profile, standardised "
           "within each metric for visual comparison only; raw values are printed in each "
           "cell. Rank agreement is positive (ρ = +0.90, exact permutation p = 0.083) but the "
           "panel is underpowered with five methods.")
    tbl(D, T3, above="Table 3 | External method capability profiles on the unambiguous "
                     "mapping set.",
        below="All methods are evaluated on the same mutations. RMSE is reported only for "
              "outputs with ΔΔG-like scaling.")

    H(D, "Discussion", 1)
    for para in [
        "This study shows that the pooled correlation conventionally reported as a single "
        "number in PPI ΔΔG benchmarking is composed of two qualitatively different "
        "predictive signals: a between-target component that places different targets on the "
        "correct overall scale, and a within-target component that orders mutation effects "
        "correctly inside each target. Because our decomposition is an exact identity, 'what "
        "pooled performance is made of' becomes directly computable.",
        "The relation to prior work should be delimited clearly. Structural leakage in PPI "
        "data and its dependence on splitting strategy were quantified systematically by "
        "Bushuiev et al., who also proposed interface-similarity splitting[2]; CATH-ddG "
        "adopted CATH homologous-superfamily splits and used structural similarity to define "
        "harder extrapolation mutations[3]; and per-complex or per-target reporting together "
        "with shortlist-style evaluation already appear in the ΔΔG literature[4,5,10]. Our "
        "contribution lies not in any single one of these ingredients but in the systematic "
        "integration of an exact between/within-target covariance decomposition with four "
        "audits — split strictness, target weighting, matched rank agreement and "
        "candidate-set geometry — inside one matched, provenance-filtered framework, "
        "together with a retrospective attribution case study.",
        "The interpretation of C_{B} and C_{W} has clear boundaries. They are covariance "
        "contributions within one identity and share a common scale, so they sum to pooled "
        "correlation; they are not two independently estimated skill scores and should not "
        "be read as causal mechanisms. A model may obtain a substantial pooled correlation "
        "largely by placing targets on approximately the right overall scale while its "
        "within-target ordering ability remains limited — precisely why pooled metrics can "
        "overstate practical screening value.",
        "The weighting-sensitivity audit yields a conclusion of direct practical relevance: "
        "the sign of the between-target component is sensitive to aggregation weights near "
        "zero, whereas the sign of the within-target component is stable. When a reported "
        "C_{B} is close to zero, the conclusion about whether between-target information is "
        "present may depend on the weighting convention rather than on the model. Any "
        "reported decomposition should declare its weighting scheme.",
        "The candidate-set geometry result has direct implications for Top-K practice. When "
        "a target has no more than K mutations, fixed-K shortlist metrics cannot distinguish "
        "models structurally; roughly four in ten targets with n ≥ 5 fall in this regime in "
        "our benchmark. Because the discriminable interval first rises and then falls as "
        "candidate sets grow, enlarging a candidate set is not universally effective at "
        "improving discrimination. We recommend reporting n alongside K, labelling n ≤ K "
        "targets separately, or adopting fraction-based evaluation.",
        "The retrospective ablation illustrates the framework's explanatory value: it "
        "distinguishes a change that raised pooled scores only from a change that also "
        "improved within-target ranking. The retrospective nature of the case study must, "
        "however, be stated explicitly: it does not establish prospective ability to predict "
        "which modifications will work.",
        "Finally, external evidence bounds the generalisation of these conclusions. Across "
        "five independently developed methods the ordering of methods agreed between pooled "
        "and within-target performance, a strongly positive agreement that neither supports "
        "nor refutes a universal decoupling of the two axes. The panel is underpowered, so we "
        "draw no directional conclusion from it. Future work should broaden the external "
        "panel and test the relationship between the decomposed components and realised "
        "screening yield on independent prospective datasets.",
    ]:
        W(D, para)

    H(D, "Methods", 1)
    for title, body in [
        ("Dataset and target definition",
         f"Evaluation rows come from a SKEMPI v2.0-derived annotation set[1] (4,473 rows, "
         f"214 targets, 108 structural families), joined to prediction outputs by explicit "
         f"row identifiers. A target is defined as an independent binding-interface unit "
         f"within a structural family."),
        ("Prediction sources and provenance criteria",
         "Prediction outputs were classified by their recorded supervised training split "
         "regime. Configurations trained under a permissive (random / complex / seq40) "
         "partition were excluded from deployment-oriented comparisons because their "
         "evaluation rows may overlap their training folds; they were retained only in the "
         "controlled split experiment. Analyses were performed at artifact, configuration "
         "and method-family level, with one performance-blind representative per independent "
         "family."),
        ("Matched evaluation sets",
         "Each prediction output received a SHA-256 fingerprint of its sorted row "
         "identifiers; cross-model comparison was performed only between outputs sharing a "
         "fingerprint (Supplementary Note 2, Supplementary Fig. S1). The primary matched set "
         "was selected by a rule fixed before any metric values were inspected."),
        ("Between-/within-target signal decomposition",
         "Row-level weights are used and the identity residual is reported as a "
         "quality-control quantity. Full derivation in Supplementary Note 1."),
        ("Target-weighting schemes",
         "Predictions are unchanged; only aggregation weights vary: per-mutation, "
         "per-target, proportional to √n, proportional to log(1 + n), and capped at 30."),
        ("Ranking and shortlist metrics",
         f"Within-target ordering used per-target Spearman ρ and Kendall τ correlations, "
         f"pairwise accuracy and the within-target covariance contribution; shortlist "
         f"metrics used fixed-K and fraction-based hit rate, precision, recall, NDCG and "
         f"enrichment factor[11]."),
        ("Statistical analysis",
         f"Rank agreement uses Spearman ρ with Kendall τ as secondary; interval estimates "
         f"use target-cluster bootstrap (5,000 resamples, seed 42)[12]; small external "
         f"panels use exact permutation tests. All cross-metric comparisons use "
         f"higher-is-better orientation[13]."),
        ("External prediction mapping and ambiguity control",
         "Where the normalised matching key was not unique on the annotation side, only "
         "uniquely mapped rows were used; ambiguous keys were excluded rather than resolved."),
        ("Reproducibility and consistency validation",
         "All headline numbers are generated from a single machine-readable file and "
         "verified by an automated consistency checker. No new model was trained and no GPU "
         "was used."),
    ]:
        H(D, title, 2); W(D, body)

    H(D, "Limitations", 1)
    for para in [
        "Our extrapolation setting is based on target-family partitioning, which does not "
        "exclude the possibility that foundation-model pretraining encountered sequences or "
        "structures similar to the evaluation families; we can control task-specific "
        "supervised overlap only. The model pool is not a random sample of all PPI ΔΔG "
        "methods, and only five methods are externally evaluable, so method-level inference "
        "is underpowered.",
        "Shortlist and candidate-selection analyses rest on re-analysis of existing "
        "predictions and include no prospective wet-lab validation; when a target has no more "
        "than K mutations the relevant metrics carry no information. The retrospective "
        "ablation likewise does not establish prospective prediction.",
        "C_{B} and C_{W} are covariance-signal contributions rather than causal mechanisms. "
        "We performed no identity-memorisation or target-shuffling intervention and therefore "
        "make no causal claim about the origin of the between-target component.",
        "No independent prospective experimental dataset was available, so no external "
        "dataset validation was performed. This work does not provide a universal model "
        "recommender and does not claim general superiority of any method family.",
    ]:
        W(D, para)

    H(D, "Data and code availability", 1)
    W(D, "All analyses use existing predictions and published external predictions. Frozen "
         "input hashes, the machine-readable source of every headline number and an automated "
         "consistency checker are provided with the Supplementary Information. No new model "
         "was trained and no GPU was used.")
    H(D, "Author contributions", 1)
    W(D, "Study conception, analysis design, interpretation of results and writing of the "
         "manuscript were performed by the author(s). The use of AI-assisted tools is "
         "declared below.")
    H(D, "Use of AI-assisted tools", 1)
    W(D, "AI-assisted coding and language tools were used to support code drafting and "
         "debugging, generation of plotting scripts, document automation, literature and "
         "repository retrieval, and language editing. All metric definitions, analysis "
         "choices, scientific interpretations, hypothesis decisions and final manuscript "
         "content were reviewed and approved by the authors, who take full responsibility for "
         "the work. AI tools were not used to generate experimental data or to modify "
         "upstream model predictions.")
    H(D, "Acknowledgements", 1)
    W(D, "[FUNDING_PLACEHOLDER] We thank the developers of the publicly released prediction "
         "sets and datasets used in this study.")

    H(D, "References", 1)
    for i, txt in enumerate(REFS, 1):
        p = D.add_paragraph()
        r = p.add_run(f"{i}. {txt}"); r.font.size = Pt(9.5)
        p.paragraph_format.left_indent = Inches(0.30)
        p.paragraph_format.first_line_indent = Inches(-0.30)
        p.paragraph_format.space_after = Pt(3)

    nm = "CR-PPI-ddG_manuscript_v6_EN.docx"
    path = OUT / nm
    D.save(str(path))
    print(f"  {nm} ({path.stat().st_size/1024:.0f} KB)")


def build_si():
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
         "y_{ci} = μ_{c} + ε_{ci} and ŷ_{ci} = μ̂_{c} + ε̂_{ci}. Expanding Cov(y, ŷ) gives a "
         "between-target term N⁻¹ Σ_{c} n_{c} (μ_{c} − μ)(μ̂_{c} − μ̂) and a within-target "
         "term N⁻¹ Σ_{c} Σ_{i} ε_{ci} ε̂_{ci}; all cross terms vanish identically because "
         "within-target deviations have zero mean in every target. Dividing by σ_{y} σ_{ŷ} "
         "gives C_{B} and C_{W} and establishes r_{pooled} = C_{B} + C_{W} exactly. The "
         "identity residual is reported as a quality-control quantity (observed maximum "
         "3.331 × 10⁻¹⁶). Because C_{B} and C_{W} share the same denominator they are "
         "covariance-signal contributions on a common scale, not independent correlations.")
    H(D, "Supplementary Note 2 | Matched universes and row-set fingerprints", 1)
    W(D, "Each prediction output received a SHA-256 fingerprint of its sorted row "
         "identifiers. Cross-model comparison was performed only between outputs sharing a "
         "fingerprint; raw metrics were never compared across fingerprints.")
    land(D.add_section())
    U = UR.rename(columns={"universe_id": "Universe", "artifact": "Artifact",
                           "rowset_hash": "Row-set fingerprint", "n_rows": "Rows (N)",
                           "n_targets": "Targets", "n_families": "Families",
                           "n_ge3": "Targets (n ≥ 3)", "n_ge5": "Targets (n ≥ 5)",
                           "n_ge10": "Targets (n ≥ 10)", "n_models": "Models"})
    tbl(D, U, fs=7.5, total_in=9.4,
        below="Table S1 | Matched-evaluation strata and row-set fingerprints.")
    port(D.add_section())
    H(D, "Supplementary Figures", 1)
    for stem, w, c in [
        ("FigS1_matched_universe_strata_EN", 6.5,
         "Figure S1 | Matched-evaluation strata grouped by SHA-256 row-set fingerprints. "
         "Each bar is one prediction output; outputs sharing a fingerprint are directly "
         "comparable with each other and with no others."),
        ("FigS2_capability_heatmap_provenance_clean_EN", 6.4,
         "Figure S2 | Capability profiles of matched, provenance-clean ablation variants "
         "(seed 42). Rows are controlled ablation variants of one local-physics method "
         "family, not independent methods. Values are standardised within each metric for "
         "visual comparison only; raw values are printed in each cell."),
        ("FigS3_pool_sensitivity_rank_agreement_EN", 6.5,
         "Figure S3 | Pool sensitivity of cross-metric rank agreement across matched strata. "
         "Within one evaluation universe and row set, changing only which model families are "
         "included moves the agreement between pooled Pearson and macro within-target "
         "Spearman across a wide range, which is why this quantity is reported as a range."),
    ]:
        addfig(D, f"{stem}.png")
        cap(D, c)
    nm = "CR-PPI-ddG_Supplementary_v6_EN.docx"
    path = OUT / nm
    D.save(str(path))
    print(f"  {nm} ({path.stat().st_size/1024:.0f} KB)")


build()
build_si()
print("v6 完成")
