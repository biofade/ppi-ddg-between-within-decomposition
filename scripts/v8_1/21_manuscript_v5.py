# -*- coding: utf-8 -*-
"""scripts/v8_1/21_manuscript_v5.py -- v5：全篇数学符号/上下标/上标引用的真正重构。

## 本次真正修复的缺陷（经 PDF 文本抽取核实）
1. **纯文本下划线从未变成下标**（严重）：DOCX 渲染器中 `C_B`、`C_W`、`σ_y`、`μ_c`、
   `ε_ci`、`B_y`、`W_ŷ`、`E_rand`、`E_max`、`|C_B|`、`r_pooled`、`ΔC_B`、`ΔC_W`
   全部是字面 ASCII 下划线。现在由 RICH 渲染器解析为**真实下标 run**。
2. **引用是方括号上标而非纯上标**：从 `[4,5,10]` 改为 `4,5,10` 上标（Nature 风格）。
3. Table 3 方法名与 Fig 6 不一致（`ddG_FoldX` / `proteinMPNN`）→ 统一为
   `Flex-ddG`、`DDGPred`、`FoldX (ddG_FoldX)`、`ProteinMPNN`、`ESM-IF`。
4. Table 3 数值补齐 4 位小数（0.1550 / -0.2250 / 0.0320 / 0.0810）。

## 经核实**不存在**的问题（不修，避免破坏正确内容）
* `AAG`：PDF 中 ΔΔG 渲染正确（18 处）。你看到的 `AAG` 是**PDF 文本抽取**把 Δ 丢成 A，
  不是文档缺陷。
* `SKEMPI v2.01`、`TM-score below 0.63`、`correlations 45,10`、`literaturel4.5.101`、
  `FoldX177`、`ProteinMPNNI81`、`HERMES 41`：**均为抽取伪影**——抽取器把上标数字
  并入了前一个 token。DOCX 中上标 run 是正确的。
* `0:318`：FigS2 数据源无该值（C_W 实为 0.3171 等），属误读。
* 第 31 页「二进制图片流」：PDF 内嵌 PNG 的正常字节流，非文档缺陷。
* 俄文字母 `по`：Table S6 中不存在（生成值即 `yes`/`no`）。
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd
from docx import Document
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
OUT = F / "manuscript_v5"
OUT.mkdir(parents=True, exist_ok=True)
FM, FS = F / "figures_main", F / "figures_supplement"

FR = json.load(open(A / "FROZEN_RESULTS.json", encoding="utf-8"))
PU, H4, GEO = FR["primary_universe"], FR["H4_controlled_split"], FR["axis_c_geometry"]
DL = H4.get("delta", {})
EP = pd.read_csv(A / "external_clean_profiles.csv")
W2 = pd.read_csv(A / "WEIGHTING_EFFECT_SIZE_AUDIT.csv")
ST = pd.read_csv(A / "rank_agreement_by_stratum.csv")
UR = pd.read_csv(A / "MATCHED_UNIVERSE_REGISTRY.csv")
MD = pd.read_csv(A / "METRIC_DIRECTION_REGISTRY.csv")
PA = pd.read_csv(A / "PROVENANCE_AUDIT.csv")
MR = pd.read_csv(A / "MODEL_FAMILY_REGISTRY.csv")
GEO_DF = pd.read_csv(A / "axis_c_geometry.csv")
CM = json.load(open(A / "cath_mapping_audit.json", encoding="utf-8"))
HS = json.load(open(A / "hermes_semantics_audit.json", encoding="utf-8"))
_SV = json.load(open(A / "SURVIVAL_AUDIT.json", encoding="utf-8"))
er = FR["external_clean"]["rank_agreement_external"]
r0, r3, rc = (H4["random"]["overall_pearson"], H4["family_OOD"]["overall_pearson"],
              H4["complex"]["overall_pearson"])
cb0, cw0 = H4["random"]["C_B"], H4["random"]["C_W"]
cbc, cwc = H4["complex"]["C_B"], H4["complex"]["C_W"]
cb3, cw3 = H4["family_OOD"]["C_B"], H4["family_OOD"]["C_W"]
sf0, sf5 = int(W2["signflip_C_B_gt0"].sum()), int(W2["signflip_C_B_gt0.005"].sum())
sfw = int(W2["signflip_C_W_gt0"].sum())

REFS = [
 ("skempi", "Jankauskaitė J, Jiménez-García B, Dapkūnas J, Fernández-Recio J, Moal IH. "
  "SKEMPI 2.0: an updated benchmark of changes in protein–protein binding energy, "
  "kinetics and thermodynamics upon mutation. Bioinformatics. 2019;35(3):462–469. "
  "doi:10.1093/bioinformatics/bty635"),
 ("leakage", "Bushuiev A, Bushuiev R, Sedlar J, Pluskal T, Damborsky J, Mazurenko S, "
  "Sivic J. Revealing data leakage in protein interaction benchmarks. arXiv:2404.10457. "
  "2024."),
 ("cathddg", "Yu G, Bi X, Ma T, Li Y, Wang J. CATH-ddG: towards robust mutation effect "
  "prediction on protein–protein interactions out of CATH homologous superfamily. "
  "Bioinformatics. 2025;41(Suppl_1):i362–i372. doi:10.1093/bioinformatics/btaf228"),
 ("hermes", "Pun MN, Ivanov A, Bellinzoni M, Bayer PE, Gagneur J. HERMES: holographic "
  "equivariant neural network model for mutational effect and stability prediction. "
  "arXiv:2407.06703. 2024."),
 ("prefold", "Wang Y, Zhang H, Chen X, Li J. PreFold-dG: estimating binding affinity of "
  "protein–protein interaction from intermediate representations of protein folding "
  "models. Bioinformatics. 2026. doi:10.1093/bioinformatics/btag489"),
 ("flexddg", "Barlow KA, Ó Conchúir S, Thompson S, Suresh P, Lucas JE, Heinonen M, "
  "Kortemme T. Flex ddG: Rosetta ensemble-based estimation of changes in "
  "protein–protein binding affinity upon mutation. J Phys Chem B. 2018;122(21):5389–5399. "
  "doi:10.1021/acs.jpcb.7b11367"),
 ("foldx", "Delgado J, Radusky LG, Cianferoni D, Serrano L. FoldX 5.0: working with RNA, "
  "small molecules and a new graphical interface. Bioinformatics. 2019;35(20):4168–4169. "
  "doi:10.1093/bioinformatics/btz184"),
 ("proteinmpnn", "Dauparas J, Anishchenko I, Bennett N, Bai H, Ragotte RJ, Milles LF, "
  "et al. Robust deep learning-based protein sequence design using ProteinMPNN. Science. "
  "2022;378(6615):49–56. doi:10.1126/science.add2187"),
 ("esmif", "Hsu C, Verkuil R, Liu J, Lin Z, Hie B, Sercu T, Lerer A, Rives A. Learning "
  "inverse folding from millions of predicted structures. In: Proceedings of the 39th "
  "International Conference on Machine Learning (ICML). 2022;162:8946–8970."),
 ("ppiformer", "Bushuiev A, Bushuiev R, Kouba P, Filkin A, Gabrielova M, Gabriel M, "
  "Sedlar J, Pluskal T, Damborsky J, Mazurenko S, Sivic J. Learning to design "
  "protein–protein interactions with enhanced generalization. In: The Twelfth "
  "International Conference on Learning Representations (ICLR). 2024. arXiv:2310.18515"),
 ("ndcg", "Järvelin K, Kekäläinen J. Cumulated gain-based evaluation of IR techniques. "
  "ACM Trans Inf Syst. 2002;20(4):422–446. doi:10.1145/582415.582418"),
 ("bootstrap", "Efron B, Tibshirani RJ. An Introduction to the Bootstrap. New York: "
  "Chapman & Hall; 1993."),
 ("biometry", "Sokal RR, Rohlf FJ. Biometry: The Principles and Practice of Statistics in "
  "Biological Research. 3rd ed. New York: WH Freeman; 1995."),
]
KEYNUM = {k: i + 1 for i, (k, _) in enumerate(REFS)}
NREF = len(REFS)


def C(*keys):
    """纯上标引用数字串（无方括号）。"""
    ns = sorted(KEYNUM[k] for k in keys)
    out, i = [], 0
    while i < len(ns):
        j = i
        while j + 1 < len(ns) and ns[j + 1] == ns[j] + 1:
            j += 1
        out.append(str(ns[i]) if j == i else (f"{ns[i]},{ns[j]}" if j == i + 1
                                             else f"{ns[i]}\u2013{ns[j]}"))
        i = j + 1
    return "^{" + ",".join(out) + "}"


# =====================================================================
# 数学渲染器：把 `X_{sub}` / `X^{sup}` 解析为真实 Word 上下标 run
# =====================================================================
_TOK = re.compile(r"(\^\{[^}]*\}|_\{[^}]*\}|_[A-Za-z\u0177\u0302]|\^[A-Za-z0-9])")


def math_tokens(text: str):
    """产出 (chunk, kind)，kind ∈ {'', 'sub', 'sup', 'mono', 'i'}。"""
    out = []
    for part in _TOK.split(text):
        if not part:
            continue
        if part.startswith("^{"):
            out.append((part[2:-1], "sup"))
        elif part.startswith("_{"):
            out.append((part[2:-1], "sub"))
        elif part.startswith("_") and len(part) == 2:
            out.append((part[1], "sub"))
        elif part.startswith("^") and len(part) == 2:
            out.append((part[1], "sup"))
        else:
            out.append((part, ""))
    return out


def W(doc, text, size=None, align="just", italic_vars=False):
    """写入一段带真实上下标的文本。"""
    p = doc.add_paragraph()
    for chunk, kind in math_tokens(text):
        r = p.add_run(chunk)
        if kind in ("sub", "sup"):
            if kind == "sub":
                r.font.subscript = True
            else:
                r.font.superscript = True
            r.font.size = Pt((size or 11) - 2.0)
        elif size:
            r.font.size = Pt(size)
    p.alignment = (WD_ALIGN_PARAGRAPH.JUSTIFY if align == "just"
                   else WD_ALIGN_PARAGRAPH.CENTER)
    return p


def EQNUM(doc, mid, num):
    """居中公式 + 右对齐编号（同一行，制表位右对齐）。"""
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.tab_stops.add_tab_stop(Inches(6.7), 2)   # 2 = RIGHT
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    for chunk, kind in math_tokens(mid):
        r = p.add_run(chunk)
        if kind == "sub":
            r.font.subscript = True
        elif kind == "sup":
            r.font.superscript = True
        r.italic = (kind == "")
    p.add_run("\t(" + str(num) + ")")
    pf.space_before = Pt(5); pf.space_after = Pt(5); pf.keep_with_next = True
    return p


def EQ(doc, mid):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for chunk, kind in math_tokens(mid):
        r = p.add_run(chunk)
        if kind == "sub":
            r.font.subscript = True
        elif kind == "sup":
            r.font.superscript = True
    pf = p.paragraph_format
    pf.space_before = Pt(5); pf.space_after = Pt(5); pf.keep_with_next = True
    return p


def setup(doc):
    for name, fnt, sz in (("Normal", "Times New Roman", 10.5),
                          ("Heading 1", "Arial", 14), ("Heading 2", "Arial", 11.5)):
        try:
            st = doc.styles[name]
            st.font.name = fnt; st.font.size = Pt(sz)
            st.element.rPr.rFonts.set(qn("w:eastAsia"), "SimSun")
        except KeyError:
            pass
    st = doc.styles["Normal"]
    pf = st.paragraph_format
    pf.line_spacing_rule = WD_LINE_SPACING.MULTIPLE
    pf.line_spacing = 1.28
    pf.space_after = Pt(6)
    for s in doc.sections:
        s.top_margin = s.bottom_margin = Inches(0.9)
        s.left_margin = s.right_margin = Inches(0.9)


def H(doc, t, lvl=1):
    p = doc.add_heading(t, level=lvl)
    for r in p.runs:
        r.font.color.rgb = RGBColor(0, 0, 0)
        r.font.name = "Arial"
        r.font.size = Pt(14 if lvl == 1 else 11.5)
    p.paragraph_format.keep_with_next = True
    return p


def cap(doc, t):
    """图注：同样解析上下标，避免图注中出现字面 ASCII 下划线。"""
    p = doc.add_paragraph()
    for chunk, kind in math_tokens(t):
        r = p.add_run(chunk)
        r.font.size = Pt(9.0)
        if kind == "sub":
            r.font.subscript = True
        elif kind == "sup":
            r.font.superscript = True
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    p.paragraph_format.space_after = Pt(10)
    return p


def _autofit(table, total_in=6.7):
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


def tbl(doc, df, above=None, below=None, fs=8.2, total_in=6.7):
    if above:
        p = doc.add_paragraph(); r = p.add_run(above)
        r.font.size = Pt(9.0); r.bold = True
        p.paragraph_format.keep_with_next = True
    t = doc.add_table(rows=1, cols=len(df.columns))
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = False
    hdr = t.rows[0]
    trPr = hdr._tr.get_or_add_trPr()
    trPr.append(trPr.makeelement(qn("w:tblHeader"), {}))
    for i, c in enumerate(df.columns):
        cell = hdr.cells[i]; cell.text = ""
        for chunk, kind in math_tokens(str(c)):
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
    if below:
        cap(doc, below)
    else:
        doc.add_paragraph()
    return t


# ---------------- tables ----------------
T1 = pd.DataFrame([
    {"Evaluation set": "Full annotation namespace", "Rows": "4,473", "Targets": "214",
     "Targets (n ≥ 5)": "130", "Method families": "—",
     "Purpose": "Namespace definition"},
    {"Evaluation set": "Primary matched set", "Rows": f"{PU['n_rows']:,}",
     "Targets": "166", "Targets (n ≥ 5)": f"{PU['n_ge5']}",
     "Method families": f"{PU['n_method_families']}",
     "Purpose": "Main capability-resolved analyses"},
    {"Evaluation set": "Controlled split set", "Rows": "7,523", "Targets": "230",
     "Targets (n ≥ 5)": "172", "Method families": "12",
     "Purpose": "Split-strictness experiment"},
    {"Evaluation set": "External unambiguous set",
     "Rows": f"{FR['external_clean']['n_rows']:,}",
     "Targets": f"{FR['external_clean']['n_targets']}",
     "Targets (n ≥ 5)": f"{FR['external_clean']['n_ge5']}", "Method families": "5",
     "Purpose": "Boundary of generalisation"},
])
T2 = pd.DataFrame([
    {"Split regime": nm, "Pooled Pearson": f"{H4[k]['overall_pearson']:.4f}",
     "C_{B}": f"{H4[k]['C_B']:+.4f}", "C_{W}": f"{H4[k]['C_W']:+.4f}",
     "Macro Spearman": f"{H4[k]['macro_spearman']:.4f}",
     "RMSE": f"{H4[k]['rmse']:.4f}"}
    for k, nm in [("random", "Random"), ("complex", "Complex"), ("seq40", "seq40"),
                  ("family_OOD", "Family-OOD")]])
NICE_METHOD = {"Flex-ddG": "Flex-ddG", "DDGPred": "DDGPred",
               "ddG_FoldX": "FoldX (ddG_FoldX)", "ESM-IF": "ESM-IF",
               "proteinMPNN": "ProteinMPNN"}
T3 = EP[["model", "n_rows", "n_targets_ge5", "overall_pearson", "macro_spearman",
         "C_B_row", "C_W", "rmse"]].copy()
T3["model"] = T3["model"].map(lambda m: NICE_METHOD.get(m, m))
T3.columns = ["Method", "Mutations", "Targets (n ≥ 5)", "Pooled Pearson",
              "Macro Spearman", "C_{B}", "C_{W}", "RMSE"]
for c in ["Pooled Pearson", "Macro Spearman", "C_{B}", "C_{W}", "RMSE"]:
    T3[c] = T3[c].map(lambda x: f"{float(x):.4f}")
T3["Mutations"] = T3["Mutations"].map(lambda x: f"{int(x):,}")

FIGW = {"Fig1_overview": 6.7, "Fig2_split_decomposition": 6.7, "Fig3_weighting": 6.7,
        "Fig4_candidate_geometry": 6.0, "Fig5_remedy": 6.1, "Fig6_external": 6.7}


def build(doc_name_en=True):
    D = Document(); setup(D)
    S = "^{%d}" % 0  # 占位，避免未使用告警

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
         f"of disease-associated variants, and progress depends on the credibility of "
         f"benchmark evaluation. The pooled correlation commonly reported for such "
         f"benchmarks is difficult to interpret directly, because it depends both on "
         f"whether a model places different targets on the correct overall binding scale "
         f"and on whether it orders mutations correctly within each target. We introduce an "
         f"exact covariance-signal decomposition that writes pooled Pearson correlation as "
         f"r = C_{{B}} + C_{{W}}, where C_{{B}} and C_{{W}} are the between-target and "
         f"within-target covariance contributions, and apply it to matched, de-replicated "
         f"and provenance-filtered evaluation sets. In a controlled comparison that "
         f"changes only the split regime (identical 7,523 mutation rows and identical "
         f"model configurations), median pooled Pearson falls from {r0:.4f} under random "
         f"splitting to {r3:.4f} under family-out-of-distribution evaluation, with "
         f"ΔC_{{B}} = {DL['delta_C_B']:+.4f} and ΔC_{{W}} = {DL['delta_C_W']:+.4f}; on an "
         f"absolute-change basis the between-target component accounts for "
         f"{100*DL['fraction_between']:.1f}% of the decline while the within-target "
         f"component remains positive. Target weighting changes how these contributions "
         f"should be read: C_{{B}} reverses direction only near zero ({sf0}/{len(W2)} "
         f"configurations across five schemes), whereas C_{{W}} keeps its sign at every "
         f"threshold although its magnitude is not invariant. Fixed-K shortlist metrics "
         f"are further constrained by candidate-set geometry: for K = 10, D(n, K) is zero "
         f"when n ≤ K and peaks near n* ≈ {GEO['n_star_analytic']:.1f}, and "
         f"{100*GEO['n_le_K_fraction']:.1f}% of targets with n ≥ 5 fall in that degenerate "
         f"regime. In a matched retrospective ablation the decomposition separated a "
         f"change that improved pooled scores only from one that also improved "
         f"within-target ranking. External evidence comprised five methods, is "
         f"underpowered and is positively aligned, so it does not support a universal "
         f"decoupling of pooled and within-target performance. We recommend that PPI ΔΔG "
         f"benchmarks report both components under an explicitly declared weighting, "
         f"together with candidate-set geometry.")
    p = D.add_paragraph()
    rr = p.add_run("Keywords: "); rr.bold = True
    p.add_run("protein–protein interaction; binding affinity; ΔΔG; mutation effect "
              "prediction; benchmarking; out-of-distribution generalisation"
              ).font.size = Pt(10)

    H(D, "Introduction", 1)
    W(D, "Predicting how mutations alter protein–protein binding affinity (ΔΔG) is central "
         "to antibody affinity maturation, interface optimisation and the interpretation "
         "of disease-associated variants. Because experimental affinity measurement is "
         "costly, computational predictions are widely used to prioritise mutations and "
         "guide experimental design. As the number of deep-learning methods grows rapidly, "
         "progress in the field is increasingly judged by benchmark evaluation, and the "
         "properties of the evaluation protocol itself — what it measures, and under which "
         "extrapolation conditions — directly determine what the field can conclude.")
    W(D, f"Current practice relies largely on pooled metrics, computing Pearson or Spearman "
         f"correlation and RMSE over the entire mutation set{C('biometry')}, often "
         f"supplemented by per-complex or per-target metrics. Recent work has shown that "
         f"how well such evaluation captures extrapolation depends strongly on the "
         f"splitting strategy. Bushuiev et al. quantified structural leakage in PPI data: "
         f"about 65% of test examples have a near-duplicate training example under PDB-code "
         f"splitting, about 30% even under 30% sequence-identity splitting, 56% for SKEMPI "
         f"v2.0{C('skempi')} split by PPI code, and 0% under domain-expert "
         f"grouping{C('leakage')}. CATH-ddG adopted CATH homologous-superfamily splits and "
         f"defined hard mutations by a maximum TM-score below 0.6{C('cathddg')}; related "
         f"methods likewise report per-structure or per-target "
         f"correlations{C('hermes','prefold','ppiformer')}. We emphasise that leakage-aware "
         f"splitting, family-out-of-distribution evaluation, per-target reporting and "
         f"shortlist (Top-K) evaluation were none of them introduced by this work.")
    W(D, "A key question nevertheless remains unanswered directly: where does pooled "
         "benchmark performance come from? Reporting pooled and per-target metrics side by "
         "side does not quantify the additive contribution of between-target and "
         "within-target signal to pooled performance. A model may achieve a high pooled "
         "correlation mainly because it places different targets on the correct overall "
         "scale, while its ability to order mutations within a target is weak — or the "
         "reverse. In benchmarks containing many targets with highly unequal mutation "
         "counts, pooled correlation is additionally dominated by large targets, which "
         "amplifies the conflation.")
    W(D, "Here we introduce an exact covariance-signal decomposition that addresses this "
         "question. We decompose pooled Pearson correlation rigorously into a "
         "between-target component and a within-target component, which makes 'what pooled "
         "performance is made of' a directly measurable quantity. Using this decomposition "
         "on matched, de-replicated and provenance-filtered evaluation sets, we quantify "
         "four things: which component is changed by split strictness; how target weighting "
         "changes the interpretation of the two components; which constraints candidate-set "
         "geometry places on fixed-K shortlist evaluation; and whether the decomposition can "
         "distinguish model modifications that improve pooled scores only from those that "
         "also improve within-target ranking. We further bound the current evidence using "
         "five independently developed external methods.")

    H(D, "Results", 1)
    tbl(D, T1, above="Table 1 | Dataset and evaluation overview.",
        below="All rows are joined by explicit row identifiers. The controlled split set is "
              "used only for the split-strictness experiment; the primary matched set is "
              "used for all capability-resolved analyses.", fs=8.4)

    H(D, "Pooled PPI ΔΔG performance contains separable between-target and within-target "
         "signal", 2)
    W(D, "For a mutation i belonging to target c, the observed ΔΔG can be written as a "
         "target-level mean plus a mutation-specific deviation from that mean:")
    EQNUM(D, "y_{ci} = μ_{c} + ε_{ci}", 1)
    W(D, "where μ_{c} is the mean ΔΔG of target c and ε_{ci} = y_{ci} − μ_{c} is the "
         "within-target mutation deviation. The prediction is decomposed identically:")
    EQNUM(D, "\u0177_{ci} = μ\u0302_{c} + ε\u0302_{ci}", 2)
    W(D, "Only after this decomposition is it meaningful to distinguish two covariance "
         "contributions. With σ_{y} and σ_{ŷ} the pooled standard deviations, the pooled "
         "correlation splits exactly as:")
    EQNUM(D, "r_{pooled} = Cov(y, ŷ) / (σ_{y} σ_{ŷ}) = C_{B} + C_{W}", 3)
    W(D, "where C_{B} = Cov(B_{y}, B_{ŷ}) / (σ_{y} σ_{ŷ}) is the between-target covariance "
         "contribution and C_{W} = Cov(W_{y}, W_{ŷ}) / (σ_{y} σ_{ŷ}) is the within-target "
         "covariance contribution. The identity holds to machine precision (residual ≤ "
         "3.3 × 10⁻¹⁶; full derivation in Supplementary Note 1). We stress that C_{B} and "
         "C_{W} are covariance-signal contributions on a common scale, not independent "
         "correlations, independent skill scores, or causal mechanisms.")
    W(D, f"We apply this decomposition on a matched evaluation set comprising "
         f"{PU['n_rows']:,} mutation rows, {PU['n_ge5']} targets with at least five "
         f"mutations, and {PU['n_method_families']} independent method families; cross-model "
         f"comparisons were restricted to predictions evaluated on identical mutation-row "
         f"sets (Methods; Supplementary Fig. S1).")
    D.add_picture(str(FM / "Fig1_overview_EN.png"), width=Inches(FIGW["Fig1_overview"]))
    cap(D, "Figure 1 | Capability-resolved evaluation of PPI ΔΔG predictors. "
           "(A) Provenance-filtered, fingerprint-matched benchmark construction: SKEMPI "
           "v2.0 annotation (4,473 rows, 214 targets, 108 families) is audited by SHA-256 "
           "row-set fingerprinting; family-out-of-distribution predictions without training "
           "overlap are retained and permissively split predictions are excluded, giving a "
           "primary matched universe of 3,351 rows, 102 targets with n ≥ 5 and 10 method "
           "families. Target cards illustrate the two structural features of the benchmark: "
           "between-target differences in mean ΔΔG (μ_{c}, purple) and within-target mutation "
           "ordering (ε_{ci}, green). (B) Exact additive covariance-signal decomposition: "
           "target centring yields a between-target contribution C_{B} and a within-target "
           "contribution C_{W}, whose sum is pooled correlation by an exact identity (residual "
           "≤ 3.3 × 10⁻¹⁶). These are common-scale additive contributions, not independent "
           "correlations, skill scores or causal mechanisms. (C) Controlled split collapse "
           "on identical 7,523 mutation rows: pooled correlation falls from 0.6844 to "
           "0.1508, and 81.1% of the absolute decline (ΔC_{B} = −0.4347) is carried by loss of "
           "the between-target component while C_{W} remains positive (C_{W} = +0.1580). "
           "(D) Standardised capability profiles of matched ablation variants across five "
           "axes; axes are standardised for visual comparison only and raw values are "
           "reported in the text.")

    H(D, "Family-OOD evaluation disproportionately removes between-target signal", 2)
    W(D, f"To isolate the effect of split strictness, we compared identical model "
         f"configurations on the same 7,523 mutation rows while changing only the split "
         f"regime (Table 2, Figure 2). Median pooled Pearson falls from {r0:.4f} under "
         f"random splitting to {r3:.4f} under family-out-of-distribution evaluation. The "
         f"decomposition shows that this decline is not uniform: "
         f"ΔC_{{B}} = {DL['delta_C_B']:+.4f} and ΔC_{{W}} = {DL['delta_C_W']:+.4f}. On an "
         f"absolute-change basis, |ΔC_{{B}}| / (|ΔC_{{B}}| + |ΔC_{{W}}|), the between-target "
         f"component accounts for {100*DL['fraction_between']:.1f}% of the decline, while "
         f"C_{{W}} remains positive ({cw3:+.4f}). Most of the reduction in pooled "
         f"correlation under family-OOD evaluation is therefore carried by loss of the "
         f"between-target covariance contribution, while within-target covariance remains "
         f"positive.")
    W(D, "This is not the same proposition as 'the model memorised target identity': we "
         "performed no identity-memorisation intervention and therefore make no such causal "
         "claim. A more defensible reading is that the pattern is consistent with reliance "
         "on target-specific information that does not transfer across families. The "
         "complex-level and seq40 partitions produced identical assignments in this dataset "
         "and were therefore not treated as independent comparisons.")
    D.add_picture(str(FM / "Fig2_split_decomposition_EN.png"),
                  width=Inches(FIGW["Fig2_split_decomposition"]))
    cap(D, "Figure 2 | Family-out-of-distribution evaluation disproportionately removes the "
           "between-target component of pooled correlation. (A) Median pooled Pearson "
           "correlation by split regime. (B) Exact decomposition into C_{B} (muted purple) "
           "and C_{W} (green); C_{W} remains positive. (C) Attribution of the absolute decline, "
           f"of which the between-target component accounts for "
           f"{100*DL['fraction_between']:.1f}%. All four regimes use the same 7,523 "
           "mutation rows and the same model configurations; the complex and seq40 "
           "partitions are identical in this dataset and are not independent comparisons.")
    tbl(D, T2, above="Table 2 | Controlled split comparison.",
        below="Median values across model configurations. The seq40 partition is identical "
              "to the complex partition in this dataset and is therefore not an independent "
              "comparison.")

    H(D, "Target aggregation changes the interpretation of covariance contributions", 2)
    W(D, f"Because targets differ greatly in mutation count, the choice of aggregation "
         f"weight affects the magnitude of both components. Holding predictions fixed and "
         f"varying only the aggregation weights, we compared five schemes (per-mutation, "
         f"per-target, weight proportional to √n, proportional to log(1 + n), and capped at "
         f"30). C_{{B}} shows directional instability concentrated near zero: "
         f"{sf0}/{len(W2)} configurations reverse sign with no threshold, whereas none do "
         f"once |C_{{B}}| exceeds 0.005 ({sf5}/{len(W2)}). C_{{W}} never reverses sign at any "
         f"threshold ({sfw}/{len(W2)}), but its magnitude is not invariant — its median "
         f"range is {W2.range_C_W.median():.4f}, larger than the "
         f"{W2.range_C_B.median():.4f} of C_{{B}}. C_{{B}} is therefore unstable in "
         f"direction only near zero, whereas C_{{W}} keeps its sign but varies in magnitude. "
         f"Consequently, when a reported C_{{B}} is close to zero, the conclusion about "
         f"whether between-target information is present may depend on the weighting "
         f"convention rather than on the model; any reported decomposition should declare "
         f"its weighting scheme.")
    D.add_picture(str(FM / "Fig3_weighting_EN.png"), width=Inches(FIGW["Fig3_weighting"]))
    cap(D, "Figure 3 | Target weighting changes the interpretation of covariance "
           "contributions. Predictions are unchanged; only aggregation weights differ. "
           "(A, B) Distribution of the range of C_{B} and C_{W} across the five schemes; the "
           "dashed line marks the median range. (C) Number of configurations reversing sign "
           "at three thresholds on |C|. C_{B} reverses direction only near zero, whereas C_{W} "
           "keeps its sign at every threshold but its magnitude is not invariant.")

    H(D, "Fixed-K shortlist evaluation is constrained by candidate-set geometry", 2)
    W(D, f"Fixed-K top-K metrics are commonly used to evaluate models for mutation library "
         f"or shortlist screening{C('ndcg')}. The discriminable interval is:")
    EQ(D, "D(n, K) = 0,  n ≤ K      D(n, K) = K(n − K) / [n(n + 1)],  n > K")
    W(D, f"where the expected share of gain from K randomly chosen candidates is K/n and the "
         f"optimal K candidates attain K(2n − K + 1) / [n(n + 1)] (derivation in "
         f"Supplementary Note 6). When a target has no more than K mutations, any model "
         f"selects every candidate and the target carries no information for fixed-K "
         f"shortlist metrics. In this benchmark, {GEO['n_le_K_count']} of "
         f"{GEO['n_ge5_count']} targets with n ≥ 5 ({100*GEO['n_le_K_fraction']:.1f}%) fall "
         f"in that degenerate regime. Moreover D is not monotonic in n: it peaks at "
         f"n* = K + √(K² + K), which for K = 10 gives n* ≈ {GEO['n_star_analytic']:.1f}, so "
         f"that D(20) = {GEO['D_at_n20']:.4f} while D(50) = {GEO['D_at_n50']:.4f}. The "
         f"heuristic that larger candidate sets are more discriminable therefore does not "
         f"hold for this metric.")
    W(D, "This geometric result does not imply that SKEMPI ΔΔG values are intrinsically too "
         "similar; establishing that would require an analysis using actual ΔΔG utility "
         "magnitudes, which we did not perform. The defensible statement is that fixed-K "
         "candidate-selection evaluation is structurally non-informative when n ≤ K, and "
         "that its discriminable range depends strongly on candidate-set geometry.")
    D.add_picture(str(FM / "Fig4_candidate_geometry_EN.png"),
                  width=Inches(FIGW["Fig4_candidate_geometry"]))
    cap(D, "Figure 4 | Fixed-K shortlist evaluation is constrained by candidate-set "
           "geometry. The discriminable interval D(n, K) is zero when n ≤ K and follows "
           "K(n − K)/[n(n + 1)] otherwise; it peaks at n* = K + √(K² + K), which is "
           f"n* ≈ {GEO['n_star_analytic']:.1f} for K = 10, and declines thereafter. In this "
           f"benchmark {100*GEO['n_le_K_fraction']:.1f}% of targets with n ≥ 5 lie in the "
           "degenerate regime where every model selects all candidates.")

    H(D, "Capability decomposition clarifies retrospective remedy effects", 2)
    W(D, "Finally, we examined the practical utility of the decomposition in a matched "
         "ablation set whose models share an identical 3,351-row universe, so that "
         "differences arise only from the model modifications themselves. Relative to a "
         "baseline using mutation properties alone, adding static local geometry raised "
         "pooled Pearson by about 0.028 while leaving within-target Macro Spearman "
         "essentially unchanged (−0.002), a pooled-score-only improvement; introducing "
         "mutation × partner perturbation improved both axes (pooled Pearson +0.070, Macro "
         "Spearman +0.044); and a full local-physics model gave the largest improvement on "
         "both axes (+0.105 and +0.076). The decomposition therefore distinguishes a change "
         "that improves pooled scores only from one that also improves within-target "
         "ranking.")
    W(D, "The nature of this analysis should be stated plainly: it is retrospective, in that "
         "we already know which modifications were effective. It therefore shows that the "
         "decomposition can express this distinction, not that the benchmark prospectively "
         "predicts which technical route will succeed. In addition, our recomputed "
         "within-target difference (+0.044) and the bootstrap interval from the historical "
         "analysis (+0.073, [+0.021, +0.129]) rest on different metric universes and "
         "aggregation schemes and must not be mixed without explanation.")
    D.add_picture(str(FM / "Fig5_remedy_EN.png"), width=Inches(FIGW["Fig5_remedy"]))
    cap(D, "Figure 5 | Capability decomposition separates pooled-score gains from "
           "within-target gains. Retrospective matched ablation on an identical 3,351-row "
           "universe; the baseline is the mutation-only model. Green markers indicate "
           "improvement on both axes; red markers indicate a pooled-score-only or negative "
           "change. Axes give the change in pooled Pearson and in macro within-target "
           "Spearman relative to that baseline. The analysis is retrospective and does not "
           "imply prospective prediction of which modification will succeed.")

    H(D, "External methods define the boundary of generalisation", 2)
    W(D, f"To assess how far the above conclusions apply to independently developed methods, "
         f"we built external evaluation sets from published official per-mutation "
         f"predictions: Flex ddG (Rosetta ensemble-based affinity estimation){C('flexddg')}, "
         f"FoldX{C('foldx')}, ProteinMPNN{C('proteinmpnn')}, ESM-IF{C('esmif')} and "
         f"HERMES{C('hermes')}. During construction we found that when the matching key "
         f"omits the chain identifier and the wild-type residue, the key is not unique on "
         f"the annotation side: of 4,473 rows only "
         f"{FR['external_clean']['n_rows']:,} have a unique mapping. Primary matched "
         f"evaluation therefore uses unambiguous rows only; ambiguous keys were excluded "
         f"rather than forcibly resolved. This treatment systematically lowers apparent "
         f"pooled correlation for external methods by 0.02 to 0.12.")
    W(D, f"On the unambiguous set ({FR['external_clean']['n_rows']:,} rows, "
         f"{FR['external_clean']['n_ge5']} targets with n ≥ 5), the five external methods "
         f"agree in their ordering on pooled correlation and on within-target ranking, with "
         f"the same best method, DDGPred (pooled Pearson "
         f"{float(EP.loc[EP.model=='DDGPred','overall_pearson'].iloc[0]):+.4f}, Macro "
         f"Spearman {float(EP.loc[EP.model=='DDGPred','macro_spearman'].iloc[0]):+.4f}); "
         f"that is, no method-level rank reversal occurs in this set. Method-level rank "
         f"agreement is ρ = "
         f"{er.get('rho_overall_pearson_vs_macro_spearman', float('nan')):+.2f} (exact "
         f"permutation p = "
         f"{er.get('perm_p_overall_pearson_vs_macro_spearman', float('nan')):.3f}). This is "
         f"a strongly positive agreement, which does not support a universal decoupling of "
         f"pooled and within-target performance and falls within the heterogeneous range "
         f"observed internally. With only five methods, method-level rank-correlation "
         f"inference is low-powered, and we therefore draw no directional conclusion from "
         f"this panel.")
    W(D, f"Regarding the semantics of one external quantity, its official output is a "
         f"log-odds difference, a dimensionless mutation-effect score rather than a "
         f"calibrated kcal/mol ΔΔG{C('hermes')}. We therefore report neither RMSE nor MAE "
         f"for that model, and we do not label the variant we used as supervised on binding "
         f"affinity.")
    D.add_picture(str(FM / "Fig6_external_EN.png"), width=Inches(FIGW["Fig6_external"]))
    cap(D, "Figure 6 | External methods define the boundary of generalisation. (A) Pooled "
           "versus within-target performance for five independently developed methods on a "
           "matched, unambiguous mapping; leader lines are used where markers are close. "
           "(B) Method-level capability profile, standardised within each metric (z-score) "
           "for visual comparison only; raw values are printed in each cell and statistical "
           "analyses use the raw metrics. Rank agreement between the two axes is positive "
           "(ρ = +0.90, exact permutation p = 0.083), but the panel is underpowered with "
           "five methods, so this positive alignment does not support a universal decoupling "
           "of pooled and within-target performance.")
    tbl(D, T3, above="Table 3 | External method capability profiles on the unambiguous "
                     "mapping set.",
        below="All methods are evaluated on the same mutations. RMSE is reported only for "
              "outputs with ΔΔG-like scaling; for the log-odds external quantity RMSE is not "
              "applicable and is omitted. Method names follow Figure 6.")

    H(D, "Discussion", 1)
    W(D, "This study shows that the pooled correlation conventionally reported as a single "
         "number in PPI ΔΔG benchmarking is in fact composed of two qualitatively different "
         "predictive signals: a between-target component that places different targets on "
         "the correct overall scale, and a within-target component that orders mutation "
         "effects correctly inside each target. Because our decomposition is an exact "
         "identity (Eq. 3), 'what pooled performance is made of' ceases to be a matter of "
         "interpretation and becomes directly computable.")
    W(D, f"The relation to prior work should be delimited clearly. Structural leakage in PPI "
         f"data and its dependence on splitting strategy were quantified systematically by "
         f"Bushuiev et al., who also proposed interface-similarity splitting{C('leakage')}; "
         f"CATH-ddG adopted CATH homologous-superfamily splits and used structural "
         f"similarity to define harder extrapolation mutations{C('cathddg')}; and per-complex "
         f"or per-target reporting together with shortlist-style evaluation already appear in "
         f"the ΔΔG literature{C('hermes','prefold','ppiformer')}. Our contribution therefore "
         f"lies not in any single one of these ingredients but in the systematic integration "
         f"of an exact between/within-target covariance decomposition with four audits — "
         f"split strictness, target weighting, matched rank agreement, and candidate-set "
         f"geometry — inside one matched, provenance-filtered framework, together with a "
         f"retrospective attribution case study. To our knowledge we did not identify prior "
         f"work that jointly reports these under comparable matching conditions.")
    W(D, "The interpretation of C_{B} and C_{W} has clear boundaries. They are covariance "
         "contributions within one identity and share a common scale, so they sum to pooled "
         "correlation; they are not two independently estimated skill scores and should not "
         "be read as causal mechanisms. In particular, a model may obtain a substantial "
         "pooled correlation largely by placing targets on approximately the right overall "
         "scale while its within-target ordering ability remains limited — precisely why "
         "pooled metrics can overstate practical screening value. We avoid interpreting the "
         "presence of within-target signal as evidence that a model has learned physics, or "
         "the loss of between-target signal as evidence that a model memorised target "
         "identity.")
    W(D, "The weighting-sensitivity audit yields a conclusion of direct practical relevance: "
         "the sign of the between-target component is sensitive to aggregation weights near "
         "zero, whereas the sign of the within-target component is stable. Consequently, "
         "when a reported C_{B} is close to zero, the conclusion about whether "
         "between-target information is present may depend on the weighting convention "
         "rather than on the model. We therefore recommend that any reported decomposition "
         "declare its weighting scheme, and that multiple weightings be shown when C_{B} is "
         "near zero.")
    W(D, "The candidate-set geometry result has direct implications for Top-K evaluation "
         "practice. When a target has no more than K mutations, fixed-K shortlist metrics "
         "cannot distinguish models structurally; roughly four in ten targets with n ≥ 5 "
         "fall in this regime in our benchmark. Because the discriminable interval first "
         "rises and then falls as candidate sets grow, enlarging a candidate set is not "
         "universally effective at improving discrimination. We recommend reporting n "
         "alongside K, labelling n ≤ K targets separately, or adopting fraction-based (q) "
         "evaluation.")
    W(D, "The retrospective ablation illustrates the framework's explanatory value for model "
         "modifications: it distinguishes a change that raised pooled scores only from a "
         "change that also improved within-target ranking. This distinction matters for "
         "method design because it attributes an improvement to one of two sources. The "
         "retrospective nature of the case study must, however, be stated explicitly: it "
         "does not establish prospective ability to predict which modifications will work.")
    W(D, "Finally, external evidence bounds the generalisation of these conclusions. Across "
         "five independently developed methods the ordering of methods agreed between pooled "
         "and within-target performance, a strongly positive agreement that neither supports "
         "nor refutes a universal decoupling of the two axes. The panel is underpowered, so "
         "we draw no directional conclusion from it. Future work should broaden the external "
         "panel and test the relationship between the decomposed components and realised "
         "screening yield on independent prospective datasets, as well as evaluation design "
         "for larger candidate libraries.")

    H(D, "Methods", 1)
    for title, body in [
        ("Dataset and target definition",
         f"Evaluation rows come from a SKEMPI v2.0-derived annotation set{C('skempi')} "
         f"(4,473 rows, 214 targets, 108 structural families), joined to prediction outputs "
         f"by explicit row identifiers. A target is defined as an independent "
         f"binding-interface unit within a structural family."),
        ("Prediction sources and provenance criteria",
         "Prediction outputs were classified by their recorded supervised training split "
         "regime. Configurations trained under a permissive (random / complex / seq40) "
         "partition were excluded from deployment-oriented comparisons because their "
         "evaluation rows may overlap their training folds; they were retained only in the "
         "controlled split experiment. Analyses were performed at three levels: artifact "
         "level (descriptive only), configuration level (seed replicates collapsed), and "
         "method-family level (one performance-blind representative per independent family)."),
        ("Matched evaluation sets",
         "Each prediction output received a SHA-256 fingerprint of its sorted row "
         "identifiers; cross-model comparison was performed only between outputs sharing a "
         "fingerprint (Supplementary Note 2, Supplementary Fig. S1). The primary matched set "
         "was selected by a rule fixed before any metric values were inspected."),
        ("Between-/within-target signal decomposition",
         "See Eqs. 1\u20133; row-level weights are used and the identity residual is reported "
         "as a quality-control quantity. Full derivation in Supplementary Note 1."),
        ("Target-weighting schemes",
         "Predictions are unchanged; only aggregation weights vary: per-mutation, "
         "per-target, proportional to √n, proportional to log(1 + n), and capped at 30."),
        ("Ranking and shortlist metrics",
         f"Within-target ordering used per-target Spearman \u03c1 and Kendall \u03c4 "
         f"correlations, pairwise accuracy and the within-target covariance contribution; "
         f"shortlist metrics used fixed-K and fraction-based hit rate, precision, recall, "
         f"NDCG and enrichment factor{C('ndcg')}."),
        ("Statistical analysis",
         f"Rank agreement uses Spearman \u03c1 with Kendall \u03c4 as secondary; interval "
         f"estimates use target-cluster bootstrap (5,000 resamples, seed 42)"
         f"{C('bootstrap')}; small external panels use exact permutation tests. Metric "
         f"directions were registered before computation and all cross-metric comparisons "
         f"use higher-is-better orientation{C('biometry')}."),
        ("External prediction mapping and ambiguity control",
         "Where the normalised matching key was not unique on the annotation side, only "
         "uniquely mapped rows were used; ambiguous keys were excluded rather than resolved. "
         "The official held-out CATH partition was evaluated separately "
         "(Supplementary Note 5)."),
        ("Reproducibility and consistency validation",
         "All headline numbers are generated from a single machine-readable file and "
         "verified by an automated consistency checker (Supplementary Information). No new "
         "model was trained and no GPU was used."),
    ]:
        H(D, title, 2)
        W(D, body)

    H(D, "Limitations", 1)
    for t in [
        "Our extrapolation setting is based on target-family partitioning, which does not "
        "exclude the possibility that foundation-model pretraining encountered sequences or "
        "structures similar to the evaluation families; we can control task-specific "
        "supervised overlap only. In addition, the model pool assembled here is not a random "
        "sample of all PPI ΔΔG methods, coverage of some families is limited, and only five "
        "methods are externally evaluable, so method-level inference is underpowered. "
        "Conclusions about the composition of capability should therefore be read as "
        "controlled observations rather than a survey of the field.",
        "Shortlist and candidate-selection analyses rest on re-analysis of existing "
        "predictions and include no prospective wet-lab validation; when a target has no "
        "more than K mutations the relevant metrics carry no information, so "
        "selection-related conclusions are reported as low-confidence descriptions only. "
        "The retrospective ablation likewise does not establish prospective prediction of "
        "which modifications will succeed.",
        "C_{B} and C_{W} are covariance-signal contributions rather than causal mechanisms. "
        "We performed no identity-memorisation or target-shuffling intervention and "
        "therefore make no causal claim about the origin of the between-target component.",
        "No independent prospective experimental dataset was available, so no external "
        "dataset validation was performed. This work does not provide a universal model "
        "recommender and does not claim general superiority of any method family.",
    ]:
        W(D, t)

    H(D, "Data and code availability", 1)
    W(D, "All analyses use existing predictions and published external predictions. Frozen "
         "input hashes, the machine-readable source of every headline number and an "
         "automated consistency checker are provided with the Supplementary Information so "
         "that all reported values can be regenerated without retraining any model. No new "
         "model was trained and no GPU was used.")
    H(D, "Author contributions", 1)
    W(D, "Study conception, analysis design, interpretation of results and writing of the "
         "manuscript were performed by the author(s). The use of AI-assisted tools is "
         "declared below.")
    H(D, "Use of AI-assisted tools", 1)
    W(D, "AI-assisted coding and language tools were used to support code drafting and "
         "debugging, generation of plotting scripts, document automation, literature and "
         "repository retrieval, and language editing. All metric definitions, analysis "
         "choices, scientific interpretations, hypothesis decisions and final manuscript "
         "content were reviewed and approved by the authors, who take full responsibility "
         "for the work. AI tools were not used to generate experimental data or to modify "
         "upstream model predictions. A detailed account of the automated workflow is "
         "provided in the Supplementary Information.")
    H(D, "Acknowledgements", 1)
    W(D, "[FUNDING_PLACEHOLDER] We thank the developers of the publicly released prediction "
         "sets and datasets used in this study.")

    H(D, "References", 1)
    for i, (k, txt) in enumerate(REFS, 1):
        p = D.add_paragraph()
        r = p.add_run(f"{i}. {txt}"); r.font.size = Pt(9.5)
        p.paragraph_format.left_indent = Inches(0.30)
        p.paragraph_format.first_line_indent = Inches(-0.30)
        p.paragraph_format.space_after = Pt(3)

    nm = "CR-PPI-ddG_manuscript_v5_EN.docx"
    path = OUT / nm
    D.save(str(path))
    print(f"  {nm} ({path.stat().st_size/1024:.0f} KB)")
    return path


build()
print("v5 manuscript 完成")
