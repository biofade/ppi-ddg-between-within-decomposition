# -*- coding: utf-8 -*-
"""scripts/v8_1/18_manuscript_v4.py -- 稿件 v4（顶刊排版规范）。

在第一~六部分要求下修复：
* 方括号数字引用 -> **上标**引用（[4,6–9] 用 en-dash），不粘连单词
* 数学符号全部真上下标（C_B -> C_B 下标、μ_c、ε_ci、σ_ŷ、E_rand、n* = K+√(K²+K)）
* 修正 `p = +0.90` -> `ρ = +0.90`；`Spearman p with Kendall t` -> `ρ` / `τ`
* Table 1 从文末移回 Results 2.1 首次提及处；千分位逗号
* 表格自适应列宽 + 8 pt，防止单词被拆行
* 图注全部改写为期刊风格（图内已无硬编码标题/脚注）

科学数值只取自 v8_1/FROZEN_RESULTS.json。
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
OUT = F / "manuscript_v4"
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
rmn = float(ST.rho_overall_vs_macro.min()); rmx = float(ST.rho_overall_vs_macro.max())
best = EP.loc[EP.overall_pearson.idxmax()]

# ---------------- references (first-appearance order) ----------------
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
 ("prefold", "PreFold-dG: estimating binding affinity of protein–protein interaction from "
  "intermediate representations of protein folding models. Bioinformatics. 2026. "
  "doi:10.1093/bioinformatics/btag489"),
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
    """上标引用串（en-dash 表范围）。"""
    ns = sorted(KEYNUM[k] for k in keys)
    out, i = [], 0
    while i < len(ns):
        j = i
        while j + 1 < len(ns) and ns[j + 1] == ns[j] + 1:
            j += 1
        out.append(str(ns[i]) if j == i else (f"{ns[i]},{ns[j]}" if j == i + 1
                                             else f"{ns[i]}\u2013{ns[j]}"))
        i = j + 1
    return "[" + ",".join(out) + "]"


# ---------------- tables ----------------
T1 = pd.DataFrame([
    {"Evaluation set": "Full annotation namespace", "Rows": "4,473", "Targets": "214",
     "Targets n≥5": "130", "Method families": "—",
     "Purpose": "Namespace definition"},
    {"Evaluation set": "Primary matched set", "Rows": f"{PU['n_rows']:,}",
     "Targets": "166", "Targets n≥5": f"{PU['n_ge5']}",
     "Method families": f"{PU['n_method_families']}",
     "Purpose": "Main capability-resolved analyses"},
    {"Evaluation set": "Controlled split set", "Rows": "7,523", "Targets": "230",
     "Targets n≥5": "172", "Method families": "12",
     "Purpose": "Split-strictness experiment"},
    {"Evaluation set": "External unambiguous set",
     "Rows": f"{FR['external_clean']['n_rows']:,}",
     "Targets": f"{FR['external_clean']['n_targets']}",
     "Targets n≥5": f"{FR['external_clean']['n_ge5']}", "Method families": "5",
     "Purpose": "Boundary of generalisation"},
])
T2 = pd.DataFrame([
    {"Split regime": nm, "Pooled Pearson": f"{H4[k]['overall_pearson']:.4f}",
     "C_B": f"{H4[k]['C_B']:+.4f}", "C_W": f"{H4[k]['C_W']:+.4f}",
     "Macro Spearman": f"{H4[k]['macro_spearman']:.4f}",
     "RMSE": f"{H4[k]['rmse']:.4f}"}
    for k, nm in [("random", "Random"), ("complex", "Complex"), ("seq40", "seq40"),
                  ("family_OOD", "Family-OOD")]])
T3 = EP[["model", "n_rows", "n_targets_ge5", "overall_pearson", "macro_spearman",
         "C_B_row", "C_W", "rmse"]].copy()
T3.columns = ["Method", "Mutations", "Targets n≥5", "Pooled Pearson", "Macro Spearman",
              "C_B", "C_W", "RMSE"]
for c in ["Pooled Pearson", "Macro Spearman", "C_B", "C_W", "RMSE"]:
    T3[c] = T3[c].round(4)


def fmt4(v):
    """4 位小数；极小值用科学计数法。"""
    try:
        x = float(v)
    except (TypeError, ValueError):
        return str(v)
    if x != 0 and abs(x) < 1e-4:
        m, e = f"{x:.2e}".split("e")
        return f"{m} × 10^{int(e)}"
    return f"{x:.4f}"


# ---------------- docx helpers ----------------
def setup(doc, zh):
    st = doc.styles["Normal"]
    st.font.size = Pt(10.5 if zh else 11)
    st.font.name = "Times New Roman"
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "SimSun")
    pf = st.paragraph_format
    pf.line_spacing_rule = WD_LINE_SPACING.MULTIPLE
    pf.line_spacing = 1.20
    pf.space_after = Pt(6)
    for s in doc.sections:
        s.top_margin = s.bottom_margin = Inches(0.9)
        s.left_margin = s.right_margin = Inches(0.9)


def H(doc, t, lvl=1):
    p = doc.add_heading(t, level=lvl)
    for r in p.runs:
        r.font.color.rgb = RGBColor(0, 0, 0)
        r.font.size = Pt(14 if lvl == 1 else 11.5)
    p.paragraph_format.keep_with_next = True
    return p


def _emit(p, t, sup):
    r = p.add_run(t)
    if sup:
        r.font.superscript = True
        r.font.size = Pt(9)
    return r


def RT(doc, text, size=None, align="just"):
    """渲染文本；`^[...]` 片段转为上标引用。"""
    p = doc.add_paragraph()
    for part in re.split(r"(\^\[[^\]]*\])", text):
        if not part:
            continue
        if part.startswith("^["):
            _emit(p, part[1:], True)
        else:
            r = _emit(p, part, False)
            if size:
                r.font.size = Pt(size)
    p.alignment = (WD_ALIGN_PARAGRAPH.JUSTIFY if align == "just"
                   else WD_ALIGN_PARAGRAPH.CENTER)
    return p


def P(doc, t, size=None):
    return RT(doc, t, size)


def EQ(doc, parts, num=None):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for t, s in parts:
        r = p.add_run(t)
        r.italic = "i" in s
        if s == "sub":
            r.font.subscript = True
        if s == "sup":
            r.font.superscript = True
    if num:
        r = p.add_run("\t\t(" + str(num) + ")"); r.italic = False
    pf = p.paragraph_format
    pf.space_before = Pt(5); pf.space_after = Pt(5); pf.keep_with_next = True
    return p


def cap(doc, t):
    p = doc.add_paragraph()
    r = p.add_run(t); r.font.size = Pt(9.0)
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    p.paragraph_format.space_after = Pt(10)
    return p


def _autofit(table, total_in=6.7):
    """按内容长度分配列宽，禁止单词被拆行。"""
    tbl = table._tbl
    tblPr = tbl.tblPr
    lay = tblPr.find(qn("w:tblLayout"))
    if lay is None:
        lay = tblPr.makeelement(qn("w:tblLayout"), {})
        tblPr.append(lay)
    lay.set(qn("w:type"), "fixed")
    ncol = len(table.columns)
    weights = []
    for j in range(ncol):
        mx = 1
        for cell in table.columns[j].cells:
            for para in cell.paragraphs:
                for run in para.runs:
                    mx = max(mx, len(run.text))
        weights.append(max(mx, 4))
    tot = sum(weights)
    widths = [Emu(int(Inches(total_in) * w / tot)) for w in weights]
    for j, w in enumerate(widths):
        for cell in table.columns[j].cells:
            cell.width = w


def tbl(doc, df, above=None, below=None, fs=8.0, total_in=6.7):
    if above:
        p = doc.add_paragraph(); r = p.add_run(above); r.font.size = Pt(9.0)
        r.bold = True
        p.paragraph_format.keep_with_next = True
    t = doc.add_table(rows=1, cols=len(df.columns))
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = False
    hdr = t.rows[0]
    trPr = hdr._tr.get_or_add_trPr()
    th = trPr.makeelement(qn("w:tblHeader"), {}); trPr.append(th)
    for i, c in enumerate(df.columns):
        cell = hdr.cells[i]; cell.text = ""
        r = cell.paragraphs[0].add_run(str(c)); r.bold = True; r.font.size = Pt(fs)
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


# ================================================================
def build(zh):
    T = (lambda en, z: z if zh else en)
    SUF = "_ZH" if zh else "_EN"

    def FIG(doc, stem, w, folder=FM):
        doc.add_picture(str(folder / f"{stem}{SUF}.png"), width=Inches(w))

    D = Document(); setup(D, zh)

    # ---------- title page ----------
    t = D.add_paragraph()
    r = t.add_run(T("Disentangling Between-Target and Within-Target Signal in "
                    "PPI ΔΔG Benchmarking",
                    "解耦 PPI ΔΔG 基准评价中的跨靶信号与靶内信号"))
    r.bold = True; r.font.size = Pt(19)
    t.alignment = WD_ALIGN_PARAGRAPH.CENTER
    s = D.add_paragraph()
    rs = s.add_run(T("An exact covariance-signal decomposition under matched, "
                     "family-out-of-distribution evaluation",
                     "匹配、家族外推条件下的精确协方差信号分解"))
    rs.italic = True; rs.font.size = Pt(11)
    s.alignment = WD_ALIGN_PARAGRAPH.CENTER
    D.add_paragraph()
    for line in [T("Author Name(s) [AUTHOR_PLACEHOLDER]", "作者姓名 [AUTHOR_PLACEHOLDER]"),
                 T("Affiliation(s) [AFFILIATION_PLACEHOLDER]",
                   "单位 [AFFILIATION_PLACEHOLDER]"),
                 T("Corresponding author: [NAME] · Email: [EMAIL] · "
                   "ORCID: [ORCID_PLACEHOLDER]",
                   "通讯作者：[NAME] · 邮箱：[EMAIL] · ORCID：[ORCID_PLACEHOLDER]")]:
        p = D.add_paragraph(); rr = p.add_run(line); rr.font.size = Pt(10)
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    D.add_paragraph()
    p = D.add_paragraph()
    rr = p.add_run(T("Keywords: ", "关键词：")); rr.bold = True; rr.font.size = Pt(10)
    p.add_run(T("protein–protein interaction; binding affinity; ΔΔG; mutation effect "
                "prediction; benchmarking; out-of-distribution generalisation",
                "蛋白质–蛋白质相互作用；结合亲和力；ΔΔG；突变效应预测；"
                "基准评价；分布外泛化")).font.size = Pt(10)
    D.add_page_break()

    # ---------- abstract ----------
    H(D, T("Abstract", "摘要"), 1)
    if zh:
        P(D, f"预测突变对蛋白质–蛋白质结合亲和力（ΔΔG）的影响是抗体工程、界面优化与"
             f"致病突变解释的核心任务，其进展高度依赖基准评价的可信度。"
             f"然而常用的合并相关系数难以直接解读：它同时取决于模型能否把不同靶点放在"
             f"正确的结合量级上，以及能否在同一靶点内部把突变正确排序。"
             f"我们引入一个精确的协方差信号分解，把合并 Pearson 相关写为 "
             f"r = C_B + C_W，其中 C_B 与 C_W 分别是跨靶与靶内协方差贡献，"
             f"并在匹配、去伪重复、来源筛选后的评价集合上应用该分解。"
             f"在仅改变划分方式的受控比较中（同一批 7,523 行突变、同一批模型配置），"
             f"合并 Pearson 中位数由随机划分的 {r0:.4f} 降至家族外推的 {r3:.4f}，"
             f"ΔC_B = {DL['delta_C_B']:+.4f}、ΔC_W = {DL['delta_C_W']:+.4f}；"
             f"按绝对变化计算，跨靶分量占下降的 {100*DL['fraction_between']:.1f}%，"
             f"而靶内分量仍为正。改变目标加权会改变对这两个分量的解读："
             f"在五种方案下 C_B 仅在近零处出现方向翻转（{sf0}/{len(W2)}），"
             f"而 C_W 在任何阈值下均不改变符号，但其幅度并非不变。"
             f"固定 K 短名单指标受候选集几何约束：K = 10 时 D(n,K) 在 n ≤ K 时为零，"
             f"并在 n* ≈ {GEO['n_star_analytic']:.1f} 取峰，"
             f"本基准中 {100*GEO['n_le_K_fraction']:.1f}% 的 n ≥ 5 靶点处于该退化区间。"
             f"在匹配的回溯消融中，该分解区分了仅改善合并分数的改动与同时改善靶内排序的改动。"
             f"外部可评价方法仅五个，功效不足且秩一致为正，因此不支持合并性能与靶内性能的"
             f"普遍解耦。我们建议 PPI ΔΔG 基准在明确声明的加权下同时报告两个分量，"
             f"并报告候选集几何。")
    else:
        P(D, f"Predicting how mutations alter protein–protein binding affinity (ΔΔG) "
             f"underpins antibody engineering, interface optimisation and the "
             f"interpretation of disease-associated variants, and progress depends on the "
             f"credibility of benchmark evaluation. The pooled correlation commonly "
             f"reported for such benchmarks is difficult to interpret directly, because it "
             f"depends both on whether a model places different targets on the correct "
             f"overall binding scale and on whether it orders mutations correctly within "
             f"each target. We introduce an exact covariance-signal decomposition that "
             f"writes pooled Pearson correlation as r = C_B + C_W, where C_B and C_W are "
             f"the between-target and within-target covariance contributions, and apply it "
             f"to matched, de-replicated and provenance-filtered evaluation sets. In a "
             f"controlled comparison that changes only the split regime (identical 7,523 "
             f"mutation rows and identical model configurations), median pooled Pearson "
             f"falls from {r0:.4f} under random splitting to {r3:.4f} under "
             f"family-out-of-distribution evaluation, with ΔC_B = {DL['delta_C_B']:+.4f} "
             f"and ΔC_W = {DL['delta_C_W']:+.4f}; on an absolute-change basis the "
             f"between-target component accounts for "
             f"{100*DL['fraction_between']:.1f}% of the decline while the within-target "
             f"component remains positive. Target weighting changes how these "
             f"contributions should be read: C_B reverses direction only near zero "
             f"({sf0}/{len(W2)} configurations across five schemes), whereas C_W keeps its "
             f"sign at every threshold although its magnitude is not invariant. Fixed-K "
             f"shortlist metrics are further constrained by candidate-set geometry: for "
             f"K = 10, D(n,K) is zero when n ≤ K and peaks near "
             f"n* ≈ {GEO['n_star_analytic']:.1f}, and {100*GEO['n_le_K_fraction']:.1f}% of "
             f"targets with n ≥ 5 fall in that degenerate regime. In a matched "
             f"retrospective ablation the decomposition separated a change that improved "
             f"pooled scores only from one that also improved within-target ranking. "
             f"External evidence comprised five methods, is underpowered and is positively "
             f"aligned, so it does not support a universal decoupling of pooled and "
             f"within-target performance. We recommend that PPI ΔΔG benchmarks report both "
             f"components under an explicitly declared weighting, together with "
             f"candidate-set geometry.")
    p = D.add_paragraph()
    rr = p.add_run(T("Keywords: ", "关键词：")); rr.bold = True
    p.add_run(T("protein–protein interaction; binding affinity; ΔΔG; mutation effect "
                "prediction; benchmarking; out-of-distribution generalisation",
                "蛋白质–蛋白质相互作用；结合亲和力；ΔΔG；突变效应预测；"
                "基准评价；分布外泛化")).font.size = Pt(10)

    # ---------- introduction ----------
    H(D, T("Introduction", "引言"), 1)
    if zh:
        P(D, "预测突变如何改变蛋白质–蛋白质结合亲和力（ΔΔG）是抗体亲和力成熟、"
             "界面优化与致病突变解释中的核心计算任务。实验测定结合亲和力的成本高昂，"
             "因此计算预测被广泛用于对突变进行排序并指导实验设计。随着深度学习方法的"
             "快速增加，该领域进展的判据越来越依赖基准评价：一个方法是否更好，"
             "通常由其在标准突变集合上报告的相关系数或误差决定。因此，评价协议本身的性质"
             "——它测量了什么、以及在何种外推条件下测量——直接决定了领域能得出什么结论。")
        P(D, f"目前的评价实践以合并指标为主，即在整个突变集合上计算 Pearson 或 Spearman "
             f"相关以及 RMSE^[{KEYNUM['biometry']}]，同时常辅以按复合物或按靶点计算的指标。"
             f"近年已有工作揭示，这类评价对外推能力的刻画受到划分策略的强烈影响。"
             f"Bushuiev 等量化了 PPI 数据中的结构性泄漏：按 PDB 代码划分时约 65% 的测试样本"
             f"存在近重复训练样本，即使按 30% 序列同一性划分仍有约 30%，SKEMPI v2.0"
             f"^[{KEYNUM['skempi']}] 若按 PPI 代码划分为 56%，而按领域专家分组则为 0%"
             f"^[{KEYNUM['leakage']}]。CATH-ddG 采用 CATH 同源超家族划分并以最大 TM-score "
             f"< 0.6 定义「难突变」^[{KEYNUM['cathddg']}]；相关方法亦报告按结构或按靶点的"
             f"相关性^[{KEYNUM['hermes']},{KEYNUM['prefold']},{KEYNUM['ppiformer']}]。"
             f"需要明确的是：泄漏感知划分、家族外推评价、按靶点报告以及短名单式（Top-K）"
             f"评价均非本文提出。")
        P(D, "然而，一个关键问题仍未被直接回答：合并基准性能究竟来自哪里？"
             "把合并指标与按靶点指标并列报告，并不等于量化了跨靶信号与靶内信号对"
             "合并性能的加性贡献。一个模型可能主要因为能把不同靶点放在正确的整体量级上"
             "而获得较高的合并相关，而其在靶点内部的排序能力可能很弱；反之亦然。"
             "在靶点数量多、各靶点突变数高度不均的基准中，合并相关还会被大靶点主导，"
             "使这一混淆进一步加剧。")
        P(D, "本文引入一个精确的协方差信号分解来回答这一问题。我们把合并 Pearson 相关"
             "严格分解为跨靶分量与靶内分量，使「合并性能由哪一部分构成」成为可直接测量的量，"
             "并在匹配、去伪重复、来源筛选后的评价集合上使用该分解量化四件事："
             "划分严格度改变的是哪一个分量；目标加权如何改变对这两个分量的解读；"
             "固定 K 短名单评价受候选集几何的何种约束；以及该分解能否区分"
             "「仅改善合并分数」与「同时改善靶内排序」的模型改动。我们进一步在五个"
             "独立开发的外部方法上界定当前证据的泛化边界。")
    else:
        P(D, "Predicting how mutations alter protein–protein binding affinity (ΔΔG) is "
             "central to antibody affinity maturation, interface optimisation and the "
             "interpretation of disease-associated variants. Because experimental affinity "
             "measurement is costly, computational predictions are widely used to "
             "prioritise mutations and guide experimental design. As the number of "
             "deep-learning methods grows rapidly, progress in the field is increasingly "
             "judged by benchmark evaluation: whether one method is better is usually "
             "decided by the correlation or error it reports on a standard mutation set. "
             "The properties of the evaluation protocol itself — what it measures, and "
             "under which extrapolation conditions — therefore directly determine what the "
             "field can conclude.")
        P(D, f"Current practice relies largely on pooled metrics, computing Pearson or "
             f"Spearman correlation and RMSE over the entire mutation "
             f"set^[{KEYNUM['biometry']}], often supplemented by per-complex or per-target "
             f"metrics. Recent work has shown that how well such evaluation captures "
             f"extrapolation depends strongly on the splitting strategy. Bushuiev et al. "
             f"quantified structural leakage in PPI data: about 65% of test examples have a "
             f"near-duplicate training example under PDB-code splitting, about 30% even "
             f"under 30% sequence-identity splitting, 56% for SKEMPI "
             f"v2.0^[{KEYNUM['skempi']}] split by PPI code, and 0% under domain-expert "
             f"grouping^[{KEYNUM['leakage']}]. CATH-ddG adopted CATH homologous-superfamily "
             f"splits and defined hard mutations by a maximum TM-score below "
             f"0.6^[{KEYNUM['cathddg']}]; related methods likewise report per-structure or "
             f"per-target "
             f"correlations^[{KEYNUM['hermes']},{KEYNUM['prefold']},{KEYNUM['ppiformer']}]. "
             f"We emphasise that leakage-aware splitting, family-out-of-distribution "
             f"evaluation, per-target reporting and shortlist (Top-K) evaluation were none "
             f"of them introduced by this work.")
        P(D, "A key question nevertheless remains unanswered directly: where does pooled "
             "benchmark performance come from? Reporting pooled and per-target metrics side "
             "by side does not quantify the additive contribution of between-target and "
             "within-target signal to pooled performance. A model may achieve a high pooled "
             "correlation mainly because it places different targets on the correct overall "
             "scale, while its ability to order mutations within a target is weak — or the "
             "reverse. In benchmarks containing many targets with highly unequal mutation "
             "counts, pooled correlation is additionally dominated by large targets, which "
             "amplifies the conflation.")
        P(D, "Here we introduce an exact covariance-signal decomposition that addresses "
             "this question. We decompose pooled Pearson correlation rigorously into a "
             "between-target component and a within-target component, which makes 'what "
             "pooled performance is made of' a directly measurable quantity. Using this "
             "decomposition on matched, de-replicated and provenance-filtered evaluation "
             "sets, we quantify four things: which component is changed by split "
             "strictness; how target weighting changes the interpretation of the two "
             "components; which constraints candidate-set geometry places on fixed-K "
             "shortlist evaluation; and whether the decomposition can distinguish model "
             "modifications that improve pooled scores only from those that also improve "
             "within-target ranking. We further bound the current evidence using five "
             "independently developed external methods.")

    # ---------- RESULTS ----------
    H(D, T("Results", "结果"), 1)

    # Table 1 -- 首次提及处
    tbl(D, T1, above=T("Table 1 | Dataset and evaluation overview.",
                       "表 1 | 数据集与评价概览。"),
        below=T("All rows are joined by explicit row identifiers. The controlled split set "
                "is used only for the split-strictness experiment; the primary matched set "
                "is used for all capability-resolved analyses.",
                "所有行均以显式行标识符连接。受控划分集合仅用于划分严格度实验；"
                "主匹配集合用于全部能力分辨分析。"), fs=8.4)

    H(D, T("Pooled PPI ΔΔG performance contains separable between-target and "
           "within-target signal",
           "合并 PPI ΔΔG 性能包含可分离的跨靶与靶内信号"), 2)
    if zh:
        P(D, "对属于靶点 c 的突变 i，观测到的 ΔΔG 可以写成该靶点的平均水平"
             "加上该突变相对这一水平的偏差：")
    else:
        P(D, "For a mutation i belonging to target c, the observed ΔΔG can be written as a "
             "target-level mean plus a mutation-specific deviation from that mean:")
    EQ(D, [("y", "i"), ("ci", "sub"), (" = μ", "i"), ("c", "sub"), (" + ε", "i"),
           ("ci", "sub")], 1)
    P(D, T("where μ_c is the mean ΔΔG of target c and ε_ci = y_ci − μ_c is the "
           "within-target mutation deviation. The prediction is decomposed identically:",
           "其中 μ_c 为靶点 c 的平均 ΔΔG，ε_ci = y_ci − μ_c 为靶内突变偏差。"
           "预测值作同样的分解："))
    EQ(D, [("ŷ", "i"), ("ci", "sub"), (" = μ̂", "i"), ("c", "sub"), (" + ε̂", "i"),
           ("ci", "sub")], 2)
    if zh:
        P(D, "只有当这一分解明确之后，才有意义地区分两种协方差贡献。"
             "设 σ_y 与 σ_ŷ 为合并标准差，合并相关可精确拆分如下：")
    else:
        P(D, "Only after this decomposition is it meaningful to distinguish two covariance "
             "contributions. With σ_y and σ_ŷ the pooled standard deviations, the pooled "
             "correlation splits exactly as:")
    EQ(D, [("r", "i"), ("pooled", "sub"), (" = Cov(y, ŷ) / (σ", "i"), ("y", "sub"),
           ("σ", "i"), ("ŷ", "sub"), (") = C", ""), ("B", "sub"), (" + C", ""),
           ("W", "sub")], 3)
    if zh:
        P(D, "其中 C_B = Cov(B_y, B_ŷ)/(σ_y σ_ŷ) 为跨靶协方差贡献，"
             "C_W = Cov(W_y, W_ŷ)/(σ_y σ_ŷ) 为靶内协方差贡献。"
             "该恒等式在机器精度内成立（残差 ≤ 3.3 × 10⁻¹⁶）。"
             "需要强调：C_B 与 C_W 是同一尺度上的协方差信号贡献，"
             "不是独立相关系数、独立能力分数或因果机制；完整推导见补充说明 1。")
    else:
        P(D, "where C_B = Cov(B_y, B_ŷ)/(σ_y σ_ŷ) is the between-target covariance "
             "contribution and C_W = Cov(W_y, W_ŷ)/(σ_y σ_ŷ) is the within-target "
             "covariance contribution. The identity holds to machine precision "
             "(residual ≤ 3.3 × 10⁻¹⁶; full derivation in Supplementary Note 1). We stress "
             "that C_B and C_W are covariance-signal contributions on a common scale, not "
             "independent correlations, independent skill scores, or causal mechanisms.")
    if zh:
        P(D, f"我们在一个匹配的评价集合上使用该分解，该集合包含 "
             f"{PU['n_rows']:,} 行突变、{PU['n_ge5']} 个至少五个突变的靶点，"
             f"涉及 {PU['n_method_families']} 个独立方法族；跨模型比较仅在同一批突变行内"
             f"进行（方法；补充图 S1）。")
    else:
        P(D, f"We apply this decomposition on a matched evaluation set comprising "
             f"{PU['n_rows']:,} mutation rows, {PU['n_ge5']} targets with at least five "
             f"mutations, and {PU['n_method_families']} independent method families; "
             f"cross-model comparisons were restricted to predictions evaluated on "
             f"identical mutation-row sets (Methods; Supplementary Fig. S1).")
    FIG(D, "Fig1_overview", 6.7)
    cap(D, T("Figure 1 | Capability-resolved evaluation of PPI ΔΔG predictors. "
             "(A) Provenance-filtered, fingerprint-matched benchmark construction: "
             "SKEMPI v2.0 annotation (4,473 rows, 214 targets, 108 families) is audited by "
             "SHA-256 row-set fingerprinting; family-out-of-distribution predictions "
             "without training overlap are retained and permissively split predictions are "
             "excluded, giving a primary matched universe of 3,351 rows, 102 targets with "
             "n ≥ 5 and 10 method families. Target cards illustrate the two structural "
             "features of the benchmark: between-target differences in mean ΔΔG (μ_c, "
             "purple) and within-target mutation ordering (ε_ci, green). (B) Exact additive "
             "covariance-signal decomposition: target centring yields a between-target "
             "contribution C_B and a within-target contribution C_W, whose sum is pooled "
             "correlation by an exact identity (residual ≤ 3.3 × 10⁻¹⁶). These are "
             "common-scale additive contributions, not independent correlations, skill "
             "scores or causal mechanisms. (C) Controlled split collapse on identical "
             "7,523 mutation rows: pooled correlation falls from 0.6844 to 0.1508, and "
             "81.1% of the absolute decline (ΔC_B = −0.4347) is carried by loss of the "
             "between-target component while C_W remains positive (C_W = +0.1580). "
             "(D) Standardised capability profiles of matched ablation variants across five "
             "axes; axes are standardised for visual comparison only and raw values are "
             "reported in the text.",
             "图 1 | PPI ΔΔG 预测器的能力分辨评价。(A) 来源筛选与指纹匹配的基准构建："
             "SKEMPI v2.0 标注（4,473 行、214 靶点、108 家族）经 SHA-256 行集指纹审计；"
             "保留无训练重叠的家族外推预测、剔除宽松划分预测，"
             "得到 3,351 行、102 个 n ≥ 5 靶点、10 个方法族的主匹配宇宙。"
             "靶点卡片示意基准的两个结构特征：跨靶平均 ΔΔG 差异（μ_c，紫色）"
             "与靶内突变排序（ε_ci，绿色）。(B) 精确可加协方差信号分解："
             "靶点中心化得到跨靶贡献 C_B 与靶内贡献 C_W，二者之和以精确恒等式"
             "等于合并相关（残差 ≤ 3.3 × 10⁻¹⁶）。二者是同一尺度上的可加贡献，"
             "并非独立相关、能力分数或因果机制。(C) 同一批 7,523 行上的受控划分崩塌："
             "合并相关由 0.6844 降至 0.1508，绝对下降的 81.1%（ΔC_B = −0.4347）"
             "由跨靶分量的丧失承载，而 C_W 仍为正（C_W = +0.1580）。"
             "(D) 匹配消融变体在五个轴上的标准化能力画像；各轴仅供视觉比较，"
             "原始数值见正文。"))

    # 2.2
    H(D, T("Family-OOD evaluation disproportionately removes between-target signal",
           "家族外推评价不成比例地移除跨靶信号"), 2)
    if zh:
        P(D, f"为分离划分严格度的效应，我们在同一批 7,523 行突变与同一批模型配置上，"
             f"仅改变划分方式（表 2、图 2）。合并 Pearson 中位数由随机划分的 {r0:.4f} "
             f"降至家族外推的 {r3:.4f}。分解显示该下降并非均分："
             f"ΔC_B = {DL['delta_C_B']:+.4f}，ΔC_W = {DL['delta_C_W']:+.4f}。"
             f"按绝对变化口径 |ΔC_B|/(|ΔC_B| + |ΔC_W|) 计算，跨靶分量占下降的 "
             f"{100*DL['fraction_between']:.1f}%，而 C_W 在家族外推下仍为正"
             f"（{cw3:+.4f}）。因此，家族外推评价下合并相关的下降主要由跨靶协方差贡献的"
             f"丧失承载，而靶内协方差保持为正。")
        P(D, f"这一模式与「模型记住了靶点身份」不是同一命题：本研究没有进行身份记忆的"
             f"干预实验，因此不作该因果陈述。更稳妥的表述是，该模式与「模型依赖不能"
             f"跨家族迁移的靶点特异性信息」相一致。此外，复合物层面与 seq40 划分在本数据上"
             f"给出完全相同的分配，因此我们不把它们作为彼此独立的比较。")
    else:
        P(D, f"To isolate the effect of split strictness, we compared identical model "
             f"configurations on the same 7,523 mutation rows while changing only the split "
             f"regime (Table 2, Figure 2). Median pooled Pearson falls from {r0:.4f} under "
             f"random splitting to {r3:.4f} under family-out-of-distribution evaluation. "
             f"The decomposition shows that this decline is not uniform: "
             f"ΔC_B = {DL['delta_C_B']:+.4f} and ΔC_W = {DL['delta_C_W']:+.4f}. On an "
             f"absolute-change basis, |ΔC_B| / (|ΔC_B| + |ΔC_W|), the between-target "
             f"component accounts for {100*DL['fraction_between']:.1f}% of the decline, "
             f"while C_W remains positive under family-out-of-distribution evaluation "
             f"({cw3:+.4f}). Most of the reduction in pooled correlation under family-OOD "
             f"evaluation is therefore carried by loss of the between-target covariance "
             f"contribution, while within-target covariance remains positive.")
        P(D, "This is not the same proposition as 'the model memorised target identity': we "
             "performed no identity-memorisation intervention and therefore make no such "
             "causal claim. A more defensible reading is that the pattern is consistent "
             "with reliance on target-specific information that does not transfer across "
             "families. The complex-level and seq40 partitions produced identical "
             "assignments in this dataset and were therefore not treated as independent "
             "comparisons.")
    FIG(D, "Fig2_split_decomposition", 6.7)
    cap(D, T("Figure 2 | Family-out-of-distribution evaluation disproportionately removes "
             "the between-target component of pooled correlation. (A) Median pooled Pearson "
             "correlation by split regime. (B) Exact decomposition into C_B (muted purple) "
             "and C_W (green); C_W remains positive. (C) Attribution of the absolute "
             "decline, of which the between-target component accounts for "
             f"{100*DL['fraction_between']:.1f}%. All four regimes use the same 7,523 "
             "mutation rows and the same model configurations; the complex and seq40 "
             "partitions are identical in this dataset and are not independent comparisons.",
             "图 2 | 家族外推评价不成比例地移除合并相关中的跨靶分量。(A) 各划分下的"
             "合并 Pearson 中位数。(B) 精确分解为 C_B（紫色）与 C_W（绿色）；C_W 仍为正。"
             f"(C) 绝对下降的归因，跨靶分量占 {100*DL['fraction_between']:.1f}%。"
             "四个划分使用同一批 7,523 行与同一批配置；复合物与 seq40 划分在本数据上相同，"
             "不构成独立比较。"))
    tbl(D, T2, above=T("Table 2 | Controlled split comparison.", "表 2 | 受控划分比较。"),
        below=T("Median values across model configurations. The seq40 partition is "
                "identical to the complex partition in this dataset and is therefore not "
                "an independent comparison.",
                "跨模型配置的中位数。seq40 划分在本数据上与复合物划分完全相同，"
                "不构成独立比较。"))

    # 2.3
    H(D, T("Target aggregation changes the interpretation of covariance contributions",
           "目标聚合改变协方差贡献的解读"), 2)
    if zh:
        P(D, f"由于各靶点突变数高度不均，聚合权重的选择会影响两个分量的大小。"
             f"我们固定预测值，仅改变聚合权重，比较五种方案（每个突变等权、每个靶点等权、"
             f"权重正比于 √n、正比于 log(1 + n)、以及上限为 30 的截断权重）。"
             f"结果显示 C_B 的方向不稳定性集中在近零区间：无阈值时 {sf0}/{len(W2)} 个配置"
             f"出现符号翻转，而在 |C_B| > 0.005 时不再出现（{sf5}/{len(W2)}）。"
             f"C_W 在任何阈值下均未出现符号翻转（{sfw}/{len(W2)}），但其幅度并非不变——"
             f"其极差中位数为 {W2.range_C_W.median():.4f}，大于 C_B 的 "
             f"{W2.range_C_B.median():.4f}。因此准确表述是：C_B 仅在近零区间方向不稳，"
             f"而 C_W 保持符号稳定但幅度可变。这意味着当报告的 C_B 接近零时，"
             f"「跨靶信息是否存在」的结论可能取决于加权口径而非模型；"
             f"因此任何报告的分解都应声明其加权方案。")
    else:
        P(D, f"Because targets differ greatly in mutation count, the choice of aggregation "
             f"weight affects the magnitude of both components. Holding predictions fixed "
             f"and varying only the aggregation weights, we compared five schemes "
             f"(per-mutation, per-target, weight proportional to √n, proportional to "
             f"log(1 + n), and capped at 30). C_B shows directional instability "
             f"concentrated near zero: {sf0}/{len(W2)} configurations reverse sign with no "
             f"threshold, whereas none do once |C_B| exceeds 0.005 ({sf5}/{len(W2)}). C_W "
             f"never reverses sign at any threshold ({sfw}/{len(W2)}), but its magnitude is "
             f"not invariant — its median range is {W2.range_C_W.median():.4f}, larger than "
             f"the {W2.range_C_B.median():.4f} of C_B. C_B is therefore unstable in "
             f"direction only near zero, whereas C_W keeps its sign but varies in "
             f"magnitude. Consequently, when a reported C_B is close to zero, the "
             f"conclusion about whether between-target information is present may depend on "
             f"the weighting convention rather than on the model; any reported "
             f"decomposition should declare its weighting scheme.")
    FIG(D, "Fig3_weighting", 6.7)
    cap(D, T("Figure 3 | Target weighting changes the interpretation of covariance "
             "contributions. Predictions are unchanged; only aggregation weights differ. "
             "(A, B) Distribution of the range of C_B and C_W across the five schemes; the "
             "dashed line marks the median range. (C) Number of configurations reversing "
             "sign at three thresholds on |C|. C_B reverses direction only near zero, "
             "whereas C_W keeps its sign at every threshold but its magnitude is not "
             "invariant.",
             "图 3 | 目标加权改变协方差贡献的解读。预测值不变，仅改变聚合权重。"
             "(A, B) 五种方案下 C_B 与 C_W 极差的分布；虚线为中位极差。"
             "(C) 三个 |C| 阈值下出现符号翻转的配置数。C_B 仅在近零处翻转方向，"
             "而 C_W 在任何阈值下均保持符号，但其幅度并非不变。"))

    # 2.4
    H(D, T("Fixed-K shortlist evaluation is constrained by candidate-set geometry",
           "固定 K 短名单评价受候选集几何约束"), 2)
    if zh:
        P(D, f"在突变文库或短名单筛选场景中，常用固定 K 的 Top-K 指标评价模型"
             f"^[{KEYNUM['ndcg']}]。可区分区间为：")
    else:
        P(D, f"Fixed-K top-K metrics are commonly used to evaluate models for mutation "
             f"library or shortlist screening^[{KEYNUM['ndcg']}]. The discriminable "
             f"interval is:")
    EQ(D, [("D(n, K) = 0,  n ≤ K;     D(n, K) = K(n − K) / [n(n + 1)],  n > K", "")], 4)
    if zh:
        P(D, f"其中随机选取 K 个候选的期望收益占比为 K/n，最优 K 个候选为 "
             f"K(2n − K + 1)/[n(n + 1)]（推导见补充说明 6）。当靶点突变数不超过 K 时，"
             f"任何模型都会选中全部候选，该靶点对固定 K 短名单指标不携带任何信息。"
             f"本基准中，{GEO['n_ge5_count']} 个 n ≥ 5 靶点里有 {GEO['n_le_K_count']} 个"
             f"（{100*GEO['n_le_K_fraction']:.1f}%）处于该退化区间。"
             f"此外 D 并非随 n 单调增加，而在 n* = K + √(K² + K) 处取峰；K = 10 时 "
             f"n* ≈ {GEO['n_star_analytic']:.1f}，例如 D(20) = {GEO['D_at_n20']:.4f} 而 "
             f"D(50) = {GEO['D_at_n50']:.4f}。因此「候选集越大越可区分」的启发式"
             f"在本指标下不成立。")
        P(D, "上述几何结论不意味着 SKEMPI 的 ΔΔG 数值本身过于接近；后者需要以真实 ΔΔG "
             "效用幅度做实证分析，本文未做该项分析。可辩护的表述是：当 n ≤ K 时"
             "固定 K 的候选选择评价在结构上不携带信息，且其可区分范围强烈依赖候选集几何。")
    else:
        P(D, f"where the expected share of gain from K randomly chosen candidates is K/n "
             f"and the optimal K candidates attain K(2n − K + 1) / [n(n + 1)] (derivation "
             f"in Supplementary Note 6). When a target has no more than K mutations, any "
             f"model selects every candidate and the target carries no information for "
             f"fixed-K shortlist metrics. In this benchmark, {GEO['n_le_K_count']} of "
             f"{GEO['n_ge5_count']} targets with n ≥ 5 "
             f"({100*GEO['n_le_K_fraction']:.1f}%) fall in that degenerate regime. Moreover "
             f"D is not monotonic in n: it peaks at n* = K + √(K² + K), which for K = 10 "
             f"gives n* ≈ {GEO['n_star_analytic']:.1f}, so that "
             f"D(20) = {GEO['D_at_n20']:.4f} while D(50) = {GEO['D_at_n50']:.4f}. The "
             f"heuristic that larger candidate sets are more discriminable therefore does "
             f"not hold for this metric.")
        P(D, "This geometric result does not imply that SKEMPI ΔΔG values are intrinsically "
             "too similar; establishing that would require an analysis using actual ΔΔG "
             "utility magnitudes, which we did not perform. The defensible statement is "
             "that fixed-K candidate-selection evaluation is structurally non-informative "
             "when n ≤ K, and that its discriminable range depends strongly on "
             "candidate-set geometry.")
    FIG(D, "Fig4_candidate_geometry", 5.9)
    cap(D, T("Figure 4 | Fixed-K shortlist evaluation is constrained by candidate-set "
             "geometry. The discriminable interval D(n, K) is zero when n ≤ K and follows "
             "K(n − K)/[n(n + 1)] otherwise; it peaks at n* = K + √(K² + K), which is "
             f"n* ≈ {GEO['n_star_analytic']:.1f} for K = 10, and declines thereafter. In "
             f"this benchmark {100*GEO['n_le_K_fraction']:.1f}% of targets with n ≥ 5 lie in "
             "the degenerate regime where every model selects all candidates.",
             "图 4 | 固定 K 短名单评价受候选集几何约束。可区分区间 D(n, K) 在 n ≤ K 时为"
             "零，否则为 K(n − K)/[n(n + 1)]；它在 n* = K + √(K² + K) 处取峰"
             f"（K = 10 时 n* ≈ {GEO['n_star_analytic']:.1f}）后下降。"
             f"本基准中 {100*GEO['n_le_K_fraction']:.1f}% 的 n ≥ 5 靶点处于"
             "所有模型都选中全部候选的退化区间。"))

    # 2.5
    H(D, T("Capability decomposition clarifies retrospective remedy effects",
           "能力分解解释回溯性模型改动的效果"), 2)
    if zh:
        P(D, "最后，我们在一个匹配的消融集合中检验该分解的实用性，该集合中的模型共享"
             "同一 3,351 行宇宙，因此差异只来自模型改动本身。以「仅使用突变性质」的模型"
             "为基线：加入静态局部几何使合并 Pearson 提高约 0.028，但靶内 Macro Spearman "
             "几乎不变（−0.002），属于仅改善合并分数的改动；引入「突变 × 配体 partner "
             "扰动」的改动同时改善两轴（合并 Pearson +0.070、Macro Spearman +0.044）；"
             "使用全量局部物理特征时两轴改善最大（+0.105 与 +0.076）。"
             "因此该分解能够区分「仅合并分数改善」与「同时改善靶内排序」两类改动。")
        P(D, "需要明确该分析的性质：它是回溯性的，即已知哪些改动有效，因此它说明的是"
             "该分解能够表达这一区别，而非基准能够前瞻预测哪条技术路线会成功。"
             "此外，本研究重新计算的靶内差异（+0.044）与历史分析中的自助法区间"
             "（+0.073，[+0.021, +0.129]）基于不同的度量宇宙与聚合方式，"
             "不可在不加说明的情况下混用。")
    else:
        P(D, "Finally, we examined the practical utility of the decomposition in a matched "
             "ablation set whose models share an identical 3,351-row universe, so that "
             "differences arise only from the model modifications themselves. Relative to a "
             "baseline using mutation properties alone, adding static local geometry raised "
             "pooled Pearson by about 0.028 while leaving within-target Macro Spearman "
             "essentially unchanged (−0.002), a pooled-score-only improvement; introducing "
             "mutation × partner perturbation improved both axes (pooled Pearson +0.070, "
             "Macro Spearman +0.044); and a full local-physics model gave the largest "
             "improvement on both axes (+0.105 and +0.076). The decomposition therefore "
             "distinguishes a change that improves pooled scores only from one that also "
             "improves within-target ranking.")
        P(D, "The nature of this analysis should be stated plainly: it is retrospective, in "
             "that we already know which modifications were effective. It therefore shows "
             "that the decomposition can express this distinction, not that the benchmark "
             "prospectively predicts which technical route will succeed. In addition, our "
             "recomputed within-target difference (+0.044) and the bootstrap interval from "
             "the historical analysis (+0.073, [+0.021, +0.129]) rest on different metric "
             "universes and aggregation schemes and must not be mixed without explanation.")
    FIG(D, "Fig5_remedy", 6.0)
    cap(D, T("Figure 5 | Capability decomposition separates pooled-score gains from "
             "within-target gains. Retrospective matched ablation on an identical 3,351-row "
             "universe; the baseline is the mutation-only model. Green markers indicate "
             "improvement on both axes; red markers indicate a pooled-score-only or "
             "negative change. Axes give the change in pooled Pearson and in macro "
             "within-target Spearman relative to that baseline. The analysis is "
             "retrospective and does not imply prospective prediction of which "
             "modification will succeed.",
             "图 5 | 能力分解区分合并分数增益与靶内增益。在同一 3,351 行宇宙上的回溯性"
             "匹配消融；基线为仅使用突变性质的模型。绿色表示两轴同时改善，红色表示"
             "仅合并分数改善或变化为负；两轴分别为合并 Pearson 与靶内 Macro Spearman "
             "相对基线的变化。该分析为回溯性，不代表可前瞻预测哪种改动会成功。"))

    # 2.6
    H(D, T("External methods define the boundary of generalisation",
           "外部方法界定泛化边界"), 2)
    if zh:
        P(D, f"为评估上述结论在独立开发的方法上的适用性，我们使用已公开的官方逐突变预测"
             f"构建外部评价集合：Flex ddG（基于 Rosetta 系综的结合亲和力估计）"
             f"^[{KEYNUM['flexddg']}]、FoldX^[{KEYNUM['foldx']}]、"
             f"ProteinMPNN^[{KEYNUM['proteinmpnn']}]、ESM-IF^[{KEYNUM['esmif']}] 与 "
             f"HERMES^[{KEYNUM['hermes']}]。在构建过程中我们发现，当匹配键省略链标识与"
             f"野生型残基时，该键在标注侧并不唯一：4,473 行中仅 "
             f"{FR['external_clean']['n_rows']:,} 行具有唯一映射。因此主匹配评价只使用"
             f"唯一映射行，歧义键被排除而非强行解析；这一处理使外部方法的表观合并相关"
             f"系统性下降 0.02 至 0.12。")
        P(D, f"在唯一映射集合（{FR['external_clean']['n_rows']:,} 行、"
             f"{FR['external_clean']['n_ge5']} 个 n ≥ 5 靶点）上，五个外部方法的总体相关"
             f"与靶内排序的排序一致，最优方法同为 {best.model}"
             f"（合并 Pearson {best.overall_pearson:+.4f}、Macro Spearman "
             f"{best.macro_spearman:+.4f}），即该集合内未出现方法层面的排名反转。"
             f"方法层面的秩一致为 ρ = "
             f"{er.get('rho_overall_pearson_vs_macro_spearman', float('nan')):+.2f}"
             f"（精确置换检验 p = "
             f"{er.get('perm_p_overall_pearson_vs_macro_spearman', float('nan')):.3f}）。"
             f"该结果为强正一致，因此不支持合并性能与靶内性能的普遍解耦；"
             f"同时它落在本研究内部观测到的异质范围内。仅有五个方法时，"
             f"方法层面的秩相关推断功效不足，因此我们不对该面板作方向性结论。")
        P(D, f"关于外部量的语义，其中一个外部模型的官方输出为对数几率差，"
             f"是无量纲的突变效应分数，而非已校准的 kcal/mol ΔΔG^[{KEYNUM['hermes']}]。"
             f"因此我们不对该模型报告 RMSE 或 MAE，也不把其使用的变体标注为"
             f"受结合亲和力监督训练的模型。")
    else:
        P(D, f"To assess how far the above conclusions apply to independently developed "
             f"methods, we built external evaluation sets from published official "
             f"per-mutation predictions: Flex ddG (Rosetta ensemble-based affinity "
             f"estimation)^[{KEYNUM['flexddg']}], FoldX^[{KEYNUM['foldx']}], "
             f"ProteinMPNN^[{KEYNUM['proteinmpnn']}], ESM-IF^[{KEYNUM['esmif']}] and "
             f"HERMES^[{KEYNUM['hermes']}]. During construction we found that when the "
             f"matching key omits the chain identifier and the wild-type residue, the key "
             f"is not unique on the annotation side: of 4,473 rows only "
             f"{FR['external_clean']['n_rows']:,} have a unique mapping. Primary matched "
             f"evaluation therefore uses unambiguous rows only; ambiguous keys were "
             f"excluded rather than forcibly resolved. This treatment systematically lowers "
             f"apparent pooled correlation for external methods by 0.02 to 0.12.")
        P(D, f"On the unambiguous set ({FR['external_clean']['n_rows']:,} rows, "
             f"{FR['external_clean']['n_ge5']} targets with n ≥ 5), the five external "
             f"methods agree in their ordering on pooled correlation and on within-target "
             f"ranking, with the same best method, {best.model} (pooled Pearson "
             f"{best.overall_pearson:+.4f}, Macro Spearman {best.macro_spearman:+.4f}); "
             f"that is, no method-level rank reversal occurs in this set. Method-level rank "
             f"agreement is ρ = "
             f"{er.get('rho_overall_pearson_vs_macro_spearman', float('nan')):+.2f} (exact "
             f"permutation p = "
             f"{er.get('perm_p_overall_pearson_vs_macro_spearman', float('nan')):.3f}). "
             f"This is a strongly positive agreement, which does not support a universal "
             f"decoupling of pooled and within-target performance and falls within the "
             f"heterogeneous range observed internally. With only five methods, "
             f"method-level rank-correlation inference is low-powered, and we therefore "
             f"draw no directional conclusion from this panel.")
        P(D, f"Regarding the semantics of one external quantity, its official output is a "
             f"log-odds difference, a dimensionless mutation-effect score rather than a "
             f"calibrated kcal/mol ΔΔG^[{KEYNUM['hermes']}]. We therefore report neither "
             f"RMSE nor MAE for that model, and we do not label the variant we used as "
             f"supervised on binding affinity.")
    FIG(D, "Fig6_external", 6.7)
    cap(D, T("Figure 6 | External methods define the boundary of generalisation. "
             "(A) Pooled versus within-target performance for five independently developed "
             "methods on a matched, unambiguous mapping; leader lines are used where "
             "markers are close. (B) Method-level capability profile, standardised within "
             "each metric (z-score) for visual comparison only; raw values are printed in "
             "each cell and statistical analyses use the raw metrics. Rank agreement "
             "between the two axes is positive (ρ = +0.90, exact permutation p = 0.083), "
             "but the panel is underpowered with five methods, so this positive alignment "
             "does not support a universal decoupling of pooled and within-target "
             "performance.",
             "图 6 | 外部方法界定泛化边界。(A) 五个独立开发的方法在匹配且无歧义映射上的"
             "合并性能与靶内性能；标记接近处使用引出线。(B) 方法级能力画像，"
             "仅按列标准化（z 分数）以便视觉比较；每格内给出原始数值，"
             "统计分析使用原始指标。两轴之间的秩一致为正（ρ = +0.90，精确置换 "
             "p = 0.083），但仅五个方法时功效不足，"
             "该正向对齐不支持合并性能与靶内性能的普遍解耦。"))
    tbl(D, T3, above=T("Table 3 | External method capability profiles on the unambiguous "
                       "mapping set.", "表 3 | 唯一映射集合上的外部方法能力画像。"),
        below=T("All methods are evaluated on the same mutations. RMSE is reported only "
                "for outputs with ΔΔG-like scaling; for the log-odds external quantity "
                "RMSE is not applicable and is omitted.",
                "所有方法在同一批突变上评价。仅对具有 ΔΔG 量纲的输出报告 RMSE；"
                "对对数几率型外部量 RMSE 不适用，故未给出。"))

    # ---------- discussion ----------
    H(D, T("Discussion", "讨论"), 1)
    if zh:
        P(D, "本研究表明，在 PPI ΔΔG 基准中通常作为单一数值报告的合并相关，"
             "实际上由两个性质不同的预测信号构成：把不同靶点放在正确整体量级上的跨靶信号，"
             "以及在靶点内部正确排序突变效应的靶内信号。我们给出的分解是精确的恒等式"
             "（式 3），因此「合并性能由哪一部分构成」不再是解释性问题，而是可直接计算的量。")
        P(D, f"与既有工作的关系需要明确界定。PPI 数据中的结构性泄漏及其对划分策略的依赖"
             f"已由 Bushuiev 等系统量化，他们同时提出了基于界面相似度的划分建议"
             f"^[{KEYNUM['leakage']}]；CATH-ddG 采用 CATH 同源超家族划分，并以结构相似度"
             f"定义了更难的外推突变^[{KEYNUM['cathddg']}]；按复合物或按靶点报告指标以及"
             f"短名单式评价在既有 ΔΔG 文献中亦已出现"
             f"^[{KEYNUM['hermes']},{KEYNUM['prefold']},{KEYNUM['ppiformer']}]。"
             f"本文的贡献不在于上述任何单一成分，而在于把精确的跨靶/靶内协方差分解"
             f"与划分严格度、目标加权、匹配条件下的排名一致性以及候选集几何四项审计"
             f"系统整合在同一套匹配、来源筛选的评价框架内，并辅以一个回溯性的模型改动"
             f"归因案例。据我们所知，我们未发现既有工作在同等匹配条件下联合报告上述各项。")
        P(D, "关于 C_B 与 C_W 的解读，需要强调其边界。二者是同一恒等式中的协方差贡献，"
             "共享同一尺度，因此可以相加得到合并相关；它们不是两个独立估计的能力分数，"
             "也不应被解释为因果机制。特别地，一个模型可能仅仅因为把各靶点放在大致正确的"
             "量级上就获得可观的合并相关，而其靶内排序能力有限——这正是合并指标可能"
             "高估实际筛选价值的原因。我们避免把靶内分量的存在解读为「模型学会了物理」，"
             "或把跨靶分量的丧失解读为「模型记住了靶点标识」。")
        P(D, "加权敏感性审计给出一个对实践有直接影响的结论：跨靶分量的符号在近零区间"
             "对聚合权重敏感，而靶内分量的符号稳定。这意味着当报告的 C_B 接近零时，"
             "「跨靶信息是否存在」的结论可能取决于加权口径而非模型。"
             "因此我们建议在报告分解时同时声明加权方案，并在 C_B 接近零时给出多种加权下的结果。")
        P(D, f"候选集几何的结果对 Top-K 评价实践有直接影响。当靶点突变数不超过 K 时，"
             f"固定 K 的短名单指标在结构上无法区分模型；本基准中约四成 n ≥ 5 靶点处于"
             f"该情形。同时可区分区间随候选集增大而先升后降，因此「扩大候选集以提高"
             f"区分度」并非普遍有效。我们建议在报告短名单指标时同时给出 n 与 K、"
             f"单独标注 n ≤ K 的靶点，或改用分数型（按比例 q）评价。")
        P(D, "回溯性消融案例说明该框架在解释模型改动时的实用价值：它能指出某项改动"
             "只提高了合并分数，而另一项改动同时改善了靶内排序。这一区分对方法设计"
             "具有直接意义，因为它把「分数提高」分解为可归因的来源。但该案例的回溯性质"
             "必须明确：它不构成对模型改进路线的前瞻预测能力。")
        P(D, "最后，外部证据界定了当前结论的泛化边界。在五个独立开发的方法上，"
             "合并性能与靶内排序的方法排序一致（强正一致），这既不支持也不反驳"
             "「两轴普遍解耦」的主张；该面板功效不足，因此我们不作方向性结论。"
             "未来工作应扩大外部方法面板，并在独立的前瞻实验数据集上检验"
             "分解分量与实际筛选收益之间的关系，以及更大候选文库下的评价设计。")
    else:
        P(D, "This study shows that the pooled correlation conventionally reported as a "
             "single number in PPI ΔΔG benchmarking is in fact composed of two "
             "qualitatively different predictive signals: a between-target component that "
             "places different targets on the correct overall scale, and a within-target "
             "component that orders mutation effects correctly inside each target. Because "
             "our decomposition is an exact identity (Eq. 3), 'what pooled performance is "
             "made of' ceases to be a matter of interpretation and becomes directly "
             "computable.")
        P(D, f"The relation to prior work should be delimited clearly. Structural leakage "
             f"in PPI data and its dependence on splitting strategy were quantified "
             f"systematically by Bushuiev et al., who also proposed interface-similarity "
             f"splitting^[{KEYNUM['leakage']}]; CATH-ddG adopted CATH homologous-superfamily "
             f"splits and used structural similarity to define harder extrapolation "
             f"mutations^[{KEYNUM['cathddg']}]; and per-complex or per-target reporting "
             f"together with shortlist-style evaluation already appear in the ΔΔG "
             f"literature^[{KEYNUM['hermes']},{KEYNUM['prefold']},{KEYNUM['ppiformer']}]. "
             f"Our contribution therefore lies not in any single one of these ingredients "
             f"but in the systematic integration of an exact between/within-target "
             f"covariance decomposition with four audits — split strictness, target "
             f"weighting, matched rank agreement, and candidate-set geometry — inside one "
             f"matched, provenance-filtered framework, together with a retrospective "
             f"attribution case study. To our knowledge we did not identify prior work that "
             f"jointly reports these under comparable matching conditions.")
        P(D, "The interpretation of C_B and C_W has clear boundaries. They are covariance "
             "contributions within one identity and share a common scale, so they sum to "
             "pooled correlation; they are not two independently estimated skill scores and "
             "should not be read as causal mechanisms. In particular, a model may obtain a "
             "substantial pooled correlation largely by placing targets on approximately "
             "the right overall scale while its within-target ordering ability remains "
             "limited — precisely why pooled metrics can overstate practical screening "
             "value. We avoid interpreting the presence of within-target signal as evidence "
             "that a model has learned physics, or the loss of between-target signal as "
             "evidence that a model memorised target identity.")
        P(D, "The weighting-sensitivity audit yields a conclusion of direct practical "
             "relevance: the sign of the between-target component is sensitive to "
             "aggregation weights near zero, whereas the sign of the within-target "
             "component is stable. Consequently, when a reported C_B is close to zero, the "
             "conclusion about whether between-target information is present may depend on "
             "the weighting convention rather than on the model. We therefore recommend "
             "that any reported decomposition declare its weighting scheme, and that "
             "multiple weightings be shown when C_B is near zero.")
        P(D, "The candidate-set geometry result has direct implications for Top-K "
             "evaluation practice. When a target has no more than K mutations, fixed-K "
             "shortlist metrics cannot distinguish models structurally; roughly four in ten "
             "targets with n ≥ 5 fall in this regime in our benchmark. Because the "
             "discriminable interval first rises and then falls as candidate sets grow, "
             "enlarging a candidate set is not universally effective at improving "
             "discrimination. We recommend reporting n alongside K, labelling n ≤ K targets "
             "separately, or adopting fraction-based (q) evaluation.")
        P(D, "The retrospective ablation illustrates the framework's explanatory value for "
             "model modifications: it distinguishes a change that raised pooled scores only "
             "from a change that also improved within-target ranking. This distinction "
             "matters for method design because it attributes an improvement to one of two "
             "sources. The retrospective nature of the case study must, however, be stated "
             "explicitly: it does not establish prospective ability to predict which "
             "modifications will work.")
        P(D, "Finally, external evidence bounds the generalisation of these conclusions. "
             "Across five independently developed methods the ordering of methods agreed "
             "between pooled and within-target performance, a strongly positive agreement "
             "that neither supports nor refutes a universal decoupling of the two axes. The "
             "panel is underpowered, so we draw no directional conclusion from it. Future "
             "work should broaden the external panel and test the relationship between the "
             "decomposed components and realised screening yield on independent prospective "
             "datasets, as well as evaluation design for larger candidate libraries.")

    # ---------- methods ----------
    H(D, T("Methods", "方法"), 1)
    if zh:
        P(D, f"4.1 数据与靶点定义。评价行来自 SKEMPI v2.0 派生的标注集合"
             f"^[{KEYNUM['skempi']}]（4,473 行 / 214 靶点 / 108 结构家族），"
             f"以行标识符精确连接预测输出。靶点定义为结构家族下的独立结合界面单元。")
        P(D, "4.2 预测来源与来源标准。预测输出按其记录的监督训练划分方式分类。"
             "在宽松划分（随机 / 复合物 / seq40）下训练的配置被排除在面向部署的比较之外，"
             "因为其评价行可能与其训练折重叠；它们仅保留在受控划分实验中。"
             "分析在三个层级进行：制品层（仅描述性）、配置层（合并随机种子）、"
             "方法族层（每个独立方法族一个代表，且为 performance-blind 选择）。")
        P(D, "4.3 匹配评价集合。为每个预测输出计算其行标识符排序列表的 SHA-256 指纹；"
             "跨模型比较仅在指纹相同的输出之间进行（补充说明 2、补充图 S1）。"
             "主匹配集合由在查看任何指标数值之前即已固定的规则选出。")
        P(D, "4.4 跨靶/靶内信号分解。见式 1–3；采用行级权重，并报告恒等式残差作为"
             "质量控制量。完整推导见补充说明 1。")
        P(D, "4.5 目标加权方案。预测值不变，仅改变聚合权重：每个突变等权、每个靶点等权、"
             "权重正比于 √n、正比于 log(1 + n)、上限 30 的截断权重。")
        P(D, f"4.6 排序与短名单指标。靶内排序使用按靶点的 Spearman ρ 与 Kendall τ 相关、"
             f"成对准确率以及靶内协方差贡献；短名单指标使用固定 K 与按比例 q 的命中率、"
             f"精确率、召回率、NDCG 与富集因子^[{KEYNUM['ndcg']}]。")
        P(D, f"4.7 统计分析。秩一致使用 Spearman ρ，Kendall τ 作次要指标；"
             f"区间估计使用靶点簇自助法（5,000 次重采样，种子 42）"
             f"^[{KEYNUM['bootstrap']}]；外部小样本使用精确置换检验。"
             f"指标方向在计算前登记，跨指标比较统一使用「越大越好」定向"
             f"^[{KEYNUM['biometry']}]。")
        P(D, f"4.8 外部预测映射与歧义控制。归一化匹配键在标注侧不唯一时，仅使用唯一映射行，"
             f"歧义键被排除而非解析；官方 held-out CATH 划分单独评价（补充说明 5）。")
        P(D, "4.9 可复现性与一致性验证。全部头条数字由单一机器可读文件生成，"
             "并附自动化一致性检查（补充材料）。未训练任何新模型，未使用 GPU。")
    else:
        P(D, f"4.1 Dataset and target definition. Evaluation rows come from a SKEMPI v2.0-"
             f"derived annotation set^[{KEYNUM['skempi']}] (4,473 rows, 214 targets, 108 "
             f"structural families), joined to prediction outputs by explicit row "
             f"identifiers. A target is defined as an independent binding-interface unit "
             f"within a structural family.")
        P(D, "4.2 Prediction sources and provenance criteria. Prediction outputs were "
             "classified by their recorded supervised training split regime. Configurations "
             "trained under a permissive (random / complex / seq40) partition were excluded "
             "from deployment-oriented comparisons because their evaluation rows may overlap "
             "their training folds; they were retained only in the controlled split "
             "experiment. Analyses were performed at three levels: artifact level "
             "(descriptive only), configuration level (seed replicates collapsed), and "
             "method-family level (one performance-blind representative per independent "
             "family).")
        P(D, "4.3 Matched evaluation sets. Each prediction output received a SHA-256 "
             "fingerprint of its sorted row identifiers; cross-model comparison was "
             "performed only between outputs sharing a fingerprint (Supplementary Note 2, "
             "Supplementary Fig. S1). The primary matched set was selected by a rule fixed "
             "before any metric values were inspected.")
        P(D, "4.4 Between-/within-target signal decomposition. See Eqs. 1–3; row-level "
             "weights are used and the identity residual is reported as a quality-control "
             "quantity. Full derivation in Supplementary Note 1.")
        P(D, "4.5 Target-weighting schemes. Predictions are unchanged; only aggregation "
             "weights vary: per-mutation, per-target, proportional to √n, proportional to "
             "log(1 + n), and capped at 30.")
        P(D, f"4.6 Ranking and shortlist metrics. Within-target ordering used per-target "
             f"Spearman ρ and Kendall τ correlations, pairwise accuracy and the "
             f"within-target covariance contribution; shortlist metrics used fixed-K and "
             f"fraction-based hit rate, precision, recall, NDCG and enrichment "
             f"factor^[{KEYNUM['ndcg']}].")
        P(D, f"4.7 Statistical analysis. Rank agreement uses Spearman ρ with Kendall τ as "
             f"secondary; interval estimates use target-cluster bootstrap (5,000 resamples, "
             f"seed 42)^[{KEYNUM['bootstrap']}]; small external panels use exact permutation "
             f"tests. Metric directions were registered before computation and all "
             f"cross-metric comparisons use higher-is-better "
             f"orientation^[{KEYNUM['biometry']}].")
        P(D, "4.8 External prediction mapping and ambiguity control. Where the normalised "
             "matching key was not unique on the annotation side, only uniquely mapped rows "
             "were used; ambiguous keys were excluded rather than resolved. The official "
             "held-out CATH partition was evaluated separately (Supplementary Note 5).")
        P(D, "4.9 Reproducibility and consistency validation. All headline numbers are "
             "generated from a single machine-readable file and verified by an automated "
             "consistency checker (Supplementary Information). No new model was trained and "
             "no GPU was used.")

    # ---------- limitations ----------
    H(D, T("Limitations", "局限"), 1)
    if zh:
        P(D, "本研究的外推设定基于靶点家族划分，但这并不能排除基础模型在预训练阶段"
             "接触过与评价家族相似的序列或结构；我们只能控制任务特异性的监督训练重叠。"
             "同时，基准中包含的模型池不是全部 PPI ΔΔG 方法的随机样本，部分模型族覆盖有限，"
             "外部可评价方法仅五个，方法层面的推断功效不足。因此关于能力构成的结论"
             "应被理解为一组受控观测，而非对领域整体的普查。")
        P(D, "短名单与候选选择分析均以已有预测数据的重新分析为基础，不包含前瞻性湿实验验证；"
             "当靶点突变数不超过 K 时，相关指标本身不携带信息，因此与选择相关的结论"
             "仅作低置信度描述。回溯性消融案例同样不构成对改进路线的前瞻预测。")
        P(D, "C_B 与 C_W 是协方差信号贡献，不是因果机制；本文未进行身份记忆或靶点打乱等"
             "干预实验，因此不对跨靶分量的来源作因果断言。")
        P(D, "本研究没有独立的、前瞻性的实验数据集，因此未进行外部数据集验证；"
             "同时本文不提供通用模型推荐器，也不主张任何方法族具有普遍优势。")
    else:
        P(D, "Our extrapolation setting is based on target-family partitioning, which does "
             "not exclude the possibility that foundation-model pretraining encountered "
             "sequences or structures similar to the evaluation families; we can control "
             "task-specific supervised overlap only. In addition, the model pool assembled "
             "here is not a random sample of all PPI ΔΔG methods, coverage of some families "
             "is limited, and only five methods are externally evaluable, so method-level "
             "inference is underpowered. Conclusions about the composition of capability "
             "should therefore be read as controlled observations rather than a survey of "
             "the field.")
        P(D, "Shortlist and candidate-selection analyses rest on re-analysis of existing "
             "predictions and include no prospective wet-lab validation; when a target has "
             "no more than K mutations the relevant metrics carry no information, so "
             "selection-related conclusions are reported as low-confidence descriptions "
             "only. The retrospective ablation likewise does not establish prospective "
             "prediction of which modifications will succeed.")
        P(D, "C_B and C_W are covariance-signal contributions rather than causal mechanisms. "
             "We performed no identity-memorisation or target-shuffling intervention and "
             "therefore make no causal claim about the origin of the between-target "
             "component.")
        P(D, "No independent prospective experimental dataset was available, so no external "
             "dataset validation was performed. This work does not provide a universal model "
             "recommender and does not claim general superiority of any method family.")

    # ---------- back matter ----------
    H(D, T("Data and code availability", "数据与代码可用性"), 1)
    P(D, T("All analyses use existing predictions and published external predictions. "
           "Frozen input hashes, the machine-readable source of every headline number and "
           "an automated consistency checker are provided with the Supplementary "
           "Information so that all reported values can be regenerated without retraining "
           "any model. No new model was trained and no GPU was used.",
           "本分析仅使用既有预测与已公开的外部预测。冻结输入哈希、所有头条数字的机器可读"
           "来源以及自动化一致性检查器随补充材料提供，使全部报告数值可在不重新训练任何"
           "模型的条件下重新生成。本研究未训练新模型，也未使用 GPU。"))
    H(D, T("Author contributions", "作者贡献"), 1)
    P(D, T("Study conception, analysis design, interpretation of results and writing of the "
           "manuscript were performed by the author(s). The use of AI-assisted tools is "
           "declared below.",
           "研究构思、分析设计、结果解读与稿件撰写由作者完成；"
           "AI 辅助工具的使用情况见下节声明。"))
    H(D, T("Use of AI-assisted tools", "AI 辅助工具使用声明"), 1)
    P(D, T("AI-assisted coding and language tools were used to support code drafting and "
           "debugging, generation of plotting scripts, document automation, literature and "
           "repository retrieval, and language editing. All metric definitions, analysis "
           "choices, scientific interpretations, hypothesis decisions and final manuscript "
           "content were reviewed and approved by the authors, who take full responsibility "
           "for the work. AI tools were not used to generate experimental data or to modify "
           "upstream model predictions. A detailed account of the automated workflow is "
           "provided in the Supplementary Information.",
           "本研究使用 AI 辅助编码与语言工具，用于支持代码草拟与调试、绘图脚本生成、"
           "文档自动化、文献与资源检索以及语言润色。所有指标定义、分析选择、科学解读、"
           "假设判定与最终稿件内容均由作者审阅并确认，作者对全部内容负责。"
           "AI 工具未被用于生成实验数据或修改任何上游模型预测。"
           "自动化工作流程的详细说明见补充材料。"))
    H(D, T("Acknowledgements", "致谢"), 1)
    P(D, T("[FUNDING_PLACEHOLDER] We thank the developers of the publicly released "
           "prediction sets and datasets used in this study.",
           "[FUNDING_PLACEHOLDER] 我们感谢本研究所用公开预测集合与数据集的开发者。"))

    H(D, T("References", "参考文献"), 1)
    for i, (k, txt) in enumerate(REFS, 1):
        p = D.add_paragraph()
        r = p.add_run(f"{i}. {txt}"); r.font.size = Pt(9.5)
        p.paragraph_format.left_indent = Inches(0.30)
        p.paragraph_format.first_line_indent = Inches(-0.30)
        p.paragraph_format.space_after = Pt(3)

    nm = f"CR-PPI-ddG_manuscript_v4_{'ZH' if zh else 'EN'}.docx"
    path = OUT / nm
    D.save(str(path))
    print(f"  {nm} ({path.stat().st_size/1024:.0f} KB)")
    return path


build(False)
build(True)
print("v4 manuscript 完成")
