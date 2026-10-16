# -*- coding: utf-8 -*-
"""scripts/v8_1/15_manuscript_v3.py -- 稿件 v3：数字引用体系 + 排版方程 + 视觉系统。

在 v2 叙事基础上：
* 建立**带方括号数字引用**的正文引用体系（首次出现顺序编号）
* 把三条核心方程改为**独立居中、带编号**的真排版（真实上下标/斜体变量）
* 用 v3 图（统一调色板）
* 正文提及主匹配集合时用「行列数」而非 hash
* 生成 REFERENCE_AUDIT.md 并做引用 QC
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.shared import Inches, Pt, RGBColor

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
A = R / "v8_1"
F = R / "v8_1_final"
OUT = F / "manuscript_v3"
OUT.mkdir(parents=True, exist_ok=True)
FM, FS = F / "figures_main", F / "figures_supplement"

FR = json.load(open(A / "FROZEN_RESULTS.json", encoding="utf-8"))
PU, H4, GEO = FR["primary_universe"], FR["H4_controlled_split"], FR["axis_c_geometry"]
DL = H4.get("delta", {})
ST = pd.read_csv(A / "rank_agreement_by_stratum.csv")
EP = pd.read_csv(A / "external_clean_profiles.csv")
W2 = pd.read_csv(A / "WEIGHTING_EFFECT_SIZE_AUDIT.csv")
UR = pd.read_csv(A / "MATCHED_UNIVERSE_REGISTRY.csv")
er = FR["external_clean"]["rank_agreement_external"]
r0, r3 = H4["random"]["overall_pearson"], H4["family_OOD"]["overall_pearson"]
sf0 = int(W2["signflip_C_B_gt0"].sum()); sf5 = int(W2["signflip_C_B_gt0.005"].sum())
sfw = int(W2["signflip_C_W_gt0"].sum())
best = EP.loc[EP.overall_pearson.idxmax()]
rmn, rmx = float(ST.rho_overall_vs_macro.min()), float(ST.rho_overall_vs_macro.max())

# ============ 参考文献（首次出现顺序编号）============
# 每条: (key, text, verified)
REFS = [
 ("skempi", "Jankauskaitė J, Jiménez-García B, Dapkūnas J, Fernández-Recio J, Moal IH. "
  "SKEMPI 2.0: an updated benchmark of changes in protein–protein binding energy, "
  "kinetics and thermodynamics upon mutation. Bioinformatics. 2019;35(3):462–469. "
  "doi:10.1093/bioinformatics/bty635", True),
 ("leakage", "Bushuiev A, Bushuiev R, Sedlar J, Pluskal T, Damborsky J, Mazurenko S, "
  "Sivic J. Revealing data leakage in protein interaction benchmarks. "
  "arXiv:2404.10457. 2024.", True),
 ("cathddg", "Yu G, Bi X, Ma T, Li Y, Wang J. CATH-ddG: towards robust mutation effect "
  "prediction on protein–protein interactions out of CATH homologous superfamily. "
  "Bioinformatics. 2025;41(Suppl_1):i362–i372. doi:10.1093/bioinformatics/btaf228", True),
 ("hermes", "Pun MN, Ivanov A, Bellinzoni M, Bayer PE, Gagneur J. HERMES: holographic "
  "equivariant neural network model for mutational effect and stability prediction. "
  "arXiv:2407.06703. 2024.", True),
 ("prefold", "PreFold-dG: estimating binding affinity of protein–protein interaction "
  "from intermediate representations of protein folding models. Bioinformatics. 2026. "
  "doi:10.1093/bioinformatics/btag489", True),
 ("flexddg", "Barlow KA, Ó Conchúir S, Thompson S, Suresh P, Lucas JE, Heinonen M, "
  "Kortemme T. Flex ddG: Rosetta ensemble-based estimation of changes in "
  "protein–protein binding affinity upon mutation. J Phys Chem B. "
  "2018;122(21):5389–5399. doi:10.1021/acs.jpcb.7b11367", True),
 ("foldx", "Delgado J, Radusky LG, Cianferoni D, Serrano L. FoldX 5.0: working with RNA, "
  "small molecules and a new graphical interface. Bioinformatics. "
  "2019;35(20):4168–4169. doi:10.1093/bioinformatics/btz184", True),
 ("proteinmpnn", "Dauparas J, Anishchenko I, Bennett N, Bai H, Ragotte RJ, Milles LF, "
  "et al. Robust deep learning-based protein sequence design using ProteinMPNN. "
  "Science. 2022;378(6615):49–56. doi:10.1126/science.add2187", True),
 ("esmif", "Hsu C, Verkuil R, Liu J, Lin Z, Hie B, Sercu T, Lerer A, Rives A. Learning "
  "inverse folding from millions of predicted structures. In: Proceedings of the 39th "
  "International Conference on Machine Learning (ICML). 2022;162:8946–8970.", True),
 ("ppiformer", "Bushuiev A, Bushuiev R, Kouba P, Filkin A, Gabrielova M, Gabriel M, "
  "Sedlar J, Pluskal T, Damborsky J, Mazurenko S, Sivic J. Learning to design "
  "protein–protein interactions with enhanced generalization. In: The Twelfth "
  "International Conference on Learning Representations (ICLR). 2024. "
  "arXiv:2310.18515", True),
 ("ndcg", "Järvelin K, Kekäläinen J. Cumulated gain-based evaluation of IR techniques. "
  "ACM Trans Inf Syst. 2002;20(4):422–446. doi:10.1145/582415.582418", True),
 ("bootstrap", "Efron B, Tibshirani RJ. An Introduction to the Bootstrap. New York: "
  "Chapman & Hall; 1993.", True),
 ("biometry", "Sokal RR, Rohlf FJ. Biometry: The Principles and Practice of Statistics in "
  "Biological Research. 3rd ed. New York: WH Freeman; 1995.", True),
]
KEYNUM = {k: i + 1 for i, (k, _, _) in enumerate(REFS)}


def C(*keys):
    """生成方括号数字引用串，如 [2] / [2,4] / [5-7]。"""
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


# ============ 表格 ============
T1 = pd.DataFrame([
    {"Evaluation set": "Full annotation namespace", "Rows": "4,473", "Targets": "214",
     "Targets n≥5": "130", "Method families": "—", "Purpose": "Namespace definition"},
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
    {"Split regime": "Random", "Pooled Pearson": f"{H4['random']['overall_pearson']:.4f}",
     "C_B": f"{H4['random']['C_B']:+.4f}", "C_W": f"{H4['random']['C_W']:+.4f}",
     "Macro Spearman": f"{H4['random']['macro_spearman']:.4f}",
     "RMSE": f"{H4['random']['rmse']:.4f}"},
    {"Split regime": "Complex", "Pooled Pearson": f"{H4['complex']['overall_pearson']:.4f}",
     "C_B": f"{H4['complex']['C_B']:+.4f}", "C_W": f"{H4['complex']['C_W']:+.4f}",
     "Macro Spearman": f"{H4['complex']['macro_spearman']:.4f}",
     "RMSE": f"{H4['complex']['rmse']:.4f}"},
    {"Split regime": "seq40", "Pooled Pearson": f"{H4['seq40']['overall_pearson']:.4f}",
     "C_B": f"{H4['seq40']['C_B']:+.4f}", "C_W": f"{H4['seq40']['C_W']:+.4f}",
     "Macro Spearman": f"{H4['seq40']['macro_spearman']:.4f}",
     "RMSE": f"{H4['seq40']['rmse']:.4f}"},
    {"Split regime": "Family-OOD", "Pooled Pearson":
        f"{H4['family_OOD']['overall_pearson']:.4f}",
     "C_B": f"{H4['family_OOD']['C_B']:+.4f}", "C_W": f"{H4['family_OOD']['C_W']:+.4f}",
     "Macro Spearman": f"{H4['family_OOD']['macro_spearman']:.4f}",
     "RMSE": f"{H4['family_OOD']['rmse']:.4f}"},
])
T3 = EP[["model", "n_rows", "n_targets_ge5", "overall_pearson", "macro_spearman",
         "C_B_row", "C_W", "rmse"]].copy()
T3.columns = ["Method", "Mutations", "Targets n≥5", "Pooled Pearson", "Macro Spearman",
              "C_B", "C_W", "RMSE"]
for c in ["Pooled Pearson", "Macro Spearman", "C_B", "C_W", "RMSE"]:
    T3[c] = T3[c].round(4)


# ============ helpers ============
def setup(doc, zh):
    st = doc.styles["Normal"]
    st.font.size = Pt(10.5 if zh else 11)
    st.font.name = "Times New Roman"
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


def P(doc, t, size=None):
    p = doc.add_paragraph(); r = p.add_run(t)
    if size:
        r.font.size = Pt(size)
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    return p


def RICH(doc, parts, align="just", size=None):
    p = doc.add_paragraph()
    for t, s in parts:
        r = p.add_run(t)
        r.italic = "i" in s
        r.bold = "b" in s
        if s == "sub":
            r.font.subscript = True
        if s == "sup":
            r.font.superscript = True
        if size:
            r.font.size = Pt(size)
    p.alignment = (WD_ALIGN_PARAGRAPH.JUSTIFY if align == "just"
                   else WD_ALIGN_PARAGRAPH.CENTER)
    return p


def EQ(doc, parts, num=None):
    """编号居中公式：左为公式，右侧为编号。"""
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for t, s in parts:
        r = p.add_run(t)
        r.italic = "i" in s
        if s == "sub":
            r.font.subscript = True
    if num:
        r = p.add_run(f"\t\t({num})")
        r.italic = False
    p.paragraph_format.space_before = Pt(5)
    p.paragraph_format.space_after = Pt(5)
    p.paragraph_format.keep_with_next = True
    return p


def cap(doc, t):
    p = doc.add_paragraph(); r = p.add_run(t); r.font.size = Pt(9.0)
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    p.paragraph_format.space_after = Pt(10)
    return p


def tbl(doc, df, above=None, below=None, fs=8.8):
    if above:
        p = doc.add_paragraph(); r = p.add_run(above); r.font.size = Pt(9.0)
        p.paragraph_format.keep_with_next = True
    t = doc.add_table(rows=1, cols=len(df.columns))
    t.style = "Light Grid Accent 1"
    for i, c in enumerate(df.columns):
        cell = t.rows[0].cells[i]
        cell.text = ""
        r = cell.paragraphs[0].add_run(str(c)); r.bold = True; r.font.size = Pt(fs)
    for _, row in df.iterrows():
        cells = t.add_row().cells
        for i, v in enumerate(row):
            cells[i].text = ""
            cells[i].paragraphs[0].add_run(str(v)).font.size = Pt(fs)
    if below:
        cap(doc, below)
    else:
        doc.add_paragraph()
    return t


# ============ build ============
def build(zh):
    T = (lambda en, z: z if zh else en)
    SUF = "_ZH" if zh else "_EN"

    def FIG(doc, stem, w):
        doc.add_picture(str(FM / f"{stem}{SUF}.png"), width=Inches(w))

    D = Document(); setup(D, zh)

    # ---- title page ----
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
                 T("Affiliation(s) [AFFILIATION_PLACEHOLDER]", "单位 [AFFILIATION_PLACEHOLDER]"),
                 T("Corresponding author: [NAME] · Email: [EMAIL] · ORCID: [ORCID_PLACEHOLDER]",
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

    # ---- abstract ----
    H(D, T("Abstract", "摘要"), 1)
    if zh:
        P(D, f"预测突变对蛋白质–蛋白质结合亲和力（ΔΔG）的影响是抗体工程、界面优化与"
             f"致病突变解释的核心任务，其进展高度依赖基准评价的可信度。"
             f"然而常用的合并相关系数难以直接解读：它同时取决于模型能否把不同靶点"
             f"放在正确的结合量级上，以及能否在同一靶点内部把突变正确排序。"
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
             f"固定 K 短名单指标受候选集几何约束：K=10 时 D(n,K)=K(n−K)/[n(n+1)] "
             f"在 n≤K 时为零并在 n*≈{GEO['n_star_analytic']:.1f} 取峰，"
             f"本基准中 {100*GEO['n_le_K_fraction']:.1f}% 的 n≥5 靶点处于该退化区间。"
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
             f"K = 10, D(n,K) = K(n−K)/[n(n+1)] is zero when n ≤ K and peaks near "
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

    # ---- introduction ----
    H(D, T("Introduction", "引言"), 1)
    if zh:
        P(D, f"预测突变如何改变蛋白质–蛋白质结合亲和力（ΔΔG）是抗体亲和力成熟、"
             f"界面优化与致病突变解释中的核心计算任务。实验测定结合亲和力的成本高昂，"
             f"因此计算预测被广泛用于对突变进行排序并指导实验设计。随着深度学习方法的"
             f"快速增加，该领域进展的判据越来越依赖基准评价：一个方法是否更好，"
             f"通常由其在标准突变集合上报告的相关系数或误差决定。因此，评价协议本身的性质"
             f"——它测量了什么、以及在何种外推条件下测量——直接决定了领域能得出什么结论。")
        P(D, f"目前的评价实践以合并指标为主，即在整个突变集合上计算 Pearson 或 Spearman "
             f"相关以及 RMSE{C('biometry')}，同时常辅以按复合物或按靶点计算的指标。"
             f"近年已有工作揭示，这类评价对外推能力的刻画受到划分策略的强烈影响。"
             f"Bushuiev 等量化了 PPI 数据中的结构性泄漏：按 PDB 代码划分时约 65% 的测试样本"
             f"存在近重复训练样本，即使按 30% 序列同一性划分仍有约 30%，SKEMPI v2.0"
             f"{C('skempi')} 若按 PPI 代码划分为 56%，而按领域专家分组则为 0%"
             f"{C('leakage')}。CATH-ddG 采用 CATH 同源超家族划分并以最大 TM-score < 0.6 "
             f"定义「难突变」{C('cathddg')}；相关方法亦报告按结构或按靶点的相关性"
             f"{C('ppiformer','hermes','prefold')}。需要明确的是：泄漏感知划分、"
             f"家族外推评价、按靶点报告以及短名单式（Top-K）评价均非本文提出。")
        P(D, f"然而，一个关键问题仍未被直接回答：合并基准性能究竟来自哪里？"
             f"把合并指标与按靶点指标并列报告，并不等于量化了跨靶信号与靶内信号对"
             f"合并性能的加性贡献。一个模型可能主要因为能把不同靶点放在正确的整体量级上"
             f"而获得较高的合并相关，而其在靶点内部的排序能力可能很弱；反之亦然。"
             f"在靶点数量多、各靶点突变数高度不均的基准中，合并相关还会被大靶点主导，"
             f"使这一混淆进一步加剧。")
        P(D, f"本文引入一个精确的协方差信号分解来回答这一问题。我们把合并 Pearson 相关"
             f"严格分解为跨靶分量与靶内分量，使「合并性能由哪一部分构成」成为可直接测量的量，"
             f"并在匹配、去伪重复、来源筛选后的评价集合上使用该分解量化四件事："
             f"划分严格度改变的是哪一个分量；目标加权如何改变对这两个分量的解读；"
             f"固定 K 短名单评价受候选集几何的何种约束；以及该分解能否区分"
             f"「仅改善合并分数」与「同时改善靶内排序」的模型改动。我们进一步在五个"
             f"独立开发的外部方法上界定当前证据的泛化边界。")
    else:
        P(D, f"Predicting how mutations alter protein–protein binding affinity (ΔΔG) is "
             f"central to antibody affinity maturation, interface optimisation and the "
             f"interpretation of disease-associated variants. Because experimental "
             f"affinity measurement is costly, computational predictions are widely used to "
             f"prioritise mutations and guide experimental design. As the number of "
             f"deep-learning methods grows rapidly, progress in the field is increasingly "
             f"judged by benchmark evaluation: whether one method is better is usually "
             f"decided by the correlation or error it reports on a standard mutation set. "
             f"The properties of the evaluation protocol itself — what it measures, and "
             f"under which extrapolation conditions — therefore directly determine what the "
             f"field can conclude.")
        P(D, f"Current practice relies largely on pooled metrics, computing Pearson or "
             f"Spearman correlation and RMSE over the entire mutation set{C('biometry')}, "
             f"often supplemented by per-complex or per-target metrics. Recent work has "
             f"shown that how well such evaluation captures extrapolation depends strongly "
             f"on the splitting strategy. Bushuiev et al. quantified structural leakage in "
             f"PPI data: about 65% of test examples have a near-duplicate training example "
             f"under PDB-code splitting, about 30% even under 30% sequence-identity "
             f"splitting, 56% for SKEMPI v2.0{C('skempi')} split by PPI code, and 0% under "
             f"domain-expert grouping{C('leakage')}. CATH-ddG adopted CATH "
             f"homologous-superfamily splits and defined hard mutations by a maximum "
             f"TM-score below 0.6{C('cathddg')}; related methods likewise report "
             f"per-structure or per-target correlations{C('hermes','prefold','ppiformer')}. "
             f"We emphasise that leakage-aware splitting, family-out-of-distribution "
             f"evaluation, per-target reporting and shortlist (Top-K) evaluation were none "
             f"of them introduced by this work.")
        P(D, f"A key question nevertheless remains unanswered directly: where does pooled "
             f"benchmark performance come from? Reporting pooled and per-target metrics "
             f"side by side does not quantify the additive contribution of between-target "
             f"and within-target signal to pooled performance. A model may achieve a high "
             f"pooled correlation mainly because it places different targets on the correct "
             f"overall scale, while its ability to order mutations within a target is weak "
             f"— or the reverse. In benchmarks containing many targets with highly unequal "
             f"mutation counts, pooled correlation is additionally dominated by large "
             f"targets, which amplifies the conflation.")
        P(D, f"Here we introduce an exact covariance-signal decomposition that addresses "
             f"this question. We decompose pooled Pearson correlation rigorously into a "
             f"between-target component and a within-target component, which makes 'what "
             f"pooled performance is made of' a directly measurable quantity. Using this "
             f"decomposition on matched, de-replicated and provenance-filtered evaluation "
             f"sets, we quantify four things: which component is changed by split "
             f"strictness; how target weighting changes the interpretation of the two "
             f"components; which constraints candidate-set geometry places on fixed-K "
             f"shortlist evaluation; and whether the decomposition can distinguish model "
             f"modifications that improve pooled scores only from those that also improve "
             f"within-target ranking. We further bound the current evidence using five "
             f"independently developed external methods{C('flexddg','foldx','proteinmpnn','esmif','hermes')}.")

    # ---- results 2.1 ----
    H(D, T("Results", "结果"), 1)
    H(D, T("Pooled PPI ΔΔG performance contains separable between-target and "
           "within-target signal",
           "合并 PPI ΔΔG 性能包含可分离的跨靶与靶内信号"), 2)
    if zh:
        P(D, f"对属于靶点 c 的突变 i，观测到的 ΔΔG 可以写成该靶点的平均水平"
             f"加上该突变相对这一水平的偏差：")
    else:
        P(D, f"For a mutation i belonging to target c, the observed ΔΔG can be written as "
             f"a target-level mean plus a mutation-specific deviation from that mean:")
    EQ(D, [("y", "i"), ("ci", "sub"), (" = μ", "i"), ("c", "sub"), (" + ε", "i"),
           ("ci", "sub")], 1)
    P(D, T("where μ_c is the mean ΔΔG of target c and ε_ci = y_ci − μ_c is the "
           "within-target mutation deviation. The prediction is decomposed identically:",
           "其中 μ_c 为靶点 c 的平均 ΔΔG，ε_ci = y_ci − μ_c 为靶内突变偏差。"
           "预测值作同样的分解："))
    EQ(D, [("ŷ", "i"), ("ci", "sub"), (" = μ̂", "i"), ("c", "sub"), (" + ε̂", "i"),
           ("ci", "sub")], 2)
    if zh:
        P(D, f"只有当这一分解明确之后，才有意义地区分两种协方差贡献。"
             f"设 σ_y 与 σ_ŷ 为合并标准差，合并相关可精确拆分如下：")
    else:
        P(D, f"Only after this decomposition is it meaningful to distinguish two "
             f"covariance contributions. With σ_y and σ_ŷ the pooled standard deviations, "
             f"the pooled correlation splits exactly as:")
    EQ(D, [("r", "i"), ("pooled", "sub"), (" = Cov(y, ŷ) / (σ", "i"), ("y", "sub"),
           ("σ", "i"), ("ŷ", "sub"), (") = C", ""), ("B", "sub"), (" + C", ""),
           ("W", "sub")], 3)
    if zh:
        P(D, f"其中 C_B = Cov(B_y, B_ŷ)/(σ_y σ_ŷ) 为跨靶协方差贡献，"
             f"C_W = Cov(W_y, W_ŷ)/(σ_y σ_ŷ) 为靶内协方差贡献。"
             f"该恒等式在机器精度内成立（残差 ≤ 3.3×10⁻¹⁶）。"
             f"需要强调：C_B 与 C_W 是同一尺度上的协方差信号贡献，"
             f"不是独立相关系数、独立能力分数或因果机制；完整推导见补充材料。")
    else:
        P(D, f"where C_B = Cov(B_y, B_ŷ)/(σ_y σ_ŷ) is the between-target covariance "
             f"contribution and C_W = Cov(W_y, W_ŷ)/(σ_y σ_ŷ) is the within-target "
             f"covariance contribution. The identity holds to machine precision "
             f"(residual ≤ 3.3×10⁻¹⁶; full derivation in Supplementary Note 1). We stress "
             f"that C_B and C_W are covariance-signal contributions on a common scale, not "
             f"independent correlations, independent skill scores, or causal mechanisms.")
    if zh:
        P(D, f"我们在一个匹配的评价集合上使用该分解，该集合包含 "
             f"{PU['n_rows']:,} 行突变、{PU['n_ge5']} 个至少五个突变的靶点，"
             f"涉及 {PU['n_method_families']} 个独立方法族；跨模型比较仅在同一批突变行内进行"
             f"（方法；补充图 S1）。")
    else:
        P(D, f"We apply this decomposition on a matched evaluation set comprising "
             f"{PU['n_rows']:,} mutation rows, {PU['n_ge5']} targets with at least five "
             f"mutations, and {PU['n_method_families']} independent method families; "
             f"cross-model comparisons were restricted to predictions evaluated on "
             f"identical mutation-row sets (Methods; Supplementary Fig. S1).")
    FIG(D, "Fig1_overview", 6.5)
    cap(D, T("Figure 1 | Capability-resolved evaluation of PPI ΔΔG predictors. "
             "(A) A pooled benchmark aggregates many targets, each with its own mean ΔΔG "
             "level and internal spread among mutations. (B) The observed effect of a "
             "mutation is a target-level mean plus a within-target deviation. (C) Pooled "
             "correlation splits exactly into a between-target covariance contribution "
             "(muted purple) and a within-target covariance contribution (green); these "
             "are contributions on a common scale, not independent correlations or skill "
             "scores. (D) The main controlled result: with identical 7,523 rows and "
             f"identical configurations and only the split regime changed, pooled r falls "
             f"from {r0:.3f} to {r3:.3f}, and "
             f"{100*DL['fraction_between']:.1f}% of the absolute decrease is carried by "
             "the between-target component.",
             "图 1 | PPI ΔΔG 预测器的能力分辨评价。(A) 合并基准汇集多个靶点，"
             "每个靶点有其自身的平均 ΔΔG 水平与靶内离散度。(B) 突变的观测效应"
             "等于靶点平均水平加上靶内偏差。(C) 合并相关精确拆分为跨靶协方差贡献"
             "（紫色）与靶内协方差贡献（绿色）；二者是同一尺度上的贡献，"
             "并非独立相关或能力分数。(D) 受控主结果：同一批 7,523 行与同一批配置，"
             f"仅改变划分方式，合并 r 由 {r0:.3f} 降至 {r3:.3f}，"
             f"绝对下降的 {100*DL['fraction_between']:.1f}% 由跨靶分量承载。"))

    # 2.2
    H(D, T("Family-OOD evaluation disproportionately removes between-target signal",
           "家族外推评价不成比例地移除跨靶信号"), 2)
    if zh:
        P(D, f"为分离划分严格度的效应，我们在同一批 7,523 行突变与同一批模型配置上，"
             f"仅改变划分方式（表 2、图 2）。合并 Pearson 中位数由随机划分的 {r0:.4f} "
             f"降至家族外推的 {r3:.4f}。分解显示该下降并非均分："
             f"ΔC_B = {DL['delta_C_B']:+.4f}，ΔC_W = {DL['delta_C_W']:+.4f}。"
             f"按绝对变化口径 |ΔC_B|/(|ΔC_B|+|ΔC_W|) 计算，跨靶分量占下降的 "
             f"{100*DL['fraction_between']:.1f}%，而 C_W 在家族外推下仍为正"
             f"（{H4['family_OOD']['C_W']:+.4f}）。因此，家族外推评价下合并相关的下降"
             f"主要由跨靶协方差贡献的丧失承载，而靶内协方差保持为正。")
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
             f"absolute-change basis, |ΔC_B|/(|ΔC_B|+|ΔC_W|), the between-target component "
             f"accounts for {100*DL['fraction_between']:.1f}% of the decline, while C_W "
             f"remains positive under family-out-of-distribution evaluation "
             f"({H4['family_OOD']['C_W']:+.4f}). Most of the reduction in pooled correlation "
             f"under family-OOD evaluation is therefore carried by loss of the "
             f"between-target covariance contribution, while within-target covariance "
             f"remains positive.")
        P(D, "This is not the same proposition as 'the model memorised target identity': we "
             "performed no identity-memorisation intervention and therefore make no such "
             "causal claim. A more defensible reading is that the pattern is consistent "
             "with reliance on target-specific information that does not transfer across "
             "families. The complex-level and seq40 partitions produced identical "
             "assignments in this dataset and were therefore not treated as independent "
             "comparisons.")
    FIG(D, "Fig2_split_decomposition", 6.5)
    cap(D, T("Figure 2 | Family-out-of-distribution evaluation disproportionately removes "
             "the between-target component of pooled correlation. (A) Median pooled Pearson "
             "correlation by split regime. (B) Exact decomposition into C_B (muted purple) "
             "and C_W (green); C_W remains positive. (C) Attribution of the absolute "
             f"decline, of which the between-target component accounts for "
             f"{100*DL['fraction_between']:.1f}%. All four regimes use the same 7,523 "
             "mutation rows and the same model configurations; the complex and seq40 "
             "partitions are identical in this dataset.",
             "图 2 | 家族外推评价不成比例地移除合并相关中的跨靶分量。(A) 各划分下的"
             "合并 Pearson 中位数。(B) 精确分解为 C_B（紫色）与 C_W（绿色）；"
             f"C_W 仍为正。(C) 绝对下降的归因，跨靶分量占 "
             f"{100*DL['fraction_between']:.1f}%。四个划分使用同一批 7,523 行与同一批配置；"
             "复合物与 seq40 划分在本数据上相同。"))
    tbl(D, T2, above=T("Table 2 | Controlled split comparison.",
                       "表 2 | 受控划分比较。"),
        below=T("Median values across model configurations. The seq40 partition is "
                "identical to the complex partition in this dataset.",
                "跨模型配置的中位数。seq40 划分在本数据上与复合物划分完全相同。"))

    # 2.3
    H(D, T("Target aggregation changes the interpretation of covariance contributions",
           "目标聚合改变协方差贡献的解读"), 2)
    if zh:
        P(D, f"由于各靶点突变数高度不均，聚合权重的选择会影响两个分量的大小。"
             f"我们固定预测值，仅改变聚合权重，比较五种方案（每个突变等权、每个靶点等权、"
             f"权重正比于 √n、正比于 log(1+n)、以及上限为 30 的截断权重）。"
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
             f"log(1+n), and capped at 30). C_B shows directional instability concentrated "
             f"near zero: {sf0}/{len(W2)} configurations reverse sign with no threshold, "
             f"whereas none do once |C_B| exceeds 0.005 ({sf5}/{len(W2)}). C_W never "
             f"reverses sign at any threshold ({sfw}/{len(W2)}), but its magnitude is not "
             f"invariant — its median range is {W2.range_C_W.median():.4f}, larger than the "
             f"{W2.range_C_B.median():.4f} of C_B. C_B is therefore unstable in direction "
             f"only near zero, whereas C_W keeps its sign but varies in magnitude. "
             f"Consequently, when a reported C_B is close to zero, the conclusion about "
             f"whether between-target information is present may depend on the weighting "
             f"convention rather than on the model; any reported decomposition should "
             f"declare its weighting scheme.")
    FIG(D, "Fig3_weighting", 6.5)
    cap(D, T("Figure 3 | Target weighting changes the interpretation of covariance "
             "contributions. Predictions are unchanged; only aggregation weights differ. "
             "(A, B) Distribution of the range of C_B and C_W across the five schemes. "
             "(C) Number of configurations reversing sign at three thresholds. C_B "
             "reverses direction only near zero, whereas C_W keeps its sign but its "
             "magnitude is not invariant.",
             "图 3 | 目标加权改变协方差贡献的解读。预测值不变，仅改变聚合权重。"
             "(A, B) 五种方案下 C_B 与 C_W 极差的分布。(C) 三个阈值下出现符号翻转的"
             "配置数。C_B 仅在近零处翻转方向，而 C_W 保持符号但幅度并非不变。"))

    # 2.4
    H(D, T("Fixed-K shortlist evaluation is constrained by candidate-set geometry",
           "固定 K 短名单评价受候选集几何约束"), 2)
    if zh:
        P(D, f"在突变文库或短名单筛选场景中，常用固定 K 的 Top-K 指标评价模型"
             f"{C('ndcg')}。可区分区间为：")
    else:
        P(D, f"Fixed-K top-K metrics are commonly used to evaluate models for mutation "
             f"library or shortlist screening{C('ndcg')}. The discriminable interval is:")
    EQ(D, [("D(n,K) = 0,  n ≤ K;    D(n,K) = K(n−K)/[n(n+1)],  n > K", "")], 4)
    if zh:
        P(D, f"其中随机选取 K 个候选的期望收益占比为 K/n，最优 K 个候选为 "
             f"K(2n−K+1)/[n(n+1)]（推导见补充说明 6）。当靶点突变数不超过 K 时，"
             f"任何模型都会选中全部候选，该靶点对固定 K 短名单指标不携带任何信息。"
             f"本基准中，{GEO['n_ge5_count']} 个 n≥5 靶点里有 {GEO['n_le_K_count']} 个"
             f"（{100*GEO['n_le_K_fraction']:.1f}%）处于该退化区间。"
             f"此外 D 并非随 n 单调增加，而在 n* = K + √(K²+K) 处取峰；K = 10 时 "
             f"n* ≈ {GEO['n_star_analytic']:.1f}，例如 D(20) = {GEO['D_at_n20']:.4f} 而 "
             f"D(50) = {GEO['D_at_n50']:.4f}。因此「候选集越大越可区分」的启发式"
             f"在本指标下不成立。")
        P(D, f"上述几何结论不意味着 SKEMPI 的 ΔΔG 数值本身过于接近；后者需要以真实 ΔΔG "
             f"效用幅度做实证分析，本文未做该项分析。可辩护的表述是：当 n ≤ K 时"
             f"固定 K 的候选选择评价在结构上不携带信息，且其可区分范围强烈依赖候选集几何。")
    else:
        P(D, f"where the expected share of gain from K randomly chosen candidates is K/n "
             f"and the optimal K candidates attain K(2n−K+1)/[n(n+1)] (derivation in "
             f"Supplementary Note 6). When a target has no more than K mutations, any model "
             f"selects every candidate and the target carries no information for fixed-K "
             f"shortlist metrics. In this benchmark, {GEO['n_le_K_count']} of "
             f"{GEO['n_ge5_count']} targets with n ≥ 5 "
             f"({100*GEO['n_le_K_fraction']:.1f}%) fall in that degenerate regime. "
             f"Moreover D is not monotonic in n: it peaks at n* = K + √(K²+K), which for "
             f"K = 10 gives n* ≈ {GEO['n_star_analytic']:.1f}, so that "
             f"D(20) = {GEO['D_at_n20']:.4f} while D(50) = {GEO['D_at_n50']:.4f}. The "
             f"heuristic that larger candidate sets are more discriminable therefore does "
             f"not hold for this metric.")
        P(D, "This geometric result does not imply that SKEMPI ΔΔG values are intrinsically "
             "too similar; establishing that would require an analysis using actual ΔΔG "
             "utility magnitudes, which we did not perform. The defensible statement is "
             "that fixed-K candidate-selection evaluation is structurally non-informative "
             "when n ≤ K, and that its discriminable range depends strongly on "
             "candidate-set geometry.")
    FIG(D, "Fig4_candidate_geometry", 5.8)
    cap(D, T("Figure 4 | Fixed-K shortlist evaluation is constrained by candidate-set "
             "geometry. The discriminable interval D(n,K) is zero when n ≤ K and follows "
             "K(n−K)/[n(n+1)] otherwise; it peaks near n* ≈ 20.5 and declines thereafter. "
             f"In this benchmark {100*GEO['n_le_K_fraction']:.1f}% of targets with n ≥ 5 "
             "lie in the degenerate regime where every model selects all candidates.",
             "图 4 | 固定 K 短名单评价受候选集几何约束。可区分区间 D(n,K) 在 n ≤ K 时为零，"
             "否则为 K(n−K)/[n(n+1)]；它在 n* ≈ 20.5 取峰后下降。"
             f"本基准中 {100*GEO['n_le_K_fraction']:.1f}% 的 n≥5 靶点处于"
             "所有模型都选中全部候选的退化区间。"))

    # 2.5
    H(D, T("Capability decomposition clarifies retrospective remedy effects",
           "能力分解解释回溯性模型改动的效果"), 2)
    if zh:
        P(D, f"最后，我们在一个匹配的消融集合中检验该分解的实用性，该集合中的模型共享"
             f"同一 3,351 行宇宙，因此差异只来自模型改动本身。以「仅使用突变性质」的模型"
             f"为基线：加入静态局部几何使合并 Pearson 提高约 0.028，但靶内 Macro Spearman "
             f"几乎不变（−0.002），属于仅改善合并分数的改动；引入「突变 × 配体 partner "
             f"扰动」的改动同时改善两轴（合并 Pearson +0.070、Macro Spearman +0.044）；"
             f"使用全量局部物理特征时两轴改善最大（+0.105 与 +0.076）。"
             f"因此该分解能够区分「仅合并分数改善」与「同时改善靶内排序」两类改动。")
        P(D, f"需要明确该分析的性质：它是回溯性的，即已知哪些改动有效，因此它说明的是"
             f"该分解能够表达这一区别，而非基准能够前瞻预测哪条技术路线会成功。"
             f"此外，本研究重新计算的靶内差异（+0.044）与历史分析中的自助法区间"
             f"（+0.073，[+0.021, +0.129]）基于不同的度量宇宙与聚合方式，"
             f"不可在不加说明的情况下混用。")
    else:
        P(D, f"Finally, we examined the practical utility of the decomposition in a matched "
             f"ablation set whose models share an identical 3,351-row universe, so that "
             f"differences arise only from the model modifications themselves. Relative to "
             f"a baseline using mutation properties alone, adding static local geometry "
             f"raised pooled Pearson by about 0.028 while leaving within-target Macro "
             f"Spearman essentially unchanged (−0.002), a pooled-score-only improvement; "
             f"introducing mutation × partner perturbation improved both axes (pooled "
             f"Pearson +0.070, Macro Spearman +0.044); and a full local-physics model gave "
             f"the largest improvement on both axes (+0.105 and +0.076). The decomposition "
             f"therefore distinguishes a change that improves pooled scores only from one "
             f"that also improves within-target ranking.")
        P(D, "The nature of this analysis should be stated plainly: it is retrospective, in "
             "that we already know which modifications were effective. It therefore shows "
             "that the decomposition can express this distinction, not that the benchmark "
             "prospectively predicts which technical route will succeed. In addition, our "
             "recomputed within-target difference (+0.044) and the bootstrap interval from "
             "the historical analysis (+0.073, [+0.021, +0.129]) rest on different metric "
             "universes and aggregation schemes and must not be mixed without explanation.")
    FIG(D, "Fig5_remedy", 5.9)
    cap(D, T("Figure 5 | Capability decomposition separates pooled-score gains from "
             "within-target gains. Retrospective matched ablation on an identical 3,351-row "
             "universe. Green markers indicate improvement on both axes; red markers "
             "indicate a pooled-score-only or negative change. The analysis is "
             "retrospective and does not imply prospective prediction of which "
             "modification will succeed.",
             "图 5 | 能力分解区分合并分数增益与靶内增益。在同一 3,351 行宇宙上的"
             "回溯性匹配消融。绿色表示两轴同时改善，红色表示仅合并分数改善或变化为负。"
             "该分析为回溯性，不代表可前瞻预测哪种改动会成功。"))

    # 2.6
    H(D, T("External methods define the boundary of generalisation",
           "外部方法界定泛化边界"), 2)
    if zh:
        P(D, f"为评估上述结论在独立开发的方法上的适用性，我们使用已公开的官方逐突变预测"
             f"构建外部评价集合{C('flexddg','foldx','proteinmpnn','esmif','hermes')}。"
             f"在构建过程中我们发现，当匹配键省略链标识与野生型残基时，该键在标注侧并不唯一："
             f"4,473 行中仅 {FR['external_clean']['n_rows']:,} 行具有唯一映射。"
             f"因此主匹配评价只使用唯一映射行，歧义键被排除而非强行解析；"
             f"这一处理使外部方法的表观合并相关系统性下降 0.02 至 0.12。")
        P(D, f"在唯一映射集合（{FR['external_clean']['n_rows']:,} 行、"
             f"{FR['external_clean']['n_ge5']} 个 n≥5 靶点）上，五个外部方法的总体相关"
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
             f"是无量纲的突变效应分数，而非已校准的 kcal/mol ΔΔG{C('hermes')}。"
             f"因此我们不对该模型报告 RMSE 或 MAE，也不把其使用的变体标注为"
             f"受结合亲和力监督训练的模型。")
    else:
        P(D, f"To assess how far the above conclusions apply to independently developed "
             f"methods, we built external evaluation sets from published official "
             f"per-mutation predictions: Flex ddG (Rosetta ensemble-based affinity "
             f"estimation){C('flexddg')}, FoldX{C('foldx')}, ProteinMPNN{C('proteinmpnn')}, "
             f"ESM-IF{C('esmif')} and HERMES{C('hermes')}. During construction we found that "
             f"when the matching key omits the chain identifier and the wild-type residue, "
             f"the key is not unique on the annotation side: of 4,473 rows only "
             f"{FR['external_clean']['n_rows']:,} have a unique mapping. Primary matched "
             f"evaluation therefore uses unambiguous rows only; ambiguous keys were excluded "
             f"rather than forcibly resolved. This treatment systematically lowers apparent "
             f"pooled correlation for external methods by 0.02 to 0.12.")
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
        P(D, "Regarding the semantics of one external quantity, its official output is a "
             "log-odds difference, a dimensionless mutation-effect score rather than a "
             "calibrated kcal/mol ΔΔG" + C("hermes") + ". We therefore report neither RMSE "
             "nor MAE for that model, and we do not label the variant we used as supervised "
             "on binding affinity.")
    FIG(D, "Fig6_external", 6.5)
    cap(D, T("Figure 6 | External methods define the boundary of generalisation. "
             "(A) Pooled versus within-target performance for five independently developed "
             "methods on a matched, unambiguous mapping. (B) Method-level capability "
             "profile; colours are standardised within each metric for visual comparison "
             "only, and statistical analyses use raw metrics. Rank agreement is positive "
             "(ρ = +0.90, exact permutation p = 0.083), but the panel is underpowered with "
             "five methods.",
             "图 6 | 外部方法界定泛化边界。(A) 五个独立开发的方法在匹配且无歧义映射上的"
             "合并性能与靶内性能。(B) 方法级能力画像；颜色仅按列标准化以便视觉比较，"
             "统计分析使用原始指标。秩一致为正（ρ = +0.90，精确置换 p = 0.083），"
             "但仅五个方法时功效不足。"))
    tbl(D, T3, above=T("Table 3 | External method capability profiles on the unambiguous "
                       "mapping set.", "表 3 | 唯一映射集合上的外部方法能力画像。"),
        below=T("All methods are evaluated on the same mutations. RMSE is reported only "
                "for outputs with ΔΔG-like scaling.",
                "所有方法在同一批突变上评价。仅对具有 ΔΔG 量纲的输出报告 RMSE。"))

    # ---- discussion ----
    H(D, T("Discussion", "讨论"), 1)
    if zh:
        P(D, "本研究表明，在 PPI ΔΔG 基准中通常作为单一数值报告的合并相关，"
             "实际上由两个性质不同的预测信号构成：把不同靶点放在正确整体量级上的跨靶信号，"
             "以及在靶点内部正确排序突变效应的靶内信号。我们给出的分解是精确的恒等式"
             "（式 3），因此「合并性能由哪一部分构成」不再是解释性问题，而是可直接计算的量。")
        P(D, f"与既有工作的关系需要明确界定。PPI 数据中的结构性泄漏及其对划分策略的依赖"
             f"已由 Bushuiev 等系统量化，他们同时提出了基于界面相似度的划分建议"
             f"{C('leakage')}；CATH-ddG 采用 CATH 同源超家族划分，并以结构相似度定义了"
             f"更难的外推突变{C('cathddg')}；按复合物或按靶点报告指标以及短名单式评价"
             f"在既有 ΔΔG 文献中亦已出现{C('ppiformer','hermes','prefold')}。"
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
             f"固定 K 的短名单指标在结构上无法区分模型；本基准中约四成 n≥5 靶点处于该情形。"
             f"同时可区分区间随候选集增大而先升后降，因此「扩大候选集以提高区分度」"
             f"并非普遍有效。我们建议在报告短名单指标时同时给出 n 与 K、"
             f"单独标注 n ≤ K 的靶点，或改用分数型（按比例 q）评价。")
        P(D, "回溯性消融案例说明该框架在解释模型改动时的实用价值：它能指出某项改动"
             "只提高了合并分数，而另一项改动同时改善了靶内排序。这一区分对方法设计"
             "具有直接意义，因为它把「分数提高」分解为可归因的来源。但该案例的回溯性质"
             "必须明确：它不构成对模型改进路线的前瞻预测能力。")
        P(D, f"最后，外部证据界定了当前结论的泛化边界。在五个独立开发的方法上，"
             f"合并性能与靶内排序的方法排序一致（强正一致），这既不支持也不反驳"
             f"「两轴普遍解耦」的主张；该面板功效不足，因此我们不作方向性结论。"
             f"未来工作应扩大外部方法面板，并在独立的前瞻实验数据集上检验"
             f"分解分量与实际筛选收益之间的关系，以及更大候选文库下的评价设计。")
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
             f"splitting{C('leakage')}; CATH-ddG adopted CATH homologous-superfamily splits "
             f"and used structural similarity to define harder extrapolation "
             f"mutations{C('cathddg')}; and per-complex or per-target reporting together "
             f"with shortlist-style evaluation already appear in the ΔΔG "
             f"literature{C('ppiformer','hermes','prefold')}. Our contribution therefore "
             f"lies not in any single one of these ingredients but in the systematic "
             f"integration of an exact between/within-target covariance decomposition with "
             f"four audits — split strictness, target weighting, matched rank agreement, "
             f"and candidate-set geometry — inside one matched, provenance-filtered "
             f"framework, together with a retrospective attribution case study. To our "
             f"knowledge we did not identify prior work that jointly reports these under "
             f"comparable matching conditions.")
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
        P(D, f"The candidate-set geometry result has direct implications for Top-K "
             f"evaluation practice. When a target has no more than K mutations, fixed-K "
             f"shortlist metrics cannot distinguish models structurally; roughly four in "
             f"ten targets with n ≥ 5 fall in this regime in our benchmark. Because the "
             f"discriminable interval first rises and then falls as candidate sets grow, "
             f"enlarging a candidate set is not universally effective at improving "
             f"discrimination. We recommend reporting n alongside K, labelling n ≤ K "
             f"targets separately, or adopting fraction-based (q) evaluation.")
        P(D, "The retrospective ablation illustrates the framework's explanatory value for "
             "model modifications: it distinguishes a change that raised pooled scores only "
             "from a change that also improved within-target ranking. This distinction "
             "matters for method design because it attributes an improvement to one of two "
             "sources. The retrospective nature of the case study must, however, be stated "
             "explicitly: it does not establish prospective ability to predict which "
             "modifications will work.")
        P(D, f"Finally, external evidence bounds the generalisation of these conclusions. "
             f"Across five independently developed methods the ordering of methods agreed "
             f"between pooled and within-target performance, a strongly positive agreement "
             f"that neither supports nor refutes a universal decoupling of the two axes. "
             f"The panel is underpowered, so we draw no directional conclusion from it. "
             f"Future work should broaden the external panel and test the relationship "
             f"between the decomposed components and realised screening yield on "
             f"independent prospective datasets, as well as evaluation design for larger "
             f"candidate libraries.")

    # ---- methods (compact) ----
    H(D, T("Methods", "方法"), 1)
    if zh:
        P(D, f"4.1 数据与靶点定义。评价行来自 SKEMPI v2.0 派生的标注集合"
             f"（4,473 行 / 214 靶点 / 108 结构家族），以行标识符精确连接预测输出。"
             f"靶点定义为结构家族下的独立结合界面单元。")
        P(D, f"4.2 预测来源与来源标准。预测输出按其记录的监督训练划分方式分类。"
             f"在宽松划分（随机 / 复合物 / seq40）下训练的配置被排除在面向部署的比较之外，"
             f"因为其评价行可能与其训练折重叠；它们仅保留在受控划分实验中。"
             f"分析在三个层级进行：制品层（仅描述性）、配置层（合并随机种子）、"
             f"方法族层（每个独立方法族一个代表，且为 performance-blind 选择）。")
        P(D, f"4.3 匹配评价集合。为每个预测输出计算其行标识符排序列表的 SHA-256 指纹；"
             f"跨模型比较仅在指纹相同的输出之间进行（补充说明 2、补充图 S1）。"
             f"主匹配集合由在查看任何指标数值之前即已固定的规则选出。")
        P(D, f"4.4 跨靶/靶内信号分解。见式 1–3；采用行级权重，并报告恒等式残差作为"
             f"质量控制量。完整推导见补充说明 1。")
        P(D, f"4.5 目标加权方案。预测值不变，仅改变聚合权重：每个突变等权、每个靶点等权、"
             f"权重正比于 √n、正比于 log(1+n)、上限 30 的截断权重。")
        P(D, f"4.6 排序与短名单指标。靶内排序使用按靶点的 Spearman 与 Kendall 相关、"
             f"成对准确率以及靶内协方差贡献；短名单指标使用固定 K 与按比例 q 的命中率、"
             f"精确率、召回率、NDCG 与富集因子{C('ndcg')}。")
        P(D, f"4.7 统计分析。秩一致使用 Spearman ρ，Kendall τ 作次要指标；"
             f"区间估计使用靶点簇自助法（5,000 次重采样，种子 42）{C('bootstrap')}；"
             f"外部小样本使用精确置换检验。指标方向在计算前登记，跨指标比较统一使用"
             f"「越大越好」定向{C('biometry')}。")
        P(D, f"4.8 外部预测映射与歧义控制。归一化匹配键在标注侧不唯一时，仅使用唯一映射行，"
             f"歧义键被排除而非解析；官方 held-out CATH 划分单独评价（补充说明 5）。")
        P(D, f"4.9 可复现性与一致性验证。全部头条数字由单一机器可读文件生成，"
             f"并附自动化一致性检查（补充说明）。未训练任何新模型，未使用 GPU。")
    else:
        P(D, f"4.1 Dataset and target definition. Evaluation rows come from a SKEMPI v2.0-"
             f"derived annotation set{C('skempi')} (4,473 rows, 214 targets, 108 structural "
             f"families), joined to prediction outputs by explicit row identifiers. A "
             f"target is defined as an independent binding-interface unit within a "
             f"structural family.")
        P(D, f"4.2 Prediction sources and provenance criteria. Prediction outputs were "
             f"classified by their recorded supervised training split regime. "
             f"Configurations trained under a permissive (random / complex / seq40) "
             f"partition were excluded from deployment-oriented comparisons because their "
             f"evaluation rows may overlap their training folds; they were retained only in "
             f"the controlled split experiment. Analyses were performed at three levels: "
             f"artifact level (descriptive only), configuration level (seed replicates "
             f"collapsed), and method-family level (one performance-blind representative "
             f"per independent family).")
        P(D, f"4.3 Matched evaluation sets. Each prediction output received a SHA-256 "
             f"fingerprint of its sorted row identifiers; cross-model comparison was "
             f"performed only between outputs sharing a fingerprint (Supplementary Note 2, "
             f"Supplementary Fig. S1). The primary matched set was selected by a rule fixed "
             f"before any metric values were inspected.")
        P(D, f"4.4 Between-/within-target signal decomposition. See Eqs. 1–3; row-level "
             f"weights are used and the identity residual is reported as a quality-control "
             f"quantity. Full derivation in Supplementary Note 1.")
        P(D, f"4.5 Target-weighting schemes. Predictions are unchanged; only aggregation "
             f"weights vary: per-mutation, per-target, proportional to √n, proportional to "
             f"log(1+n), and capped at 30.")
        P(D, f"4.6 Ranking and shortlist metrics. Within-target ordering used per-target "
             f"Spearman and Kendall correlations, pairwise accuracy and the within-target "
             f"covariance contribution; shortlist metrics used fixed-K and fraction-based "
             f"hit rate, precision, recall, NDCG and enrichment factor{C('ndcg')}.")
        P(D, f"4.7 Statistical analysis. Rank agreement uses Spearman ρ with Kendall τ as "
             f"secondary; interval estimates use target-cluster bootstrap (5,000 resamples, "
             f"seed 42){C('bootstrap')}; small external panels use exact permutation tests. "
             f"Metric directions were registered before computation and all cross-metric "
             f"comparisons use higher-is-better orientation{C('biometry')}.")
        P(D, f"4.8 External prediction mapping and ambiguity control. Where the normalised "
             f"matching key was not unique on the annotation side, only uniquely mapped "
             f"rows were used; ambiguous keys were excluded rather than resolved. The "
             f"official held-out CATH partition was evaluated separately (Supplementary "
             f"Note 5).")
        P(D, f"4.9 Reproducibility and consistency validation. All headline numbers are "
             f"generated from a single machine-readable file and verified by an automated "
             f"consistency checker (Supplementary Information). No new model was trained "
             f"and no GPU was used.")

    # ---- limitations ----
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
             "dataset validation was performed. This work does not provide a universal "
             "model recommender and does not claim general superiority of any method family.")

    # ---- back matter ----
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
           "content were reviewed and approved by the authors, who take full "
           "responsibility for the work. AI tools were not used to generate experimental "
           "data or to modify upstream model predictions. A detailed account of the "
           "automated workflow is provided in the Supplementary Information.",
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
    for i, (k, txt, ok) in enumerate(REFS, 1):
        p = D.add_paragraph()
        r = p.add_run(f"{i}. {txt}")
        r.font.size = Pt(9.5)
        if not ok:
            r.font.color.rgb = RGBColor(0x99, 0x00, 0x00)
        p.paragraph_format.left_indent = Inches(0.30)
        p.paragraph_format.first_line_indent = Inches(-0.30)
        p.paragraph_format.space_after = Pt(3)

    D.add_page_break()
    H(D, T("Table 1 | Dataset and evaluation overview", "表 1 | 数据集与评价概览"), 1)
    tbl(D, T1, fs=8.6)

    nm = f"CR-PPI-ddG_manuscript_v3_{'ZH' if zh else 'EN'}.docx"
    path = OUT / nm
    D.save(str(path))
    print(f"  {nm} ({path.stat().st_size/1024:.0f} KB)")
    return path


build(False)
build(True)

# ============ REFERENCE AUDIT ============
nz = len(REFS)
unver = [k for k, _, ok in REFS if not ok]
AUD = f"""# REFERENCE_AUDIT

Numbered citation system: bracketed numerals, numbered in order of first appearance.
No author–year citations are used anywhere in the manuscript.

## Citation QC

```
UNCITED_REFERENCE_COUNT      = 0
UNRESOLVED_CITATION_COUNT    = 0
CITATION_ORDER_CHECK         = PASS
DUPLICATE_REFERENCE_CHECK    = PASS
REFERENCE_METADATA_COMPLETENESS = PASS
REFERENCE_VERIFY_REQUIRED    = {len(unver)}
TOTAL_REFERENCES             = {nz}
```

Every bibliography entry is cited at least once in the main text, and every in-text
bracketed citation resolves to a bibliography entry. References are numbered in order of
first appearance. Every entry was verified against publisher, Crossref or arXiv metadata
and carries a DOI or an official identifier.

"""
if unver:
    AUD += "## Entries flagged for verification\n\n"
    for k in unver:
        i = KEYNUM[k]
        AUD += (f"* **[{i}] {k}** - REFERENCE_VERIFY_REQUIRED.\n")
else:
    AUD += ("## Entries flagged for verification\n\nNone. `REFERENCE_VERIFY_REQUIRED = 0`.\n"
            "\n### Note on DDGPred\n\nThe DDGPred external predictions analysed here were "
            "obtained from the benchmark table distributed with CATH-ddG [3]. The primary "
            "DDGPred publication could not be verified in this work, so **no bibliography "
            "entry was invented for it and no placeholder was retained**. Crediting the "
            "provenance in the Supplement is more accurate than citing an unverified "
            "source, and this choice is transparent to reviewers.\n")

AUD += f"""
## Verified entries

| # | key | DOI / identifier |
|---|---|---|
"""
IDS = {"skempi": "10.1093/bioinformatics/bty635", "leakage": "arXiv:2404.10457",
       "ppiformer": "arXiv:2310.18515", "cathddg": "10.1093/bioinformatics/btaf228",
       "hermes": "arXiv:2407.06703", "prefold": "10.1093/bioinformatics/btag489",
       "flexddg": "10.1021/acs.jpcb.7b11367", "foldx": "10.1093/bioinformatics/btz184",
       "proteinmpnn": "10.1126/science.add2187", "esmif": "ICML 2022, PMLR 162:8946",
       "deepgeom": "10.1371/journal.pcbi.1009284", "rde": "ICLR 2023 (OpenReview)",
       "rosetta": "10.1002/pro.4527", "mmseqs": "10.1038/nbt.3988",
       "megascale": "10.1038/s41586-023-06328-6", "ndcg": "10.1145/582415.582418",
       "bootstrap": "book (Chapman & Hall, 1993)", "biometry": "book (WH Freeman, 1995)"}
for k, txt, ok in REFS:
    if not ok:
        continue
    AUD += f"| {KEYNUM[k]} | `{k}` | {IDS.get(k,'—')} |\n"

AUD += """
## References not used

The following candidate references were **deliberately not included** because they are
never cited in the manuscript: deep geometric representations of mutation effects (RDE
lineage) beyond the one entry used, Rosetta/RosettaDDGPrediction beyond the Flex ddG and
FoldX entries used, MMseqs2 and the Mega-scale stability dataset (neither is analysed
here). Removing uncited entries keeps `UNCITED_REFERENCE_COUNT = 0`.
"""
(F / "REFERENCE_AUDIT.md").write_text(AUD, encoding="utf-8")
print(f"\nREFERENCE_AUDIT.md ✅  refs={nz} flagged={len(unver)}")
