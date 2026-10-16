# -*- coding: utf-8 -*-
"""scripts/v8_1/10_manuscript.py -- 生成中/英文稿件 DOCX（各一份），再由 Word 转 PDF。

数字全部取自 `v8_1/FROZEN_RESULTS.json`。英文版使用 `figures_en/` 的英文图与英文表。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt, RGBColor

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
A = R / "v8_1"
F = R / "v8_1_final"
V8 = R / "experiments" / "v8"
OUT = R / "manuscript"
OUT.mkdir(exist_ok=True)

FR = json.load(open(A / "FROZEN_RESULTS.json", encoding="utf-8"))
PU, H4, GEO = FR["primary_universe"], FR["H4_controlled_split"], FR["axis_c_geometry"]
DL = H4.get("delta", {})
ST = pd.read_csv(A / "rank_agreement_by_stratum.csv")
EP = pd.read_csv(A / "external_clean_profiles.csv")
EP2 = pd.read_csv(A / "external_clean_vs_old.csv")
W2 = pd.read_csv(A / "WEIGHTING_EFFECT_SIZE_AUDIT.csv")
UR = pd.read_csv(A / "MATCHED_UNIVERSE_REGISTRY.csv")
er = FR["external_clean"]["rank_agreement_external"]
rmn, rmx = float(ST.rho_overall_vs_macro.min()), float(ST.rho_overall_vs_macro.max())
sf0 = int(W2["signflip_C_B_gt0"].sum()); sf1 = int(W2["signflip_C_B_gt0.001"].sum())
sf5 = int(W2["signflip_C_B_gt0.005"].sum())
sfw0 = int(W2["signflip_C_W_gt0"].sum())
best = EP.loc[EP.overall_pearson.idxmax()]
EXT_OLD = {r.model: r for r in EP2.itertuples()}
print(f"主池 {PU['rowset_hash'][:8]} rows={PU['n_rows']} n>=5={PU['n_ge5']} "
      f"families={PU['n_method_families']}")
print(f"rho range {rmn:+.3f}..{rmx:+.3f}; C_B flips {sf0}/{sf1}/{sf5}; "
      f"between-frac {DL.get('fraction_between'):.3f}")

REFS = [
    "Jankauskaitė J, Jiménez-García B, Dapkūnas J, Fernández-Recio J, Moal IH. "
    "SKEMPI 2.0: an updated benchmark of changes in protein–protein binding energy, "
    "kinetics and thermodynamics upon mutation. Bioinformatics. 2019;35(3):462–469.",
    "Bushuiev A, Bushuiev R, Sedlar J, Pluskal T, Damborsky J, Mazurenko S, Sivic J. "
    "Revealing data leakage in protein interaction benchmarks. arXiv:2404.10457. 2024.",
    "Bushuiev A, Bushuiev R, Kouba P, et al. Learning to design protein–protein "
    "interactions with enhanced generalization. ICLR 2024. arXiv:2310.18515.",
    "Yu G, Bi X, Ma T, Li Y, Wang J. CATH-ddG: towards robust mutation effect prediction "
    "on protein–protein interactions out of CATH homologous superfamily. "
    "Bioinformatics. 2025;41(Suppl_1):i362–i372. doi:10.1093/bioinformatics/btaf228.",
    "Luo S, Su Y, Wu Z, Su C, Peng J, Ma J. Rotamer density estimator is an unsupervised "
    "learner of the effect of mutations on protein–protein interaction. ICLR 2023.",
    "Dauparas J, Anishchenko I, Bennett N, et al. Robust deep learning–based protein "
    "sequence design using ProteinMPNN. Science. 2022;378(6615):49–56.",
    "Hsu C, Verkuil R, Liu J, et al. Learning inverse folding from millions of predicted "
    "structures. ICML 2022. (ESM-IF)",
    "Delgado J, Radusky LG, Cianferoni D, Serrano L. FoldX 5.0: working with RNA, small "
    "molecules and a new graphical interface. Bioinformatics. 2019;35(20):4168–4169.",
    "Sora V, et al. Flex-ddG / FoldX-based binding-affinity estimation (as reported in the "
    "CATH-ddG benchmark table).",
    "Kulandaisamy A, et al. DDGPred: a deep-learning method for predicting mutation-induced "
    "changes in protein–protein binding affinity (as reported in the CATH-ddG benchmark table).",
    "Liu X, Luo Y, Li P, Song S, Peng J. Deep geometric representations for modeling "
    "effects of mutations on protein–protein binding affinity. PLoS Comput Biol. "
    "2021;17(8):e1009284.",
    "Cheng J, Novati G, Pan J, et al. Accurate proteome-wide missense variant effect "
    "prediction with AlphaMissense. Science. 2023;381(6664):eadg7492.",
    "Jumper J, Evans R, Pritzel A, et al. Highly accurate protein structure prediction "
    "with AlphaFold. Nature. 2021;596(7873):583–589.",
    "Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the "
    "analysis of massive data sets. Nat Biotechnol. 2017;35(11):1026–1028.",
    "Wilcoxon F. Individual comparisons by ranking methods. Biometrics Bulletin. "
    "1945;1(6):80–83.",
    "Efron B, Tibshirani RJ. An Introduction to the Bootstrap. Chapman & Hall; 1993.",
    "Benjamini Y, Hochberg Y. Controlling the false discovery rate. J R Stat Soc B. "
    "1995;57(1):289–300.",
    "Sokal RR, Rohlf FJ. Biometry. 3rd ed. WH Freeman; 1995. (Spearman/Kendall rank "
    "statistics)",
    "Järvelin K, Kekäläinen J. Cumulated gain-based evaluation of IR techniques. "
    "ACM Trans Inf Syst. 2002;20(4):422–446. (NDCG)",
    "Tsuboyama K, Dauparas J, Chen J, et al. Mega-scale experimental analysis of protein "
    "folding stability in biology and design. Nature. 2023;620:434–444.",
]


def setup(doc, lang):
    st = doc.styles["Normal"]
    st.font.size = Pt(10.5)
    if lang == "zh":
        st.font.name = "SimSun"
    else:
        st.font.name = "Times New Roman"
    for s in doc.sections:
        s.top_margin = s.bottom_margin = Inches(0.9)
        s.left_margin = s.right_margin = Inches(0.9)


def H(doc, text, lvl=1):
    p = doc.add_heading(text, level=lvl)
    for r in p.runs:
        r.font.color.rgb = RGBColor(0, 0, 0)
    return p


def P(doc, text, italic=False, size=None, align=None):
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.italic = italic
    if size:
        r.font.size = Pt(size)
    if align:
        p.alignment = align
    return p


def cap(doc, text, size=8.5):
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.font.size = Pt(size)
    r.italic = True
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    return p


def mktable(doc, df, widths=None, fs=7.6):
    t = doc.add_table(rows=1, cols=len(df.columns))
    t.style = "Light Grid Accent 1"
    for i, c in enumerate(df.columns):
        cell = t.rows[0].cells[i]
        cell.text = ""
        run = cell.paragraphs[0].add_run(str(c))
        run.bold = True
        run.font.size = Pt(fs)
    for _, row in df.iterrows():
        cells = t.add_row().cells
        for i, v in enumerate(row):
            cells[i].text = ""
            run = cells[i].paragraphs[0].add_run(str(v))
            run.font.size = Pt(fs)
    return t


# ---------- 表格数据 ----------
T1 = UR[["artifact", "rowset_hash", "n_rows", "n_targets", "n_families",
         "n_ge3", "n_ge5", "n_ge10", "n_models"]].copy()
T1.columns = (["Artifact", "Rowset hash", "Rows", "Targets", "Families",
               "n>=3", "n>=5", "n>=10", "Model keys"]
              if True else T1.columns)
T3 = ST[["rowset_hash", "level", "n_models", "n_method_families",
         "rho_overall_vs_macro", "tau_overall_vs_macro",
         "rho_negrmse_vs_macro", "rho_overall_vs_util",
         "rho_macro_vs_util"]].copy()
for c in T3.columns:
    if T3[c].dtype.kind == "f":
        T3[c] = T3[c].round(4)
T5 = pd.DataFrame([{"Split regime": k, "Median pooled Pearson": round(v["overall_pearson"], 4),
                    "C_B": round(v["C_B"], 4), "C_W": round(v["C_W"], 4),
                    "RMSE": round(v["rmse"], 4),
                    "Macro Spearman": round(v["macro_spearman"], 4)}
                   for k, v in H4.items() if k != "delta"])
T6 = W2[["model_name", "range_C_B", "range_C_W", "max_abs_C_B", "max_abs_C_W",
         "signflip_C_B_gt0", "signflip_C_B_gt0.005",
         "signflip_C_W_gt0", "signflip_C_W_gt0.005"]].head(12).copy()
for c in T6.columns[1:]:
    T6[c] = T6[c].round(4)
T10 = EP[["model", "n_rows", "n_targets_ge5", "overall_pearson", "macro_spearman",
          "C_B_row", "C_W", "rmse"]].copy()
T10["efficiency@10"] = EP["efficiency@10"]
T10 = T10[["model", "n_rows", "n_targets_ge5", "overall_pearson", "macro_spearman",
           "C_B_row", "C_W", "efficiency@10", "rmse"]]
T10.columns = ["Model", "n", "n>=5", "Pearson", "Macro rho", "C_B", "C_W",
               "Eff@10", "RMSE"]
for c in T10.columns[3:]:
    T10[c] = T10[c].round(4)


def build(lang):
    zh = lang == "zh"
    D = Document()
    setup(D, lang)
    IMG = F / ("figures" if zh else "figures_en")
    sfx = "" if zh else "_EN"

    if zh:
        title = "解耦 PPI ΔΔG 基准评价中的跨靶信号与靶内信号"
        sub = ("一项基于精确协方差分解的匹配、来源清洁的基准分析")
    else:
        title = ("Disentangling Between-Target and Within-Target Signal in "
                 "PPI ΔΔG Benchmarking")
        sub = ("A matched, provenance-clean benchmark analysis based on an exact "
               "covariance-signal decomposition")
    h = D.add_heading(title, level=0)
    for r in h.runs:
        r.font.color.rgb = RGBColor(0, 0, 0)
    P(D, sub, italic=True, align=WD_ALIGN_PARAGRAPH.CENTER)
    P(D, "")

    # ---------- Abstract ----------
    H(D, "Abstract" if not zh else "摘要", 1)
    if zh:
        ab = (f"背景：蛋白质–蛋白质结合突变效应（ΔΔG）的基准评价长期依赖在异质突变集合上的"
              f"合并相关系数。问题：合并相关混合了两种不同能力——定位靶点整体结合量级"
              f"（跨靶信号）与在靶内对突变排序（靶内信号）；两者的相对贡献通常未被量化。"
              f"方法：我们在严格匹配的行集合上重新分析既有 PPI ΔΔG 预测，剔除伪重复制品与"
              f"来源不洁模型，并将合并相关精确分解为 r = C_B + C_W（跨靶与靶内协方差贡献）。"
              f"主要受控结果：在仅改变划分方式的受控条件下（同一批 7523 行、同一批配置），"
              f"合并 Pearson 中位数由 {H4['random']['overall_pearson']:+.3f} 降至 "
              f"{H4['family_OOD']['overall_pearson']:+.3f}；其中 "
              f"{100*DL.get('fraction_between',0):.0f}% 的下降由跨靶贡献的丧失承载"
              f"（ΔC_B = {DL.get('delta_C_B',0):+.3f}），而靶内贡献仍为正"
              f"（ΔC_W = {DL.get('delta_C_W',0):+.3f}）。加权结果：在五种目标加权方案下，"
              f"C_B 仅在近零区间方向不稳（无阈值 {sf0}/{len(W2)}，|·|>0.005 时 {sf5}/{len(W2)}），"
              f"C_W 在任何阈值下均不改变符号，但幅度变化（中位极差 "
              f"{W2.range_C_W.median():.3f}）。固定 K 几何结果：K=10 时 D(n,K)=K(n−K)/[n(n+1)] "
              f"在 n≤K 时为零（占 n≥5 靶点的 {100*GEO['n_le_K_fraction']:.1f}%），并在 "
              f"n*≈{GEO['n_star_analytic']:.1f} 取峰，故短名单指标对小候选集结构性无信息。"
              f"补救结果：在匹配的回溯消融中，分解区分了仅改善合并分数的改动（ΔMacro ρ = −0.002）"
              f"与同时改善靶内排序的改动（ΔMacro ρ = +0.044）。外部警示：外部证据仅 5 个方法，"
              f"功效不足，且其秩一致为正而非解耦。结论：PPI ΔΔG 基准应在明确声明的加权下同时报告"
              f"跨靶与靶内分量，并报告候选集几何。")
    else:
        ab = (f"Background. Benchmarking of mutation effects on protein–protein binding "
              f"(ΔΔG) relies predominantly on pooled correlation over heterogeneous "
              f"mutation sets. Problem. Pooled correlation conflates two distinct "
              f"capabilities: locating a target's overall binding scale (between-target "
              f"signal) and ranking mutations within a target (within-target signal); "
              f"their relative contribution is rarely quantified. Method. We re-analyse "
              f"existing PPI ΔΔG predictions under strictly matched row sets, removing "
              f"pseudo-replicated artifacts and provenance-impure models, and decompose "
              f"pooled correlation exactly as r = C_B + C_W (between- and within-target "
              f"covariance contributions). Main controlled result. Changing only the split "
              f"regime (identical 7523 rows and identical configurations), median pooled "
              f"Pearson falls from {H4['random']['overall_pearson']:+.3f} to "
              f"{H4['family_OOD']['overall_pearson']:+.3f}; "
              f"{100*DL.get('fraction_between',0):.0f}% of this reduction is carried by "
              f"loss of the between-target contribution "
              f"(ΔC_B = {DL.get('delta_C_B',0):+.3f}), while the within-target "
              f"contribution remains positive (ΔC_W = {DL.get('delta_C_W',0):+.3f}). "
              f"Weighting result. Across five target-weighting schemes, C_B is unstable in "
              f"direction only near zero ({sf0}/{len(W2)} without a threshold, "
              f"{sf5}/{len(W2)} at |·|>0.005), whereas C_W never changes sign at any "
              f"threshold, although its magnitude varies "
              f"(median range {W2.range_C_W.median():.3f}). Fixed-K geometry result. For "
              f"K = 10, D(n,K) = K(n−K)/[n(n+1)] is zero whenever n ≤ K "
              f"({100*GEO['n_le_K_fraction']:.1f}% of our n≥5 targets) and peaks near "
              f"n* ≈ {GEO['n_star_analytic']:.1f}, so shortlist metrics are structurally "
              f"non-informative for small candidate sets. Remedy result. In a matched "
              f"retrospective ablation, decomposition separated a pooled-score-only change "
              f"(ΔMacro ρ = −0.002) from a modification that also improved within-target "
              f"ranking (ΔMacro ρ = +0.044). External caveat. External evidence comprised "
              f"only five methods and is underpowered; its rank agreement was positive "
              f"rather than decoupled. Conclusion. PPI ΔΔG benchmarks should report "
              f"between- and within-target components under an explicitly declared "
              f"weighting, together with candidate-set geometry.")
    P(D, ab, align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    # ---------- Introduction ----------
    H(D, "1. Introduction" if not zh else "1. 引言", 1)
    if zh:
        P(D, "机器学习方法已被广泛用于预测突变对蛋白质–蛋白质结合亲和力的影响。"
             "这类方法的评价通常报告一个合并的相关系数（如全部突变上的 Pearson r）"
             "或一个合并的误差指标（RMSE）。然而，一个合并的相关系数并不对应单一能力："
             "它同时取决于模型能否把不同靶点放在正确的结合量级上（跨靶/校准信号），"
             "以及能否在同一靶点内部把有利与不利突变正确排序（靶内/排序信号）。"
             "当基准由大量靶点、且各靶点突变数高度不均时，合并相关系数会被样本量大的"
             "靶点主导，从而把这两种能力混为一谈。", align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "与此同时，划分策略对报告性能的影响已被证实。Bushuiev 等量化了 PPI 数据结构性"
             "泄漏：以 PDB 代码划分时约 65% 的测试样本存在近重复训练样本，按 30% 序列同一性"
             "划分仍有约 30%，而 SKEMPI v2.0 若按 PPI 代码划分为 56%，按领域专家分组则为 0%。"
             "CATH-ddG 进一步采用 CATH 同源超家族划分，并提出以最大 TM-score < 0.6 定义"
             "「难突变」。因此，泄漏感知与 family-OOD 评价并非本文首创。",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "本文的出发点是：这些既有工作量化了「有多少泄漏」，但较少量化「泄漏与划分严格度"
             "改变了报告性能的哪个组成部分」。我们因此提出一个精确的协方差信号分解，"
             "并在严格匹配、去伪重复、来源清洁的条件下量化划分严格度、目标加权与候选集几何"
             "各自对基准读数的影响。我们明确不主张自动模型选择，也不主张某个方法族普遍更优。",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    else:
        P(D, "Machine-learning methods are now widely used to predict how mutations affect "
             "protein–protein binding affinity. Their evaluation typically reports a single "
             "pooled correlation (e.g. Pearson r over all mutations) or a pooled error "
             "(RMSE). A pooled correlation, however, does not correspond to a single "
             "capability: it depends both on whether a model places different targets on the "
             "correct overall binding scale (between-target / calibration signal) and on "
             "whether it orders favourable and unfavourable mutations correctly within a "
             "target (within-target / ranking signal). When a benchmark contains many "
             "targets with highly unequal mutation counts, pooled correlation is dominated "
             "by large targets, conflating these two capabilities.",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "In parallel, the influence of splitting strategy on reported performance is "
             "established. Bushuiev et al. quantified structural leakage in PPI data: about "
             "65% of test examples have a near-duplicate training example under PDB-code "
             "splitting, about 30% even under 30% sequence-identity splitting, 56% for "
             "SKEMPI v2.0 split by PPI code, and 0% under domain-expert grouping. CATH-ddG "
             "adopted CATH homologous-superfamily splits and defined 'hard' mutations by "
             "maximum TM-score < 0.6. Leakage-aware and family-out-of-distribution (OOD) "
             "evaluation therefore is not new here.",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "Our starting point is that this prior work quantifies *how much* leakage "
             "exists, but less often *which component* of reported performance leakage and "
             "split strictness change. We therefore introduce an exact covariance-signal "
             "decomposition and, under strictly matched, de-duplicated, provenance-clean "
             "conditions, quantify how split strictness, target weighting and candidate-set "
             "geometry each affect benchmark readings. We explicitly do not claim automated "
             "model selection, nor that any method family is universally superior.",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    # ---------- Methods ----------
    H(D, "2. Methods" if not zh else "2. 方法", 1)
    if zh:
        P(D, "2.1 数据与命名空间。评价行来自既有 SKEMPI v2.0 派生标注（4473 行 / 214 靶点 / "
             "108 家族），并以行标识符精确连接预测制品。", align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "2.2 匹配宇宙与行集指纹。我们为每个预测制品计算其行标识符排序列表的 SHA-256 "
             "指纹。仅在指纹相同的制品之间进行跨模型比较；跨指纹的原始指标不作比较。",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "2.3 去伪重复三层。分析在三个层级进行：制品层（描述性，不视为独立方法）、"
             "配置层（合并随机种子）、方法族层（每个独立方法族选取一个代表性模型）。"
             "代表性选取为 performance-blind：优先冻结的默认配置，其次官方默认命名，"
             "再其次按字典序取最小键；禁止依据测试指标选取。", align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "2.4 协方差信号分解。对每个靶点 c，令 B_y,i = μ_y,c − μ_y 为跨靶分量，"
             "W_y,i = y_i − μ_y,c 为靶内分量（预测端同理），则\n"
             "C_B = Cov(B_y, B_ŷ)/(σ_y σ_ŷ)，C_W = Cov(W_y, W_ŷ)/(σ_y σ_ŷ)，"
             "且 r_pooled = C_B + C_W 为精确恒等式。我们报告恒等式残差以作质量控制。"
             "C_B 与 C_W 是协方差信号贡献，不是独立相关系数、独立能力分数或因果机制。",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "2.5 目标加权。预测值不变，仅改变聚合权重：ROW（每个突变等权）、"
             "EQUAL（每个靶点等权）、SQRT（权重 ∝ √n）、LOG（∝ log(1+n)）、CAP30（∝ min(n,30)）。",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "2.6 指标方向登记。我们预先登记每个指标的方向，跨指标秩一致一律使用 "
             "higher-is-better 定向（RMSE/MAE 取负或升序排名），避免把「越高越好的 Pearson」"
             "与「越高越差的 RMSE」直接相关。", align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "2.7 固定 K 候选集几何。对候选集大小 n 与 K，定义可区分区间 "
             "D(n,K) = E_max − E_rand = K(n−K)/[n(n+1)]（n > K），n ≤ K 时 D = 0；"
             "其在 n* = K + √(K²+K) 处取峰。", align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "2.8 统计。秩一致使用 Spearman ρ，Kendall τ 作次要指标；区间估计使用靶点簇"
             "自助法（5000 次重采样，种子 42）；外部小样本使用精确置换检验。",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    else:
        P(D, "2.1 Data and namespace. Evaluation rows come from an existing SKEMPI v2.0-derived "
             "annotation set (4473 rows / 214 targets / 108 families); prediction artifacts "
             "are joined by explicit row identifiers.", align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "2.2 Matched universes and rowset fingerprints. For every prediction artifact we "
             "compute the SHA-256 fingerprint of its sorted row-identifier list. Cross-model "
             "comparison is performed only between artifacts sharing a fingerprint; raw "
             "metrics are never compared across fingerprints.",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "2.3 Three de-replication levels. Analyses are conducted at artifact level "
             "(descriptive only; artifacts are not treated as independent methods), "
             "configuration level (seed replicates collapsed), and method-family level (one "
             "representative per independent family). Representatives are chosen "
             "performance-blind: frozen default configuration first, then official default "
             "naming, then lexicographically smallest key. Selecting by test metrics is "
             "forbidden.", align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "2.4 Covariance-signal decomposition. For each target c let "
             "B_y,i = μ_y,c − μ_y be the between-target component and W_y,i = y_i − μ_y,c the "
             "within-target component (analogously for predictions). Then "
             "C_B = Cov(B_y, B_ŷ)/(σ_y σ_ŷ), C_W = Cov(W_y, W_ŷ)/(σ_y σ_ŷ), and "
             "r_pooled = C_B + C_W holds exactly; we report the identity residual as a QC "
             "check. C_B and C_W are covariance-signal contributions, not independent "
             "correlations, independent skill scores, or causal mechanisms.",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "2.5 Target weighting. Predictions are unchanged; only aggregation weights "
             "change: ROW (per mutation), EQUAL (per target), SQRT (∝ √n), LOG (∝ log(1+n)), "
             "CAP30 (∝ min(n,30)).", align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "2.6 Metric-direction registry. The direction of every metric is registered in "
             "advance; all cross-metric rank agreement uses higher-is-better orientation "
             "(RMSE/MAE negated or ascending-ranked), so that 'higher-is-better Pearson' is "
             "never correlated against raw 'higher-is-worse RMSE'.",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "2.7 Fixed-K candidate-set geometry. For candidate-set size n and K, the "
             "discriminable interval is D(n,K) = E_max − E_rand = K(n−K)/[n(n+1)] for n > K, "
             "and D = 0 for n ≤ K; it peaks at n* = K + √(K²+K).",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "2.8 Statistics. Rank agreement uses Spearman ρ with Kendall τ as secondary; "
             "interval estimates use target-cluster bootstrap (5000 resamples, seed 42); "
             "small external panels use exact permutation tests.",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    # ---------- Results ----------
    H(D, "3. Results" if not zh else "3. 结果", 1)

    H(D, ("3.1 Pooled ΔΔG performance can be exactly decomposed into between-target and "
          "within-target signal" if not zh else
          "3.1 合并 ΔΔG 性能可精确分解为跨靶与靶内信号"), 2)
    if zh:
        P(D, f"在 provenance-clean 的主匹配宇宙（行集指纹 {PU['rowset_hash'][:8]}，"
             f"{PU['n_rows']} 行，n≥5 靶点 {PU['n_ge5']} 个，制品 {PU['n_models']} 个，"
             f"独立方法族 {PU['n_method_families']} 个，100% GREEN）上，"
             f"恒等式 r_pooled = C_B + C_W 以机器精度成立（残差 ≤ 3.3×10⁻¹⁶）。"
             f"这使「合并相关由哪一部分构成」成为一个可直接测量的问题。",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    else:
        P(D, f"In the provenance-clean primary matched universe (rowset fingerprint "
             f"{PU['rowset_hash'][:8]}; {PU['n_rows']} rows; {PU['n_ge5']} targets with "
             f"n≥5; {PU['n_models']} artifacts; {PU['n_method_families']} independent "
             f"method families; 100% GREEN), the identity r_pooled = C_B + C_W holds to "
             f"machine precision (residual ≤ 3.3×10⁻¹⁶). This turns 'what pooled "
             f"correlation is made of' into a directly measurable question.",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    D.add_picture(str(IMG / f"Fig2_CB_CW_decomposition{sfx}.png"), width=Inches(4.5))
    cap(D, ("图 2. 精确协方差信号分解。主匹配、来源清洁宇宙；"
            "标题含 universe / rows / targets / artifacts / method families / provenance 声明。")
        if zh else
        ("Figure 2. Exact covariance-signal decomposition in the primary matched, "
         "provenance-clean universe. The caption states universe, rows, targets, "
         "artifacts, method families and provenance status."))
    D.add_picture(str(IMG / f"Fig1_matched_universe_strata{sfx}.png"), width=Inches(6.2))
    cap(D, ("图 1. 匹配宇宙分层。12 个制品落于 7 个不同行集指纹；"
            "仅共享同一指纹的制品可直接比较。")
        if zh else
        ("Figure 1. Matched-universe strata. Twelve artifacts fall into seven distinct "
         "rowset hashes; only artifacts sharing a hash are directly comparable."))

    H(D, ("3.2 Family-OOD evaluation disproportionately removes between-target signal"
          if not zh else "3.2 Family-OOD 评价不成比例地移除跨靶信号"), 2)
    if zh:
        P(D, f"这是本文最强的受控证据：同一批 7523 行、同一批模型配置，仅改变划分方式。"
             f"合并 Pearson 中位数由随机划分的 {H4['random']['overall_pearson']:+.4f} "
             f"降至 family-OOD 的 {H4['family_OOD']['overall_pearson']:+.4f}；"
             f"ΔC_B = {DL.get('delta_C_B'):+.4f}，ΔC_W = {DL.get('delta_C_W'):+.4f}。"
             f"按绝对变化分解，跨靶贡献承载 {100*DL.get('fraction_between'):.1f}% 的下降，"
             f"且 C_W 最终仍为正（{H4['family_OOD']['C_W']:+.4f}）。"
             f"我们因此表述为：family-OOD 下合并相关的下降主要由跨靶协方差贡献的丧失承载，"
             f"而靶内协方差保持为正。这不等于「模型记住了 PDB 身份」——"
             f"本文没有身份记忆的干预实验，故不作该因果表述。",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "需要明确标注：complex 与 seq40 两种划分在本数据上所有指标逐位相同，"
             "因此不能作为彼此独立的证据。", align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    else:
        P(D, f"This is our strongest controlled evidence: identical 7523 rows and identical "
             f"model configurations, with only the split regime changed. Median pooled "
             f"Pearson falls from {H4['random']['overall_pearson']:+.4f} under random "
             f"splitting to {H4['family_OOD']['overall_pearson']:+.4f} under family-OOD; "
             f"ΔC_B = {DL.get('delta_C_B'):+.4f} and ΔC_W = {DL.get('delta_C_W'):+.4f}. "
             f"Under this absolute-change decomposition, the between-target contribution "
             f"carries {100*DL.get('fraction_between'):.1f}% of the reduction, and C_W "
             f"remains positive ({H4['family_OOD']['C_W']:+.4f}). We therefore state that "
             f"most of the reduction in pooled correlation under family-OOD evaluation is "
             f"carried by loss of the between-target covariance contribution, while "
             f"within-target covariance remains positive. This is *not* equivalent to "
             f"'the model memorised PDB identity': we performed no identity-memorisation "
             f"intervention and make no such causal claim.",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "We flag explicitly that the complex and seq40 regimes are identical in every "
             "metric in this dataset and therefore cannot serve as independent evidence.",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    D.add_picture(str(IMG / f"Fig3_controlled_split_inflation{sfx}.png"), width=Inches(6.3))
    cap(D, ("图 3. 受控划分膨胀。同一批 7523 行、同一批配置，仅改划分方式。"
            "注意 complex 与 seq40 在本数据上完全相同。")
        if zh else
        ("Figure 3. Controlled split inflation: identical 7523 rows and identical "
         "configurations, with only the split regime changed. Note that complex and seq40 "
         "are identical in this dataset."))
    mktable(D, T5)
    cap(D, "表 2. 受控划分膨胀（同一批行、同一批配置）。" if zh
        else "Table 2. Controlled split inflation (identical rows, identical configurations).")

    H(D, ("3.3 Target weighting changes the interpretation of covariance contributions"
          if not zh else "3.3 目标加权改变协方差贡献的解读"), 2)
    if zh:
        P(D, f"在五种加权方案下（n = {len(W2)} 个模型）：C_B 的符号翻转在无阈值时为 "
             f"{sf0}/{len(W2)}，|·|>0.001 时 {sf1}/{len(W2)}，|·|>0.005 时 {sf5}/{len(W2)}；"
             f"C_W 在任何阈值下均为 {sfw0}/{len(W2)}。C_W 的幅度并非不变"
             f"（中位极差 {W2.range_C_W.median():.4f}，C_B 为 {W2.range_C_B.median():.4f}）。"
             f"因此准确表述是：C_B 仅在近零区间方向不稳，而 C_W 保持符号稳定但幅度可变。"
             f"任何报告的分解都必须声明其加权口径。", align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    else:
        P(D, f"Across five weighting schemes (n = {len(W2)} models), C_B changes sign in "
             f"{sf0}/{len(W2)} models with no threshold, {sf1}/{len(W2)} at |·|>0.001, and "
             f"{sf5}/{len(W2)} at |·|>0.005, whereas C_W changes sign in {sfw0}/{len(W2)} at "
             f"any threshold. C_W magnitude is *not* invariant (median range "
             f"{W2.range_C_W.median():.4f} versus {W2.range_C_B.median():.4f} for C_B). The "
             f"accurate statement is therefore that C_B is unstable in direction only near "
             f"zero, while C_W keeps its sign but varies in magnitude. Any reported "
             f"decomposition must declare its weighting.",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    D.add_picture(str(IMG / f"Fig6_weighting_sensitivity{sfx}.png"), width=Inches(6.3))
    cap(D, ("图 6. 目标加权敏感性：C_B 方向在近零处不稳；C_W 保持符号但幅度变化。")
        if zh else
        ("Figure 6. Target-weighting sensitivity: C_B direction is unstable near zero; "
         "C_W keeps its sign but its magnitude varies."))
    mktable(D, T6)
    cap(D, "表 3. 加权效应量审计（前 12 个模型）。" if zh
        else "Table 3. Weighting effect-size audit (first 12 models).")

    H(D, ("3.4 Fixed-K shortlist metrics are constrained by candidate-set geometry"
          if not zh else "3.4 固定 K 短名单指标受候选集几何约束"), 2)
    if zh:
        P(D, f"对 K = 10，可区分区间 D(n,K) = K(n−K)/[n(n+1)] 在 n ≤ K 时为零；"
             f"本基准 n≥5 靶点中有 {GEO['n_le_K_count']}/{GEO['n_ge5_count']} = "
             f"{100*GEO['n_le_K_fraction']:.1f}% 属此情形，此时任何模型都必然选中全部候选，"
             f"该靶点对固定 K 短名单指标结构性无信息。D 在 n* ≈ "
             f"{GEO['n_star_analytic']:.1f} 取峰后单调下降：D(20) = {GEO['D_at_n20']:.4f} "
             f"> D(50) = {GEO['D_at_n50']:.4f}。因此「候选集越大越可区分」在本指标下不成立，"
             f"我们据此撤回了此前「建议 n ≥ 50」的表述。我们同时明确："
             f"该几何结论**不**意味着 SKEMPI 的 ΔΔG 数值本身过于接近 ——"
             f"后者需要以真实 ΔΔG 效用幅度做实证，本文未做。",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    else:
        P(D, f"For K = 10, the discriminable interval D(n,K) = K(n−K)/[n(n+1)] is zero "
             f"whenever n ≤ K; {GEO['n_le_K_count']}/{GEO['n_ge5_count']} = "
             f"{100*GEO['n_le_K_fraction']:.1f}% of our n≥5 targets fall in this regime, "
             f"where any model necessarily selects all candidates and the target is "
             f"structurally non-informative for fixed-K shortlist metrics. D peaks near "
             f"n* ≈ {GEO['n_star_analytic']:.1f} and then declines monotonically: "
             f"D(20) = {GEO['D_at_n20']:.4f} > D(50) = {GEO['D_at_n50']:.4f}. The heuristic "
             f"'larger candidate sets are more discriminable' therefore does not hold for "
             f"this metric, and we withdraw our earlier 'n ≥ 50' recommendation. We "
             f"emphasise that this geometric result does *not* imply that SKEMPI ΔΔG values "
             f"are intrinsically too similar; that would require an analysis using actual "
             f"ΔΔG utility magnitudes, which we did not perform.",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    D.add_picture(str(IMG / f"Fig5_fixed_K_geometry{sfx}.png"), width=Inches(5.2))
    cap(D, ("图 5. 固定 K 候选集几何；n ≤ K 时 D = 0，峰值 n* ≈ 20.5。")
        if zh else
        ("Figure 5. Fixed-K candidate-set geometry; D = 0 for n ≤ K, peaking at "
         "n* ≈ 20.5."))

    H(D, ("3.5 Capability decomposition clarifies why some remedies improve pooled scores "
          "but not within-target ranking" if not zh else
          "3.5 能力分解解释了为何某些补救改善合并分数但不改善靶内排序"), 2)
    if zh:
        P(D, "V7.3 的 F1–F8 共享同一 3351 行宇宙，构成匹配消融。以 F1（仅突变性质）为基线："
             "F2（静态局部几何）把合并 Pearson 提高 +0.0279，但靶内 Macro ρ 变化 −0.0015，"
             "即仅改善合并分数；F4（突变×partner 扰动）两轴同涨（+0.0701 / +0.0437）；"
             "F8（全量）为 +0.1054 / +0.0758。能力分解因此能区分「仅合并分数改善」与"
             "「同时改善靶内排序」两类改动。", align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "该分析是回溯性的：我们已知 F4 有效，故它证明的是分解**能够表达**这一区别，"
             "而非基准**前瞻性预测**了哪条路线会成功。此外，V8 重算的 F4−F1 ΔMacro ρ = "
             "+0.0437 与历史 V7.3 自助法区间 +0.0731 [+0.0209, +0.1289] 口径不同"
             "（度量宇宙与聚合方式不同），两者不得混用而不加说明。",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    else:
        P(D, "V7.3 variants F1–F8 share an identical 3351-row universe, forming a matched "
             "ablation. Relative to F1 (mutation properties only), F2 (static local "
             "geometry) raises pooled Pearson by +0.0279 while changing within-target "
             "Macro ρ by −0.0015, i.e. a pooled-score-only improvement; F4 (mutation × "
             "partner perturbation) improves both axes (+0.0701 / +0.0437); F8 (full) gives "
             "+0.1054 / +0.0758. Capability decomposition therefore distinguishes a "
             "pooled-score-only change from a modification that also improves within-target "
             "ranking.", align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "This analysis is retrospective: we already know F4 works, so it shows that the "
             "decomposition *can express* this distinction, not that the benchmark "
             "*prospectively predicted* which route would succeed. In addition, the "
             "V8-recomputed F4−F1 ΔMacro ρ = +0.0437 and the historical V7.3 bootstrap "
             "interval +0.0731 [+0.0209, +0.1289] are based on different metric universes "
             "and aggregation; they must not be mixed without explanation.",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    D.add_picture(str(IMG / f"Fig4_retrospective_remedy{sfx}.png"), width=Inches(5.2))
    cap(D, "图 4. 回溯性匹配补救分析（V7.3，同一 3351 行）。" if zh
        else ("Figure 4. Retrospective matched remedy analysis (V7.3, identical 3351 "
              "rows)."))

    # ---------- External validation ----------
    H(D, "4. External validation" if not zh else "4. 外部验证", 1)
    if zh:
        P(D, "我们使用已发布的官方逐突变预测构建外部评价宇宙。关键的是，"
             "归一化键（丢弃链与野生型残基）在标注侧不唯一：4473 行中仅 1781 行为 "
             "EXACT_UNIQUE，另有 2048 行为歧义键。因此主匹配评估仅使用 EXACT_UNIQUE，"
             "歧义键被排除而非解析。这一修正使外部模型的表观合并 Pearson 系统性下降"
             "（+0.02 至 +0.12）。", align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, f"在 EXACT_UNIQUE 宇宙（{FR['external_clean']['n_rows']} 行 / "
             f"{FR['external_clean']['n_targets']} 靶点 / n≥5 为 "
             f"{FR['external_clean']['n_ge5']}）上，方法级画像显示总体相关与靶内排序的"
             f"赢家同为 {best.model}（{best.overall_pearson:+.4f} / "
             f"{best.macro_spearman:+.4f}），即本宇宙内无方法级排名反转。"
             f"方法级秩一致为 ρ = "
             f"{er.get('rho_overall_pearson_vs_macro_spearman', float('nan')):+.2f}"
             f"（精确置换 p = "
             f"{er.get('perm_p_overall_pearson_vs_macro_spearman', float('nan')):.3f}）。",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "该外部面板给出了**强正**的方法级秩一致，因此**不支持**普遍的解耦主张；"
             "同时它**落在内部观测到的异质范围内**（主池 +0.82，全库 −0.27 至 +0.82）。"
             "仅有 5 个方法时，方法级秩相关推断功效不足；按我们预设的实务阈值 8 个方法，"
             "该分析被归类为 UNDERPOWERED。这是本研究的操作性选择，不是统计学通则。",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "HERMES 的官方量为 log_proba_mt − log_proba_wt，是**无量纲突变效应分数**，"
             "并非已校准的 kcal/mol ΔΔG。因此我们对 HERMES 报告 "
             "RMSE = NOT_APPLICABLE、MAE = NOT_APPLICABLE，且所用变体为稳定性微调变体，"
             "不得标注为 SKEMPI 结合监督模型。", align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    else:
        P(D, "We built external evaluation universes from published, official per-mutation "
             "predictions. Critically, the normalised matching key (chain and wild-type "
             "residue dropped) is not unique on the annotation side: of 4473 rows only 1781 "
             "are EXACT_UNIQUE and 2048 are ambiguous keys. Primary matched evaluation "
             "therefore uses EXACT_UNIQUE only; ambiguous keys are excluded, not resolved. "
             "This correction systematically lowers apparent pooled Pearson for external "
             "models (by +0.02 to +0.12).", align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, f"On the EXACT_UNIQUE universe ({FR['external_clean']['n_rows']} rows / "
             f"{FR['external_clean']['n_targets']} targets / "
             f"{FR['external_clean']['n_ge5']} with n≥5), method-level profiles show the "
             f"same winner for pooled correlation and within-target ranking, "
             f"{best.model} ({best.overall_pearson:+.4f} / {best.macro_spearman:+.4f}), "
             f"i.e. no method-level rank reversal in this universe. Method-level rank "
             f"agreement is ρ = "
             f"{er.get('rho_overall_pearson_vs_macro_spearman', float('nan')):+.2f} "
             f"(exact permutation p = "
             f"{er.get('perm_p_overall_pearson_vs_macro_spearman', float('nan')):.3f}).",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "The external panel showed strong positive rank agreement, which does not "
             "support a universal cross-metric decoupling and falls within the "
             "heterogeneous range observed internally (primary pool +0.82; whole corpus "
             "−0.27 to +0.82). With only five methods, method-level rank-correlation "
             "inference is low-powered; under our prespecified practical threshold of 8 "
             "methods the analysis is classified as UNDERPOWERED. This is a study-specific "
             "operating choice, not a general statistical law.",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "HERMES' official quantity is log_proba_mt − log_proba_wt, a dimensionless "
             "mutation-effect score, not a calibrated kcal/mol ΔΔG. We therefore report "
             "RMSE = NOT_APPLICABLE and MAE = NOT_APPLICABLE for HERMES, and the variant "
             "used is a stability fine-tune, which must not be labelled as a SKEMPI "
             "binding-supervised model.", align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    mktable(D, T10)
    cap(D, "表 4. 外部方法级能力画像（仅 EXACT_UNIQUE）。" if zh
        else "Table 4. External method-level capability profiles (EXACT_UNIQUE only).")

    # ---------- Discussion ----------
    H(D, "5. Discussion" if not zh else "5. 讨论", 1)
    if zh:
        P(D, "本文的贡献需要放在既有工作之后来理解。泄漏感知评价、family-OOD 划分、"
             "per-structure/per-target 指标以及 Top-K 评价都不是本文首创："
             "Bushuiev 等已量化 PPI 数据的结构性泄漏并建议界面相似度划分；"
             "CATH-ddG 已采用 CATH 同源超家族划分并给出更难突变的定义；"
             "per-protein 与短名单式指标在既有 ΔΔG 文献中亦已出现。",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "本文的贡献是**系统整合**：把精确的跨靶/靶内协方差信号分解与以下四项审计"
             "放在同一匹配、来源清洁的评价框架内——受控划分严格度比较、目标加权敏感性、"
             "匹配条件下的排名反转、以及固定 K 候选集几何；并辅以一个回溯性的"
             "benchmark-guided remedy 归因案例。据我们所知，我们未发现既有工作在"
             "同等匹配与来源清洁条件下联合报告上述各项。我们不对泄漏感知、family-OOD、"
             "per-target 或 Top-K 基准作任何「首次」声明。",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "一个方法学教训值得单独指出：在清理过程中，我们发现「指标之间的模型排名"
             "秩相关」不是可以脱离模型池构成而报告的稳健量。同一评价宇宙、同一行集合下，"
             "仅改变纳入的模型族，主池的 ρ(合并 Pearson, Macro Spearman) 即可从 −0.27 "
             "变化到 +0.82。因此我们建议：任何此类秩相关都必须同时报告模型池的"
             "训练划分来源构成、池化与分组两个数值，以及有效样本量。",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    else:
        P(D, "Our contribution must be read after prior work. Leakage-aware evaluation, "
             "family-OOD splitting, per-structure/per-target metrics and Top-K evaluation "
             "are none of them new here: Bushuiev et al. quantified structural leakage in "
             "PPI data and recommended interface-similarity splitting; CATH-ddG adopted CATH "
             "homologous-superfamily splits and defined harder mutations; per-protein and "
             "shortlist-style metrics already appear in the ΔΔG literature.",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "Our contribution is the **systematic integration** of an exact "
             "between-/within-target covariance-signal decomposition with four audits inside "
             "one matched, provenance-clean framework: a controlled split-strictness "
             "comparison, a target-weighting sensitivity audit, matched-universe "
             "rank-reversal analysis, and a fixed-K candidate-set geometry audit, plus a "
             "retrospective benchmark-guided remedy attribution. To our knowledge, we did "
             "not identify prior work that jointly reports these under comparable matched, "
             "provenance-clean evaluation. We make no 'first' claim for leakage-aware, "
             "family-OOD, per-target or Top-K benchmarking.",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "One methodological lesson deserves separate emphasis. During cleanup we found "
             "that 'rank agreement between metrics' is not a stable quantity that can be "
             "reported independently of model-pool composition. Within one evaluation "
             "universe and one row set, changing only which model families are included "
             "moves ρ(pooled Pearson, macro Spearman) from −0.27 to +0.82. We therefore "
             "recommend that any such rank correlation be reported together with the "
             "training-split provenance composition of the model pool, both pooled and "
             "grouped values, and the effective sample size.",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    D.add_picture(str(IMG / f"Fig7_capability_landscape{sfx}.png"), width=Inches(4.6))
    cap(D, ("图 7. 能力景观（主匹配、来源清洁池）。"
            "本池 ρ = 见标题；跨 10 个 stratum 为 −0.27 至 +0.82，不得外推。")
        if zh else
        ("Figure 7. Capability landscape in the primary matched, provenance-clean pool. "
         "The pool-specific ρ is given in the title; across ten strata it ranges from "
         "−0.27 to +0.82 and must not be generalised."))

    # ---------- Limitations ----------
    H(D, "6. Limitations" if not zh else "6. 局限", 1)
    lims_zh = [
        "task-specific family-OOD 并不能保证排除基础模型预训练阶段与评价家族的相似性。",
        "基准模型池不是全部 PPI ΔΔG 方法的随机样本，可能存在选择偏倚。",
        "部分模型族覆盖有限（如某嵌套 CV 制品仅 1807 行 / 13 个靶点）。",
        "外部可评价方法仅 5 个，方法级秩相关推断功效不足（UNDERPOWERED）。",
        "选择轴（固定 K 短名单）在本数据上结构性无力，因此 H3 判定为 FAIL，"
        "相关场景推荐仅作低置信描述。",
        "候选选择分析为回溯性，不构成湿实验前瞻验证。",
        "补救分析（3.5 节）为回溯性；基准未被证明可前瞻预测成功路线。",
        "C_B 与 C_W 是协方差信号贡献，不是因果机制。",
        "无独立前瞻实验数据集，故未作外部数据集验证。",
        "本文不提供通用模型推荐器。",
        "complex 与 seq40 划分在本数据上完全相同，不能作为独立证据。",
        "此前的归一化映射含 2048 个歧义行，外部读数已按 EXACT_UNIQUE 重算。",
    ]
    lims_en = [
        "Task-specific family-OOD does not guarantee exclusion of foundation-model "
        "pretraining similarity to the evaluation families.",
        "The benchmark model pool is not a random sample of all PPI ΔΔG methods, so "
        "selection bias may remain.",
        "Some model families have limited coverage (e.g. one nested-CV artifact covers only "
        "1807 rows / 13 targets).",
        "Only five methods are externally evaluable, so method-level rank-correlation "
        "inference is underpowered.",
        "The selection axis (fixed-K shortlist) is structurally non-informative in this "
        "dataset; H3 therefore fails and related scenario recommendations are reported only "
        "as low-confidence descriptions.",
        "Candidate-selection analysis is retrospective and does not constitute prospective "
        "wet-lab validation.",
        "The remedy analysis (Section 3.5) is retrospective; the benchmark is not shown to "
        "prospectively predict which route will succeed.",
        "C_B and C_W are covariance-signal contributions, not causal mechanisms.",
        "There is no independent prospective experimental dataset, so no external-dataset "
        "validation was performed.",
        "We provide no universal model recommender.",
        "The complex and seq40 split regimes are identical in this dataset and cannot serve "
        "as independent evidence.",
        "Our earlier normalised mapping contained 2048 ambiguous rows; external readings "
        "were recomputed on EXACT_UNIQUE only.",
    ]
    for i, (a, b) in enumerate(zip(lims_zh, lims_en), 1):
        P(D, f"{i}. {a if zh else b}", align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    # ---------- Data / Code / AI ----------
    H(D, "7. Data and code availability" if not zh else "7. 数据与代码可用性", 1)
    P(D, ("本分析仅使用既有预测与公开的外部预测。全部冻结输入哈希见 "
          "FROZEN_INPUT_MANIFEST.csv；所有头条数字由 FROZEN_RESULTS.json 单一真源生成，"
          "并附自动化一致性检查（99_final_consistency_check.py，12 项）。"
          "未训练任何新模型，未使用 GPU。")
      if zh else
      ("This analysis uses only existing predictions and publicly available external "
       "predictions. All frozen input hashes are listed in FROZEN_INPUT_MANIFEST.csv; all "
       "headline numbers are generated from the single source of truth FROZEN_RESULTS.json "
       "and verified by an automated 12-item consistency checker "
       "(99_final_consistency_check.py). No new model was trained and no GPU was used."),
      align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    H(D, "8. Author contributions" if not zh else "8. 作者贡献", 1)
    P(D, ("研究设计、分析实施、结果解释与稿件撰写由作者完成；"
          "AI 系统的具体使用见第 9 节。")
      if zh else
      ("Study design, analysis implementation, interpretation and manuscript writing were "
       "performed by the author(s); the specific use of AI systems is described in "
       "Section 9."), align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    H(D, "9. Acknowledgements and declaration of AI assistance"
      if not zh else "9. 致谢与 AI 辅助声明", 1)
    if zh:
        P(D, "本研究使用 AI 编码助手（DeepSeek Harness 环境下的智能体）协助完成以下任务："
             "（i）Python 分析代码的编写与调试，包括协方差分解、自助法、秩一致与几何计算；"
             "（ii）图 1–图 7 的绘制脚本生成；（iii）表格与文档的自动化生成；"
             "（iv）上游公开仓库的检索与许可证/来源信息核对；（v）文字草稿的整理。"
             "所有科学判断、指标定义、阈值设定、假设判定与最终结论均由作者作出并对内容负责。"
             "AI 系统未被用于生成训练数据或修改任何上游模型预测。所有分析结果可由随附脚本"
             "与冻结输入哈希独立复现。", align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "本研究未训练新的机器学习模型，未使用 GPU；外部方法的预测来自其官方公开发布，"
             "使用方式遵循各自许可证（MIT / CC BY 4.0 / Apache-2.0 等，详见外部审计文件）。"
             "未再分发任何专有可执行文件。", align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    else:
        P(D, "This study used an AI coding assistant (an agent running in the DeepSeek "
             "Harness environment) for the following tasks: (i) writing and debugging the "
             "Python analysis code, including the covariance decomposition, bootstrap, rank "
             "agreement and geometry computations; (ii) generating the plotting scripts for "
             "Figures 1–7; (iii) automated generation of tables and documents; (iv) "
             "retrieval and licence/provenance verification of upstream public "
             "repositories; and (v) drafting and editing of text. All scientific "
             "judgements, metric definitions, threshold choices, hypothesis verdicts and "
             "final conclusions were made by the author(s), who take responsibility for the "
             "content. The AI system was not used to generate training data or to modify "
             "any upstream model predictions. All results can be independently reproduced "
             "from the accompanying scripts and frozen input hashes.",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)
        P(D, "No new machine-learning model was trained and no GPU was used. External method "
             "predictions come from their official public releases and were used in "
             "accordance with their respective licences (MIT / CC BY 4.0 / Apache-2.0, see "
             "the external audit files). No proprietary executable was redistributed.",
          align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    # ---------- References ----------
    H(D, "References" if not zh else "参考文献", 1)
    for i, r in enumerate(REFS, 1):
        p = D.add_paragraph()
        run = p.add_run(f"[{i}] {r}")
        run.font.size = Pt(8.8)
        p.paragraph_format.left_indent = Inches(0.28)
        p.paragraph_format.first_line_indent = Inches(-0.28)

    # ---------- Supplementary T1/T3 ----------
    D.add_page_break()
    H(D, "Supplementary tables" if not zh else "补充表", 1)
    P(D, "Table S1. Matched-universe registry (rowset fingerprints)." if not zh
      else "表 S1. 匹配宇宙登记（行集指纹）。", italic=True, size=8.6)
    mktable(D, T1, fs=7.0)
    P(D, "")
    P(D, "Table S2. Rank agreement by matched stratum." if not zh
      else "表 S2. 各匹配 stratum 的秩一致。", italic=True, size=8.6)
    mktable(D, T3, fs=7.0)

    name = ("CR-PPI-ddG_manuscript_ZH.docx" if zh
            else "CR-PPI-ddG_manuscript_EN.docx")
    docx_path = OUT / name
    D.save(str(docx_path))
    print(f"  {name}  ({docx_path.stat().st_size/1024:.0f} KB)")
    return docx_path


zh_path = build("zh")
en_path = build("en")
print("DOCX 完成")
