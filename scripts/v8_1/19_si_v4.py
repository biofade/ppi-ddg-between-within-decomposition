# -*- coding: utf-8 -*-
"""scripts/v8_1/19_si_v4.py -- 补充材料 v4：表格排版/精度/中文残留修复 + 新图。

修复（第六部分）：
* Note 10 的 disposition 中文残留 -> 英文；状态码补空格
* Table S3 标题 "first 30 of 29" -> "(all 29 metrics)"
* Table S4/S5 浮点统一 4 位小数（极小值科学计数法）
* Table S2/S3/S5/S6 列宽自适应 + 8 pt，避免单词被拆行
* 图改为 v4 版（图内无硬编码标题/脚注，编号 S1/S2/S3 正确）
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd
from docx import Document
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
OUT = F / "manuscript_v4"
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

# Note 10：disposition 全部英文化 + 补空格
DISP = {
    "S1": "Single-value near-orthogonality (ρ = −0.26) does not hold; the cross-pool "
          "range includes strong positive agreement (up to +0.82), so near-orthogonality "
          "is not a universal conclusion.",
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
                    "Remains defensible": "No" if not v.get("remains_defensible")
                    else "Yes",
                    "Disposition": DISP.get(k, "—")}
                   for k, v in _SV.items()])


def f4(v):
    """统一 4 位小数；极小值用科学计数法。"""
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
    st.font.size = Pt(10.5); st.font.name = "Times New Roman"
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "SimSun")
    st.paragraph_format.space_after = Pt(6)
    st.paragraph_format.line_spacing = 1.15
    for s in doc.sections:
        s.top_margin = s.bottom_margin = Inches(0.9)
        s.left_margin = s.right_margin = Inches(0.9)


def H(doc, t, lvl=1):
    p = doc.add_heading(t, level=lvl)
    for r in p.runs:
        r.font.color.rgb = RGBColor(0, 0, 0)
        r.font.size = Pt(13.5 if lvl == 1 else 11.5)
    p.paragraph_format.keep_with_next = True
    return p


def P(doc, t, size=None):
    p = doc.add_paragraph(); r = p.add_run(t)
    if size:
        r.font.size = Pt(size)
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    return p


def cap(doc, t):
    p = doc.add_paragraph(); r = p.add_run(t); r.font.size = Pt(8.8)
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    return p


def _autofit(table, total_in=6.7):
    tbl = table._tbl
    tblPr = tbl.tblPr
    lay = tblPr.find(qn("w:tblLayout"))
    if lay is None:
        lay = tblPr.makeelement(qn("w:tblLayout"), {}); tblPr.append(lay)
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
    for j, w in enumerate(weights):
        ew = Emu(int(Inches(total_in) * w / tot))
        for cell in table.columns[j].cells:
            cell.width = ew


def tbl(doc, df, fs=8.0, note=None, total_in=6.7):
    t = doc.add_table(rows=1, cols=len(df.columns))
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = False
    hdr = t.rows[0]
    trPr = hdr._tr.get_or_add_trPr()
    trPr.append(trPr.makeelement(qn("w:tblHeader"), {}))
    for i, c in enumerate(df.columns):
        cell = hdr.cells[i]; cell.text = ""
        r = cell.paragraphs[0].add_run(str(c)); r.bold = True; r.font.size = Pt(fs)
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


# ---------- 准备表格数据（全部修精度/列名） ----------
UR_t = UR.copy()
for c in [c for c in UR_t.columns if UR_t[c].dtype.kind in "if"]:
    UR_t[c] = UR_t[c].map(lambda x: f"{int(x)}" if float(x).is_integer() else f4(x))

PA_t = PA.head(40).copy()

MD_t = MD.copy()
for c in ["raw_direction", "ranking_direction", "transform_for_rank_comparison",
          "eligible_for_calibration"]:
    if c in MD_t.columns:
        MD_t[c] = MD_t[c].astype(str).str.replace("_", " ", regex=False)

GEO_t = GEO_DF.iloc[::40].head(11).copy()
for c in ["E_rand", "E_max", "D"]:
    GEO_t[c] = GEO_t[c].map(f4)

W2_t = W2.head(30).copy()
for c in W2_t.columns:
    if c != "model_name":
        W2_t[c] = W2_t[c].map(f4)
W2_t = W2_t.rename(columns={"model_name": "Model configuration",
                            "range_C_B": "Range C_B", "range_C_W": "Range C_W",
                            "max_abs_C_B": "Max |C_B|", "max_abs_C_W": "Max |C_W|",
                            "signflip_C_B_gt0": "C_B flips (>0)",
                            "signflip_C_W_gt0": "C_W flips (>0)",
                            "signflip_C_B_gt0.001": "C_B flips (>0.001)",
                            "signflip_C_W_gt0.001": "C_W flips (>0.001)",
                            "signflip_C_B_gt0.005": "C_B flips (>0.005)",
                            "signflip_C_W_gt0.005": "C_W flips (>0.005)",
                            "signflip_C_B_gt0.01": "C_B flips (>0.01)",
                            "signflip_C_W_gt0.01": "C_W flips (>0.01)"})

keep = [c for c in ["artifact_id", "method_family", "seed", "split_regime",
                    "rowset_hash", "is_ablation", "is_primary_representative"]
        if c in MR.columns]
MR_t = MR[keep].head(60).copy()
for c in ["is_ablation", "is_primary_representative"]:
    if c in MR_t.columns:
        MR_t[c] = MR_t[c].map(lambda x: "yes" if bool(x) else "no")
MR_t = MR_t.rename(columns={"artifact_id": "Artifact", "method_family": "Method family",
                            "seed": "Seed", "split_regime": "Split regime",
                            "rowset_hash": "Row-set fingerprint",
                            "is_ablation": "Ablation",
                            "is_primary_representative": "Primary rep."})


def build(zh):
    T = (lambda en, z: z if zh else en)
    SUF = "_ZH" if zh else "_EN"
    D = Document(); setup(D)
    p = D.add_paragraph(); r = p.add_run(T("Supplementary Information", "补充材料"))
    r.bold = True; r.font.size = Pt(20)
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p = D.add_paragraph()
    r = p.add_run(T("Disentangling Between-Target and Within-Target Signal in "
                    "PPI ΔΔG Benchmarking",
                    "解耦 PPI ΔΔG 基准评价中的跨靶信号与靶内信号"))
    r.italic = True; r.font.size = Pt(11)
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    D.add_page_break()

    H(D, T("Supplementary Note 1 | Derivation of the covariance decomposition",
           "补充说明 1 | 协方差分解的推导"), 1)
    P(D, T("Write the true value and prediction for mutation i in target c as "
           "y_ci = μ_c + ε_ci and ŷ_ci = μ̂_c + ε̂_ci, where μ_c and μ̂_c are target means "
           "and ε, ε̂ are within-target deviations whose target-wise means vanish by "
           "construction. Expanding Cov(y, ŷ) gives a between-target term "
           "Σ_c n_c (μ_c − μ̄)(μ̂_c − μ̄̂)/N and a within-target term Σ_c Σ_i ε_ci ε̂_ci/N; "
           "all cross terms vanish identically because within-target deviations have zero "
           "mean in every target. Dividing by σ_y σ_ŷ gives C_B and C_W and establishes "
           "r_pooled = C_B + C_W exactly. Row-level weights are used throughout and the "
           "identity residual is reported as a quality-control quantity (observed maximum "
           "3.331 × 10⁻¹⁶). Because C_B and C_W share the same denominator they are "
           "covariance-signal contributions on a common scale, not independent "
           "correlations or independently estimated skill scores.",
           "把靶点 c 中突变 i 的真值与预测写为 y_ci = μ_c + ε_ci 与 ŷ_ci = μ̂_c + ε̂_ci，"
           "其中 μ_c 与 μ̂_c 为靶点均值，ε 与 ε̂ 为靶内偏差（按构造其靶内均值为零）。"
           "展开 Cov(y, ŷ) 得到一个跨靶项与一个靶内项；由于靶内偏差在每个靶点内均值为零，"
           "交叉项恒为零。除以 σ_y σ_ŷ 即得 C_B 与 C_W，故 r_pooled = C_B + C_W 精确成立。"
           "全文使用行级权重，并报告恒等式残差作为质量控制量（实测最大 3.331 × 10⁻¹⁶）。"
           "由于 C_B 与 C_W 共享同一分母，它们是同一尺度上的协方差信号贡献，"
           "而非独立相关系数或独立能力分数。"), size=10)

    H(D, T("Supplementary Note 2 | Matched universes and row-set fingerprints",
           "补充说明 2 | 匹配宇宙与行集指纹"), 1)
    P(D, T("Each prediction output received a row-set fingerprint, the SHA-256 hash of its "
           "sorted row identifiers. Cross-model comparison was performed only between "
           "outputs sharing a fingerprint; raw metrics were never compared across "
           "fingerprints. The primary matched set was selected by a rule fixed before any "
           "metric values were inspected: provenance cleanliness, then number of "
           "independent method families, then number of targets with n ≥ 5, then number of "
           "rows.",
           "每个预测输出都被赋予一个行集指纹，即其行标识符排序列表的 SHA-256 哈希。"
           "跨模型比较仅在同一指纹内进行；绝不在不同指纹之间比较原始指标。"
           "主匹配集合由在查看任何指标数值之前即已固定的规则选出：来源清洁度、"
           "独立方法族数、n ≥ 5 靶点数、行数。"), size=10)
    tbl(D, UR_t, note=T("Table S1 | Matched-evaluation strata and row-set fingerprints. "
                        "Universe identifiers are given for cross-reference to the "
                        "machine-readable results file.",
                        "表 S1 | 匹配评价分层与行集指纹。"))

    H(D, T("Supplementary Note 3 | Provenance filtering and de-replication",
           "补充说明 3 | 来源筛选与去伪重复"), 1)
    P(D, T("Provenance was assigned from the recorded training split regime rather than "
           "inferred from names. Configurations trained under a permissive (random / "
           "complex / seq40) partition were excluded from deployment-oriented comparisons "
           "and retained only in the controlled split experiment, where the split regime "
           "is itself the manipulated variable. Analyses were performed at artifact, "
           "configuration and method-family level; representatives were chosen "
           "performance-blind.",
           "来源状态依据记录的训练划分方式判定，而非依据名称推断。在宽松划分"
           "（随机/复合物/seq40）下训练的配置被排除在面向部署的比较之外，"
           "仅保留在受控划分实验中。分析在制品层、配置层与方法族层进行；"
           "代表的选择为 performance-blind。"), size=10)
    tbl(D, PA_t, note=T(f"Table S2 | Provenance audit, showing the first 40 of "
                        f"{len(PA)} entries. Training split regime determines provenance "
                        f"status.",
                        f"表 S2 | 来源审计，共 {len(PA)} 条，示前 40 条。"))

    H(D, T("Supplementary Note 4 | Metric-direction registry and RMSE-direction correction",
           "补充说明 4 | 指标方向登记与 RMSE 方向修正"), 1)
    P(D, T("A metric-direction registry was established before any rank-agreement "
           "computation. RMSE and MAE are lower-is-better; all cross-metric rank agreement "
           "uses negated RMSE (or ascending ranks) so that a higher-is-better correlation "
           "is never correlated against a raw higher-is-worse error. An earlier version of "
           "the calibration-axis summary violated this rule and has been superseded; the "
           "corrected values are used throughout. The registry comprises all 29 metrics "
           "used in this work.",
           "在任何秩一致计算之前即建立了指标方向登记表。RMSE 与 MAE 为越小越好；"
           "所有跨指标秩一致均使用取负后的 RMSE（或升序名次），以避免把「越大越好」的"
           "相关系数与「越大越差」的原始误差直接相关。早期版本的校准轴汇总违反了该规则"
           "并已被取代；本稿全部使用修正后的数值。登记表涵盖本研究使用的全部 29 个指标。"),
      size=10)
    tbl(D, MD_t, fs=7.4,
        note=T(f"Table S3 | Metric-direction registry (all {len(MD)} metrics).",
               f"表 S3 | 指标方向登记（全部 {len(MD)} 个指标）。"))

    H(D, T("Supplementary Note 5 | External mapping ambiguity",
           "补充说明 5 | 外部映射歧义"), 1)
    P(D, T(f"The normalised matching key, formed by dropping chain identity and the "
           f"wild-type residue, is not unique on the annotation side: of 4,473 rows "
           f"{CM['n_exact_unique']:,} are uniquely mapped, {CM['n_ambiguous_key']:,} are "
           f"ambiguous and {CM['n_unmatched']:,} are unmatched. Primary matched evaluation "
           f"used uniquely mapped rows only; ambiguous keys were excluded rather than "
           f"resolved. The official held-out CATH partition carries its own unambiguous "
           f"rows and was evaluated separately without canonical re-mapping. For the "
           f"external quantity log_proba_mt − log_proba_wt, RMSE and MAE are "
           f"NOT_APPLICABLE because the quantity is a dimensionless mutation-effect score "
           f"rather than a calibrated kcal/mol ΔΔG.",
           f"由丢弃链标识与野生型残基得到的归一化匹配键在标注侧并不唯一：4,473 行中 "
           f"{CM['n_exact_unique']:,} 行为唯一映射，{CM['n_ambiguous_key']:,} 行为歧义键，"
           f"{CM['n_unmatched']:,} 行未匹配。主匹配评价仅使用唯一映射行；"
           f"歧义键被排除而非解析。官方 held-out CATH 划分自带无歧义行，"
           f"已单独评价且未做 canonical 重映射。"), size=10)
    cap(D, T(f"Mapping QC: canonical rows before merge = "
             f"{HS['qc']['canonical_rows_before_merge']}, after merge = "
             f"{HS['qc']['canonical_rows_after_merge']}, duplicate row identifiers = "
             f"{HS['qc']['duplicate_canonical_row_id']}, exact matches = "
             f"{HS['qc']['exact_match_count']}, ambiguous = {HS['qc']['ambiguous_count']}, "
             f"unmatched = {HS['qc']['unmatched_count']}.",
             f"映射 QC：合并前 canonical 行 = "
             f"{HS['qc']['canonical_rows_before_merge']}，合并后 = "
             f"{HS['qc']['canonical_rows_after_merge']}，重复行标识符 = "
             f"{HS['qc']['duplicate_canonical_row_id']}，精确匹配 = "
             f"{HS['qc']['exact_match_count']}，歧义 = {HS['qc']['ambiguous_count']}，"
             f"未匹配 = {HS['qc']['unmatched_count']}。"))

    H(D, T("Supplementary Note 6 | Fixed-K geometry derivation",
           "补充说明 6 | 固定 K 几何推导"), 1)
    P(D, T("For a target with n mutations and fixed K, let the rank-derived gain of the K "
           "most favourable mutations be the available utility. Choosing K candidates "
           "uniformly at random without replacement yields an expected gain share "
           "E_rand = K/n by symmetry, while the optimal K candidates attain "
           "E_max = K(2n − K + 1)/[n(n + 1)]. The discriminable interval is therefore "
           "D(n, K) = E_max − E_rand = K(n − K)/[n(n + 1)] for n > K, and D = 0 for "
           "n ≤ K, since every candidate is then selected. D is not monotonic in n: "
           "dD/dn = 0 at n* = K + √(K² + K), giving n* ≈ 20.5 for K = 10, after which D "
           "declines; for example D(20) = 0.2381 and D(50) = 0.1569.",
           "对突变数为 n、固定 K 的靶点，设最优 K 个突变提供可用收益。"
           "从 n 个候选中无放回随机抽取 K 个，由对称性得期望收益占比 E_rand = K/n；"
           "最优 K 个候选为 E_max = K(2n − K + 1)/[n(n + 1)]。因此可区分区间为 "
           "D(n, K) = E_max − E_rand = K(n − K)/[n(n + 1)]（n > K）；n ≤ K 时所有候选"
           "都被选中，故 D = 0。D 对 n 并非单调：在 n* = K + √(K² + K) 处导数为零，"
           "K = 10 时 n* ≈ 20.5，此后 D 下降；例如 D(20) = 0.2381 而 D(50) = 0.1569。"),
      size=10)
    tbl(D, GEO_t, note=T("Table S4 | Discriminable interval versus candidate-set size "
                         "(K = 10), every fortieth value.",
                         "表 S4 | 可区分区间与候选集大小的关系（K = 10），每 40 个取一个。"))

    H(D, T("Supplementary Note 7 | Full weighting audit",
           "补充说明 7 | 完整加权审计"), 1)
    P(D, T("Values are rounded to four decimal places; ranges smaller than 1 × 10⁻⁴ are "
           "given in scientific notation. Sign-flip columns count the weighting schemes in "
           "which the component changes sign, at the stated threshold on its absolute value.",
           "数值四舍五入至 4 位小数；小于 1 × 10⁻⁴ 的极差以科学计数法给出。"
           "符号翻转列统计在给定绝对值阈值下该分量发生符号改变的加权方案数。"), size=10)
    tbl(D, W2_t, fs=7.4,
        note=T(f"Table S5 | Weighting effect-size audit (first 30 of {len(W2)} "
               f"configurations).",
               f"表 S5 | 加权效应量审计（共 {len(W2)} 个配置，示前 30）。"))

    H(D, T("Supplementary Note 8 | Full model registry",
           "补充说明 8 | 完整模型注册表"), 1)
    tbl(D, MR_t, fs=7.4,
        note=T("Table S6 | Model family registry (first 60 entries). Row-set fingerprints "
               "allow outputs to be grouped into comparable matched universes.",
               "表 S6 | 模型族注册表（示前 60 条）。行集指纹用于把输出归入可比较的匹配宇宙。"))

    H(D, T("Supplementary Note 9 | Automated workflow and reproducibility",
           "补充说明 9 | 自动化工作流与可复现性"), 1)
    P(D, T("AI-assisted coding tools were used to draft and debug the analysis and "
           "plotting scripts, and to automate document generation, as declared in the main "
           "text. All headline numbers are generated from a single machine-readable results "
           "file; an automated consistency checker verifies that the documents, tables and "
           "figures agree with it. Metric definitions, analysis choices, threshold "
           "settings, hypothesis verdicts and final conclusions were made by the authors. "
           "The AI system was not used to generate data or to modify any upstream model "
           "prediction, and all results can be regenerated from the accompanying scripts "
           "and frozen input hashes. No new model was trained and no GPU was used.",
           "如正文所声明，本研究使用 AI 辅助编码工具草拟与调试分析与绘图脚本，"
           "并自动化生成文档。所有头条数字由单一机器可读结果文件生成；"
           "自动化一致性检查器验证文档、表格与图与该文件一致。"
           "指标定义、分析选择、阈值设定、假设判定与最终结论均由作者作出。"
           "AI 系统未被用于生成数据或修改任何上游模型预测，"
           "全部结果可由随附脚本与冻结输入哈希重新生成。未训练新模型，未使用 GPU。"),
      size=10)

    H(D, T("Supplementary Note 10 | Historical claim audit",
           "补充说明 10 | 历史主张审计"), 1)
    P(D, T("For transparency we record that four headline statements from an earlier "
           "version of this analysis did not survive the matched, de-duplicated, "
           "provenance-clean re-analysis and were withdrawn. Their dispositions are "
           "tabulated below.",
           "为透明起见，我们记录：本分析早期版本的四个头条陈述在匹配、去伪重复、"
           "来源清洁的重新分析后未能成立并已撤回。其处置见下表。"), size=10)
    tbl(D, SV[["Historical claim", "Remains defensible", "Disposition"]], fs=8.4)

    H(D, T("Supplementary Figures", "补充图"), 1)
    for stem, folder, w, c_en, c_zh in [
        ("FigS1_matched_universes", FS, 6.4,
         "Figure S1 | Matched-evaluation strata and row-set fingerprints. Each bar is one "
         "prediction output; outputs sharing a row-set fingerprint are directly comparable "
         "with each other and with no others.",
         "图 S1 | 匹配评价分层与行集指纹。每根条为一个预测输出；"
         "共享同一行集指纹的输出之间可直接比较。"),
        ("Capability_heatmap_provenance_clean", FS, 6.4,
         "Figure S2 | Capability profiles of matched, provenance-clean method families. "
         "Rows are the controlled ablation variants of a single local-physics method "
         "family (seed 42), not independent methods. Values are standardised within each "
         "metric (z-score) for visual comparison only; raw values are printed in each cell "
         "and statistical analyses use the raw metrics. This display does not imply any "
         "cross-metric correlation.",
         "图 S2 | 匹配且来源清洁的方法族能力画像。行为同一局部物理方法族的受控消融变体"
         "（seed 42），并非独立方法。数值仅按列标准化（z 分数）以便视觉比较；"
         "每格内给出原始数值，统计分析使用原始指标。该展示不隐含任何跨指标相关。"),
        ("FigS3_capability_pool_sensitivity", FS, 6.4,
         "Figure S3 | Pool sensitivity of cross-metric rank agreement. Within one "
         "evaluation universe and row set, changing only which model families are included "
         "moves the agreement between pooled Pearson and macro within-target Spearman "
         "across a wide range, which is why this quantity is reported as a range rather "
         "than as a single property of a benchmark.",
         "图 S3 | 跨指标秩一致的池敏感性。在同一评价宇宙与同一行集合内，"
         "仅改变纳入的模型族即可使合并 Pearson 与靶内 Macro Spearman 的一致性"
         "在很宽范围内变化，因此该量以范围而非单一性质报告。"),
    ]:
        D.add_picture(str(folder / f"{stem}{SUF}.png"), width=Inches(w))
        cap(D, c_zh if zh else c_en)

    nm = f"CR-PPI-ddG_Supplementary_v4_{'ZH' if zh else 'EN'}.docx"
    path = OUT / nm
    D.save(str(path))
    print(f"  {nm} ({path.stat().st_size/1024:.0f} KB)")


build(False)
build(True)
print("v4 supplementary 完成")
