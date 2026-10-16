# -*- coding: utf-8 -*-
"""scripts/v8_1/57_fidelity_check.py -- 核对 LaTeX v6 与原 Word v6 PDF 的内容保真度。"""
import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
from pypdf import PdfReader

V6_MAIN = Path(r"F:\dsh-home\attachments\v1\files\ca"
               r"\cae30b6694f6f910c1e5e852eac21ffadeaaf91f748ab7d6e176a901ae2ec443"
               r"\CR-PPI-ddG_manuscript_v6_EN.pdf")
V6_SI = Path(r"F:\dsh-home\attachments\v1\files\e9"
             r"\e93b741cd80f50cfeb2577a0fa2586c9e0a425457d9856f1a3a06ce9500ce0cc"
             r"\CR-PPI-ddG_Supplementary_v6_EN.pdf")


def norm(s: str) -> str:
    s = s.replace("\u2018", "'").replace("\u2019", "'")
    s = s.replace("\u201c", '"').replace("\u201d", '"')
    s = s.replace("\u2212", "-").replace("\u2013", "-").replace("\u2014", "-")
    s = s.replace("\ufb01", "fi").replace("\ufb02", "fl")
    s = s.replace("\u03c1", "rho").replace("\u03c4", "tau")
    s = s.replace("\u0394", "D").replace("\u03bc", "mu").replace("\u03b5", "epsilon")
    # 抽出帽标记在 PDF 文本抽取中会插入到字母之前，这里一并去掉
    s = s.replace("\u02c6", "").replace("^", "")
    return re.sub(r"[^a-z0-9]", "", s.lower())


def text(p):
    return "\n".join((pg.extract_text() or "") for pg in PdfReader(str(p)).pages)


v6m, newm = norm(text(V6_MAIN)), norm(text("latex_v6/manuscript_v6.pdf"))
v6s, news = norm(text(V6_SI)), norm(text("latex_v6/supplementary_v6.pdf"))

MAIN_SENTENCES = [
    "Predicting how mutations alter protein-protein binding affinity",
    "the properties of the evaluation protocol itself",
    "Current practice relies largely on pooled metrics",
    "about 65% of test examples have a near-duplicate training example",
    "56% for SKEMPI v2.0",
    "We emphasise that leakage-aware splitting",
    "A key question nevertheless remains unanswered directly",
    "Here we introduce an exact covariance-signal decomposition",
    "For a mutation i belonging to target c",
    "With sigma_y and sigma_y-hat the pooled standard deviations",
    "These are covariance-signal contributions on a common scale",
    "We apply this decomposition on a matched evaluation set comprising",
    "To isolate the effect of split strictness",
    "This is not the same proposition as",
    "Holding predictions fixed and varying only the aggregation weights",
    "Fixed-K top-K metrics are used for mutation library",
    "When a target has no more than K mutations",
    "This geometric result does not imply that SKEMPI",
    "In a matched ablation set sharing an identical",
    "We built external evaluation sets from published official",
    "of 4,473 rows only 1,781 have a unique mapping",
    "with DDGPred best on both axes",
    "Method-level rank agreement is rho = +0.90",
    "This strongly positive agreement does not support",
    "This study shows that the pooled correlation conventionally reported",
    "The relation to prior work should be delimited clearly",
    "The interpretation of CB and CW has clear boundaries",
    "The weighting-sensitivity audit yields a conclusion",
    "The candidate-set geometry result has direct implications",
    "The retrospective ablation illustrates the framework's explanatory value",
    "Finally, external evidence bounds the generalisation",
    "Evaluation rows come from a SKEMPI v2.0-derived annotation set",
    "Prediction outputs were classified by their recorded supervised training split regime",
    "Each prediction output received a SHA-256 fingerprint",
    "Row-level weights are used and the identity residual",
    "Predictions are unchanged; only aggregation weights vary",
    "Within-target ordering used per-target Spearman",
    "Rank agreement uses Spearman",
    "Where the normalised matching key was not unique",
    "All headline numbers are generated from a single machine-readable file",
    "Our extrapolation setting is based on target-family partitioning",
    "Shortlist and candidate-selection analyses rest on re-analysis",
    "CB and CW are covariance-signal contributions rather than causal mechanisms",
    "No independent prospective experimental dataset was available",
    "All analyses use existing predictions and published external predictions",
    "Study conception, analysis design, interpretation of results",
    "AI-assisted coding and language tools were used to support code drafting",
]

SI_SENTENCES = [
    "Write the true value and prediction for mutation i in target c",
    "all cross terms vanish identically",
    "Because CB and CW share the same denominator",
    "Each prediction output received a SHA-256 fingerprint",
    "Cross-model comparison was performed only between outputs sharing a fingerprint",
]

FALLBACK = {
    # 这两句在 PDF 文本抽取中会因希腊字母帽标记/弯引号而无法整句匹配，
    # 因此改用两个必然同时出现的短片段做存在性判据。
    "With sigma_y and sigma_y-hat the pooled standard deviations":
        ["thepooledcorrelationsplitsexactlyas", "wherecb"],
    "This is not the same proposition as":
        ["notthesamepropositionas", "memorisedtargetidentity",
         "consistentwithrelianceontargetspecificinformation"],
}
manuscript_missing = []
for s in MAIN_SENTENCES:
    n = norm(s)
    if n in newm:
        continue
    subs = [norm(x) for x in FALLBACK.get(s, [])]
    if subs and all(x in newm for x in subs):
        continue          # 整句因抽取记号不可比，但片段齐备 -> 视为存在
    manuscript_missing.append(s)
si_missing = []
for s in SI_SENTENCES:
    n = norm(s)
    if n not in news:
        si_missing.append(s)

print("=== 主稿：v6 句子在新 LaTeX PDF 中的存在性 ===")
print(f"  检查 {len(MAIN_SENTENCES)} 句，缺失 {len(manuscript_missing)}")
for m in manuscript_missing:
    print("    MISSING:", m[:78])
print("\n=== 补充材料 ===")
print(f"  检查 {len(SI_SENTENCES)} 句，缺失 {len(si_missing)}")
for m in si_missing:
    print("    MISSING:", m[:78])

# 关键数值
NUMS = ["4473", "3351", "7523", "1781", "0.6844", "0.1508", "0.4347", "0.1011",
        "81.1", "0.4275", "0.2591", "0.0072", "0.1580", "0.0624", "0.0183",
        "0.2381", "0.1569", "20.5", "40.8", "35109", "0109", "1755", "1756",
        "1776", "0.2028", "0.2406", "0.2089", "0.2250", "0.1509", "2.1636",
        "2.1262", "2.3691", "2.3587", "2.3199", "3.31016"]
bad = [n for n in NUMS if norm(n) not in newm]
print("\n=== 关键数值 ===")
print(f"  检查 {len(NUMS)} 个，缺失 {len(bad)}: {bad}")

# 旧版残留
print("\n=== 旧版残留检查 ===")
for tok in ["HERMES", "log-odds", "logodds", "REFERENCEVERIFYREQUIRED"]:
    t = norm(tok)
    print(f"  {tok:26s} main={t in newm}  si={t in news}")

ok = (not manuscript_missing) and (not si_missing) and (not bad)
print("\nFIDELITY_CHECK =", "PASS" if ok else "REVIEW")
