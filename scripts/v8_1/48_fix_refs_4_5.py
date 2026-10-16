# -*- coding: utf-8 -*-
"""scripts/v8_1/48_fix_refs_4_5.py -- 修正参考文献 [4] 与 [5]，并给 [2] 补 DOI。

用户核查结论（Google Scholar / PubMed / Oxford Academic / OpenReview）：
* [4] DDGPred 占位符 -> 真实原始文献 Shan et al. 2022, PNAS
* [5] PreFold-dG 作者为 AI 幻觉名 -> 真实作者 Park S 等（LG AI Research/KAIST）
      标题末尾 model（单数）；卷期 2026;42(Suppl_2)
* [2] Bushuiev 泄露论文补 DOI:10.48550/arXiv.2404.10457
* [10] 保持不动（改为 [4]=DDGPred 后 [10] 仍在正文出现，无孤立文献）

Source of truth: 用户提供的逐条核查结果。不新增任何未核实的书目字段。
"""
import io
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "scripts/v8_1/31_manuscript_v6.py"
s = io.open(P, encoding="utf-8").read()
orig = s

# ---------------- [4] DDGPred: 占位符 -> Shan et al. 2022 PNAS ----------------
OLD4 = (' "DDGPred: primary method publication. REFERENCE_VERIFY_REQUIRED. The per-mutation "\n'
        ' "DDGPred values analysed here were obtained from the benchmark table distributed with "\n'
        ' "CATH-ddG [3]; the primary DDGPred publication could not be verified in this work and no "\n'
        ' "bibliographic fields have been invented.",')
NEW4 = (' "Shan S, Luo S, Yang Z, Hong J, Su Y, Ding F, et al. Deep learning guided optimization "\n'
        ' "of human antibody against SARS-CoV-2 variants with broad neutralization. Proc Natl Acad "\n'
        ' "Sci U S A. 2022;119(11):e2122954119. doi:10.1073/pnas.2122954119",')

# ---------------- [5] PreFold-dG: 修正作者、单数 model、卷期 ----------------
OLD5 = (' "Wang Y, Zhang H, Chen X, Li J. PreFold-dG: estimating binding affinity of protein–protein "\n'
        ' "interaction from intermediate representations of protein folding models. Bioinformatics. "\n'
        ' "2026. doi:10.1093/bioinformatics/btag489",')
NEW5 = (' "Park S, Yim S, Kim D, Yoo K, Hwang D, Lee K, Jang J, Kim K. PreFold-dG: estimating "\n'
        ' "binding affinity of protein–protein interaction from intermediate representations of "\n'
        ' "protein folding model. Bioinformatics. 2026;42(Suppl_2):btag489. "\n'
        ' "doi:10.1093/bioinformatics/btag489",')

# ---------------- [2] Bushuiev 泄露论文补 DOI ----------------
OLD2 = (' "data leakage in protein interaction benchmarks. arXiv:2404.10457. 2024.",')
NEW2 = (' "data leakage in protein interaction benchmarks. arXiv:2404.10457. 2024. "\n'
        ' "doi:10.48550/arXiv.2404.10457",')

applied = []
for old, new, tag in ((OLD4, NEW4, "[4] DDGPred -> Shan et al. 2022 PNAS"),
                      (OLD5, NEW5, "[5] PreFold-dG authors corrected"),
                      (OLD2, NEW2, "[2] Bushuiev DOI added")):
    n = s.count(old)
    if n:
        s = s.replace(old, new)
        applied.append(f"{tag} (x{n})")
    else:
        applied.append(f"{tag} -- NOT FOUND")

io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print("syntax OK | changed:", s != orig)
for a in applied:
    print("  ", a)
print("REFERENCE_VERIFY_REQUIRED left:", s.count("REFERENCE_VERIFY_REQUIRED"))
print("hallucinated authors left:", s.count("Wang Y, Zhang H, Chen X, Li J"))
print("i362 lowercase preserved:", "i362–i372" in s)
