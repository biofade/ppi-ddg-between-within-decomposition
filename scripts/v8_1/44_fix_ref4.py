# -*- coding: utf-8 -*-
"""scripts/v8_1/44_fix_ref4.py -- 把参考文献 [4] 从 HERMES 改为 DDGPred。

理由：正文 external panel 现在列出 Flex-ddG[6], FoldX[7], ProteinMPNN[8], ESM-IF[9]
与 DDGPred[4]。原 [4] 是 HERMES 条目 —— 而 HERMES 既不在冻结的 external panel 中，
也不是正文提到的方法，因此该条目为旧版残留。

DDGPred 的原始论文无法核实，故按项目一贯做法**不编造**作者/期刊，而是标注
REFERENCE_VERIFY_REQUIRED，并写明其逐突变预测取自 CATH-ddG 附带的基准表 [3]。
"""
import io
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "scripts/v8_1/31_manuscript_v6.py"
s = io.open(P, encoding="utf-8").read()

OLD = (' "Pun MN, Ivanov A, Bellinzoni M, Bayer PE, Gagneur J. HERMES: holographic equivariant "\n'
       ' "neural network model for mutational effect and stability prediction. arXiv:2407.06703. 2024.",')
NEW = (' "DDGPred: primary method publication. REFERENCE_VERIFY_REQUIRED. The per-mutation "\n'
       ' "DDGPred values analysed here were obtained from the benchmark table distributed with "\n'
       ' "CATH-ddG [3]; the primary DDGPred publication could not be verified in this work and no "\n'
       ' "bibliographic fields have been invented.",')

n = s.count(OLD)
s = s.replace(OLD, NEW)
io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
try:
    ast.parse(s)
    print("syntax OK")
except SyntaxError as e:
    print("SyntaxError", e.lineno, e.text)
    raise SystemExit(1)
print("reference [4] repointed to DDGPred (x%d); HERMES left: %d" % (n, s.count("HERMES")))
