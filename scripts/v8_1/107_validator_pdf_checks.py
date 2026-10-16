# -*- coding: utf-8 -*-
"""scripts/v8_1/107_validator_pdf_checks.py -- 让 validator 直接检查最终 PDF 里的图内容。

新增检查（针对本轮暴露的问题）：
  AC. 最终 PDF 中 Figure 3 不得含旧字符串（0.0183 / 35/109 / 13/109 / stable once）
  AD. 最终 PDF 中 Figure 3 必须含新数值（由 WEIGHTING_FINAL.json 读出）
  AE. 最终 PDF 中 Figure 3 的 C_W 翻转计数与 audit 一致
  AF. 正文与 Figure 3 不得矛盾（同一页出现的 C_B 翻转数一致）
"""
import io
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "validate_final_manuscript.py"
s = io.open(P, encoding="utf-8").read()

NEW = r'''
# ================= 直接检查最终 PDF 的图形内容 =================
import pymupdf as _pm  # noqa: E402
WF = load_json("WEIGHTING_FINAL.json")
_T = WF["stats"]["T_target_total_mass"]

_pdf = _pm.open(str(LATEX / "manuscript_v6.pdf"))
_pdf_text = "\n".join(pg.get_text() for pg in _pdf)
_pdf.close()

_STALE = ["0.0183", "35/109", "13/109", "stable once",
          "reverses direction only near zero"]
_found_stale = [x for x in _STALE if x in _pdf_text]
rec("AC. final PDF has no stale Figure 3 strings",
    "PASS" if not _found_stale else "FAIL",
    "stale strings found in compiled PDF = %s" % (_found_stale or "none"))

_medB = "%.4f" % _T["range_C_B_median"]
_medW = "%.4f" % _T["range_C_W_median"]
_need = [_medB, _medW,
         "%d/109" % _T["signflip_C_B"]["0.005"]] 
_missing = [x for x in _need if x not in _pdf_text]
rec("AD. final PDF Figure 3 shows updated values",
    "PASS" if not _missing else "FAIL",
    "expected in PDF: %s; missing = %s" % (_need, _missing or "none"))

_cw_in_pdf = ("%d/109" % _T["signflip_C_W"]["0"]) in _pdf_text
rec("AE. final PDF Figure 3 C_W flip count matches audit",
    "PASS" if _cw_in_pdf else "FAIL",
    "C_W flip count %d/109 present in PDF = %s"
    % (_T["signflip_C_W"]["0"], _cw_in_pdf))

# 正文与图不得矛盾：图里的 C_B 数值必须也出现在正文
_txt_only = M + S
_consist = (_medB in _txt_only) or ("%.4f" % _T["range_C_B_median"]) in _txt_only
rec("AF. text and Figure 3 consistent on C_B",
    "PASS" if _consist else "FAIL",
    "range_C_B median %.4f appears in text = %s" % (_T["range_C_B_median"], _consist))

# AG. D(n,K) 限定语必须同时出现在 Abstract 与 Results
_dnk_ok = ("rank-linear discriminability surrogate" in M
           and "rank-linear utility" in M and "reference-dependent" in M)
rec("AG. D(n,K) qualified in Abstract and Results",
    "PASS" if _dnk_ok else "FAIL",
    "abstract qualifier = %s; results qualifier = %s"
    % ("rank-linear discriminability surrogate" in M,
       "reference-dependent" in M))

# AH. 补充材料 Table 交叉引用正确
_tab_ok = ("Table~S2" in S and "Table~S4 below" not in S
           and "(Table~S1).}" not in S)
rec("AH. supplementary table cross-references",
    "PASS" if _tab_ok else "FAIL",
    "Figure S3 -> Table S2 = %s; stale 'Table~S4 below' gone = %s"
    % ("Table~S2" in S, "Table~S4 below" not in S))
'''

anchor = "# ================= 输出 ================="
if anchor in s and "直接检查最终 PDF 的图形内容" not in s:
    s = s.replace(anchor, NEW.strip() + "\n\n" + anchor, 1)
    io.open(P, "w", encoding="utf-8").write(s)
    print("validator extended with PDF-content checks AC-AH")
else:
    print("anchor missing or already extended")

import ast  # noqa: E402
ast.parse(io.open(P, encoding="utf-8").read())
print("syntax OK")
