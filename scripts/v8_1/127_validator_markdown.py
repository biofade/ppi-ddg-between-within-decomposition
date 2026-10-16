# -*- coding: utf-8 -*-
"""scripts/v8_1/127_validator_markdown.py -- 加入 Markdown 残留 + 符号约定的最终-PDF 检查。"""
import io
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "validate_final_manuscript.py"
s = io.open(P, encoding="utf-8").read()

NEW = r'''
# ================= 第六轮：Markdown 残留 / 符号约定 =================
_dm = pymupdf.open(str(LATEX / "manuscript_v6.pdf"))
_tm = "\n".join(pg.get_text() for pg in _dm)
_dm.close()
_ds = pymupdf.open(str(LATEX / "supplementary_v6.pdf"))
_ts = "\n".join(pg.get_text() for pg in _ds)
_ds.close()

# AS. 最终 PDF 不得有 Markdown 强调残留
import re as _re2
_md_bold = _re2.findall(r"\*\*[^*\n]{1,60}\*\*", _tm + "\n" + _ts)
_md_ital = _re2.findall(r"(?<!\*)\*[^*\n]{1,60}\*(?!\*)", _tm + "\n" + _ts)
rec("AS. no Markdown emphasis residue in final PDFs",
    "PASS" if (not _md_bold and not _md_ital) else "FAIL",
    "bold residue = %s; italic residue = %s"
    % (_md_bold[:3] or "none", _md_ital[:3] or "none"))

# AT. ΔΔG 符号约定必须明确定义
_sign_ok = ("weaker binding" in _tm
            and ("DeltaDeltaG" in _tm or "\u0394\u0394G" in _tm)
            and "no sign inversion" in _tm)
rec("AT. DeltaDeltaG sign convention defined in Methods",
    "PASS" if _sign_ok else "FAIL",
    "sign convention stated = %s; 'no sign inversion' present = %s"
    % ("weaker binding" in _tm, "no sign inversion" in _tm))

# AU. positive affine rescaling 措辞
rec("AU. 'positive affine rescaling' wording",
    "PASS" if ("positive affine rescaling" in _tm
               and "invariant to affine rescaling" not in _tm) else "FAIL",
    "positive-affine wording = %s; unqualified wording gone = %s"
    % ("positive affine rescaling" in _tm,
       "invariant to affine rescaling" not in _tm))

# AV. |DeltaDeltaG|-like 不得出现（绝对值会丢方向）
rec("AV. no absolute-value DDG notation for signed predictors",
    "PASS" if ("|\u0394\u0394G|-like" not in _ts
               and "|\\Delta\\Delta G|-like" not in S) else "FAIL",
    "absolute-value notation present = %s"
    % ("|\u0394\u0394G|-like" in _ts))

# AW. external 来源明确为 CATH-ddG skempi_v2.csv 且无 "official"
rec("AW. external source stated precisely",
    "PASS" if ("CATH-ddG" in _tm and "official per-mutation" not in _tm
               and "skempi" in _ts.lower()) else "FAIL",
    "CATH-ddG named = %s; 'official per-mutation' removed = %s; "
    "skempi_v2.csv named in supplement = %s"
    % ("CATH-ddG" in _tm, "official per-mutation" not in _tm,
       "skempi" in _ts.lower()))

# AX. Table S3 caption 不得自重复
rec("AX. Table S3 caption not self-duplicated",
    "PASS" if "columns (Table S3)." not in _ts else "FAIL",
    "self-duplicating caption present = %s" % ("columns (Table S3)." in _ts))
'''

anchor = "# ================= 输出 ================="
if anchor in s and "Markdown 残留 / 符号约定" not in s:
    s = s.replace(anchor, NEW.strip() + "\n\n" + anchor, 1)
    io.open(P, "w", encoding="utf-8").write(s)
    print("validator extended with AS-AX")
else:
    print("anchor missing or already present")

import ast  # noqa: E402
ast.parse(io.open(P, encoding="utf-8").read())
print("syntax OK")
