# -*- coding: utf-8 -*-
"""scripts/v8_1/119_validator_fig6.py -- 加入 Figure 6 标签的最终-PDF 硬检查 + 其余第四轮检查。"""
import io
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "validate_final_manuscript.py"
s = io.open(P, encoding="utf-8").read()

NEW = r'''
# ================= 第四轮：Figure 6 标签与外部覆盖 =================
_d = pymupdf.open(str(LATEX / "manuscript_v6.pdf"))
_pdf_main = "\n".join(pg.get_text() for pg in _d)
_d.close()

_F6_BAD = ["Matched unambiguous mapping", "1,781 mutations; 67 targets"]
_F6_GOOD = ["Unambiguous annotation universe", "method-specific coverage", "1,755-1,776"]
_f6_bad = [x for x in _F6_BAD if x in _pdf_main]
_f6_good = [x for x in _F6_GOOD if x in _pdf_main]
rec("AJ. final PDF Figure 6 labels updated",
    "PASS" if (not _f6_bad and len(_f6_good) == len(_F6_GOOD)) else "FAIL",
    "stale labels = %s; required labels present = %d/%d"
    % (_f6_bad or "none", len(_f6_good), len(_F6_GOOD)))

# AK. 图内归因措辞统一为 summed absolute component change
_attr_bad = ("absolute component decline" in _pdf_main)
rec("AK. in-figure attribution wording consistent",
    "PASS" if (not _attr_bad and "summed absolute" in _pdf_main) else "FAIL",
    "stale 'absolute component decline' in PDF = %s; 'summed absolute' present = %s"
    % (_attr_bad, "summed absolute" in _pdf_main))

# AL. external common intersection 已做且有 Note
XI = load_json("EXTERNAL_COMMON_INTERSECTION.json")
_ci_ok = ("Supplementary Note 7" in S and XI["n_intersection_rows"] > 0
          and XI["rank_variable"]["ranking_pooled"] == XI["rank_intersection"]["ranking_pooled"]
          and XI["rank_variable"]["ranking_macro"] == XI["rank_intersection"]["ranking_macro"])
rec("AL. external common-intersection sensitivity performed",
    "PASS" if _ci_ok else "FAIL",
    "Note 7 present = %s; intersection rows = %d; pooled ranking identical = %s; "
    "macro ranking identical = %s"
    % ("Supplementary Note 7" in S, XI["n_intersection_rows"],
       XI["rank_variable"]["ranking_pooled"] == XI["rank_intersection"]["ranking_pooled"],
       XI["rank_variable"]["ranking_macro"] == XI["rank_intersection"]["ranking_macro"]))

# AM. 不得宣称未做的验证
_claim_bad = ("we verified that\nit does not drive the ranking conclusion" in M
              or "de-replicated" in M)
rec("AM. no unverified coverage claim, no 'de-replicated'",
    "PASS" if not _claim_bad else "FAIL",
    "stale claim present = %s; 'de-replicated' present = %s"
    % ("we verified that" in M and "not drive the ranking conclusion" in M,
       "de-replicated" in M))

# AN. Discussion 不得把 surrogate 推成一般规律
_disc_bad = "enlarging a candidate set is not universally effective" in M
_disc_ok = ("Under the canonical rank-linear reference used here" in M
            and "reference-free conclusion is narrower" in M)
rec("AN. Discussion does not generalise the surrogate",
    "PASS" if (not _disc_bad and _disc_ok) else "FAIL",
    "over-general sentence present = %s; qualified wording present = %s"
    % (_disc_bad, _disc_ok))
'''

anchor = "# ================= 输出 ================="
if anchor in s and "Figure 6 标签与外部覆盖" not in s:
    s = s.replace(anchor, NEW.strip() + "\n\n" + anchor, 1)
    io.open(P, "w", encoding="utf-8").write(s)
    print("validator extended with AJ-AN")
else:
    print("anchor missing or already present")

import ast  # noqa: E402
ast.parse(io.open(P, encoding="utf-8").read())
print("syntax OK")
