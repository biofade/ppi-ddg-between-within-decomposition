# -*- coding: utf-8 -*-
"""scripts/v8_1/113_fix_figcheck.py -- 让 AI 检查真正读取编译后的 PDF（而非 .tex 源）。"""
import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "validate_final_manuscript.py"
s = io.open(P, encoding="utf-8").read()

# 替换 AI 检查块
m = re.search(r"# =+ 最终 PDF 必须内嵌最新图形.*?rec\(\"AI\..*?\n\n", s, re.S)
if not m:
    print("AI block not found"); sys.exit(1)

NEW = r'''# ================= 最终 PDF 必须内嵌最新图形（AI 检查） =================
_FIGMAP = {
    "Fig1_overview_EN": "manuscript", "Fig2_controlled_split_collapse_EN": "manuscript",
    "Fig3_target_weighting_sensitivity_EN": "manuscript",
    "Fig4_fixed_K_geometry_EN": "manuscript", "Fig5_retrospective_remedy_EN": "manuscript",
    "Fig6_external_methods_EN": "manuscript",
    "FigS1_matched_universe_strata_EN": "supplementary",
    "FigS2_capability_heatmap_provenance_clean_EN": "supplementary",
    "FigS3_pool_sensitivity_rank_agreement_EN": "supplementary",
}
_compiled = {}
for _k, _fn in (("manuscript", "manuscript_v6.pdf"),
                ("supplementary", "supplementary_v6.pdf")):
    _d = pymupdf.open(str(LATEX / _fn))
    _compiled[_k] = re.sub(r"\s+", " ", "\n".join(pg.get_text() for pg in _d)).lower()
    _d.close()

_fig_bad = []
for _n, _w in _FIGMAP.items():
    _srcp = ROOT / "figures_v5" / (_n + ".pdf")
    if not _srcp.exists():
        _fig_bad.append(_n + ":src-missing"); continue
    _t = re.sub(r"\s+", " ", pymupdf.open(str(_srcp))[0].get_text()).lower()
    _fr, _seen = [], set()
    for _m in re.finditer(r"[a-z][a-z0-9 \.\-\+%]{7,60}", _t):
        _f = _m.group(0).strip()
        if len(_f) >= 10 and not re.fullmatch(r"[\d\.\s%\-]+", _f) and _f not in _seen:
            _seen.add(_f); _fr.append(_f)
    _fr = _fr[:12]
    if not _fr:
        continue
    _hit = sum(1 for _f in _fr if _f in _compiled[_w])
    if _hit / len(_fr) < 0.6:
        _fig_bad.append("%s(%d/%d)" % (_n, _hit, len(_fr)))
rec("AI. final PDF embeds current figure sources",
    "PASS" if not _fig_bad else "FAIL",
    "9 figures compared against figures_v5 vector sources using the compiled PDF text; "
    "failures = %s" % (_fig_bad or "none"))

'''
s = s[:m.start()] + NEW + s[m.end():]
io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print("AI check now reads compiled PDFs; syntax OK")
