# -*- coding: utf-8 -*-
"""scripts/v8_1/112_validator_figcheck.py -- 把「最终 PDF 必须内嵌最新图形」加入 validator。"""
import io
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "validate_final_manuscript.py"
s = io.open(P, encoding="utf-8").read()

NEW = r'''
# ================= 最终 PDF 必须内嵌最新图形（AI 检查） =================
_FIGMAP = {
    "Fig1_overview_EN": 0, "Fig2_controlled_split_collapse_EN": 0,
    "Fig3_target_weighting_sensitivity_EN": 0, "Fig4_fixed_K_geometry_EN": 0,
    "Fig5_retrospective_remedy_EN": 0, "Fig6_external_methods_EN": 0,
    "FigS1_matched_universe_strata_EN": 1, "FigS2_capability_heatmap_provenance_clean_EN": 1,
    "FigS3_pool_sensitivity_rank_agreement_EN": 1,
}
_blob = {0: re.sub(r"\s+", " ", M).lower(), 1: re.sub(r"\s+", " ", S).lower()}
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
    _hit = sum(1 for _f in _fr if _f in _blob[_w])
    if _hit / len(_fr) < 0.6:
        _fig_bad.append("%s(%d/%d)" % (_n, _hit, len(_fr)))
rec("AI. final PDF embeds current figure sources",
    "PASS" if not _fig_bad else "FAIL",
    "9 figures compared against figures_v5 vector sources by text fingerprint; "
    "failures = %s" % (_fig_bad or "none"))
'''

anchor = "# ================= 输出 ================="
if anchor in s and "最终 PDF 必须内嵌最新图形" not in s:
    s = s.replace(anchor, NEW.strip() + "\n\n" + anchor, 1)
    io.open(P, "w", encoding="utf-8").write(s)
    print("validator extended with check AI (figure embedding)")
else:
    print("anchor missing or already present")

import ast  # noqa: E402
ast.parse(io.open(P, encoding="utf-8").read())
print("syntax OK")
