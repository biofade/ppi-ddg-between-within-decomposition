# -*- coding: utf-8 -*-
"""scripts/v8_1/40_fig1_layout.py -- 把 Figure 1 的安放方式改为横向页面 + 2 列并排。"""
import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "scripts/v8_1/31_manuscript_v6.py"
s = io.open(P, encoding="utf-8").read()

start = s.find("    # Figure 1 原图宽高比")
end = s.find("    cap(D, \"Figure 1 |")
if start == -1 or end == -1:
    print("markers not found", start, end)
    raise SystemExit(1)

new = '''    # ---- Figure 1 放大安放 ----------------------------------------------
    # 原图 1916x821（宽高比 2.33）。放 6.5 in 正文宽时高度仅 2.79 in，图内文字过小。
    # 改为：整页横向 + 上下两半左右并排（每半 4.55 in 宽），
    # 使单位面积文字相对单张整图放大约 1.4 倍。
    _sec = D.add_section()
    _w, _h = _sec.page_width, _sec.page_height
    _sec.orientation = WD_ORIENT.LANDSCAPE
    _sec.page_width, _sec.page_height = _h, _w
    _sec.left_margin = _sec.right_margin = Inches(0.45)
    _sec.top_margin = _sec.bottom_margin = Inches(0.5)
    _tb = D.add_table(rows=1, cols=2)
    _tb.autofit = False
    _tb.alignment = WD_TABLE_ALIGNMENT.CENTER
    for _i, _f in enumerate(("Fig1_half1_EN.png", "Fig1_half2_EN.png")):
        _cell = _tb.rows[0].cells[_i]
        _cell.width = Inches(4.95)
        _para = _cell.paragraphs[0]
        _para.alignment = WD_ALIGN_PARAGRAPH.CENTER
        _run = _para.add_run()
        _run.add_picture(str(FIGS / _f), width=Inches(4.85))
    _sec2 = D.add_section()
    _w2, _h2 = _sec2.page_width, _sec2.page_height
    _sec2.orientation = WD_ORIENT.PORTRAIT
    _sec2.page_width, _sec2.page_height = _h2, _w2
    _sec2.left_margin = _sec2.right_margin = Inches(0.9)
    _sec2.top_margin = _sec2.bottom_margin = Inches(0.9)
'''

s = s[:start] + new + s[end:]
io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print("Fig1 layout replaced; syntax OK")
print("TABLE_ALIGNMENT imported:", "WD_TABLE_ALIGNMENT" in s.split("\n")[0:60].__str__()
      or "WD_TABLE_ALIGNMENT" in s)
