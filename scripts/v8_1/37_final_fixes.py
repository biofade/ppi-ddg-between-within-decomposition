# -*- coding: utf-8 -*-
"""scripts/v8_1/37_final_fixes.py -- 三项收尾修改。

1. **引用改上标**：正文 `[13]` / `[4,5,10]` 等方括号引用 -> 真正的 Word 上标 run。
2. **Figure 4 换图**：安装用户新提供的 figure4（1619x971）。
3. **Figure 1 放大**：Figure 1 原图 1916x821（宽高比 2.33），放在 6.5 in 正文宽下
   渲染高度仅约 2.8 in，字因此极小。改为放在**整页横向（landscape）页面**上，
   可用宽度 9.4 in（+44%），并把插页比例写清楚，使图内文字实际放大。
"""
from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

from PIL import Image

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(__file__).resolve().parent.parent.parent
FIGS = R / "figures_v5"

# ---------------------------------------------------------------- 1) Fig4 换图
NEW4 = Path(r"F:\dsh-home\attachments\v1\files\ff"
            r"\ff942b146670da55db52fddf8b94bf35b42f7e76624618aa43e90d52ed9f7196"
            r"\figure4")
dst4 = FIGS / "Fig4_fixed_K_geometry_EN.png"
if NEW4.exists():
    im = Image.open(NEW4)
    if im.mode not in ("RGB", "RGBA"):
        im = im.convert("RGB")
    im.save(dst4, "PNG")
    print(f"[1] Fig4 replaced -> {im.size[0]}x{im.size[1]}  ({dst4.name})")
else:
    print("[1] WARNING: new figure4 not found; Fig4 left unchanged")
    raise SystemExit(1)

# ------------------------------------------------- 2) 在主稿脚本中启用上标引用
P = R / "scripts" / "v8_1" / "31_manuscript_v6.py"
s = P.read_text(encoding="utf-8")

# 2a) W() 渲染器：把 [n] / [n,m] / [n-m] 解析为上文标 run
old_hdr = '_TOK = re.compile(r"(\\^\\{[^}]*\\}|_\\{[^}]*\\}|_[A-Za-z\\u0177\\u0302])")'
new_hdr = ('_TOK = re.compile(r"(\\^\\{[^}]*\\}|_\\{[^}]*\\}|_[A-Za-z\\u0177\\u0302]'
           '|\\[\\d+(?:[,\\u2013]\\d+)*\\])")')
if old_hdr in s:
    s = s.replace(old_hdr, new_hdr)
    print("[2a] citation token added to tokenizer")
else:
    print("[2a] WARNING: tokenizer pattern not found")

# 2b) toks(): 处理 [..] -> sup（去掉方括号）
old_toks = '''        elif part.startswith("_") and len(part) == 2:
            out.append((part[1], "sub"))
        else:
            out.append((part, ""))'''
new_toks = '''        elif part.startswith("_") and len(part) == 2:
            out.append((part[1], "sub"))
        elif part.startswith("[") and part.endswith("]"):
            inner = part[1:-1].replace("\\u2013", "\\u2013")
            out.append((inner, "sup"))
        else:
            out.append((part, ""))'''
if old_toks in s:
    s = s.replace(old_toks, new_toks)
    print("[2b] toks() now emits superscript citations")
else:
    print("[2b] WARNING: toks() branch not found")

# 2c) 确保 W() 在 sup 分支使用上标字号
s = s.replace(
    '''    for c, k in toks(text):
        r = p.add_run(c)
        if size:
            r.font.size = Pt(size)
        if k == "sub":
            r.font.subscript = True
        elif k == "sup":
            r.font.superscript = True
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    return p''',
    '''    for c, k in toks(text):
        r = p.add_run(c)
        if k == "sup":
            r.font.superscript = True
        elif k == "sub":
            r.font.subscript = True
        if size and k == "":
            r.font.size = Pt(size)
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    return p''')

# ------------------------------------------------- 3) Figure 1 放到整页横向
s = s.replace(
    'from docx.enum.section import WD_ORIENT',
    'from docx.enum.section import WD_ORIENT')

old_fig1 = '    addfig(D, "Fig1_overview_EN.png")'
new_fig1 = '''    # Figure 1 原图宽高比 2.33，放在 6.5 in 正文宽下渲染高度仅约 2.8 in，
    # 图内文字过小。改用整页横向（landscape）安放，可用宽度 9.4 in（+44%）。
    _sec = D.add_section()
    _w, _h = _sec.page_width, _sec.page_height
    _sec.orientation = WD_ORIENT.LANDSCAPE
    _sec.page_width, _sec.page_height = _h, _w
    _sec.left_margin = _sec.right_margin = Inches(0.55)
    _sec.top_margin = _sec.bottom_margin = Inches(0.55)
    D.add_picture(str(FIGS / "Fig1_overview_EN.png"), width=Inches(9.4))
    _sec2 = D.add_section()
    _w2, _h2 = _sec2.page_width, _sec2.page_height
    _sec2.orientation = WD_ORIENT.PORTRAIT
    _sec2.page_width, _sec2.page_height = _h2, _w2
    _sec2.left_margin = _sec2.right_margin = Inches(0.9)
    _sec2.top_margin = _sec2.bottom_margin = Inches(0.9)'''
if old_fig1 in s:
    s = s.replace(old_fig1, new_fig1, 1)
    print("[3] Figure 1 now on a dedicated landscape page at 9.4 in")
else:
    print("[3] WARNING: Fig1 call not found")

P.write_text(s, encoding="utf-8")
import ast  # noqa: E402
ast.parse(s)
print("    manuscript script syntax OK")
