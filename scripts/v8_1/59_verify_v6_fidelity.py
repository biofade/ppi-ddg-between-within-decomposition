# -*- coding: utf-8 -*-
"""scripts/v8_1/59_verify_v6_fidelity.py -- 严格核对 LaTeX 版是否忠实于用户提供的 v6。

核查三项，全部机器判定：
  A. 图片：新 PDF 内嵌图形 vs v6 原件对应图形，逐像素比对
  B. 图注：v6 与新版每个 Figure/Table caption 文本一致
  C. 正文：v6 全文（去掉栏目/页码噪声）与新版全文的字符多重集差异
"""
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pymupdf

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

V6M = Path(r"F:\dsh-home\attachments\v1\files\ca"
           r"\cae30b6694f6f910c1e5e852eac21ffadeaaf91f748ab7d6e176a901ae2ec443"
           r"\CR-PPI-ddG_manuscript_v6_EN.pdf")
V6S = Path(r"F:\dsh-home\attachments\v1\files\e9"
           r"\e93b741cd80f50cfeb2577a0fa2586c9e0a425457d9856f1a3a06ce9500ce0cc"
           r"\CR-PPI-ddG_Supplementary_v6_EN.pdf")
NEWM = Path("latex_v6/manuscript_v6.pdf")
NEWS = Path("latex_v6/supplementary_v6.pdf")

print("=" * 74)
print("A. 图片：v6 内嵌图形 vs 新版内嵌图形（逐像素）")
print("=" * 74)


def page_graphics(doc):
    """返回每页的 (drawing_count, list_of_raster_pixmaps)。"""
    info = []
    for pg in doc:
        pixes = []
        for im in pg.get_images(full=True):
            try:
                pix = pymupdf.Pixmap(doc, im[0])
                if pix.n > 4:
                    pix = pymupdf.Pixmap(pymupdf.csRGB, pix)
                pixes.append(pix)
            except Exception:
                pass
        info.append((len(pg.get_drawings()), pixes))
    return info


d6 = pymupdf.open(str(V6M))
dn = pymupdf.open(str(NEWM))
g6, gn = page_graphics(d6), page_graphics(dn)
print(f"  v6   : 总 drawings = {sum(a for a, _ in g6)}")
print(f"  新版 : 总 drawings = {sum(a for a, _ in gn)}")
r6 = [p for _, ps in g6 for p in ps]
rn = [p for _, ps in gn for p in ps]
print(f"  v6   : 位图图形 {len(r6)} 个   新版: 位图图形 {len(rn)} 个")
for i, p in enumerate(r6[:6]):
    print(f"    v6  #{i+1}  {p.width}x{p.height}")
for i, p in enumerate(rn[:6]):
    print(f"    new #{i+1}  {p.width}x{p.height}")


def pix_bytes(p):
    return p.samples


match_n = 0
for i, p in enumerate(r6):
    if any(pix_bytes(p) == pix_bytes(q) for q in rn):
        match_n += 1
print(f"  位图逐像素一致: {match_n}/{len(r6)}")
d6.close()
dn.close()

print()
print("=" * 74)
print("B. 图注文本一致性")
print("=" * 74)


def captions(path):
    d = pymupdf.open(str(path))
    t = "\n".join(pg.get_text() for pg in d)
    d.close()
    out = {}
    for m in re.finditer(r"((?:Figure|Table)\s*S?\d+)\s*[:|.]\s*(.{0,240}?)"
                         r"(?=(?:Figure|Table)\s*S?\d+\s*[:|.]|\n\n|$)", t, re.S):
        out[m.group(1).replace(" ", "")] = re.sub(r"\s+", " ", m.group(2)).strip()
    return out


c6, cn = captions(V6M), captions(NEWM)
keys = sorted(set(c6) | set(cn))
for k in keys:
    a, b = c6.get(k, ""), cn.get(k, "")
    na, nb = re.sub(r"[^a-z0-9]", "", a.lower()), re.sub(r"[^a-z0-9]", "", b.lower())
    same = (na[:110] == nb[:110]) if na or nb else False
    tag = "一致" if same else ("仅v6有" if not b else ("仅新版有" if not a else "不同"))
    print(f"  {k:10s} {tag}")
    if not same and a and b:
        print(f"      v6 : {a[:110]}")
        print(f"      new: {b[:110]}")

print()
print("=" * 74)
print("C. 正文文本差异（字符多重集）")
print("=" * 74)


def body_counts(path, drop_numbers=True):
    d = pymupdf.open(str(path))
    t = "\n".join(pg.get_text() for pg in d)
    d.close()
    t = re.sub(r"\s+", "", t)
    if drop_numbers:
        t = re.sub(r"[^A-Za-z]", "", t)      # 只比字母，规避编号/上标差异
    return Counter(t.lower())


b6, bn = body_counts(V6M), body_counts(NEWM)
missing = b6 - bn          # v6 有而新版没有
extra = bn - b6            # 新版多出来的
print(f"  v6 字母总数   : {sum(b6.values())}")
print(f"  新版字母总数  : {sum(bn.values())}")
print(f"  缺失字母数    : {sum(missing.values())}")
print(f"  多出字母数    : {sum(extra.values())}")
if missing:
    print("  缺失字符 (top):", missing.most_common(12))
if extra:
    print("  多出字符 (top):", extra.most_common(12))
print()
print("结论: 字母多重集差异越小，说明正文改动越少；",
      "完全一致则为 0。")
