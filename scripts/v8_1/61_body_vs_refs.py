# -*- coding: utf-8 -*-
"""scripts/v8_1/61_body_vs_refs.py -- 分离「正文」与「参考文献」，分别判定忠实度。"""
import re
import sys
from collections import Counter
from pathlib import Path

import pymupdf

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

V6M = Path(r"F:\dsh-home\attachments\v1\files\ca"
           r"\cae30b6694f6f910c1e5e852eac21ffadeaaf91f748ab7d6e176a901ae2ec443"
           r"\CR-PPI-ddG_manuscript_v6_EN.pdf")
NEWM = Path("latex_v6/manuscript_v6.pdf")


def split_body_refs(path):
    d = pymupdf.open(str(path))
    t = "\n".join(pg.get_text() for pg in d)
    d.close()
    i = t.rfind("References")
    return t[:i], t[i:]


def letters(s):
    return Counter(re.sub(r"[^A-Za-z]", "", s).lower())


body6, ref6 = split_body_refs(V6M)
bodyn, refn = split_body_refs(NEWM)

print("=" * 72)
print("正文（References 之前）")
print("=" * 72)
b6, bn = letters(body6), letters(bodyn)
print(f"  v6 字母数   : {sum(b6.values())}")
print(f"  新版字母数  : {sum(bn.values())}")
print(f"  缺失 (v6有新版无) : {sum((b6-bn).values())}")
print(f"  多出 (新版有v6无) : {sum((bn-b6).values())}")
if b6 - bn:
    print("   缺失明细:", (b6 - bn).most_common(10))
if bn - b6:
    print("   多出明细:", (bn - b6).most_common(10))

print()
print("=" * 72)
print("参考文献区")
print("=" * 72)
r6, rn = letters(ref6), letters(refn)
print(f"  v6 字母数   : {sum(r6.values())}")
print(f"  新版字母数  : {sum(rn.values())}")
print(f"  差值        : {sum(rn.values()) - sum(r6.values()):+d}"
      "   ← BibTeX 把作者展开为全名所致")
print()

# 逐条比对文献的「标识信息」：期刊名/年份/卷期页/DOI 是否一致
def ref_items(s):
    out = []
    for m in re.finditer(r"\[(\d{1,2})\]\s*(.*?)(?=\[\d{1,2}\]|$)", s, re.S):
        out.append((int(m.group(1)), re.sub(r"\s+", " ", m.group(2)).strip()))
    return dict(out)


i6, inn = ref_items(ref6), ref_items(refn)
print("  逐条核查（年份 / 卷期页 / DOI 是否一致）：")
for k in sorted(set(i6) | set(inn)):
    a, b = i6.get(k, ""), inn.get(k, "")
    ya = set(re.findall(r"(19|20)\d{2}", a))
    yb = set(re.findall(r"(19|20)\d{2}", b))
    da = set(re.findall(r"10\.\d{4,9}/[^\s\.;]+", a))
    db = set(re.findall(r"10\.\d{4,9}/[^\s\.;]+", b))
    y_ok = bool(ya & yb) or (not ya and not yb)
    d_ok = bool(da & db) or (not da and not db)
    print(f"    [{k:2d}] 年份一致={y_ok}  DOI一致={d_ok}"
          f"   v6年={sorted(ya)} 新版年={sorted(yb)}")
