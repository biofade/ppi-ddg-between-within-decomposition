# -*- coding: utf-8 -*-
"""scripts/v8_1/60_diff_v6_text.py -- 定位新版相对 v6 多出的文本在哪里。"""
import re
import sys
from pathlib import Path

import pymupdf

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

V6M = Path(r"F:\dsh-home\attachments\v1\files\ca"
           r"\cae30b6694f6f910c1e5e852eac21ffadeaaf91f748ab7d6e176a901ae2ec443"
           r"\CR-PPI-ddG_manuscript_v6_EN.pdf")
NEWM = Path("latex_v6/manuscript_v6.pdf")


def paras(path):
    d = pymupdf.open(str(path))
    t = "\n".join(pg.get_text() for pg in d)
    d.close()
    t = t.replace("\u00ad", "")
    # 按行聚合为段：行尾无句读则视为续行
    lines = [l.strip() for l in t.split("\n")]
    out, cur = [], ""
    for l in lines:
        if not l:
            if cur:
                out.append(cur); cur = ""
            continue
        cur = (cur + " " + l).strip() if cur else l
        if l.endswith((".", ":", "?", "!", "。")) and len(cur) > 80:
            out.append(cur); cur = ""
    if cur:
        out.append(cur)
    return [p for p in out if len(p) > 40]


def key(s):
    return re.sub(r"[^a-z0-9]", "", s.lower())


p6, pn = paras(V6M), paras(NEWM)
kn = [key(p) for p in pn]
kn_all = "".join(kn)

print(f"v6 段落 {len(p6)} 个，新版段落 {len(pn)} 个\n")
print("=== 新版中「v6 没有」的段落 ===")
n_new = 0
for p in pn:
    k = key(p)
    if len(k) < 40:
        continue
    # 若该段有 >=120 字符的连续子串在 v6 中出现，视为同源
    hit = False
    for L in (200, 150, 120, 90, 60):
        for i in range(0, max(1, len(k) - L), 40):
            if k[i:i + L] in "".join(key(x) for x in p6):
                hit = True; break
        if hit:
            break
    if not hit and len(k) > 60:
        n_new += 1
        print(f"  [{len(p):4d}] {p[:150]}")
if n_new == 0:
    print("  （无：新版没有引入 v6 之外的新段落）")

print()
print("=== 参考文献区对比 ===")
for name, path in (("v6", V6M), ("new", NEWM)):
    d = pymupdf.open(str(path))
    t = "\n".join(pg.get_text() for pg in d)
    d.close()
    i = t.rfind("References")
    seg = t[i:]
    print(f"  {name}: 文献区字符数 = {len(re.sub(r'[^A-Za-z]','',seg))}")
    nums = re.findall(r"(?m)^\s*\[?(\d{1,2})[\].]", seg)[:16]
    print(f"       条目起始数字: {nums}")
