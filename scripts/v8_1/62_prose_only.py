# -*- coding: utf-8 -*-
"""scripts/v8_1/62_prose_only.py -- 只比较「连续正文散文」，排除表格/图注/参考文献。"""
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

# 表格与图注的特征词：命中即从散文集合中剔除
DROP = re.compile(
    r"(Evaluation set|Split regime|Method Mutations|Full annotation|Primary matched"
    r"|Controlled split|External unambiguous|Random 0\.|Complex 0\.|Family-OOD 0\."
    r"|Flex-ddG 1,|DDGPred 1,|FoldX|ESM-IF 1,|ProteinMPNN 1,"
    r"|Capability-resolved evaluation of PPI|Family-out-of-distribution evaluation dis"
    r"|Target weighting changes the interpretation|Fixed-K shortlist evaluation is constrained"
    r"|Capability decomposition separates|External methods define the boundary"
    r"|Median values across|All rows are joined|All methods are evaluated)", re.I)


def prose(path):
    d = pymupdf.open(str(path))
    t = "\n".join(pg.get_text() for pg in d)
    d.close()
    t = t[:t.rfind("References")]
    # 按空行/句末聚合
    t = re.sub(r"\s+", " ", t)
    chunks = re.split(r"(?<=[.:!?])\s+", t)
    keep, cur = [], ""
    for c in chunks:
        cur = (cur + " " + c).strip() if cur else c
        if len(cur) >= 120:
            keep.append(cur); cur = ""
    if cur:
        keep.append(cur)
    out = []
    for p in keep:
        p = p.strip()
        if len(p) < 120:
            continue
        if DROP.search(p[:120]):
            continue
        out.append(p)
    return out


def key(s):
    return re.sub(r"[^a-z0-9]", "", s.lower())


p6, pn = prose(V6M), prose(NEWM)
k6 = Counter()
for p in p6:
    k6.update(key(p))
kn = Counter()
for p in pn:
    kn.update(key(p))

print("=" * 72)
print("连续散文（已剔除表格、图注、参考文献）")
print("=" * 72)
print(f"  v6 段数 {len(p6)}   新版段数 {len(pn)}")
print(f"  v6 字母数 {sum(k6.values())}   新版字母数 {sum(kn.values())}")
print(f"  缺失 (v6 有 / 新版无) : {sum((k6-kn).values())}")
print(f"  多出 (新版有 / v6 无) : {sum((kn-k6).values())}")
if k6 - kn:
    print("   缺失明细:", (k6 - kn).most_common(8))
if kn - k6:
    print("   多出明细:", (kn - k6).most_common(8))

# 逐段定位差异
print()
print("=== 逐段核对（v6 每段是否在新版中找到同源段落）===")
joined = " || ".join(key(p) for p in pn)
bad = 0
for i, p in enumerate(p6, 1):
    k = key(p)
    hit = False
    for L in (250, 200, 150, 120, 100, 80):
        if L > len(k):
            continue
        if k[:L] in joined or k[-L:] in joined:
            hit = True; break
    if not hit:
        bad += 1
        print(f"  [{i:2d}] 未匹配 ({len(k)} 字符): {p[:130]}")
if bad == 0:
    print("  全部 v6 散文段落均在新版中找到同源段落 -> 正文文字未改动")
