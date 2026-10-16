# -*- coding: utf-8 -*-
"""scripts/v8_1/42_external_text_fix.py -- 清理正文的旧版 external-method 残留。

机器可读源 `v8_1/external_clean_profiles.csv` 判定：
    最终 external panel = Flex-ddG, DDGPred, ddG_FoldX, ESM-IF, proteinMPNN
    HERMES 不在其中；五个方法 RMSE 全部有值（缺失 0）

因此正文两处为旧版残留，予以修正：
1. 方法清单 `... ESM-IF[9] and HERMES[4]` -> `... ESM-IF[9] and DDGPred[4]`
   （[4] 原为 HERMES 引用，现改指 DDGPred，条目仍标注 REFERENCE_VERIFY_REQUIRED）
2. 删除 HERMES 专属的 log-odds / 「不报告 RMSE 或 MAE」两句
   （Table 3 五行 RMSE 均有值，该句与表格矛盾）
"""
import io
import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = Path("scripts/v8_1/31_manuscript_v6.py")
s = io.open(P, encoding="utf-8").read()
orig = s

# ---- 1) 方法清单：HERMES -> DDGPred ----
s = s.replace(
    'f"predictions: Flex ddG[6], FoldX[7], ProteinMPNN[8], ESM-IF[9] and HERMES[4]. "',
    'f"predictions: Flex ddG[6], FoldX[7], ProteinMPNN[8], ESM-IF[9] and DDGPred[4]. "')

# ---- 2) 删除 log-odds / 不报告 RMSE 的句子 ----
pat = re.compile(
    r'\s*f"We also note that one "\s*\n\s*f"external quantity is a log-odds difference[\s\S]{0,420}?'
    r'f"MAE for that model\."\)', re.M)
s2, n = pat.subn(')', s)
if n:
    s = s2
    print(f"[2] removed log-odds / no-RMSE sentence (x{n})")
else:
    print("[2] WARNING: log-odds sentence pattern not found")

io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print("changed:", s != orig, "| syntax OK")
print("HERMES occurrences left in script:", s.count("HERMES"))
print("log-odds occurrences left:", s.count("log-odds"))
for i, ln in enumerate(s.split("\n"), 1):
    if "DDGPred[4]" in ln or "HERMES" in ln or "log-odds" in ln:
        print(f"{i:5d}| {ln.strip()[:120]}")
