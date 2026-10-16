# -*- coding: utf-8 -*-
"""scripts/v8_1/150_fig1_title_fit.py -- 修 Figure 1 Panel A 标题被裁切 + 画布加宽满页。

问题：Panel A 标题 "Provenance-filtered, fingerprint-matched benchmark construction"
      在 Panel B 表头覆盖层下方被裁（lw 覆盖），且画布较窄时标题过长。
处理：1) 缩短 Panel A 标题；2) 画布加宽，使满页时面板更宽、标题有空间。
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = Path(r"F:\dsh\cache\tmp\figuregen\make_all_figures.py")
s = P.read_text(encoding="utf-8")
done, miss = [], []

# 1) 标题缩短（保留含义）
OLD_T = '"Provenance-filtered, fingerprint-matched benchmark construction"'
NEW_T = '"Provenance-filtered benchmark construction"'
if OLD_T in s:
    s = s.replace(OLD_T, NEW_T, 1); done.append("fig1-panelA-title")
else:
    miss.append("fig1-panelA-title")

# 2) 画布加宽（高度略增），使 textwidth 缩放后图更高
OLD_F = "    fig = plt.figure(figsize=(10.6, 12.9))"
NEW_F = "    fig = plt.figure(figsize=(13.0, 13.6))"
if OLD_F in s:
    s = s.replace(OLD_F, NEW_F, 1); done.append("fig1-canvas")
else:
    miss.append("fig1-canvas")

P.write_text(s, encoding="utf-8")
print("fixed:", done)
if miss:
    print("MISSED:", miss)
