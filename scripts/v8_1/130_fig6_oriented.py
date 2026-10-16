# -*- coding: utf-8 -*-
"""scripts/v8_1/130_fig6_oriented.py -- Figure 6 改为使用 sign-oriented scores。

依据 v8_1/EXTERNAL_ORIENTATION.json：
  raw      : ESM-IF = -0.2250 pooled / +0.0213 macro
  oriented : ESM-IF = +0.2250 pooled / -0.0213 macro（其余方法不变）

Figure 6 的 (A) 面板用于**方法性能比较**，必须使用同一 performance orientation，
因此用 oriented 值；Panel B 的热图同样改用 oriented，并在 caption 说明。
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = Path(r"F:\dsh\cache\tmp\figuregen\make_all_figures.py")
s = P.read_text(encoding="utf-8")

OLD = '''    methods=["Flex-ddG","DDGPred","FoldX\\n(ddG_FoldX)","ESM-IF","ProteinMPNN"]
    pooled=np.array([0.2028,0.2406,0.2089,-0.2250,0.1509])
    macro=np.array([0.1092,0.1287,0.1088,0.0213,0.0320])
    cb=np.array([0.0478,0.0715,0.0657,-0.1149,0.0699])
    cw=np.array([0.1550,0.1690,0.1433,-0.1101,0.0810])'''

NEW = '''    methods=["Flex-ddG","DDGPred","FoldX\\n(ddG_FoldX)","ESM-IF","ProteinMPNN"]
    # sign-oriented to a common performance orientation (y = ddG, larger = weaker binding).
    # Only ESM-IF needs a sign flip (its published score is larger = stronger binding);
    # the other four columns are already oriented like ddG. See EXTERNAL_ORIENTATION.json.
    pooled=np.array([0.2028,0.2406,0.2089,0.2250,0.1509])
    macro=np.array([0.1092,0.1287,0.1088,-0.0213,0.0320])
    cb=np.array([0.0478,0.0715,0.0657,0.1149,0.0699])
    cw=np.array([0.1550,0.1690,0.1433,0.1101,0.0810])'''

if OLD in s:
    s = s.replace(OLD, NEW, 1)
    print("Figure 6 data -> sign-oriented")
else:
    print("MISS: Figure 6 data anchor")

# 图内说明文字也更新
OLD2 = ('            "Unambiguous annotation universe\\n"\n'
        '            "1,781 potentially mappable rows\\n"\n'
        '            "method-specific coverage = 1,755-1,776\\n"\n'
        '            "rank agreement \\u03c1 = +0.90",')
NEW2 = ('            "Unambiguous annotation universe\\n"\n'
        '            "1,781 potentially mappable rows\\n"\n'
        '            "method-specific coverage = 1,755-1,776\\n"\n'
        '            "sign-oriented scores; rank agreement \\u03c1 = +0.30",')
if OLD2 in s:
    s = s.replace(OLD2, NEW2, 1)
    print("Figure 6 in-figure text -> oriented rho")
else:
    print("MISS: Figure 6 text anchor (checking alternates)")
    s2 = s.replace("rank agreement \\u03c1 = +0.90", "sign-oriented scores; rank agreement \\u03c1 = +0.30")
    if s2 != s:
        s = s2; print("  text updated via fallback")
    else:
        print("  text NOT updated")

P.write_text(s, encoding="utf-8")
print("wrote", P)
