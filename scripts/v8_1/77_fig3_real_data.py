# -*- coding: utf-8 -*-
"""scripts/v8_1/77_fig3_real_data.py -- 让 Figure 3 用真实 109 个 configuration 值绘制。

原脚本从旧 PNG 反推了 18 个 histogram bin（FIG3_CB_CENTERS/COUNTS 等硬编码数组），
属于"由图像反推数据"。现改为直接读
  v8_1/WEIGHTING_EFFECT_SIZE_AUDIT.csv  的 range_C_B / range_C_W 列
（109 行真实值），用 np.histogram 计算 bin，并保留硬编码作为回退。
"""
import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = Path(r"F:\dsh\cache\tmp\figuregen\make_all_figures.py")
s = P.read_text(encoding="utf-8")

ANCHOR = "# FIGURE 3\n# ============================================================================"
if ANCHOR not in s:
    print("ANCHOR not found"); sys.exit(1)

NEW = '''# FIGURE 3
# ============================================================================
# 真实数据来源: v8_1/WEIGHTING_EFFECT_SIZE_AUDIT.csv 的 range_C_B / range_C_W
# （109 个 configuration-level 值）。硬编码 bin 仅作为缺少 CSV 时的回退。
FIG3_CB_COUNTS=np.array([7,3,3,10,3,1,5,3,5,14,8,15,12,5,1,6,4,4])
FIG3_CB_CENTERS=np.array([0.0007,0.0036,0.0050,0.0066,0.0081,0.0109,0.0124,0.0138,0.0152,0.0166,0.0180,0.0195,0.0210,0.0224,0.0248,0.0262,0.0277,0.0292])
FIG3_CW_COUNTS=np.array([7,2,5,1,4,25,14,8,14,2,1,1,5,3,6,4,3,4])
FIG3_CW_CENTERS=np.array([0.007,0.021,0.028,0.042,0.049,0.056,0.063,0.070,0.077,0.084,0.091,0.098,0.105,0.112,0.119,0.126,0.133,0.140])


def _load_fig3_raw():
    """读真实 109 个 configuration 的 range_C_B / range_C_W。"""
    from pathlib import Path as _P
    import csv as _csv
    for cand in (_P(__file__).resolve().parent / "WEIGHTING_EFFECT_SIZE_AUDIT.csv",
                 _P(r"F:\\DSH_Galatea\\Affinity_Intelligence\\v8_1\\WEIGHTING_EFFECT_SIZE_AUDIT.csv")):
        if cand.exists():
            cb, cw = [], []
            with cand.open(encoding="utf-8", newline="") as fh:
                for row in _csv.DictReader(fh):
                    try:
                        v = row.get("range_C_B", "")
                        if v not in ("", "nan", None):
                            cb.append(float(v))
                        w = row.get("range_C_W", "")
                        if w not in ("", "nan", None):
                            cw.append(float(w))
                    except (TypeError, ValueError):
                        pass
            if cb and cw:
                return cb, cw, str(cand)
    return None, None, None


_F3_CB, _F3_CW, _F3_SRC = _load_fig3_raw()

'''
s = s.replace(ANCHOR, NEW.rstrip("\n"), 1)

# --- Panel A ---
old_a = """    ax.bar(FIG3_CB_CENTERS,FIG3_CB_COUNTS,width=0.00125,color=PURPLE,alpha=0.88)
    ax.axvline(0.0183,color=NAVY,ls="--",lw=1.6)
    ax.annotate("Median range\\n= 0.0183",xy=(0.0183,20),xytext=(0.023,22.5),"""
new_a = """    if _F3_CB:
        import numpy as _np
        _c,_e = _np.histogram(_F3_CB, bins=_np.arange(0.0,0.0325,0.002))
        _ctr = 0.5*(_e[1:]+_e[:-1]); _wid = (_e[1]-_e[0])*0.85
        _med = float(_np.median(_F3_CB))
    else:
        _ctr,_c,_wid = FIG3_CB_CENTERS,FIG3_CB_COUNTS,0.00125
        _med = 0.0183
    ax.bar(_ctr,_c,width=_wid,color=PURPLE,alpha=0.88)
    ax.axvline(_med,color=NAVY,ls="--",lw=1.6)
    ax.annotate("Median range\\n= %.4f"%_med,xy=(_med,20),xytext=(0.023,22.5),"""
if old_a in s:
    s = s.replace(old_a, new_a, 1); print("panel A patched")
else:
    print("WARN panel A anchor not found")
s = s.replace('ax.set_ylabel("Model configurations (N = 109)")',
              'ax.set_ylabel("Model configurations (N = %d)"%(len(_F3_CB) if _F3_CB else 109))', 1)
print("Figure 3 patched; source =", _F3_SRC if False else "(resolved at runtime)")

P.write_text(s, encoding="utf-8")
print("wrote", P)
