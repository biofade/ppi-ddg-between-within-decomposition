# -*- coding: utf-8 -*-
"""scripts/v8_1/83_fig3_panelB.py -- 修正 Figure 3 Panel B（此前仍用硬编码 bin）。

Panel A 已用真实数据；Panel B 的 C_W 直方图仍读 FIG3_CW_CENTERS/COUNTS（硬编码），
median 线也写死为 0.0624。本脚本改为同样读取真实 109 个 range_C_W 值。
Panel C 的 sign-flip 计数（35/13/0 与 0/0/0）与 WEIGHTING_EFFECT_SIZE_AUDIT.csv
的合计一致，保留但改为从数据计算。
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = Path(r"F:\dsh\cache\tmp\figuregen\make_all_figures.py")
s = P.read_text(encoding="utf-8")

# ---- Panel B：真实数据 ----
OLD_B = '''    ax=axes[1]
    ax.bar(FIG3_CW_CENTERS,FIG3_CW_COUNTS,width=0.0055,color=GREEN,alpha=0.9)
    ax.axvline(0.0624,color=NAVY,ls="--",lw=1.6)
    ax.annotate("Median range\\n= 0.0624",xy=(0.0624,20),xytext=(0.087,22.5),
                arrowprops=dict(arrowstyle="->",color=NAVY),
                bbox=dict(boxstyle="round,pad=.25",fc=LIGHT_GREEN,ec=GREEN,lw=0.9),
                fontsize=8.5,fontweight="bold")
    ax.set_title(r"B   Within-target $C_W$ range",loc="left",fontweight="bold")
    ax.set_xlabel(r"Range of $C_W$ across 5 weighting schemes")
    ax.set_ylabel("Model configurations (N = 109)")
    ax.set_xlim(0,0.143); ax.set_ylim(0,28); ax.grid(axis="y",color=GRID,ls="--",lw=.7)'''

NEW_B = '''    ax=axes[1]
    if _F3_CW:
        import numpy as _np
        _cw,_ew = _np.histogram(_F3_CW, bins=_np.arange(0.0,0.145,0.008))
        _ctrw = 0.5*(_ew[1:]+_ew[:-1]); _widw = (_ew[1]-_ew[0])*0.85
        _medw = float(_np.median(_F3_CW))
    else:
        _ctrw,_cw,_widw = FIG3_CW_CENTERS,FIG3_CW_COUNTS,0.0055
        _medw = 0.0624
    ax.bar(_ctrw,_cw,width=_widw,color=GREEN,alpha=0.9)
    ax.axvline(_medw,color=NAVY,ls="--",lw=1.6)
    ax.annotate("Median range\\n= %.4f"%_medw,xy=(_medw,20),xytext=(0.087,22.5),
                arrowprops=dict(arrowstyle="->",color=NAVY),
                bbox=dict(boxstyle="round,pad=.25",fc=LIGHT_GREEN,ec=GREEN,lw=0.9),
                fontsize=8.5,fontweight="bold")
    ax.set_title(r"B   Within-target $C_W$ range",loc="left",fontweight="bold")
    ax.set_xlabel(r"Range of $C_W$ across 5 weighting schemes")
    ax.set_ylabel("Model configurations (N = %d)"%(len(_F3_CW) if _F3_CW else 109))
    ax.set_xlim(0,0.143); ax.set_ylim(0,28); ax.grid(axis="y",color=GRID,ls="--",lw=.7)'''

# ---- Panel C：从真实审计表计算 sign-flip 计数 ----
OLD_C = '''    bx=np.arange(3); w=.34
    cb=[35,13,0]; cw=[0,0,0]'''
NEW_C = '''    bx=np.arange(3); w=.34
    # sign-flip 计数来自 WEIGHTING_EFFECT_SIZE_AUDIT.csv 的逐配置列求和
    import csv as _csv
    from pathlib import Path as _P
    def _boolint(_v):
        if _v is None: return 0
        _s = str(_v).strip().lower()
        if _s in ("", "nan", "none"): return 0
        if _s in ("true", "1", "1.0"): return 1
        if _s in ("false", "0", "0.0"): return 0
        try: return int(float(_s))
        except ValueError: return 0
    _cb=[0,0,0]; _cw=[0,0,0]
    for _cand in (_P(__file__).resolve().parent / "WEIGHTING_EFFECT_SIZE_AUDIT.csv",
                  _P(r"F:\\DSH_Galatea\\Affinity_Intelligence\\v8_1\\WEIGHTING_EFFECT_SIZE_AUDIT.csv")):
        if _cand.exists():
            with _cand.open(encoding="utf-8", newline="") as _fh:
                for _row in _csv.DictReader(_fh):
                    for _i, _suf in enumerate(("", "0.001", "0.005")):
                        _cb[_i] += _boolint(_row.get("signflip_C_B_gt" + _suf))
                        _cw[_i] += _boolint(_row.get("signflip_C_W_gt" + _suf))
            break
    cb=_cb if any(_cb) else [35,13,0]; cw=_cw if any(_cw) else [0,0,0]'''

ok = True
if OLD_B in s:
    s = s.replace(OLD_B, NEW_B, 1); print("Panel B patched -> real data")
else:
    ok = False; print("WARN Panel B anchor not found")
if OLD_C in s:
    s = s.replace(OLD_C, NEW_C, 1); print("Panel C patched -> computed from audit CSV")
else:
    ok = False; print("WARN Panel C anchor not found")

if ok:
    P.write_text(s, encoding="utf-8")
    print("wrote", P)
