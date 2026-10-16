# -*- coding: utf-8 -*-
"""scripts/v8_1/148_figs_seq40_keyfix.py -- Figure 1 Panel C 与 Figure 2 反映 seq40 键级修正。

新值（来自 FROZEN_RESULTS.H4_controlled_split，9 配置中位）：
    Random  0.6844 | Complex 0.1931 | seq40 0.1961 | Family-OOD 0.1508
原先 Figure 1 Panel C 把 Complex 与 seq40 合并成一张卡（0.1931），
Figure 2 也把 seq40 画成与 complex 相同的点并标注 "identical partition"。
现在四档各自独立呈现。
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = Path(r"F:\dsh\cache\tmp\figuregen\make_all_figures.py")
s = P.read_text(encoding="utf-8")
done, miss = [], []

# ---------- Figure 1 Panel C：三卡 -> 四卡 ----------
OLD_C = '''    split_cards=[
        (0.05,"Random","permutation",0.6844,BLUE,LIGHT_BLUE),
        (0.36,"Complex / seq40","interface / sequence",0.1931,CORAL,LIGHT_CORAL),
        (0.67,"Family-OOD","between-family",0.1508,CORAL,LIGHT_CORAL),
    ]
    for x,title,sub,val,c,fc in split_cards:
        rounded(ax,(x,0.61),(0.26,0.20),fc=fc,ec="#DCE1E9")
        ax.text(x+0.13,0.775,title,ha="center",fontsize=10,fontweight="bold",color=NAVY)
        ax.text(x+0.13,0.743,f"({sub})",ha="center",fontsize=7.2,color=GRAY)
        ax.text(x+0.13,0.66,f"{val:.4f}",ha="center",fontsize=18,fontweight="bold",color=c)
        ax.text(x+0.13,0.625,"pooled Pearson",ha="center",fontsize=7.3,color=NAVY)'''
NEW_C = '''    split_cards=[
        (0.030,"Random","permutation",0.6844,BLUE,LIGHT_BLUE),
        (0.268,"Complex","interface",0.1931,CORAL,LIGHT_CORAL),
        (0.506,"seq40","sequence",0.1961,CORAL,LIGHT_CORAL),
        (0.744,"Family-OOD","between-family",0.1508,CORAL,LIGHT_CORAL),
    ]
    for x,title,sub,val,c,fc in split_cards:
        rounded(ax,(x,0.61),(0.226,0.20),fc=fc,ec="#DCE1E9")
        ax.text(x+0.113,0.775,title,ha="center",fontsize=10,fontweight="bold",color=NAVY)
        ax.text(x+0.113,0.743,f"({sub})",ha="center",fontsize=7.2,color=GRAY)
        ax.text(x+0.113,0.66,f"{val:.4f}",ha="center",fontsize=17,fontweight="bold",color=c)
        ax.text(x+0.113,0.625,"pooled Pearson",ha="center",fontsize=7.0,color=NAVY)'''
if OLD_C in s:
    s = s.replace(OLD_C, NEW_C, 1); done.append("fig1-panelC-4cards")
else:
    miss.append("fig1-panelC-4cards")

# ---------- Figure 2：数据 + 去掉 identical 标注 ----------
OLD_D = '''    pooled=np.array([0.6844,0.1931,0.1931,0.1508])
    cb=np.array([0.4275,0.0233,0.0233,-0.0072])
    cw=np.array([0.2591,0.1685,0.1685,0.1580])'''
NEW_D = '''    pooled=np.array([0.6844,0.1931,0.1961,0.1508])
    cb=np.array([0.4275,0.0233,0.0239,-0.0072])
    cw=np.array([0.2591,0.1685,0.1721,0.1580])'''
if OLD_D in s:
    s = s.replace(OLD_D, NEW_D, 1); done.append("fig2-data")
else:
    miss.append("fig2-data")

OLD_A = ('    ax.axvspan(0.75,2.25,color="#EEF0F3",alpha=0.9,zorder=-2)\n'
         '    ax.text(1.5,0.40,"identical partition\\nin this dataset",ha="center",'
         'fontsize=8,color=GRAY)\n')
if OLD_A in s:
    s = s.replace(OLD_A, "", 1); done.append("fig2-remove-identical-annotation")
else:
    miss.append("fig2-remove-identical-annotation")

P.write_text(s, encoding="utf-8")
print("fixed:", done)
if miss:
    print("MISSED:", miss)
