# -*- coding: utf-8 -*-
"""scripts/v8_1/30_fig_tune_final.py -- 消除 Fig1 与 Fig4 的剩余文字碰撞（只调坐标/字号）。"""
import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "generate_publication_figures_v5.py"
s = io.open(P, encoding="utf-8").read()
n = 0


def rep(a, b):
    global s, n
    if a in s:
        s = s.replace(a, b)
        n += 1
    else:
        print("MISS:", a[:78].replace("\n", " "))


# ---- Fig1 Panel D: 雷达图整体下移，标签 pad 减小，图例更靠下 ----
rep("ax_polar = fig.add_axes([0.625, 0.225, 0.215, 0.160], polar=True)",
    "ax_polar = fig.add_axes([0.640, 0.190, 0.200, 0.145], polar=True)")
rep("ax_polar.tick_params(pad=11)", "ax_polar.tick_params(pad=7)")
rep("ax_polar.set_xticklabels(categories, fontsize=8.8",
    "ax_polar.set_xticklabels(categories, fontsize=8.2")
rep("bbox_to_anchor=(0.5, -0.24), ncol=1,", "bbox_to_anchor=(0.5, -0.30), ncol=1,")

# ---- Fig1 Panel D: 说明文字下移，避开图例 ----
rep('ax_d.text(0.50, 0.048, "Axes standardised',
    'ax_d.text(0.50, 0.030, "Axes standardised')

# ---- Fig1 Panel A: 右侧标签字号减小，避免超出目标卡 ----
rep('ax_a.text(tx + 0.054, mu_y, "$\\\\mu_c$", ha=\'left\', va=\'center\', fontsize=9.5,',
    'ax_a.text(tx + 0.046, mu_y, "$\\\\mu_c$", ha=\'left\', va=\'center\', fontsize=8.2,')
rep("fontsize=8.5, fontweight='bold', color=PALETTE['cw_green'])",
    "fontsize=7.8, fontweight='bold', color=PALETTE['cw_green'])")

# ---- Fig1 Panel B: 棱镜图标右移，给说明文字让位 ----
rep("prism = Polygon([(0.11, 0.71), (0.17, 0.71), (0.14, 0.82)], closed=True,",
    "prism = Polygon([(0.075, 0.715), (0.128, 0.715), (0.10, 0.812)], closed=True,")
rep("ax_b.plot([0.06, 0.125], [0.765, 0.765], color='#1F2937', lw=2.2, zorder=5)",
    "ax_b.plot([0.035, 0.088], [0.765, 0.765], color='#1F2937', lw=2.2, zorder=5)")
rep("ax_b.plot([0.155, 0.21], [0.785, 0.81], color=PALETTE['cb_purple'], lw=2.2, zorder=5)",
    "ax_b.plot([0.115, 0.165], [0.785, 0.81], color=PALETTE['cb_purple'], lw=2.2, zorder=5)")
rep("ax_b.plot([0.155, 0.21], [0.745, 0.72], color=PALETTE['cw_green'], lw=2.2, zorder=5)",
    "ax_b.plot([0.115, 0.165], [0.745, 0.72], color=PALETTE['cw_green'], lw=2.2, zorder=5)")
rep("ax_b.text(0.60, 0.722,", "ax_b.text(0.63, 0.722,")

# ---- Fig1 Panel C: 柱状图 x 轴标签与审计条分离 ----
rep("ax_c_bar = fig.add_axes([0.085, 0.175, 0.255, 0.185])",
    "ax_c_bar = fig.add_axes([0.085, 0.205, 0.255, 0.160])")
rep('ax_c_bar.set_xlabel("Covariance contribution ($C_{\\\\mathrm{B}} + C_{\\\\mathrm{W}}$)",\n'
    '                    fontsize=10, fontweight=\'bold\')',
    'ax_c_bar.set_xlabel("Covariance contribution", fontsize=9.2, fontweight="bold")')

# ---- Fig4: 内嵌图例与强调框分离 ----
rep("ax.legend(loc='upper right', frameon=True, facecolor='white', fontsize=10.5)",
    "ax.legend(loc='upper right', frameon=True, facecolor='white', fontsize=9.6)")
rep("xytext=(33, 0.228),", "xytext=(34, 0.238),")

io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print(f"applied {n} edits; syntax OK")
