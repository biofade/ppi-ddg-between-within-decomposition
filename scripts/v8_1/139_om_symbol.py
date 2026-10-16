# -*- coding: utf-8 -*-
"""scripts/v8_1/139_om_symbol.py -- orientation multiplier mu_c -> o_m（避免与 decomposition 的 mu_c 冲突）。

注意：只改 **方向乘子** 的符号，不动 decomposition 里表示 target 均值的 mu_c。
decomposition 的用法出现在 y_ci = mu_c + eps_ci 一带，方向乘子出现在
s^{aligned} = o_m s^{raw} 与 Table S3 表头。
"""
import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
done, miss = [], []

# ---------- 补充材料 ----------
p = Path("latex_v6/supplementary_v6.tex")
s = p.read_text(encoding="utf-8")

REPL = [
    (r"sign-only orientation, \$s\^\{\\mathrm\{aligned\}\} = \\mu_c\\, s\^\{\\mathrm\{raw\}\}",
     r"sign-only orientation, $s^{\mathrm{aligned}} = o_m\, s^{\mathrm{raw}}"),
    (r"\$\\mu_\{\\mathrm\{ESM\\text\{-\}IF\}\} = -1\$ and \$\\mu_c = \+1\$ for the other four",
     r"$o_{\mathrm{ESM\text{-}IF}} = -1$ and $o_m = +1$ for the other four"),
    (r"orientation and \$\\mu_c\$, so the unmodified provenance and the aligned",
     r"orientation and $o_m$, so the unmodified provenance and the aligned"),
    (r"Table~S4 records each raw column with its raw\\norientation and \$\\mu_c\$",
     r"Table~S3 records each raw column with its raw\\norientation and $o_m$"),
    (r"Method & Raw column & Units & Raw orientation & \$\\mu_c\$ & RMSE meaning",
     r"Method & Raw column & Units & Raw orientation & $o_m$ & RMSE meaning"),
]
for pat, rep in REPL:
    s2, n = re.subn(pat, lambda m: rep, s, count=1)
    if n:
        s = s2; done.append("supp:" + pat[:40])
    else:
        miss.append("supp:" + pat[:40])

# Table S3 数据行的 +1/-1 列保持（已是 $+1$/$-1$）
p.write_text(s, encoding="utf-8")

# ---------- 主稿 ----------
mp = Path("latex_v6/manuscript_v6.tex")
m = mp.read_text(encoding="utf-8")
MR = [
    (r"s_m\^\{\\mathrm\{aligned\}\} = o_m \\, s_m\^\{\\mathrm\{raw\}\}",
     r"s_m^{\mathrm{aligned}} = o_m \, s_m^{\mathrm{raw}}"),
    (r"Supplementary Table~S3, which also records \$o_m\$",
     r"Supplementary Table~S3, which also records $o_m$"),
]
for pat, rep in MR:
    m2, n = re.subn(pat, lambda mm: rep, m, count=1)
    if n:
        m = m2; done.append("main:" + pat[:36])
    else:
        miss.append("main:" + pat[:36])
mp.write_text(m, encoding="utf-8")

print("renamed:", len(done))
if miss:
    print("MISSED:", miss)
