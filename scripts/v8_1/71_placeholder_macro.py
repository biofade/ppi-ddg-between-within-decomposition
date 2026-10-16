# -*- coding: utf-8 -*-
"""scripts/v8_1/71_placeholder_macro.py -- 用占位符代替易错的反斜杠。

问题：SI 模板里 figfit 宏的反斜杠数量反复出错（需要 2 个，实际写成 6 个），
      而 read 工具显示与真实字节不一致，肉眼校对不可靠。

做法：源码里只写 %%FIGFIT_DEF%% / %%FIGFIT_CALL%% 这类纯 ASCII 占位符，
      反斜杠只在生成脚本里以显式字节构造，彻底绕开计数问题。
"""
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "scripts/v8_1/53_build_latex_v6.py"
raw = open(P, "rb").read()
B = 0x5C  # 反斜杠
BS = bytes([B, B])          # 模板里需要 2 个 -> 输出 1 个
BS4 = bytes([B, B, B, B])   # 模板里需要 4 个 -> 输出 2 个

# 1) 把 SI 模板里那两行（6 个/4 个反斜杠）整行换成占位符
lines = raw.split(b"\n")
fixed = 0
for i, ln in enumerate(lines):
    if ln.count(B) >= 4 and b"figfit" in ln:
        lines[i] = b"%%FIGFIT_DEF%%"
        fixed += 1
    elif ln.count(B) >= 4 and b"includegraphics" in ln and b"keepaspectratio" in ln:
        lines[i] = b"%%FIGFIT_DEF2%%"
        fixed += 1
raw = b"\n".join(lines)
print(f"replaced {fixed} macro line(s) with placeholders")

# 2) 生成时做替换：占位符 -> 正确字节
#    在 .replace("__TS1__", ...) 那一串之后追加
MARK = b'SI = SI.replace("__TS1__", tS1)'
if MARK in raw:
    raw = raw.replace(
        MARK,
        b'SI = SI.replace("%%FIGFIT_DEF%%", "\\\\newcommand{\\\\figfit}[2][\\\\textwidth]{%")\n'
        b'SI = SI.replace("%%FIGFIT_DEF2%%", "  \\\\includegraphics[width=#1,height=0.82'
        b'\\\\textheight,keepaspectratio]{#2}}")\n'
        + MARK, 1)
    print("substitution hook added")
else:
    print("WARNING: hook anchor not found")

open(P, "wb").write(raw)
import ast  # noqa: E402
ast.parse(open(P, encoding="utf-8").read())
print("syntax OK")

# 3) 报告 SI 模板里剩余的 figfit 行
for i, ln in enumerate(open(P, "rb").read().split(b"\n"), 1):
    if b"figfit" in ln:
        print(f"  L{i} {ln[:52]}  bs={ln.count(B)}")
