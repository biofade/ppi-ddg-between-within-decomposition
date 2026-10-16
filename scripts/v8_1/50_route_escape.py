# -*- coding: utf-8 -*-
"""scripts/v8_1/50_route_escape.py -- 把所有表格行的拼接统一走 esc()。"""
import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "scripts/v8_1/49_build_latex.py"
s = io.open(P, encoding="utf-8").read()
orig = s

# 所有 " & ".join(...) + r" \\"  统一替换为 row(...)
s = re.sub(r'" & "\.join\((r|rows\d+|c)?\)?[^)]*\)\s*\+ r" \\\\"',
           lambda m: m.group(0), s)   # no-op guard

# 手工列出全部构造点，逐个改成 row()
pairs = [
    ('t1 = "\\n".join(" & ".join(r) + r" \\\\" for r in rows1)',
     't1 = "\\n".join(row(r) for r in rows1)'),
    ('t2 = "\\n".join(" & ".join(r) + r" \\\\" for r in rows2)',
     't2 = "\\n".join(row(r) for r in rows2)'),
    ('t3 = "\\n".join(" & ".join(r) + r" \\\\" for r in rows3)',
     't3 = "\\n".join(row(r) for r in rows3)'),
    ('tS1 = "\\n".join(" & ".join(r) + r" \\\\" for r in rowsS1)',
     'tS1 = "\\n".join(row(r) for r in rowsS1)'),
    ('tS4 = "\\n".join(" & ".join(r) + r" \\\\" for r in rowsS4)',
     'tS4 = "\\n".join(row(r) for r in rowsS4)'),
    ('tS5 = "\\n".join(" & ".join(r) + r" \\\\" for r in rowsS5)',
     'tS5 = "\\n".join(row(r) for r in rowsS5)'),
    ('tS6 = "\\n".join(" & ".join(r) + r" \\\\" for r in rowsS6)',
     'tS6 = "\\n".join(row(r) for r in rowsS6)'),
    ('tS3 = "\\n".join(" & ".join(r) + r" \\\\" for r in rowsS3)',
     'tS3 = "\\n".join(row(r) for r in rowsS3)'),
]
applied = 0
for a, b in pairs:
    if a in s:
        s = s.replace(a, b)
        applied += 1
    else:
        print("MISS:", a[:70])

# 取消各处手工的 .replace("_", r"\_")，交给 esc()
s = s.replace('.replace("_", r"\\_")', '')
s = s.replace('str(r.artifact_id).replace("_", r"\\_")', 'str(r.artifact_id)')
s = s.replace('str(r.method_family).replace("_", r"\\_")', 'str(r.method_family)')
s = s.replace('str(r.split_regime).replace("_", r"\\_")', 'str(r.split_regime)')
s = s.replace('str(r.metric).replace("_", r"\\_")', 'str(r.metric)')

io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print(f"routed {applied}/8 row builders through row(); syntax OK")
print("remaining raw joins:", s.count('" & ".join'))
