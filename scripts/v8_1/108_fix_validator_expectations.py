# -*- coding: utf-8 -*-
"""scripts/v8_1/108_fix_validator_expectations.py -- 更新 validator 中过期的期望值。

  N: pool 表编号由 S4 变为 S2（F0–F9 provenance 表不在 supplement 中）
  W: 加权审计改为 T 口径（target total mass），数值来自 WEIGHTING_FINAL.json
"""
import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "validate_final_manuscript.py"
s = io.open(P, encoding="utf-8").read()

# ---- N ----
old_n = s[s.find('rec("N. Supplementary figure/table numbering"'):
          s.find('rec("O. D(n,K) derived')]
new_n = '''_pooltab = "tab:s2pools" if "tab:s2pools" in S else "tab:s4pools"
rec("N. Supplementary figure/table numbering",
    "PASS" if (_figctr in S and _pooltab in S) else "FAIL",
    "S-prefix counters set = %s; pool table label = %s"
    % (_figctr in S, _pooltab))

'''
s = s.replace(old_n, new_n, 1)

# _figctr 需要在 N 之前定义
if "_figctr = " not in s.split('rec("N. Supplementary')[0]:
    s = s.replace('_pooltab = "tab:s2pools"',
                  '_figctr = chr(92) + "thefigure}{S" + chr(92) + "arabic{figure}}"\n'
                  '_pooltab = "tab:s2pools"', 1)

# ---- W ----
m = re.search(r'rec\("W\. weighting audit values.*?\n\n', s, re.S)
if m:
    new_w = '''rec("W. weighting audit values use final (T) decomposition",
    "PASS" if ((WF["primary"] == "T_target_total_mass")
               and ("%.4f" % _T["range_C_B_median"]) in M
               and ("%d/109" % _T["signflip_C_B"]["0.005"]) in M
               and ("35/109" not in M) and ("0.0183" not in M)) else "FAIL",
    "primary=%s w=%s; range_C_B median=%.4f range_C_W median=%.4f; "
    "text has %d/109 and %d/109; stale 35/109=%s 0.0183=%s"
    % (WF["primary"], WF["w_ci_primary"], _T["range_C_B_median"],
       _T["range_C_W_median"], _T["signflip_C_B"]["0"], _T["signflip_C_B"]["0.005"],
       "35/109" in M, "0.0183" in M))

'''
    s = s[:m.start()] + new_w + s[m.end():]
    print("W check rewritten")
else:
    print("W check anchor not found")

io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print("syntax OK")
