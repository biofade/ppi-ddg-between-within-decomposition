# -*- coding: utf-8 -*-
"""scripts/v8_1/43_apply_external_fix.py -- 精确替换 external-method 两处旧版残留。"""
import io
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "scripts/v8_1/31_manuscript_v6.py"
s = io.open(P, encoding="utf-8").read()
orig = s

OLD1 = 'f"predictions: Flex ddG[6], FoldX[7], ProteinMPNN[8], ESM-IF[9] and HERMES[4]. "'
NEW1 = 'f"predictions: Flex ddG[6], FoldX[7], ProteinMPNN[8], ESM-IF[9] and DDGPred[4]. "'

OLD2 = ('         f"underpowered, so we draw no directional conclusion from it. We also note that one "\n'
        '         f"external quantity is a log-odds difference, a dimensionless mutation-effect score "\n'
        '         f"rather than a calibrated kcal/mol ΔΔG[4]; we therefore report neither RMSE nor "\n'
        '         f"MAE for that model.")')
NEW2 = ('         f"underpowered, so we draw no directional conclusion from it.")')

n1 = s.count(OLD1)
n2 = s.count(OLD2)
s = s.replace(OLD1, NEW1).replace(OLD2, NEW2)
io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print(f"replaced list x{n1}, log-odds sentence x{n2}")
print("changed:", s != orig, "| syntax OK")
print("HERMES left:", s.count("HERMES"), "| log-odds left:", s.count("log-odds"))
