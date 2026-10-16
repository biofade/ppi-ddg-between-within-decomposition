# -*- coding: utf-8 -*-
"""scripts/v8_1/134_fix_checkF.py -- 更新 check F 以适配方向对齐后的 Table 3。"""
import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "validate_final_manuscript.py"
s = io.open(P, encoding="utf-8").read()

m = re.search(r"# =+ F\. External panel consistency =+.*?(?=# =+ F2\.)", s, re.S)
if not m:
    print("F block not found"); sys.exit(1)

NEW = r'''# ================= F. External panel consistency =================
ext = AC["external"]
prof = {p["model"]: p for p in ext["profiles"]}
EO0 = load_json("EXTERNAL_ORIENTATION.json")
MULT0 = EO0["multipliers"]
okF, detF = True, []
# 允许行尾的 $^{\dagger}$ 标记
tbl3 = re.findall(r"^(Flex-ddG|DDGPred|FoldX \(ddG\\_FoldX\)|ESM-IF|ProteinMPNN)"
                  r"(?:\$\^\{\\dagger\}\$)?\s*&\s*([\d,]+)\s*&\s*(\d+)\s*&\s*"
                  r"([+\-]?[\d\.]+)\s*&\s*([+\-]?[\d\.]+)", M, re.M)
names = {"Flex-ddG": "Flex-ddG", "DDGPred": "DDGPred",
         "FoldX (ddG\\_FoldX)": "ddG_FoldX", "ESM-IF": "ESM-IF",
         "ProteinMPNN": "proteinMPNN"}
for nm, rows, tg, pe, ms in tbl3:
    k = names[nm]
    mu = MULT0.get(k, 1)
    if int(rows.replace(",", "")) != int(prof[k]["n_rows"]):
        okF = False; detF.append(f"{k} rows {rows} != {prof[k]['n_rows']}")
    if int(tg) != int(prof[k]["n_targets_ge5"]):
        okF = False; detF.append(f"{k} targets {tg} != {prof[k]['n_targets_ge5']}")
    # 表中数值应为方向对齐后的值
    if abs(float(pe) - mu * float(prof[k]["overall_pearson"])) > 5e-4:
        okF = False; detF.append(f"{k} pearson {pe} != {mu}*{prof[k]['overall_pearson']:.4f}")
    if abs(float(ms) - mu * float(prof[k]["macro_spearman"])) > 5e-4:
        okF = False; detF.append(f"{k} macro {ms} != {mu}*{prof[k]['macro_spearman']:.4f}")
rec("F. Table 3 external values (sign-oriented)",
    "PASS" if okF and len(tbl3) == 5 else "FAIL",
    f"parsed {len(tbl3)} method rows; all match mu_c x raw from external_clean_profiles.csv"
    if okF else "; ".join(detF))

'''
s = s[:m.start()] + NEW + s[m.end():]
io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print("check F updated; syntax OK")
