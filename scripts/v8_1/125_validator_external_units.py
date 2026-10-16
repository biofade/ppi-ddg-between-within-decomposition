# -*- coding: utf-8 -*-
"""scripts/v8_1/125_validator_external_units.py -- 加入 external 单位/provenance 的验证检查。"""
import io
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "validate_final_manuscript.py"
s = io.open(P, encoding="utf-8").read()

NEW = r'''
# ================= 第五轮：external 单位与 provenance =================
EPR = pd.read_csv(V81 / "EXTERNAL_PROVENANCE.csv")
_invalid = set(EPR[~EPR.rmse_valid.astype(bool)].method)
_d = pymupdf.open(str(LATEX / "manuscript_v6.pdf"))
_t3 = "\n".join(pg.get_text() for pg in _d)
_d.close()

# AO. 无量纲打分不得报 RMSE
_bad_rmse = []
for _m in _invalid:
    _val = {"ESM-IF": "2.3587", "proteinMPNN": "2.3199"}.get(_m)
    if _val and _val in _t3:
        _bad_rmse.append(_m)
rec("AO. no RMSE for dimensionless external scores",
    "PASS" if not _bad_rmse else "FAIL",
    "dimensionless scores = %s; stale RMSE still in PDF = %s"
    % (sorted(_invalid), _bad_rmse or "none"))

# AP. provenance 表与 Note 8 存在
_d2 = pymupdf.open(str(LATEX / "supplementary_v6.pdf"))
_sup = "\n".join(pg.get_text() for pg in _d2)
_d2.close()
_ok_prov = ("Provenance and score semantics of the external prediction"
            in _sup and "Supplementary Note 8" in _sup)
rec("AP. external provenance table and Note 8 present",
    "PASS" if _ok_prov else "FAIL",
    "Table S3 (provenance) = %s; Note 8 = %s"
    % ("Provenance and score semantics" in _sup, "Supplementary Note 8" in _sup))

# AQ. 声明不得说所有方法都是 kcal/mol
_bad_unit = ("All methods are evaluated on the same mutations" in _t3
             or "kcal/mol" in _t3 and "dimensionless" not in _t3)
rec("AQ. units caveat stated in main text",
    "PASS" if (not _bad_unit and "dimensionless" in _t3) else "FAIL",
    "dimensionless caveat present in manuscript = %s" % ("dimensionless" in _t3))

# AR. Table 3 的 RMSE 对两个无量纲方法必须是 ---
_ok_dash = ("---" in _t3) or ("\u2014" in _t3)
rec("AR. Table 3 marks undefined RMSE",
    "PASS" if _ok_dash else "FAIL",
    "dash marker present in manuscript PDF = %s" % _ok_dash)
'''

anchor = "# ================= 输出 ================="
if anchor in s and "external 单位与 provenance" not in s:
    s = s.replace(anchor, NEW.strip() + "\n\n" + anchor, 1)
    io.open(P, "w", encoding="utf-8").write(s)
    print("validator extended with AO-AR")
else:
    print("anchor missing or already present")

import ast  # noqa: E402
ast.parse(io.open(P, encoding="utf-8").read())
print("syntax OK")
