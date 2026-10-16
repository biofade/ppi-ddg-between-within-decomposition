# -*- coding: utf-8 -*-
"""scripts/v8_1/133_validator_orientation.py -- 加入方向对齐相关的验证检查。"""
import io
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "validate_final_manuscript.py"
s = io.open(P, encoding="utf-8").read()

NEW = r'''
# ================= 第七轮：方向对齐 =================
EO = load_json("EXTERNAL_ORIENTATION.json")
_RA, _RB = EO["raw"], EO["oriented"]

_dm2 = pymupdf.open(str(LATEX / "manuscript_v6.pdf"))
_tm2 = "\n".join(pg.get_text() for pg in _dm2)
_dm2.close()
_ds2 = pymupdf.open(str(LATEX / "supplementary_v6.pdf"))
_ts2 = "\n".join(pg.get_text() for pg in _ds2)
_ds2.close()

# AY. ESM-IF 必须已方向对齐（Table 3 显示 +0.2250 / -0.0213）
_ok_orient = ("+0.2250" in _tm2 and "-0.0213" in _tm2)
rec("AY. ESM-IF sign-oriented in main table",
    "PASS" if _ok_orient else "FAIL",
    "oriented ESM-IF values (+0.2250, -0.0213) in manuscript PDF = %s" % _ok_orient)

# AZ. rank agreement 必须报 oriented 值，且同时说明 raw 值
_rho_ok = (f"{_RB['rho']:.2f}" in _tm2 or "0.30" in _tm2)
_p_ok = f"{_RB['perm_p']:.3f}"[:4] in _tm2      # 0.683
_raw_ok = ("+0.90" in _tm2 or "0.90" in _tm2)
rec("AZ. rank agreement reported with orientation stated",
    "PASS" if (_rho_ok and _p_ok and _raw_ok) else "FAIL",
    "oriented rho=%.2f present=%s; p=%.3f present=%s; raw alternative stated=%s"
    % (_RB["rho"], _rho_ok, _RB["perm_p"], _p_ok, _raw_ok))

# BA. 方向对齐不得使用 label / 校准
_ok_nocal = ("no evaluation labels" in _tm2 or "uses no evaluation labels" in _tm2)
rec("BA. orientation alignment is label-free",
    "PASS" if _ok_nocal else "FAIL",
    "'uses no evaluation labels' stated in Methods = %s" % _ok_nocal)

# BB. 补充材料所有 ESM-IF / ProteinMPNN 的 RMSE 必须为破折号
_bad_rmse = [v for v in ("2.3587", "2.3199", "2.3496", "2.3111", "2.1620", "2.1265",
                         "2.3623") if v in _ts2]
rec("BB. no RMSE for dimensionless scores anywhere in supplement",
    "PASS" if not _bad_rmse else "FAIL",
    "stale RMSE values found in supplementary PDF = %s" % (_bad_rmse or "none"))

# BC. Table S4 必须含方向列
_ok_dir = ("orientation" in _ts2.lower() and ("mu" in _ts2.lower() or "μ" in _ts2))
rec("BC. provenance table records raw orientation and multiplier",
    "PASS" if _ok_dir else "FAIL",
    "orientation column present = %s" % _ok_dir)

# BD. "same targets" 措辞已改
rec("BD. no incorrect 'same targets' claim",
    "PASS" if "evaluated against the same targets" not in _ts2 else "FAIL",
    "stale phrase present = %s" % ("evaluated against the same targets" in _ts2))
'''

anchor = "# ================= 输出 ================="
if anchor in s and "第七轮：方向对齐" not in s:
    s = s.replace(anchor, NEW.strip() + "\n\n" + anchor, 1)
    io.open(P, "w", encoding="utf-8").write(s)
    print("validator extended with AY-BD")
else:
    print("anchor missing or already present")

import ast  # noqa: E402
ast.parse(io.open(P, encoding="utf-8").read())
print("syntax OK")
