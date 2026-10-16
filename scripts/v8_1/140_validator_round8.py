# -*- coding: utf-8 -*-
"""scripts/v8_1/140_validator_round8.py -- 第八轮检查 + 修正旧检查的过期期望。"""
import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "validate_final_manuscript.py"
s = io.open(P, encoding="utf-8").read()

NEW = r'''
# ================= 第八轮：方向选择的诚实性 =================
_d8 = pymupdf.open(str(LATEX / "manuscript_v6.pdf"))
_t8 = "\n".join(pg.get_text() for pg in _d8)
_d8.close()
_s8 = pymupdf.open(str(LATEX / "supplementary_v6.pdf"))
_u8 = "\n".join(pg.get_text() for pg in _s8)
_s8.close()

# BE. 不得声称 orientation 是 label-free
_bad_labelfree = ("it uses no evaluation labels" in _t8
                  or "uses no evaluation labels" in _t8)
_ok_labelaware = ("selected using the evaluation labels" in _t8
                  or "label-informed" in _t8)
rec("BE. orientation is declared label-informed, not label-free",
    "PASS" if (not _bad_labelfree and _ok_labelaware) else "FAIL",
    "stale 'uses no evaluation labels' = %s; label-informed declaration = %s"
    % (_bad_labelfree, _ok_labelaware))

# BF. raw 为 primary，oriented 为 sensitivity，两者都报
_ok_dual = ("raw" in _t8.lower() and "sign-orient" in _t8.lower()
            and "0.90" in _t8 and "0.30" in _t8)
rec("BF. both raw and oriented results reported",
    "PASS" if _ok_dual else "FAIL",
    "raw rho present = %s; oriented rho present = %s"
    % ("0.90" in _t8, "0.30" in _t8))

# BG. 不得再出现 "no sign flip"
rec("BG. no contradictory 'no sign flip' claim",
    "PASS" if "no sign flip" not in _t8 and "no sign flip" not in _u8 else "FAIL",
    "main=%s supplement=%s"
    % ("no sign flip" in _t8, "no sign flip" in _u8))

# BH. Table S4 引用必须为 S3（且 Supplement 无 Table S4）
rec("BH. supplementary table reference is S3 not S4",
    "PASS" if ("Table S3" in _u8 and "Table S4" not in _u8
               and "Table~S4" not in M) else "FAIL",
    "Table S3 present = %s; 'Table S4' in supplement = %s"
    % ("Table S3" in _u8, "Table S4" in _u8))

# BI. Figure 6 Panel B 标题必须是 sign-oriented
rec("BI. Figure 6 Panel B title says sign-oriented",
    "PASS" if ("sign-oriented values printed" in _t8
               and "raw values printed" not in _t8) else "FAIL",
    "sign-oriented title = %s; stale raw title = %s"
    % ("sign-oriented values printed" in _t8, "raw values printed" in _t8))

# BJ. common-intersection 能量型方法的 RMSE 必须保留
_ok_rmse = all(v in _u8 for v in ("2.1620", "2.1265", "2.3623"))
rec("BJ. energy-type RMSE retained in common-intersection table",
    "PASS" if _ok_rmse else "FAIL",
    "Flex/FoldX/DDGPred RMSE present = %s" % _ok_rmse)

# BK. orientation multiplier 不得再用 mu_c
_bad_mu = ("\\mu_c\\, s^{\\mathrm{raw}}" in S or "$\\mu_c$ & RMSE" in S)
rec("BK. orientation multiplier uses o_m, not mu_c",
    "PASS" if (not _bad_mu and "o_m" in S and "o_m" in M) else "FAIL",
    "o_m present in both = %s; mu_c reused for orientation = %s"
    % ("o_m" in S and "o_m" in M, _bad_mu))

# BL. Note 7 必须报 oriented 排名
_ok_note7 = ("exploratory sign-oriented" in _u8
             and "no significant agreement" in _u8)
rec("BL. Note 7 updated to oriented analysis",
    "PASS" if _ok_note7 else "FAIL",
    "oriented wording in Note 7 = %s" % _ok_note7)
'''

anchor = "# ================= 输出 ================="
if anchor in s and "第八轮：方向选择的诚实性" not in s:
    s = s.replace(anchor, NEW.strip() + "\n\n" + anchor, 1)
    io.open(P, "w", encoding="utf-8").write(s)
    print("validator extended with BE-BL")
else:
    print("anchor missing or already present")

# ---- 修正过期检查：AZ 现在应接受 raw+oriented 并存 ----
s = io.open(P, encoding="utf-8").read()
old_az = s[s.find('rec("AZ. rank agreement reported with orientation stated"'):
           s.find('rec("BA. orientation alignment is label-free"')]
new_az = '''rec("AZ. rank agreement reported with orientation stated",
    "PASS" if (("0.90" in _tm2) and ("0.30" in _tm2) and ("0.683" in _tm2)) else "FAIL",
    "raw rho 0.90 present=%s; oriented rho 0.30 present=%s; p=0.683 present=%s"
    % ("0.90" in _tm2, "0.30" in _tm2, "0.683" in _tm2))

'''
s = s.replace(old_az, new_az, 1)
io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print("AZ expectation updated; syntax OK")
