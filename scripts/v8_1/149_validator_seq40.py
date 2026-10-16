# -*- coding: utf-8 -*-
"""scripts/v8_1/149_validator_seq40.py -- 加入 seq40 键级修正的防回归检查。"""
import io
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "validate_final_manuscript.py"
s = io.open(P, encoding="utf-8").read()

NEW = r'''
# ================= 第九轮：seq40 键级修正防回归 =================
# BM. seq40 不得再与 complex 相同（键级泄漏修复后它们是不同划分）
_h4 = FR["H4_controlled_split"]
_c, _s40 = _h4["complex"], _h4["seq40"]
_same = all(abs(_c[k] - _s40[k]) < 1e-12
            for k in ("overall_pearson", "C_B", "C_W", "rmse", "macro_spearman"))
rec("BM. seq40 is no longer identical to complex",
    "PASS" if not _same else "FAIL",
    "complex r=%.4f vs seq40 r=%.4f (identical=%s)"
    % (_c["overall_pearson"], _s40["overall_pearson"], _same))

# BN. 正文/图不得再声称两档 identical
_dm9 = pymupdf.open(str(LATEX / "manuscript_v6.pdf"))
_tm9 = "\n".join(pg.get_text() for pg in _dm9)
_dm9.close()
_bad9 = [x for x in ("identical assignments", "identical partition",
                     "Complex / seq40") if x in _tm9]
rec("BN. no 'identical' claim for seq40 in final PDF",
    "PASS" if not _bad9 else "FAIL",
    "stale wording found = %s" % (_bad9 or "none"))

# BO. 键级泄漏核验：四个 level 均无跨 fold 的 obs_key
import json as _j9
_kx = _j9.load(open(str(V81 / "KEY_LEAKAGE_EXACT.json"), encoding="utf-8"))
_leak = {k: v["keys_crossing_folds"] for k, v in _kx["detail"].items()
         if v["keys_crossing_folds"] > 0}
rec("BO. key-level leakage eliminated across all levels",
    "PASS" if not _leak else "FAIL",
    "levels with keys crossing folds = %s" % (_leak or "none"))
'''

anchor = "# ================= 输出 ================="
if anchor in s and "seq40 键级修正防回归" not in s:
    s = s.replace(anchor, NEW.strip() + "\n\n" + anchor, 1)
    io.open(P, "w", encoding="utf-8").write(s)
    print("validator extended with BM-BO")
else:
    print("anchor missing or already present")

import ast  # noqa: E402
ast.parse(io.open(P, encoding="utf-8").read())
print("syntax OK")
