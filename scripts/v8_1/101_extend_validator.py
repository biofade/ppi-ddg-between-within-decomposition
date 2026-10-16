# -*- coding: utf-8 -*-
"""scripts/v8_1/101_extend_validator.py -- 把第二轮审计的检查项加入 validate_final_manuscript.py。"""
import io
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "validate_final_manuscript.py"
s = io.open(P, encoding="utf-8").read()

NEW = r'''
# ================= 第二轮审计新增检查 =================
W = load_json("WEIGHTING_CONSISTENT.json")
DK = load_json("DUPLICATE_KEY_SENSITIVITY.json")

# U. 加权分解恒等式（逐 scheme，109 配置）
res_w = W["identity_residual_by_scheme"]
worst = max(res_w.values())
rec("U. weighted identity r_w = C_B,w + C_W,w",
    "PASS" if worst <= 1e-12 and W["identity_holds_all"] else "FAIL",
    f"109 configs x 5 schemes; max|residual| = {worst:.3e}; "
    f"by scheme {res_w}")

# V. Note 5 不得声称 C_W 跨 scheme 不变
bad_inv = ("identical across all five schemes" in S
           or "identical across all five" in S)
rec("V. Note 5 self-consistency (C_W not scheme-invariant)",
    "PASS" if not bad_inv else "FAIL",
    f"forbidden 'identical across all five schemes' present = {bad_inv}; "
    f"C_W cross-scheme range median = {W['range_C_W_median_new']:.4f}")

# W. 加权审计数值与正文/Figure 3 一致
old_vals = ("35/109" in M and "0.0183" in M)
new_vals = ("44/109" in M and "32/109" in M)
rec("W. weighting audit values updated to consistent decomposition",
    "PASS" if new_vals and not old_vals else "FAIL",
    f"new (44/109, 32/109) present = {new_vals}; stale (35/109, 0.0183) present = {old_vals}; "
    f"range_C_B median new = {W['range_C_B_median_new']:.4f} "
    f"(old {W['range_C_B_median_old']:.4f})")

# X. D(n,K) 必须声明 rank-linear reference
rank_ok = ("rank-linear" in S and "u_{(j)} = n - j + 1" in S)
rec("X. D(n,K) states rank-linear utility reference",
    "PASS" if rank_ok else "FAIL",
    f"rank-linear reference declared = {rank_ok}; "
    f"removed universal-utility claim = "
    f"{'any monotone transform' not in S}")

# Y. 重复键敏感性
rec("Y. duplicate-key sensitivity analysis",
    "PASS" if ("Duplicate-key sensitivity" in S
               and DK["decomp_verified_max_abs_diff"] < 1e-9
               and DK["direction_preserved"]) else "FAIL",
    f"Note 6 present = {'Duplicate-key sensitivity' in S}; "
    f"decomposition verified to {DK['decomp_verified_max_abs_diff']:.1e}; "
    f"attribution {100*min(v['attribution'] for v in DK['by_mode'].values()):.1f}"
    f"-{100*max(v['attribution'] for v in DK['by_mode'].values()):.1f}% "
    f"(direction preserved = {DK['direction_preserved']})")

# Z. 补充材料 caption 不得重复前缀
import re as _re
bad_cap = _re.findall(r"\\caption\{(?:Figure S\d|Table S\d)\s*[:.]\s*(?:Figure|Table)",
                      S) + _re.findall(r"\\caption\{\\textbf\{(?:Figure|Table) S\d", S)
rec("Z. supplementary caption numbering",
    "PASS" if not bad_cap else "FAIL",
    f"doubled prefixes found = {len(bad_cap)}; "
    f"S-counter set = {'\\thefigure}{S\\arabic{figure}}' in S}; "
    f"Table S4 caption clean = "
    f"{'analysis.}' in S and 'analysis Table~S4' not in S}")

# AA. Figure 6 caption 不得暗示完全 matched
fig6_bad = ("on a matched,\nunambiguous mapping" in M
            or "on a matched, unambiguous mapping" in M)
rec("AA. Figure 6 caption wording",
    "PASS" if not fig6_bad and "uniquely mapped subset" in M else "FAIL",
    f"stale 'matched unambiguous mapping' wording present = {fig6_bad}")

# AB. Abstract 措辞
abs_ok = ("We formulate an exact between-/within-target covariance decomposition" in M
          and "summed\nabsolute component change" in M)
rec("AB. Abstract novelty and attribution wording",
    "PASS" if abs_ok else "FAIL",
    f"formulate-wording present = {'We formulate an exact' in M}; "
    f"'summed absolute component change' in Abstract/S = "
    f"{'summed' in M and 'absolute component change' in M}")
'''

anchor = "# ================= 输出 ================="
if anchor in s and "第二轮审计新增检查" not in s:
    s = s.replace(anchor, NEW.strip() + "\n\n" + anchor, 1)
    io.open(P, "w", encoding="utf-8").write(s)
    print("validator extended with checks U-AB")
else:
    print("anchor missing or already extended")

import ast  # noqa: E402
ast.parse(io.open(P, encoding="utf-8").read())
print("syntax OK")
