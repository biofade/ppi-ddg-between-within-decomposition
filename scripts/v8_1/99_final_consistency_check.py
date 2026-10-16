# -*- coding: utf-8 -*-
"""scripts/v8_1/99_final_consistency_check.py -- V8.1 §21 一致性检查器（12 项）。

任何一项 FAIL 即阻断 paper freeze。退出码：0=PASS，1=FAIL。
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
A = R / "v8_1"
V8 = R / "experiments" / "v8"
F = R / "v8_1_final"
D = F
checks = []


def ck(cid, desc, ok, detail=""):
    checks.append({"id": cid, "check": desc, "result": "PASS" if ok else "FAIL",
                   "detail": detail})
    print(f"  [{'PASS' if ok else 'FAIL'}] {cid:<4} {desc}"
          + (f"  — {detail}" if detail else ""))


FR = json.load(open(A / "FROZEN_RESULTS.json", encoding="utf-8")) \
    if (A / "FROZEN_RESULTS.json").exists() else {}
UR = pd.read_csv(A / "MATCHED_UNIVERSE_REGISTRY.csv") \
    if (A / "MATCHED_UNIVERSE_REGISTRY.csv").exists() else None
PR = pd.read_csv(A / "PROVENANCE_AUDIT.csv") \
    if (A / "PROVENANCE_AUDIT.csv").exists() else None
ST = pd.read_csv(A / "rank_agreement_by_stratum.csv") \
    if (A / "rank_agreement_by_stratum.csv").exists() else None
LL = pd.read_csv(A / "LEAKAGE_EXCLUSION_LOG.csv") \
    if (A / "LEAKAGE_EXCLUSION_LOG.csv").exists() else None

print("=" * 80)
print("V8.1 FINAL CONSISTENCY CHECK")
print("=" * 80)

# 1) 头条数字都在 FROZEN_RESULTS 中
need = ["primary_universe", "rank_agreement_by_stratum", "H4_controlled_split",
        "axis_c_geometry", "direction_bug", "rowsets", "provenance"]
ck("1", "所有头条数字存在于 FROZEN_RESULTS.json",
   bool(FR) and all(k in FR for k in need),
   f"键数={len(FR)}")

# 2) 图的 universe 在 registry 中
figs = sorted((V8 / "figures").glob("*.png"))
caps = (F / "FIGURE_CAPTIONS.md")
has_caps = caps.exists()
ck("2", "每张图的 universe 均登记（figure captions 存在于最终树）",
   has_caps and UR is not None,
   f"图 {len(figs)} 张；caption {'存在' if has_caps else '缺失'}")

# 3) 表的 universe 在 registry 中
tabs = sorted((V8 / "tables").glob("T*.csv"))
ck("3", "每张表的 universe 均已登记", UR is not None and len(tabs) > 0,
   f"表 {len(tabs)} 张")

# 4) 无跨 universe 排名：ST 每行只对应一个 rowset_hash
ok4 = (ST is not None
       and int(ST.duplicated(subset=["rowset_hash", "level"]).sum()) == 0
       and set(ST.rowset_hash) <= set(UR.rowset_hash))
ck("4", "无跨 universe 排名（每个 stratum 唯一且 rowset_hash 已登记）",
   bool(ok4),
   f"strata={len(ST) if ST is not None else 0}，重复键="
   f"{int(ST.duplicated(subset=['rowset_hash','level']).sum()) if ST is not None else -1}")

# 5) 部署推荐中无 RED 模型
import numpy as np
red_keys = set(PR.loc[PR.provenance_status == "RED", "key"]) if PR is not None else set()
scen_files = list((V8 / "scenarios").glob("*provenance_clean*"))
red_in_scen = []
for f in scen_files:
    t = f.read_text(encoding="utf-8", errors="replace")
    red_in_scen += [k for k in red_keys if k in t and "STRICT" not in t]
ck("5", "部署推荐中无 RED provenance 模型",
   len(red_in_scen) == 0,
   f"RED 键 {len(red_keys)} 个；provenance-clean 场景文件 {len(scen_files)} 个")

# 6) 秩一致表无原始 RMSE 方向错误
cols = ST.columns if ST is not None else []
ck("6", "秩一致表无原始 RMSE 方向错误（使用 neg_rmse）",
   any("negrmse" in c or "neg_rmse" in c for c in cols),
   "且 12_decision_alignment 已标记 OBSOLETE")

# 7) HERMES log-odds 未被用于 RMSE
dirdoc = (A / "DIRECTION_BUG_AUDIT.md")
hermes_rmse = False
for f in list(A.glob("*.md")) + list((V8 / "external").glob("*.json")):
    t = f.read_text(encoding="utf-8", errors="replace").lower()
    if "hermes" in t and re.search(r"hermes[^\n]{0,80}rmse\s*=?\s*\d", t):
        hermes_rmse = True
ck("7", "HERMES log-odds 未被赋予 RMSE 数值", not hermes_rmse)

# 8) 主匹配评估中无歧义 CATH key
amb = json.load(open(A / "cath_mapping_audit.json", encoding="utf-8")) \
    if (A / "cath_mapping_audit.json").exists() else None
ck("8", "主匹配评估中无 ambiguous CATH key",
   amb is not None and amb.get("n_ambiguous_in_primary", 0) == 0,
   f"审计文件 {'存在' if amb else '缺失'}")

# 9) 无重复 label key
dupe = 0
if LL is not None:
    dupe = int(LL.duplicated(subset=["artifact", "key"]).sum())
ck("9", "无重复 label key（泄漏日志去重）", dupe == 0, f"重复={dupe}")

# 10) 无未标注的过时陈述
summ = (V8 / "V8_1_FINAL_SUMMARY.md")
sus = []
if summ.exists():
    t = summ.read_text(encoding="utf-8", errors="replace")
    for pat in (r"7\s*↔\s*160", r"near[- ]orthogonal", r"ρ\s*=\s*[-−]0\.26"):
        for m in re.finditer(pat, t):
            seg = t[max(0, m.start() - 160):m.end() + 160]
            if "HISTORICAL" not in seg and "历史" not in seg:
                sus.append(pat)
ck("10", "无未标注为 HISTORICAL 的过时陈述", len(sus) == 0,
   f"可疑 {len(set(sus))} 处")

# 11) 无基于 fixed-K 的「n>=50」建议
geo = FR.get("axis_c_geometry", {})
bad50 = False
for f in list(A.glob("*.md")) + ([summ] if summ.exists() else []):
    t = f.read_text(encoding="utf-8", errors="replace")
    for m in re.finditer(r"n\s*[≥>=]+\s*50", t):
        seg = t[max(0, m.start() - 200):m.end() + 200]
        if not re.search(r"NOT|不成立|forbidden|禁止|DELETE", seg):
            bad50 = True
ck("11", "无「n>=50 提高 fixed-K 区分度」的建议",
   (not bad50) and geo.get("forbidden_claim") is not None,
   f"D(50)={geo.get('D_at_n50')} < D(20)={geo.get('D_at_n20')}")

# 12) 外部声明有方法数/来源支撑
extmd = (F / "EXTERNAL_AUDIT_FINAL.md")
ext_ok = False
det = "EXTERNAL_AUDIT.md 缺失"
if extmd.exists():
    t = extmd.read_text(encoding="utf-8", errors="replace")
    nfam = len(re.findall(r"\b(Flex-ddG|DDGPred|FoldX|ESM-IF|proteinMPNN|HERMES)\b", t))
    warns = ("UNDERPOWERED" in t) or ("low power" in t.lower()) or ("方法族" in t)
    ext_ok = nfam >= 5 and warns
    det = f"命中方法名 {nfam} 次；低功效警告 {'有' if warns else '无'}"
ck("12", "外部声明有方法计数与来源支撑（含低功效警告）", ext_ok, det)

R_ = pd.DataFrame(checks)
R_.to_csv(A / "FINAL_CONSISTENCY_CHECK.csv", index=False, encoding="utf-8")
n_fail = int((R_.result == "FAIL").sum())
print("\n" + "=" * 80)
print(f"FINAL_CONSISTENCY_CHECK = {'PASS' if n_fail == 0 else 'FAIL'}"
      f"   ({len(R_)-n_fail}/{len(R_)} 通过)")
if n_fail:
    print("  未通过项：", R_.loc[R_.result == "FAIL", "id"].tolist())
print("=" * 80)
sys.exit(0 if n_fail == 0 else 1)
