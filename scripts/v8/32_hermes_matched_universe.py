# -*- coding: utf-8 -*-
"""scripts/v8/32_hermes_matched_universe.py -- HERMES × V8 canonical 的**精确行级匹配**（§21/§50）。

HERMES 官方预测列：
    PDB_filename | mutant_chain | mutant | Mutation(s)_cleaned | #Pdb | log_proba_mt__minus__log_proba_wt

**关键方向问题**：HERMES 输出的是 **log odds（突变体 − 野生型）**，不是 ΔΔG。
  - `log_proba_mt − log_proba_wt` **越大 = 越有利**
  - V8 约定 `ΔΔG > 0 = 结合变弱` ⇒ **越大 = 越不利**
  ⇒ 两者**方向相反**；构造 `hermes_ddg_proxy = −log_proba_mt__minus__log_proba_wt`
    以对齐 V8 的符号约定。**此符号约定必须在报告中显式声明**（§19）。

匹配键：canonical 的 `(target_id_canonical, mutation)` 对 HERMES 的
`(PDB_filename 前4位, Mutation(s)_cleaned)`。

输出：`external/hermes_matched.csv`（matched universe 的行级表）
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
DL = Path("F:/dsh/downloads/hermes")
OUT = R / "experiments" / "v8" / "external"
OUT.mkdir(parents=True, exist_ok=True)

lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)
print(f"canonical: {len(lab)} 行 / {lab.target_id_canonical.nunique()} targets")
print(f"canonical 全部列: {list(lab.columns)}")

# ---- canonical 侧：找突变列 ----
cand_mut = [c for c in lab.columns
            if any(k in c.lower() for k in ("mut", "variant", "aa_change",
                                            "mutation"))]
print(f"\n候选突变列: {cand_mut}")
for c in cand_mut:
    print(f"  {c!r}: 样例 {lab[c].astype(str).head(4).tolist()}")

mcol = "mutation" if "mutation" in lab.columns else (cand_mut[0] if cand_mut else None)
print(f"\n使用突变列 = {mcol!r}")

def norm_mut(s: str) -> str:
    """归一化突变串：去空白/逗号排序/大写。"""
    s = re.sub(r"\s+", "", str(s)).upper()
    if not s or s == "NAN":
        return ""
    parts = [p for p in re.split(r"[,;+]", s) if p]
    return ",".join(sorted(parts))

# ---- HERMES 侧 ----
H = {}
for tag in ("hermes_py000_cdna117k", "hermes_py050_cdna117k"):
    p = DL / f"{tag}.csv"
    if p.exists():
        H[tag] = pd.read_csv(p, low_memory=False)
print(f"\nHERMES 文件: {list(H.keys())}")
if not H:
    sys.exit("无 HERMES 文件")

h0 = H[list(H.keys())[0]]
LP = "log_proba_mt__minus__log_proba_wt"
h0 = h0.copy()
h0["pdb4"] = h0["PDB_filename"].astype(str).str.upper().str.strip().str[:4]
h0["mut_n"] = h0["Mutation(s)_cleaned"].map(norm_mut)
h0["key"] = h0["pdb4"] + "|" + h0["mut_n"]
# 方向对齐：HERMES log-odds 越大越有利；V8 ddG 越大越不利 ⇒ 取负
h0["hermes_proxy"] = -pd.to_numeric(h0[LP], errors="coerce")

print(f"\nHERMES 行 = {len(h0)}, 唯一 key = {h0.key.nunique()}")
dup = h0.key.duplicated().sum()
print(f"重复 key = {dup}  ⇒ 需按 key 聚合（取均值）")

# ---- canonical 侧 key ----
c = lab.copy()
c["pdb4"] = c.target_id_canonical.astype(str).str.upper().str.strip().str[:4]
c["mut_n"] = c[mcol].map(norm_mut)
c["key"] = c["pdb4"] + "|" + c["mut_n"]
print(f"canonical 唯一 key = {c.key.nunique()} / {len(c)}")

# ---- 合并 ----
hagg = h0.groupby("key", as_index=False).agg(
    hermes_proxy=("hermes_proxy", "mean"), n_hermes=("hermes_proxy", "size"))
M = c.merge(hagg, on="key", how="left")
hit = M.hermes_proxy.notna()
print("\n" + "=" * 70)
print("匹配结果（精确 (PDB4, 归一化突变串)）")
print("=" * 70)
print(f"  canonical 行        = {len(M)}")
print(f"  命中 HERMES 的行    = {int(hit.sum())}  ({hit.mean():.4f})")
print(f"  未命中              = {int((~hit).sum())}")

# 只用单点突变的子集（HERMES 侧重单点）
c["n_mut"] = c[mcol].astype(str).map(
    lambda s: len([x for x in re.split(r"[,;+]", s) if x.strip()]))
single = c.n_mut == 1
Ms = M[single]
print(f"\n  仅单点突变: canonical = {int(single.sum())}, "
      f"命中 = {int(Ms.hermes_proxy.notna().sum())} "
      f"({Ms.hermes_proxy.notna().mean():.4f})")

if hit.sum() > 0:
    y = pd.to_numeric(M.loc[hit, "ddG"] if "ddG" in M.columns else None,
                      errors="coerce")
    hp = M.loc[hit, "hermes_proxy"]
    print(f"\n  匹配子集: HERMES proxy n={hp.notna().sum()} "
          f"mean={hp.mean():.4f} sd={hp.std():.4f}")
    # 找 canonical 的真值列
    ycol = next((x for x in M.columns
                 if x.lower() in ("ddg", "ddg_exp", "y", "label")), None)
    print(f"  canonical 真值列 = {ycol!r}")
    if ycol:
        yy = pd.to_numeric(M.loc[hit, ycol], errors="coerce")
        m2 = yy.notna() & hp.notna()
        if m2.sum() > 10:
            r = np.corrcoef(yy[m2], hp[m2])[0, 1]
            print(f"  粗相关 r(canonical {ycol}, HERMES proxy) = {r:+.4f} "
                  f"(n={int(m2.sum())})")
            print("  ⇒ 若为**负**，说明方向对齐错误；若为正且量级合理，方向正确。")

M.to_csv(OUT / "hermes_matched.csv", index=False, encoding="utf-8")
meta = {
    "hermes_files": list(H.keys()),
    "hermes_raw_column": LP,
    "sign_convention_applied": ("hermes_proxy = -(log_proba_mt - log_proba_wt)"
                                " 以对齐 V8 的 ΔΔG>0=结合变弱"),
    "canonical_rows": int(len(M)),
    "matched_rows": int(hit.sum()),
    "match_rate": float(hit.mean()),
    "matched_single_point": int(Ms.hermes_proxy.notna().sum()),
    "evidence_tier": ("TIER 2 (secondary) — 依 §50，HERMES 在 SKEMPI v2.0 上微调，"
                      "不能自动进入 Tier 1；须先审计其训练家族与 "
                      "OUTER_FOLDS_V7_2 的重叠"),
}
(OUT / "hermes_matched_meta.json").write_text(
    json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"\n写出 hermes_matched.csv / hermes_matched_meta.json")
