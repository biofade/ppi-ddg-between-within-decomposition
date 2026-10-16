# -*- coding: utf-8 -*-
"""scripts/v8/40_external_heldout_cath.py -- 第二个外部宇宙：CATH-ddG 的 **held-out CATH 超家族测试集**。

## 为什么这是 V8 最贴近的**外部**设置

CATH-ddG 的 `cath_fold` 列给出官方划分：
    train = 5893 行,  val(=held-out CATH 测试) = **813 行**
论文描述该测试集为「54 个 CATH 超家族 / 53 个 PPI / 813 个突变」，
序列同一性 < 30%，且**不含 easy mutation**（最大 TM-score < 0.6）。

⇒ 这是**独立于 V8 的 family-OOD 构建方式**（CATH 超家族 vs 本项目手工家族）
的严格外推测试集，且**自带 6 个基线方法的逐突变预测**（零推理）。

**局限**：本表**不含 CATH-ddG 自身**的逐突变预测（未发布），
故只能评价其**基线**——但这正好是 V8 要比较的「常规方法」。

## §50 污染审计

基线列（Flex-ddG / DDGPred / FoldX / ESM-IF / proteinMPNN）在其本表语境下
是被**复算**的结果；它们的原始训练集与本测试集的重叠**未核实**。
故本宇宙记为 **TIER 2**。
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(R / "scripts" / "v8"))
import v8metrics as M  # noqa: E402

OUT = R / "experiments" / "v8"
EXT = OUT / "external"
SRC = Path("F:/dsh/downloads/cathddg/skempi_v2.csv")
MIN_N = 5
MODELS = ["Flex-ddG", "DDGPred", "ddG_FoldX", "ESM-IF", "proteinMPNN"]

E = pd.read_csv(SRC, low_memory=False)
print(f"源 = {E.shape}")
print(f"cath_fold 取值: {E.cath_fold.value_counts().to_dict()}")

T = E[E.cath_fold.astype(str).str.strip().str.lower() == "val"].copy()
print(f"\nheld-out CATH 测试集 = {len(T)} 行")
print(f"  PPI (complex_PPI) 数 = {T.complex_PPI.nunique()}")
print(f"  PDB 数 = {T.complex.nunique()}")
print(f"  num_muts 分布 = {T.num_muts.value_counts().to_dict()}")

# 逐模型数值化
for c in MODELS + ["ddG"]:
    T[c] = pd.to_numeric(T[c], errors="coerce")
run = ~T.Label.astype(str).str.contains("train", case=False, na=False)
print(f"\n  非 train 行 = {int(run.sum())} / {len(T)}")
print(f"  各模型可用数：")
for c in MODELS:
    print(f"    {c:<12} {int(T[c].notna().sum()):>4} "
          f"(非 train 中 {int(T.loc[run, c].notna().sum()):>4})")

# ---------- 三轴（target = complex_PPI，即 PPI 单位）----------
print("\n" + "=" * 78)
print("held-out CATH 三轴读数（target 单位 = complex_PPI）")
print("=" * 78)
TGT = T.complex_PPI.astype(str).to_numpy()
rows = []
print(f"{'model':<13}{'n':>6}{'targets':>9}{'tgt>=5':>8}{'Overall r':>11}"
      f"{'RMSE':>9}{'Macro rho':>11}{'C_W':>9}{'Eff@10':>9}")
for c in MODELS:
    m = T[c].notna().to_numpy() & np.isfinite(T.ddG.to_numpy(float))
    if m.sum() < 100:
        print(f"{c:<13} 可用行不足（{int(m.sum())}）")
        continue
    U = M.prep_universe(T.ddG.to_numpy(float)[m], T[c].to_numpy(float)[m],
                        TGT[m])
    r = M.compute_all(U)
    r["model"] = c
    rows.append(r)
    print(f"{c:<13}{r['n_rows']:>6}{r['n_targets']:>9}{r['n_targets_ge5']:>8}"
          f"{r['overall_pearson']:>+11.4f}{r['rmse']:>9.4f}"
          f"{r['macro_spearman']:>+11.4f}{r['C_W']:>+9.4f}"
          f"{r['efficiency@10']:>+9.4f}")

D = pd.DataFrame(rows)
if len(D):
    D.to_csv(EXT / "external_heldout_cath.csv", index=False, encoding="utf-8")
    sel = D.dropna(subset=["overall_pearson", "macro_spearman", "efficiency@10"])
    print(f"\n  可评价外部模型 n = {len(sel)}")
    if len(sel) >= 5:
        r1 = spearmanr(sel.overall_pearson, sel.macro_spearman).correlation
        r2 = spearmanr(sel.overall_pearson, sel["efficiency@10"]).correlation
        r3 = spearmanr(sel.macro_spearman, sel["efficiency@10"]).correlation
        print(f"  rho(Overall r, Macro rho) = {r1:+.4f}")
        print(f"  rho(Overall r, Eff@10)    = {r2:+.4f}")
        print(f"  rho(Macro rho, Eff@10)    = {r3:+.4f}")
        print(f"\n  ⚠ n={len(sel)} 的秩相关无判别力（n=5 时 95% CI 约 [−0.4, +0.99]）")
        print("     ⇒ 仍然 NOT_TESTABLE")
    else:
        print(f"  ⚠ 仅 {len(sel)} 个 ⇒ 秩相关不可计算")

meta = {"universe": "U_EXT_HELDOUT_CATH",
        "source": "https://github.com/ak422/CATH-ddG data/SKEMPI2/skempi_v2.csv",
        "split_column": "cath_fold == 'val'",
        "n_rows": int(len(T)), "n_ppi": int(T.complex_PPI.nunique()),
        "models": MODELS, "n_evaluable": int(len(D)),
        "evidence_tier": "TIER 2（基线原始训练集与本测试集的重叠未核实，§50）",
        "not_published": "CATH-ddG 自身的逐突变预测未发布（补充材料仅 per-protein 表）",
        "verdict": ("EXTERNAL_MODEL_POOL_TOO_SMALL: n=%d < 8 ⇒ NOT_TESTABLE"
                    % len(D))}
(EXT / "external_heldout_cath_meta.json").write_text(
    json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"\n写出 external_heldout_cath.csv / _meta.json")
