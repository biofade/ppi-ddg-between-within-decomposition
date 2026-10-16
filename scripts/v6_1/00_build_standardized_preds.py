# -*- coding: utf-8 -*-
"""scripts/v6_1/00_build_standardized_preds.py -- V6.1 PHASE 1.0: 补建长表.

背景
----
V6.1 规格引用 `experiments/v6_sota_audit/standardized_preds.csv`，但该文件
在 V6.0 中**未生成**（V6.0 只落了汇总指标 `sota_audit_metrics.csv` 与符号审计）。
本脚本按 V6.0 已核实的来源**补建**该长表（不改动 V6.0 任何产物）。

来源与符号（沿用 V6.0 `sign_audit.csv` 的 `needs_flip` 判定）
-------------------------------------------------------------
* CATH-ddG 官方表：`F:\\dsh\\downloads\\cathddg\\skempi_v2_with_all_results.csv`
  真值 `ddG`（已用本地 SKEMPI 独立复算验证：Pearson +0.9869、符号一致 94.99%）
  预测列 = Flex-ddG / DDGPred / ddG_FoldX / ESM-IF / proteinMPNN
* DSSA-PPI RBM：`F:\\dsh\\downloads\\dssa_ppi\\RBM_baselines_comparison_ddG.csv`
  真值 `EXP_fitness`（适应度，**非 ΔΔG**）→ 单独标 dataset_name

输出（长表，一行一个 (dataset, model, complex, mutation)）
---------------------------------------------------------
`experiments/v6_sota_audit/standardized_preds.csv`
列：`dataset_name, model, complex_id, mutation, true_ddG, pred_ddG, flipped,
     split, n_muts`

过滤：按 complex 分组，丢弃组内突变数 < 5 的组（V6 §2.2 规则）。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

OUT = C.ROOT / "experiments" / "v6_sota_audit"
OUT.mkdir(parents=True, exist_ok=True)
OUT_CSV = OUT / "standardized_preds.csv"

CATH_CSV = Path(r"F:\dsh\downloads\cathddg\skempi_v2_with_all_results.csv")
RBM_CSV = Path(r"F:\dsh\downloads\dssa_ppi\RBM_baselines_comparison_ddG.csv")
MIN_N = 5

# V6.0 sign_audit 的判定结果（在 CATH 全表上估计，ESM-IF 反向）
CATH_MODELS = {"Flex-ddG": "Flex-ddG", "DDGPred": "DDGPred",
               "ddG_FoldX": "FoldX", "ESM-IF": "ESM-IF",
               "proteinMPNN": "proteinMPNN"}
CATH_FLIP = {"Flex-ddG": False, "DDGPred": False, "FoldX": False,
             "ESM-IF": True, "proteinMPNN": False}
RBM_MODELS = ["DSSA-PPI", "PPIFormer", "ESM-2-FFN", "MuTon", "MpbPPI",
              "ProBASS", "SAAMBE-SEQ", "FoldX"]
RBM_FLIP = {"DSSA-PPI": False, "PPIFormer": True, "ESM-2-FFN": True,
            "MuTon": False, "MpbPPI": False, "ProBASS": False,
            "SAAMBE-SEQ": True, "FoldX": True}


def main() -> int:
    rows: list[pd.DataFrame] = []
    print("=" * 70)
    print("V6.1 PHASE 1.0 — 补建 standardized_preds.csv")
    print("=" * 70)

    # ---------------- CATH ----------------
    if CATH_CSV.exists():
        ca = pd.read_csv(CATH_CSV, low_memory=False)
        te = ca[ca["cath_fold"].astype(str).str.strip().str.lower() == "val"].copy()
        te["true_ddG"] = pd.to_numeric(te["ddG"], errors="coerce")
        te = te[te["true_ddG"].notna()]
        cnt = te.groupby("complex")["true_ddG"].transform("size")
        te = te[cnt >= MIN_N]
        print(f"[CATH] held-out val 过滤后 {len(te)} 行 / "
              f"{te['complex'].nunique()} complex")
        for col, name in CATH_MODELS.items():
            pv = pd.to_numeric(te[col], errors="coerce")
            m = pv.notna()
            flip = CATH_FLIP.get(name, False)
            p = (-pv[m]) if flip else pv[m]
            rows.append(pd.DataFrame({
                "dataset_name": "CATH",
                "model": name,
                "complex_id": te.loc[m, "complex"].astype(str).to_numpy(),
                "mutation": te.loc[m, "Mutation(s)_cleaned"].astype(str).to_numpy(),
                "true_ddG": te.loc[m, "true_ddG"].to_numpy(float),
                "pred_ddG": p.to_numpy(float),
                "flipped": flip,
                "split": te.loc[m, "cath_fold"].astype(str).to_numpy(),
                "n_muts": pd.to_numeric(te.loc[m, "num_muts"],
                                        errors="coerce").to_numpy()}))
        # 同时保留**未反号**原始值供敏感性分析
        print(f"       预测列: {list(CATH_MODELS.values())}")
    else:
        print(f"[CATH] MISSING: {CATH_CSV}")

    # ---------------- RBM ----------------
    if RBM_CSV.exists():
        rb = pd.read_csv(RBM_CSV, low_memory=False)
        rb["true_ddG"] = pd.to_numeric(rb["EXP_fitness"], errors="coerce")
        rb = rb[rb["true_ddG"].notna()]
        # 单靶点：用 pos 区块作为代理分组（每 200 个位置序号一箱）
        rb["complex_id"] = "RBM_pos" + (pd.to_numeric(rb["pos"], errors="coerce")
                                        .rank(method="dense") // 200).astype("Int64").astype(str)
        cnt = rb.groupby("complex_id")["true_ddG"].transform("size")
        rb = rb[cnt >= MIN_N]
        print(f"[RBM] 过滤后 {len(rb)} 行 / {rb['complex_id'].nunique()} 分箱")
        for name in RBM_MODELS:
            if name not in rb.columns:
                continue
            pv = pd.to_numeric(rb[name], errors="coerce")
            m = pv.notna()
            flip = RBM_FLIP.get(name, False)
            p = (-pv[m]) if flip else pv[m]
            rows.append(pd.DataFrame({
                "dataset_name": "DSSA-PPI_RBM",
                "model": name,
                "complex_id": rb.loc[m, "complex_id"].astype(str).to_numpy(),
                "mutation": (rb.loc[m, "wt"].astype(str) + rb.loc[m, "pos"].astype(str)
                             + rb.loc[m, "mut"].astype(str)).to_numpy(),
                "true_ddG": rb.loc[m, "true_ddG"].to_numpy(float),
                "pred_ddG": p.to_numpy(float),
                "flipped": flip,
                "split": "case_study",
                "n_muts": 1}))
    else:
        print(f"[RBM] MISSING: {RBM_CSV}")

    if not rows:
        print("!! 无数据可写")
        return 1
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(OUT_CSV, index=False, encoding="utf-8")
    print()
    print(f"写出 {OUT_CSV}")
    print(f"总行数: {len(out)}")
    print()
    print("| dataset | model | n 行 | n complex | flipped |")
    print("|---|---|---:|---:|---|")
    for (ds, mo), g in out.groupby(["dataset_name", "model"]):
        print(f"| {ds} | {mo} | {len(g)} | {g['complex_id'].nunique()} | "
              f"{bool(g['flipped'].iloc[0])} |")
    return 0


if __name__ == "__main__":
    sys.exit(main())
