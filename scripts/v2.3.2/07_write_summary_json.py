# -*- coding: utf-8 -*-
"""scripts/v2.3.2/07_write_summary_json.py -- 补写 V2.3.2 LightGBM 消融的汇总 JSON.

背景
----
`03_merge_esm2_v2.py` 的 JSON 汇总段在早期版本里引用了未定义的旧变量
(`use_iface`), 导致**结果 CSV / 配对检验 / 最佳模型都已正确落盘, 但 JSON 未写出**。
修复变量名后需要重跑 18 分钟才能补这一个文件 —— 而 CSV 里的数字已经是最终结果,
无需重算。因此本脚本**只读已落盘的产物**重新组装同一个 JSON,
并显式标注它是"从产物重建", 而不是假装是那次运行直接写的。

**不重新计算任何指标。**
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
FEAT = ROOT / "features"
OUT = FEAT / "v2.3.2_lgbm_results.json"
CSV = FEAT / "v2.3.2_lgbm_results.csv"
PAIRS = FEAT / "v2.3.2_lgbm_pairs.csv"
PREDS = FEAT / "v2.3.2_test_predictions.npz"
JOIN_MD = FEAT / "v2.3.2_join_diagnosis.md"

V22_BASELINE = 0.23747781344414348


def main() -> int:
    if not CSV.exists():
        print(f"!! 缺少 {CSV}")
        return 2
    res = pd.read_csv(CSV)
    pairs = pd.read_csv(PAIRS) if PAIRS.exists() else pd.DataFrame()
    base = res[res["feature_set"] == "physics_only"]
    best = res.sort_values("spearman", ascending=False).iloc[0]
    base_rho = float(base["spearman"].iloc[0]) if len(base) else float("nan")

    n_sig = int((pairs["bh_q"] < 0.05).sum()) if len(pairs) else 0
    n_better = int(((pairs["bh_q"] < 0.05) & (pairs["delta"] > 0)).sum()) if len(pairs) else 0

    payload = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "scheme": "v2.3.2",
        "rebuild_note": ("本文件由 `07_write_summary_json.py` 从已落盘的 "
                         "v2.3.2_lgbm_results.csv / _pairs.csv **重建** —— "
                         "03_merge_esm2_v2.py 的 JSON 段曾因 `use_iface` 未定义而"
                         "未写出 (结果 CSV 已正确落盘)。此处不重算任何指标。"),
        "join_key": "pdb_id",
        "join_note": ("V2.3.2: join 键由 interface_key 改为 pdb_id "
                      "(V2.3 的 interface_key 语义不一致 -> 0% 命中)。"
                      "一对多时缓存侧取平均, 特征表侧共享 -> 结构级退化。"),
        "esm2_run": True,
        "esm2_cache_present": PREDS.exists(),
        "n_configs": int(len(res)),
        "n_pair_tests": int(len(pairs)),
        "n_pair_tests_bh_sig": n_sig,
        "n_pair_tests_bh_sig_and_better": n_better,
        "session": {
            "train_rows": 6482, "val_rows": 965, "test_rows": 752,
            "source": "V2.2 Proteinea split (reused, asserted)",
        },
        "v22_recorded_test_spearman": V22_BASELINE,
        "baseline_reproduced_spearman": base_rho,
        "baseline_abs_diff": abs(base_rho - V22_BASELINE),
        "baseline_reproduced_bit_exact": bool(abs(base_rho - V22_BASELINE) < 1e-12),
        "best_config": {
            "feature_set": str(best["feature_set"]),
            "embedding_variant": str(best.get("embedding_variant", "-")),
            "pca_dim": int(best["pca_dim"]),
            "n_features": int(best["n_features"]),
            "spearman": float(best["spearman"]),
            "spearman_ci": [float(best["spearman_ci_low"]),
                            float(best["spearman_ci_high"])],
            "ndcg@10": float(best["ndcg@10"]),
            "ef@1%": float(best["ef@1%"]),
            "precision@10": float(best["precision@10"]),
            "delta_vs_baseline": float(best["spearman"] - base_rho),
        },
        "best_top_enrichment_config": None,
        "notes": [
            "PCA 只在训练集 fit (硬约束)。",
            "配对 bootstrap 5000 次 + BH-FDR; 每个模型自身 CI 用 600 次 (算力预算)。",
            "ESM-2 为 mean-pool over residues (非界面残基加权)。",
            "超参完全沿用 V2.2, 未针对 ESM-2 调参。",
        ],
    }
    te = res.sort_values(["ef@1%", "precision@10"], ascending=False).iloc[0]
    payload["best_top_enrichment_config"] = {
        "feature_set": str(te["feature_set"]),
        "embedding_variant": str(te.get("embedding_variant", "-")),
        "pca_dim": int(te["pca_dim"]),
        "spearman": float(te["spearman"]),
        "ef@1%": float(te["ef@1%"]),
        "precision@10": float(te["precision@10"]),
        "ndcg@10": float(te["ndcg@10"]),
        "note": ("全局 Spearman 最优 != 头部富集最优。该配置 Spearman 较低但"
                 "EF@1%/P@10 明显更高 —— 对只做得起 10 个克隆的湿实验可能更有用。"),
    }
    if pairs is not None and len(pairs):
        b = pairs[pairs["model"] ==
                  f"{best['feature_set']}@{best.get('embedding_variant','-')}"
                  f"@{int(best['pca_dim'])}"]
        if len(b):
            r = b.iloc[0]
            payload["best_config"]["paired_test_vs_baseline"] = {
                "delta": float(r["delta"]), "ci": [float(r["lo"]), float(r["hi"])],
                "p_two_sided": float(r["p_two_sided"]), "bh_q": float(r["bh_q"]),
                "n_boot": 5000, "n_distinct_deltas": int(r.get("n_distinct_deltas", -1)),
            }

    OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"写出 {OUT} ({OUT.stat().st_size} B)")
    print(f"  baseline reproduced: {base_rho:+.6f} (V2.2 {V22_BASELINE:.6f}, "
          f"abs diff {abs(base_rho - V22_BASELINE):.2e})")
    print(f"  best: {payload['best_config']['feature_set']} "
          f"{payload['best_config']['spearman']:+.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
