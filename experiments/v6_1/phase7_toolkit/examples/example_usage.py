# -*- coding: utf-8 -*-
"""affirigid_metrics 使用示例（合成数据自检，可直接运行）.

运行:
    python example_usage.py

预期输出:
    oracle 模型的 efficiency_at_10 ≈ +1.0, macro Spearman ≈ +1.0
    anti 模型的 efficiency_at_10 < 0,   macro Spearman ≈ -1.0
    random 模型接近 0 / 0
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from affirigid_metrics import bootstrap, topk  # noqa: E402


def make_data(n_complex: int = 30, per_complex: int = 25, seed: int = 0):
    rs = np.random.RandomState(seed)
    rows = []
    for c in range(n_complex):
        true = rs.normal(size=per_complex) * 2.0 + 1.0   # ΔΔG，越小越好
        rows.append(pd.DataFrame({
            "complex_id": f"CX{c:03d}",
            "true_ddG": true,
            # oracle: 预测 == 真值
            "pred_oracle": true,
            # anti: 完全反向
            "pred_anti": -true,
            # noisy: 部分相关
            "pred_noisy": true * 0.6 + rs.normal(size=per_complex) * 1.5,
            # random: 无关
            "pred_random": rs.normal(size=per_complex),
        }))
    return pd.concat(rows, ignore_index=True)


def main() -> int:
    df = make_data()
    print(f"合成数据: {len(df)} 行 / {df.complex_id.nunique()} 复合物")
    print(f"Base rate (ΔΔG<0) = {(df.true_ddG < 0).mean():.4f}")
    print()

    print("=" * 74)
    print("1) topk_metrics —— 单次评估")
    print("=" * 74)
    print(f"{'模型':<14} {'macro_pdb_ρ':>12} {'p_at_10':>9} {'ef_1pct':>9} "
          f"{'hit10(<0)':>10} {'Efficiency@10':>14}")
    for col in ("pred_oracle", "pred_noisy", "pred_anti", "pred_random"):
        m = topk.topk_metrics(df.true_ddG, df[col], df.complex_id)
        print(f"{col:<14} {m['macro_pdb_spearman']:>12.4f} {m['p_at_10']:>9.4f} "
              f"{m['ef_1pct']:>9.4f} {m['hit10_beneficial']:>10.4f} "
              f"{m['efficiency_at_10']:>14.4f}")
    print()
    print("读法: oracle 应 ≈1.0 / anti 应 <0 / random 应 ≈0。")
    print("     Efficiency@10 = 0 表示与随机无异，1 表示达到 oracle。")
    print()

    print("=" * 74)
    print("2) cluster_bootstrap_ci —— 把「复合物」当作抽样单位")
    print("=" * 74)
    for col in ("pred_oracle", "pred_noisy", "pred_anti"):
        ci = bootstrap.cluster_bootstrap_ci(
            topk.per_structure_spearman,
            df.true_ddG.values, df[col].values, df.complex_id.values,
            n_boot=2000, seed=42, min_groups=5)
        print(f"{col:<14} point={ci['point']:+.4f}  "
              f"95% CI=[{ci['ci_low']:+.4f}, {ci['ci_high']:+.4f}]  "
              f"n_clusters={ci['n_clusters']}  valid={ci['n_valid']}")
    print()
    print("为什么必须 cluster bootstrap：同一复合物的多个突变体不独立；")
    print("把 N 行当 N 个独立样本会低估 CI、产生假显著。")
    print()

    print("=" * 74)
    print("3) 配对比较（同一批复合物上的两个模型）")
    print("=" * 74)
    # 逐复合物取 per-complex Spearman，再做配对检验
    per_cx = {}
    for col in ("pred_oracle", "pred_noisy"):
        vals = []
        for cx, gg in df.groupby("complex_id"):
            vals.append(topk.per_structure_spearman(
                gg.true_ddG.values, gg[col].values, np.array([cx] * len(gg)),
                min_n=5))
        per_cx[col] = np.array(vals, float)
    w = bootstrap.paired_wilcoxon(per_cx["pred_oracle"], per_cx["pred_noisy"])
    print(f"配对 Wilcoxon (oracle vs noisy): W={w['W']:.1f}  p={w['p']:.4g}  "
          f"n={w['n']}  Δmean={w['delta_mean']:+.4f}")
    print()
    print("注意: 此时 **不能用 Mann-Whitney U** —— 两个模型作用于同一批突变体，")
    print("      两组不独立，必须用配对检验（Wilcoxon 符号秩）。")
    print()

    # 断言自检
    ok = True
    mo = topk.topk_metrics(df.true_ddG, df.pred_oracle, df.complex_id)
    ma = topk.topk_metrics(df.true_ddG, df.pred_anti, df.complex_id)
    ok &= abs(mo["efficiency_at_10"] - 1.0) < 0.05
    ok &= ma["efficiency_at_10"] < 0
    ok &= mo["macro_pdb_spearman"] > 0.95
    print("=" * 74)
    print(f"自检: {'全部通过 ✅' if ok else '失败 ❌'}")
    print("=" * 74)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
