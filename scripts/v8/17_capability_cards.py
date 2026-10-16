# -*- coding: utf-8 -*-
"""scripts/v8/17_capability_cards.py -- V8 PHASE 16 / §63：模型能力卡 + §65 Pareto。

对每个**主要**模型键生成 `experiments/v8/cards/<model_id>.md`（§63 模板）。
同时算 Pareto 前沿（§65）：RMSE 越低、Macro ρ 越高、Efficiency@10 越高。

限制：默认只为主要模型族各取代表（避免生成 202 张卡），
可用 `--all` 生成全部。
"""

from __future__ import annotations

import argparse
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
OUT = R / "experiments" / "v8"
CARDS = OUT / "cards"
PARETO = OUT / "pareto"
CARDS.mkdir(parents=True, exist_ok=True)
PARETO.mkdir(parents=True, exist_ok=True)

S = pd.read_csv(OUT / "BENCHMARK_SCOREBOARD.csv", low_memory=False)
B = pd.read_csv(OUT / "bootstrap" / "rank_agreement_bootstrap.csv") \
    if (OUT / "bootstrap" / "rank_agreement_bootstrap.csv").exists() else None


def safe(s: str) -> str:
    return re.sub(r"[^\w\-.@]+", "_", str(s))[:90]


def pareto_front(df, cols, senses):
    """cols/senses: 目标列与方向(+1 越大越好 / -1 越小越好)。返回布尔掩码。"""
    v = df[cols].to_numpy(float) * np.asarray(senses, float)
    n = len(v)
    keep = np.ones(n, bool)
    for i in range(n):
        if not np.all(np.isfinite(v[i])):
            keep[i] = False; continue
        dom = np.all(v >= v[i], axis=1) & np.any(v > v[i], axis=1)
        if dom.any():
            keep[i] = False
    return keep


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()

    for uni, g in S.groupby("universe"):
        g = g.dropna(subset=["rmse", "macro_spearman", "efficiency@10"])
        if len(g) < 3:
            continue
        m = pareto_front(g, ["rmse", "macro_spearman", "efficiency@10"],
                         [-1, +1, +1])
        P = g[m].copy()
        P["universe"] = uni
        P.to_csv(PARETO / f"pareto_{uni}.csv", index=False, encoding="utf-8")
        print(f"[{uni}] Pareto 前沿模型数 = {int(m.sum())} / {len(g)}")

    lines = ["# V8 Pareto-Optimal Models（§65）\n",
             "> 避免过早使用复合评分；此处只给出**非支配集**。\n",
             "> 目标：`RMSE` 越低、`Macro Spearman` 越高、`Efficiency@10` 越高。\n"]
    for f in sorted(PARETO.glob("pareto_*.csv")):
        P = pd.read_csv(f)
        uni = f.stem.replace("pareto_", "")
        lines.append(f"\n## universe `{uni}`（{len(P)} 个非支配模型）\n")
        lines.append("| model | family | RMSE | Macro ρ | Eff@10 | Overall r |")
        lines.append("|---|---|---:|---:|---:|---:|")
        for _, r in P.sort_values("macro_spearman", ascending=False).iterrows():
            lines.append(f"| `{r['model_name']}` | {r['model_family']} | {r['rmse']:.4f} | "
                         f"{r['macro_spearman']:+.4f} | {r['efficiency@10']:+.4f} | "
                         f"{r['overall_pearson']:+.4f} |")
    (OUT / "PARETO_MODELS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"写出 PARETO_MODELS.md")

    # ---- 能力卡 ----
    picks = []
    if a.all:
        picks = [r for _, r in S.iterrows()]
    else:
        for uni, g in S.groupby("universe"):
            g = g.dropna(subset=["macro_spearman"])
            if not len(g):
                continue
            for col, asc in (("macro_spearman", False), ("overall_pearson", False),
                             ("rmse", True), ("efficiency@10", False), ("C_W", False)):
                idx = g[col].idxmin() if asc else g[col].idxmax()
                picks.append(S.loc[idx])
    seen, n = set(), 0
    for rr in picks:
        mid = str(rr["model_id"])
        if mid in seen:
            continue
        seen.add(mid)

        def gv(k, default="NA"):
            v = rr.get(k, default) if hasattr(rr, "get") else default
            return v
        outp = CARDS / (safe(mid) + ".md")
        txt = f"""# Model Capability Card

**Model**: `{gv('model_name')}`
**Family**: {gv('model_family')}
**Provenance status**: {gv('provenance_status')}
**Evidence tier**: {gv('resource_tier')}
**Prediction file**: `{gv('prediction_file')}`
**Join mode**: {gv('join_mode')}

## Universe

Universe: `{gv('universe')}`
Rows: {gv('n_rows')} | Targets: {gv('n_targets')} | n>=5: {gv('n_targets_ge5')} | min_n: {gv('min_n')}

## Calibration

Overall Pearson: {gv('overall_pearson'):.4f}
Overall Spearman: {gv('overall_spearman'):.4f}
RMSE: {gv('rmse'):.4f}
MAE: {gv('mae'):.4f}
Target-mean Pearson: {gv('target_mean_pearson'):.4f}
Target-mean Spearman: {gv('target_mean_spearman'):.4f}
Attenuation slope: {gv('attenuation_slope'):.4f}
C_B row: {gv('C_B_row'):.4f}
C_B equal: {gv('C_B_equal'):.4f}

## Ranking

Macro Spearman: {gv('macro_spearman'):.4f}
Macro Kendall tau-b: {gv('macro_kendall'):.4f}
Pairwise δ=0: {gv('pairwise_d0'):.4f} | δ=0.2: {gv('pairwise_d02'):.4f} | δ=0.5: {gv('pairwise_d05'):.4f}
C_W: {gv('C_W'):.4f}

## Selection

Hit@1: {gv('hit@1'):.4f} | Hit@5: {gv('hit@5'):.4f} | Hit@10: {gv('hit@10'):.4f}
Precision@10: {gv('precision@10'):.4f} | Recall@10: {gv('recall@10'):.4f}
Efficiency@5: {gv('efficiency@5'):.4f} | Efficiency@10: {gv('efficiency@10'):.4f}
NDCG@5: {gv('ndcg@5'):.4f} | NDCG@10: {gv('ndcg@10'):.4f}
EF@1%: {gv('ef@1'):.4f} | EF@5%: {gv('ef@5'):.4f}

## Split sensitivity

本项目对该模型的**划分敏感性**记录于 `SPLIT_INFLATION_REPORT.md`（仅 V6.2 族有跨划分预测）。

## Target-weighting sensitivity

见 `WEIGHTING_SENSITIVITY_REPORT.md`；C_B 对加权敏感，C_W 稳定。

## Main strength / weakness

（依三轴读数人工填写；本卡只提供数据）

## Appropriate use case

Quantitative ΔΔG / Library ranking / Wet-lab shortlist / Other —— **需按三轴分别判断**。

## Caveats

* 本卡为**自动生成**，未做人工解读。
* `C_B`/`C_W` 是**协方差信号贡献**，不是独立标准化能力（§29）。
* 选择效用轴（Eff@10 等）在本 universe 上**区分度有限**（见 `DECISION_ALIGNMENT_REPORT.md`）。
"""
        outp.write_text(txt, encoding="utf-8")
        n += 1
    print(f"生成能力卡 {n} 张 -> {CARDS}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
