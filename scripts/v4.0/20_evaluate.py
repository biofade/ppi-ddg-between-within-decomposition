# -*- coding: utf-8 -*-
"""scripts/v4.0/20_evaluate.py -- V4.0 任务 5: ΔΔG 标准化评估脚本 + benchmark 打包.

用法
----
    python scripts/v4.0/20_evaluate.py --predictions <你的预测.npz> \
        --name <模型名> [--split-file data/v4.0/ddg_split.json]

或（自检）:
    python scripts/v4.0/20_evaluate.py --self-test

`--predictions` 指向的 npz 必须含:
* `test_y`      : (n_test,) 真值 ΔΔG
* `<name>`      : (n_test,) 预测值（顺序与 `test_y` 一致）
* 可选 `test_group` : 复合物 ID，用于 per-structure Pearson

评估协议（本 benchmark 的**唯一**口径，供后续工作引用）
-------------------------------------------------------
* **任务**: `ΔΔG = ΔG_mutant − ΔG_wild-type`（kcal/mol）
* **split**: 复合物级 + 40% 序列同一性 + 家族级三重隔离
  （见 `data/v4.0/ddg_split.json` 与 `features/v4.0_ddg_split_report.md`）
* **指标**: Pearson（主，与文献可比）/ Spearman / RMSE / MAE /
  **per-structure Pearson**（SAMAffinity / Pythia-PPI 口径）/
  AUROC（ΔΔG>0 = 削弱结合 为正类）
* **统计**: 配对 bootstrap 5000 次 + BH-FDR（多个模型对比时）

**不能声称"比文献好"**：文献（Pythia-PPI / SAMAffinity / DSSA-PPI）用的划分都比
本 benchmark **松**（见 `features/v4.0_ddg_baseline_report.md` 第 4 节的逐条对照）。
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v40_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

BENCH = C.V40 / "ddg_benchmark"


def evaluate(y: np.ndarray, p: np.ndarray, groups: np.ndarray | None = None) -> dict:
    return C.all_metrics(y, p, groups)


def build_benchmark() -> None:
    """打包可复现 benchmark: 数据 + split + README。"""
    BENCH.mkdir(parents=True, exist_ok=True)
    fwd = pd.read_csv(C.PAIRS, low_memory=False)
    fwd = fwd[fwd["direction"] == "forward"].reset_index(drop=True)
    split = json.loads(C.SPLIT.read_text(encoding="utf-8"))
    # 同步 04/10 的过滤（v2.2 映射失败的 3 对）—— 以 split.json 的 pair_index 为准
    keep = sorted(set(sum(split["pair_index"].values(), [])))
    sub = fwd.loc[keep].reset_index(drop=True)
    # 重新编号 split 索引
    remap = {old: new for new, old in enumerate(keep)}
    new_split = {s: sorted(remap[i] for i in idx)
                 for s, idx in split["pair_index"].items()}
    sub.to_csv(BENCH / "ddg_benchmark_pairs.csv", index=False, encoding="utf-8")
    (BENCH / "split.json").write_text(json.dumps(
        {"rules": split["rules"], "counts": {s: len(v) for s, v in new_split.items()},
         "pair_index": new_split,
         "family_of_complex": split["family_of_complex"]},
        ensure_ascii=False, indent=1), encoding="utf-8")
    readme = f"""# V4.0 ΔΔG Benchmark

## 任务
预测 `ΔΔG = ΔG_mutant − ΔG_wild-type`（kcal/mol）。
`ΔΔG > 0` 表示突变**削弱**结合；`< 0` 表示**增强**。

## 数据
* `ddg_benchmark_pairs.csv` — **{len(sub)}** 个 WT×mutant 配对，**{sub['pdb_id'].nunique()}** 个复合物。
  列含 `pdb_id, wt_row, mut_row, split, dg_wt, dg_mut, ddg, direction, n_mut,
  mutation, ligand_seq_changed, receptor_seq_changed, temperature_wt,
  temperature_mut, affinity_method_wt, affinity_method_mut, source`。
  全部来自 **SKEMPI v2.0**（`{C.RAW.name}`）。
* 原始来源: PPB-Affinity / Proteinea `filtered.csv`（8207 行）中的 `Mutations` 字段。

## Split（`split.json`）
`pair_index` 给出 train/val/test 的行号（指向 `ddg_benchmark_pairs.csv`）。

三重隔离规则（全部已用精确判据验证，见 `features/v4.0_ddg_split_report.md`）:
1. **复合物级** — 同一 `pdb_id` 的全部配对在同一侧（交集 0）。
2. **序列同一性 ≤40%** — 单链接聚类，3-mer Jaccard 为同源代理；
   **完全相同序列总是合并**。
   > 3-mer Jaccard 会**高估**真实同一性（V3.1 结论），故本判据是**保守方向**：
   > 可能过度合并（split 更小），但不会低估泄漏。
3. **家族级** — 家族 = 「共享任一链簇」的并查集连通分量（**不是**简单的
   (lig_cl, rec_cl) 组合；后者会把共享受体簇的复合物误判为不同家族，
   实测导致 6 条完全相同的序列跨 train/test）。家族不跨侧。

规模: {new_split['train'].__len__()} / {new_split['val'].__len__()} / {new_split['test'].__len__()}
（train/val/test）。

## 指标
* **Pearson**（主指标，与文献可比）
* Spearman、RMSE、MAE
* **per-structure Pearson** — 每个复合物内单独算 Pearson 再平均
  （SAMAffinity / Pythia-PPI 口径；默认丢弃对/结构 < {C.MIN_PS_N} 的复合物并报告数量）
* **AUROC** — `ΔΔG > 0`（削弱结合）为正类

## 评估
```bash
python scripts/v4.0/20_evaluate.py --predictions <pred.npz> --name <模型名>
```
预测 npz 需含 `test_y` 与同名预测数组（顺序一致），可选 `test_group`。

## 与文献的可比性（重要）
本 benchmark 的划分比 Pythia-PPI（per-structure 五折，**含训练折内复合物**）、
SAMAffinity（mutation-level，作者自承乐观）、DSSA-PPI（MpbPPI split）**都更严**。
因此**不得**声称"比文献好/差"；可比区间是 SAMAffinity complex-level
（Pearson 0.54–0.59）与 DSSA-PPI complex-level（PCC 0.58–0.71），且本 benchmark 更严。

## 已知局限
* 仅 **231 个复合物**（SKEMPI 子集），test 仅 61 个复合物 -> per-structure 统计噪声大。
* 59% 的配对来自**丙氨酸扫描**，`ΔΔG > 0` 占 74.5%（生物学上合理，但分布高度偏斜）。
* 条件混杂: **仅 54.6% 的配对**温度一致、57.2% 测量方法一致；
  `Temperature(K)` 列混有数值 `298` 与字符串 `"298(assumed)"`。
* 未做真实序列比对（无 mmseqs2/BLAST），序列隔离依赖 3-mer Jaccard 代理。
"""
    (BENCH / "README.md").write_text(readme, encoding="utf-8")
    print(f"benchmark -> {BENCH}")


def main() -> int:
    ap = argparse.ArgumentParser(description="V4.0 ΔΔG 标准化评估")
    ap.add_argument("--predictions", type=Path)
    ap.add_argument("--name", default="model")
    ap.add_argument("--build-benchmark", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.build_benchmark:
        build_benchmark()
        return 0

    if args.self_test:
        # 自检: 用已知的相关结构验证指标实现
        rs = np.random.RandomState(0)
        y = rs.normal(size=200)
        ok = True
        for noise, lo, hi in ((0.0, 0.99, 1.0), (1.0, 0.6, 0.8)):
            p = y + rs.normal(scale=noise, size=200)
            r = C.pearson(y, p)
            good = lo <= r <= hi
            ok &= good
            print(f"  noise={noise}: pearson={r:+.4f} (expect {lo}-{hi}) "
                  f"{'OK' if good else 'FAIL'}")
        p = -y
        print(f"  anti-correlated: pearson={C.pearson(y, p):+.4f} (expect -1.0)")
        ok &= abs(C.pearson(y, p) + 1.0) < 1e-9
        print(f"  AUROC(perfect)={C.auroc(y, y):.4f} (expect 1.0)")
        ok &= abs(C.auroc(y, y) - 1.0) < 1e-9
        print(f"  AUROC(reversed)={C.auroc(y, -y):.4f} (expect 0.0)")
        ok &= abs(C.auroc(y, -y)) < 1e-9
        # per-structure: 设计成组内完美相关
        g = np.repeat(np.arange(20), 10)
        yy = rs.normal(size=200)
        pp = np.array([yy[g == k] * (k + 1) + k for k in range(20)]).ravel()
        ps = C.per_structure_pearson(yy, pp, g, min_n=5)
        print(f"  per-structure pearson={ps['per_structure_pearson']:+.4f} "
              f"(expect 1.0, n={ps['per_structure_n']})")
        ok &= abs(ps["per_structure_pearson"] - 1.0) < 1e-9
        print(f"\nSELF-TEST {'PASSED' if ok else 'FAILED'}")
        return 0 if ok else 1

    if not args.predictions:
        ap.print_help()
        return 2
    z = np.load(args.predictions, allow_pickle=True)
    if "test_y" not in z.files:
        print("!! npz 缺少 `test_y`")
        return 3
    y = z["test_y"].astype(float)
    if args.name not in z.files:
        cands = [k for k in z.files if k != "test_y"]
        print(f"!! npz 中没有 `{args.name}`；可用: {cands}")
        return 3
    p = z[args.name].astype(float)
    g = z["test_group"] if "test_group" in z.files else None
    m = evaluate(y, p, g)
    print(f"### {args.name}")
    for k in ("n", "pearson", "spearman", "rmse", "mae", "auroc",
              "per_structure_pearson", "per_structure_n", "per_structure_skipped"):
        if k in m:
            v = m[k]
            print(f"- {k}: {v:+.4f}" if isinstance(v, float) else f"- {k}: {v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
