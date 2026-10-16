# -*- coding: utf-8 -*-
"""scripts/v2.3.2/nn_mlp.py -- V2.3.2 任务 3 / N1: MLP 基线 (判决实验).

规格 (任务书逐字)
-----------------
* 输入: ``physics + esm2_both``
* 结构: ``Linear(D,256) -> ReLU -> Dropout(0.3) -> Linear(256,64) -> ReLU
  -> Dropout(0.3) -> Linear(64,1)``
* Adam ``lr=1e-3``, ``weight_decay=1e-4``; MSE; ``batch=64``; ``max 300 epoch``;
  早停准则: 验证集 **Spearman** 连续 20 epoch 不提升
* 种子: 42 / 43 / 44 / 45
* **CPU only** (硬约束 7), 禁用 `.cuda()`

判决目标
--------
验证「**只有 pooled 嵌入时, NN 是否能打过 LightGBM**」。
不跑花哨方案 —— 只有这一个 MLP。

输入形态 (两种都跑, 让数据说话)
--------------------------------
* ``--pca-dims 32 64 128``: 输入 = 物理 60 + PCA(both) 32/64/128
* ``--raw640``: 输入 = 物理 60 + **原始 640x2 = 1280 维** (无 PCA)

输出
----
* ``features/v2.3.2_nn_mlp_results.csv``
* ``features/v2.3.2_nn_mlp_predictions.npz``
"""

from __future__ import annotations

import argparse
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nn_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass


def main() -> int:
    ap = argparse.ArgumentParser(description="V2.3.2 N1: MLP 判决实验 (CPU only)")
    ap.add_argument("--esm2-cache", type=Path, default=C.ESM2_CACHE)
    ap.add_argument("--v22-features", type=Path, default=C.V22_FEATURES)
    ap.add_argument("--pca-dims", type=int, nargs="*", default=list(C.PCA_DIMS))
    ap.add_argument("--raw640", action="store_true",
                    help="额外跑一次不做 PCA 的 1280 维输入")
    ap.add_argument("--variant", default="raw", choices=("raw", "l2"))
    ap.add_argument("--seeds", type=int, nargs="*", default=list(C.SEEDS))
    ap.add_argument("--out-csv", type=Path,
                    default=C.FEAT_DIR / "v2.3.2_nn_mlp_results.csv")
    ap.add_argument("--preds-npz", type=Path,
                    default=C.FEAT_DIR / "v2.3.2_nn_mlp_predictions.npz")
    args = ap.parse_args()

    C.LOG_DIR.mkdir(parents=True, exist_ok=True)
    import torch
    import torch.nn as nn

    print("=" * 78)
    print("V2.3.2 N1: MLP 判决实验")
    print("=" * 78)
    print(f"torch {torch.__version__} | CPU only (硬约束 7)")

    t0 = time.time()
    data = C.load_data(args.esm2_cache, args.v22_features)
    print(f"数据: {len(data['df'])} 行 | 划分 {data['split_counts']}")
    mr = C.get_metrics_module()
    mfns = C.build_metric_fns(mr, C.TOP_FRAC)

    dim_tags: list[str] = [str(d) for d in args.pca_dims]
    if args.raw640:
        dim_tags.append("raw640")
    blocks = C.pca_blocks(data, dims=tuple(args.pca_dims), seed=C.SEEDS[0],
                          include_raw640=args.raw640)

    tr, va, te = data["tr_m"], data["va_m"], data["te_m"]
    y = data["y"]
    Xp = data["Xp"].to_numpy(np.float32)

    class MLP(nn.Module):
        """任务书逐字规格。"""

        def __init__(self, d_in: int):
            super().__init__()
            self.net = nn.Sequential(
                nn.Linear(d_in, 256), nn.ReLU(), nn.Dropout(C.DROPOUT),
                nn.Linear(256, 64), nn.ReLU(), nn.Dropout(C.DROPOUT),
                nn.Linear(64, 1),
            )

        def forward(self, batch):
            return self.net(batch["x"])

    rows: list[dict] = []
    preds: dict[str, np.ndarray] = {}

    for dtag in dim_tags:
        key = f"{args.variant}|both|{dtag}"
        if key not in blocks:
            print(f"!! 跳过缺失的块 {key}")
            continue
        Z = blocks[key]
        X = np.hstack([Xp, Z]).astype(np.float32)
        d_in = X.shape[1]
        print(f"\n--- 输入 {d_in} 维 (physics 60 + both_{args.variant} {Z.shape[1]}) ---")
        for seed in args.seeds:
            try:
                r = C.train_nn(lambda n_in: MLP(n_in["d_in"]),
                               inputs={}, n_in={"d_in": d_in},
                               y_tr=y[tr], y_va=y[va],
                               X_tr_parts={"x": X[tr]}, X_va_parts={"x": X[va]},
                               X_te_parts={"x": X[te]}, seed=seed)
                y_s, p_s = -y[te], -r["test_pred"]
                met = {}
                for k, fn in mfns.items():
                    try:
                        met[k] = float(fn(y_s, p_s))
                    except Exception:  # noqa: BLE001
                        met[k] = float("nan")
                lo, hi, _ = C.bootstrap_ci(mfns["spearman"], y_s, p_s, 600, seed)
                rows.append({
                    "model": "MLP", "pca_dim": dtag, "input_dim": d_in,
                    "embedding": f"physics+both_{args.variant}", "seed": seed,
                    "best_epoch": r["best_epoch"], "epochs_ran": r["epochs_ran"],
                    "train_seconds": round(r["train_seconds"], 2),
                    "val_spearman": float(r["best_val_spearman"]),
                    "spearman_ci_low": lo, "spearman_ci_high": hi,
                    **met,
                })
                preds[f"MLP|{dtag}|seed{seed}"] = r["test_pred"].astype(np.float32)
                print(f"  seed {seed}: val_rho={r['best_val_spearman']:+.4f}  "
                      f"test_rho={met['spearman']:+.4f}  "
                      f"({r['epochs_ran']} ep, {r['train_seconds']:.1f}s)")
            except Exception as exc:  # noqa: BLE001
                print(f"  seed {seed}: FAIL {type(exc).__name__}: {exc}")
                with open(C.FAIL_LOG, "a", encoding="utf-8") as fh:
                    fh.write(f"[nn_mlp] {dtag} seed{seed}: {exc}\n"
                             f"{traceback.format_exc()}\n")

    if rows:
        pd.DataFrame(rows).to_csv(args.out_csv, index=False, encoding="utf-8")
        np.savez_compressed(args.preds_npz, test_y=y[te].astype(np.float32), **preds)
        print(f"\n写出 {args.out_csv.name} ({len(rows)} 行)")
        print(f"写出 {args.preds_npz.name}")
    print(f"总耗时 {time.time() - t0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
