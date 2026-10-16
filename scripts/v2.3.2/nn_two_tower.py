# -*- coding: utf-8 -*-
"""scripts/v2.3.2/nn_two_tower.py -- V2.3.2 任务 3 / N2: 双塔 MLP (判决实验).

规格 (任务书逐字)
-----------------
* 配体塔: ``Linear(640,256) -> ReLU -> Dropout -> Linear(256,64)``
* 受体塔: 同上
* 物理塔: ``Linear(D_phys,64) -> ReLU``
* 拼接 192 维 -> ``Linear(192,64) -> ReLU -> Linear(64,1)``
* 其余训练设置与 N1 相同 (Adam 1e-3 / wd 1e-4 / MSE / batch 64 /
  max 300 epoch / 早停 val Spearman 20 epoch / 种子 42-45)
* **CPU only** (硬约束 7)

关于塔的输入维度 (重要说明)
---------------------------
任务书写的是 ``Linear(640,256)``, 即塔直接吃**原始 640 维** ESM-2 嵌入。
但 640 维 x 6500 训练行属于**过参数化**, 因此本脚本两种都跑:

* ``--pca-dims 32 64 128``: 塔吃 **PCA 后的 32/64/128 维** (与 N1 口径一致, 可比);
* ``--raw640``: 塔吃**原始 640 维** (任务书字面规格)。

> ``physics+both`` 时的 ``D_phys``: 物理塔固定只吃 60 个物理特征 (与任务书的
> "物理塔"语义一致), ESM-2 的 ligand/receptor 分别进各自塔, 不重复进物理塔。

输出
----
* ``features/v2.3.2_nn_twotower_results.csv``
* ``features/v2.3.2_nn_twotower_predictions.npz``
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
    ap = argparse.ArgumentParser(description="V2.3.2 N2: 双塔 MLP (CPU only)")
    ap.add_argument("--esm2-cache", type=Path, default=C.ESM2_CACHE)
    ap.add_argument("--v22-features", type=Path, default=C.V22_FEATURES)
    ap.add_argument("--pca-dims", type=int, nargs="*", default=list(C.PCA_DIMS))
    ap.add_argument("--raw640", action="store_true",
                    help="额外跑一次塔吃原始 640 维 (任务书字面规格)")
    ap.add_argument("--variant", default="raw", choices=("raw", "l2"))
    ap.add_argument("--seeds", type=int, nargs="*", default=list(C.SEEDS))
    ap.add_argument("--out-csv", type=Path,
                    default=C.FEAT_DIR / "v2.3.2_nn_twotower_results.csv")
    ap.add_argument("--preds-npz", type=Path,
                    default=C.FEAT_DIR / "v2.3.2_nn_twotower_predictions.npz")
    args = ap.parse_args()

    C.LOG_DIR.mkdir(parents=True, exist_ok=True)
    import torch
    import torch.nn as nn

    print("=" * 78)
    print("V2.3.2 N2: 双塔 MLP 判决实验")
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
    d_phys = Xp.shape[1]

    class TwoTower(nn.Module):
        """任务书逐字规格: 配体塔 + 受体塔 + 物理塔 -> 192 -> 64 -> 1。"""

        def __init__(self, d_lig: int, d_rec: int, d_phys_in: int):
            super().__init__()
            self.lig = nn.Sequential(
                nn.Linear(d_lig, 256), nn.ReLU(), nn.Dropout(C.DROPOUT),
                nn.Linear(256, 64),
            )
            self.rec = nn.Sequential(
                nn.Linear(d_rec, 256), nn.ReLU(), nn.Dropout(C.DROPOUT),
                nn.Linear(256, 64),
            )
            self.phys = nn.Sequential(nn.Linear(d_phys_in, 64), nn.ReLU())
            self.head = nn.Sequential(
                nn.Linear(192, 64), nn.ReLU(), nn.Linear(64, 1),
            )

        def forward(self, batch):
            h = torch.cat([self.lig(batch["lig"]), self.rec(batch["rec"]),
                           self.phys(batch["phys"])], dim=1)
            return self.head(h)

    rows: list[dict] = []
    preds: dict[str, np.ndarray] = {}

    for dtag in dim_tags:
        klig = f"{args.variant}|ligand|{dtag}"
        krec = f"{args.variant}|receptor|{dtag}"
        if klig not in blocks or krec not in blocks:
            print(f"!! 跳过缺失的块 {klig} / {krec}")
            continue
        L, R = blocks[klig], blocks[krec]
        d_lig, d_rec = L.shape[1], R.shape[1]
        print(f"\n--- 双塔输入: lig {d_lig} + rec {d_rec} + phys {d_phys} "
              f"-> 拼接 192 -> 64 -> 1 ---")
        for seed in args.seeds:
            try:
                r = C.train_nn(
                    lambda n_in: TwoTower(n_in["d_lig"], n_in["d_rec"], n_in["d_phys"]),
                    inputs={},
                    n_in={"d_lig": d_lig, "d_rec": d_rec, "d_phys": d_phys},
                    y_tr=y[tr], y_va=y[va],
                    X_tr_parts={"lig": L[tr], "rec": R[tr], "phys": Xp[tr]},
                    X_va_parts={"lig": L[va], "rec": R[va], "phys": Xp[va]},
                    X_te_parts={"lig": L[te], "rec": R[te], "phys": Xp[te]},
                    seed=seed)
                y_s, p_s = -y[te], -r["test_pred"]
                met = {}
                for k, fn in mfns.items():
                    try:
                        met[k] = float(fn(y_s, p_s))
                    except Exception:  # noqa: BLE001
                        met[k] = float("nan")
                lo, hi, _ = C.bootstrap_ci(mfns["spearman"], y_s, p_s, 600, seed)
                rows.append({
                    "model": "TwoTowerMLP", "pca_dim": dtag,
                    "input_dim": d_lig + d_rec + d_phys,
                    "d_lig": d_lig, "d_rec": d_rec, "d_phys": d_phys,
                    "embedding": f"two_tower_{args.variant}", "seed": seed,
                    "best_epoch": r["best_epoch"], "epochs_ran": r["epochs_ran"],
                    "train_seconds": round(r["train_seconds"], 2),
                    "val_spearman": float(r["best_val_spearman"]),
                    "spearman_ci_low": lo, "spearman_ci_high": hi,
                    **met,
                })
                preds[f"TwoTowerMLP|{dtag}|seed{seed}"] = r["test_pred"].astype(np.float32)
                print(f"  seed {seed}: val_rho={r['best_val_spearman']:+.4f}  "
                      f"test_rho={met['spearman']:+.4f}  "
                      f"({r['epochs_ran']} ep, {r['train_seconds']:.1f}s)")
            except Exception as exc:  # noqa: BLE001
                print(f"  seed {seed}: FAIL {type(exc).__name__}: {exc}")
                with open(C.FAIL_LOG, "a", encoding="utf-8") as fh:
                    fh.write(f"[nn_two_tower] {dtag} seed{seed}: {exc}\n"
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
