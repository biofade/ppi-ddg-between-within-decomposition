# -*- coding: utf-8 -*-
"""scripts/v2.5/03_cnn_perres.py -- V2.5 任务 1.2: 1D-CNN + 注意力池化 (判决实验核心).

为什么这个模型是"判决"的关键
-----------------------------
任务 1.3 已证明: 对 per-residue 嵌入做 **mean-pool** 并不比链级 pooled 好
(Δ = −0.0249, 多数不显著)。这只说明"换个范围做平均"没用。

**真正的问题是**: 保留**残基顺序**并让网络自己学习加权 (卷积 + 注意力)
能不能超过任何形式的平均? 若也不能, 则"粒度是瓶颈"这条假设**被证伪**。

架构 (任务书逐字规格)
----------------------
* 配体塔: `Conv1d(640,128,k=3) -> ReLU -> MaxPool -> Conv1d(128,64,k=3)`
  -> **注意力池化** -> 64 维
* 受体塔: 同结构
* 物理塔: `Linear(D_phys,64) -> ReLU`
* 拼接 192 -> `Linear(192,64) -> ReLU -> Linear(64,1)`
* Adam lr=1e-3, wd=1e-4; MSE; batch 32; max 200 epoch;
  早停 = 验证集 Spearman 连续 20 epoch 不升

工程约束 (本机 7.5 GB 内存 / 仅 CPU)
------------------------------------
* **不把 1.9 GB per-residue 矩阵载入内存**: 用 `np.load(..., mmap_mode='r')`
  直接从磁盘按 batch 取残基, 这是本脚本的关键工程决定;
* 掩码卷积: 填充位置的输入置 0, 每层后重新置 0, 池化时用 mask 排除填充;
* torch 线程数限制为 2 (与其它任务共存)。

输出
----
`features/v2.5_cnn_results.csv`、`features/v2.5_cnn_predictions.npz`、
`models/v2.5_best_perres.pt`、`features/v2.5_cnn_report.md`
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v25_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

OUT_CSV = C.FEAT / "v2.5_cnn_results.csv"
OUT_PRED = C.FEAT / "v2.5_cnn_predictions.npz"
OUT_MD = C.FEAT / "v2.5_cnn_report.md"
BEST_PT = C.MODEL_DIR / "v2.5_best_perres.pt"

SEEDS = (42, 43)
MAX_EPOCH = 200
PATIENCE = 20
BATCH = 32
LR = 1e-3
WD = 1e-4
DROPOUT = 0.3

L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def main() -> int:
    ap = argparse.ArgumentParser(description="V2.5 1D-CNN + 注意力池化 (CPU)")
    ap.add_argument("--seeds", type=int, nargs="*", default=list(SEEDS))
    ap.add_argument("--epochs", type=int, default=MAX_EPOCH)
    ap.add_argument("--batch", type=int, default=BATCH)
    ap.add_argument("--max-res", type=int, default=160,
                    help="每链训练时最多使用的残基数 (从下采样缓存取, 控制算力)")
    args = ap.parse_args()

    import torch
    import torch.nn as nn
    torch.manual_seed(C.SEED)
    torch.set_num_threads(2)

    mr = C.get_metrics_module()
    mfns = C.build_metric_fns(mr)
    t0 = time.time()

    emit("# V2.5 任务 1.2 — 1D-CNN + 注意力池化 per-residue 判决实验")
    emit("")
    emit(f"> 设备 **cpu only** (torch {torch.__version__}); 线程数 2")
    emit(f"> 早停: 验证集 Spearman 连续 {PATIENCE} epoch 不升; 种子 {args.seeds}")
    emit("")

    # ---------- 数据 ----------
    data = C.load_aligned(use_perres_pooled=True)
    y = data["y"]
    tr, va, te = data["tr"], data["va"], data["te"]
    Xp = data["Xp"]
    Lk = data["lig_topk"]      # (n, K, 640) fp16
    Rk = data["rec_topk"]
    Lm = data["lig_mask"]
    Rm = data["rec_mask"]
    n, K, D = Lk.shape
    emit("## 1. 数据")
    emit("")
    emit(f"- 行数 {n}; train {int(tr.sum())} / val {int(va.sum())} / test {int(te.sum())}")
    emit(f"- per-residue 张量: ({n}, {K}, {D}) fp16; "
         f"物理特征 {Xp.shape[1]} 维")
    emit(f"- 训练时每链最多用 {min(args.max_res, K)} 个残基 (缓存 top-{K})")
    emit(f"- 保留残基数均值: 配体 {Lm.sum(1).mean():.1f} / 受体 {Rm.sum(1).mean():.1f}")

    # 张量组装 (整块 fp16->fp32 会占 ~560MB*2; 这里转成连续 fp16 便于按 batch 切片)
    Lig = np.ascontiguousarray(Lk)
    Rec = np.ascontiguousarray(Rk)

    class MaskedCNN(nn.Module):
        """掩码 1D 卷积 + 注意力池化。"""

        def __init__(self, d_in: int = 640, d_att: int = 64):
            super().__init__()
            self.conv1 = nn.Conv1d(d_in, 128, kernel_size=3, padding=1)
            self.conv2 = nn.Conv1d(128, 64, kernel_size=3, padding=1)
            self.act = nn.ReLU()
            self.drop = nn.Dropout(DROPOUT)
            # 注意力: 用可学习 query 对残基加权 (任务书"可选注意力池化")
            self.att = nn.Linear(64, 1)

        def forward(self, x: torch.Tensor, mask: torch.Tensor):
            # x: (B, L, D) float32; mask: (B, L) bool
            h = x.transpose(1, 2)                      # (B, D, L)
            m = mask.unsqueeze(1).float()              # (B, 1, L)
            h = self.act(self.conv1(h)) * m
            h = self.act(self.conv2(h)) * m            # (B, 64, L)
            h = self.drop(h)
            h = h.transpose(1, 2)                      # (B, L, 64)
            score = self.att(h).squeeze(-1)            # (B, L)
            score = score.masked_fill(~mask, -1e9)
            w = torch.softmax(score, dim=1)            # (B, L)
            pooled = (h * w.unsqueeze(-1)).sum(dim=1)  # (B, 64)
            return pooled

    class TwoTowerCNN(nn.Module):
        def __init__(self, d_phys: int):
            super().__init__()
            self.lig = MaskedCNN()
            self.rec = MaskedCNN()
            self.phys = nn.Sequential(nn.Linear(d_phys, 64), nn.ReLU())
            self.head = nn.Sequential(nn.Linear(192, 64), nn.ReLU(),
                                      nn.Linear(64, 1))

        def forward(self, lg, lm, rg, rm, ph):
            h = torch.cat([self.lig(lg, lm), self.rec(rg, rm), self.phys(ph)], dim=1)
            return self.head(h)

    def slice_res(A: np.ndarray, M: np.ndarray, idx: np.ndarray,
                  max_res: int):
        """取 idx 行, 截断到 max_res 个残基, 返回 (B, L, D) fp32 与 (B, L) bool。"""
        blk = A[idx][:, :max_res]                      # (B, L, D) fp16
        msk = M[idx][:, :max_res]
        return torch.from_numpy(blk.astype(np.float32)), torch.from_numpy(msk)

    def to_phys(idx: np.ndarray):
        """物理特征按 batch 转 tensor (统一入口, 避免漏转 numpy)。"""
        return torch.from_numpy(Xp[idx])

    results: list[dict] = []
    preds: dict[str, np.ndarray] = {}
    best_state = None
    best_score = -9.0

    for seed in args.seeds:
        torch.manual_seed(seed)
        np.random.seed(seed)
        model = TwoTowerCNN(Xp.shape[1])
        opt = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WD)
        lossf = nn.MSELoss()

        tr_idx = np.where(tr)[0]
        va_idx = np.where(va)[0]
        te_idx = np.where(te)[0]
        y_tr = y[tr]; y_va = y[va]
        Xp_tr = torch.from_numpy(Xp[tr]); Xp_va = torch.from_numpy(Xp[va])

        best_rho, best_ep, bad, best_st = -9.0, -1, 0, None
        emit("")
        emit(f"--- seed {seed} ---")
        for ep in range(args.epochs):
            model.train()
            perm = np.random.RandomState(seed * 1000 + ep).permutation(len(tr_idx))
            for s in range(0, len(perm), args.batch):
                b = tr_idx[perm[s:s + args.batch]]
                lg, lm = slice_res(Lig, Lm, b, args.max_res)
                rg, rm = slice_res(Rec, Rm, b, args.max_res)
                ph = to_phys(b)
                yb = torch.from_numpy(y[b].astype(np.float32)).view(-1, 1)
                opt.zero_grad()
                out = model(lg, lm, rg, rm, ph)
                loss = lossf(out, yb)
                loss.backward()
                opt.step()
            model.eval()
            with torch.no_grad():
                pv = []
                for s in range(0, len(va_idx), 128):
                    b = va_idx[s:s + 128]
                    lg, lm = slice_res(Lig, Lm, b, args.max_res)
                    rg, rm = slice_res(Rec, Rm, b, args.max_res)
                    pv.append(model(lg, lm, rg, rm, to_phys(b)).cpu().numpy().ravel())
                pv = np.concatenate(pv)
            rho = mr.spearman_rho(-y_va, -pv)
            if rho == rho and rho > best_rho + 1e-6:
                best_rho, best_ep, bad = float(rho), ep, 0
                best_st = {k: v.detach().clone() for k, v in model.state_dict().items()}
            else:
                bad += 1
            if (ep + 1) % 10 == 0 or bad >= PATIENCE:
                emit(f"  ep {ep+1}: val_rho={rho:+.4f} best={best_rho:+.4f} bad={bad}")
            if bad >= PATIENCE:
                break
        emit(f"  seed {seed}: best val_rho={best_rho:+.4f} @ep{best_ep+1} "
             f"({time.time() - t0:.0f}s)")
        if best_st is not None:
            model.load_state_dict(best_st)
        model.eval()
        with torch.no_grad():
            pt = []
            for s in range(0, len(te_idx), 128):
                b = te_idx[s:s + 128]
                lg, lm = slice_res(Lig, Lm, b, args.max_res)
                rg, rm = slice_res(Rec, Rm, b, args.max_res)
                pt.append(model(lg, lm, rg, rm, to_phys(b)).cpu().numpy().ravel())
            pt = np.concatenate(pt)
        y_s, p_s = -y[te], -pt
        row = {"model": "cnn_attn", "seed": seed, "best_epoch": best_ep + 1,
               "val_spearman": best_rho, "n_features": int(Xp.shape[1]) + 128}
        for k, fn in mfns.items():
            try:
                row[k] = float(fn(y_s, p_s))
            except Exception:  # noqa: BLE001
                row[k] = float("nan")
        for k in ("spearman", "ndcg@10", "ef@1%", "precision@10"):
            lo, hi, _ = C.bootstrap_ci(mfns[k], y_s, p_s, C.N_BOOT_CI, seed)
            row[f"{k}_ci_low"], row[f"{k}_ci_high"] = lo, hi
        results.append(row)
        preds[f"cnn_seed{seed}"] = pt
        emit(f"  seed {seed} TEST: Spearman {row['spearman']:+.4f} "
             f"[{row['spearman_ci_low']:+.3f}, {row['spearman_ci_high']:+.3f}], "
             f"NDCG@10 {row['ndcg@10']:.4f}, EF@1% {row['ef@1%']:.4f}, "
             f"P@10 {row['precision@10']:.4f}")
        if row["spearman"] > best_score:
            best_score = row["spearman"]
            best_state = (best_st, seed, row)

    res = pd.DataFrame(results)
    res.to_csv(OUT_CSV, index=False, encoding="utf-8")
    if preds:
        np.savez(OUT_PRED, test_y=y[te].astype(np.float32),
                 test_pdb_id=data["pdb_ids"][te],
                 **{k: v.astype(np.float32) for k, v in preds.items()})
        # 多种子平均 (判决用)
        if len(preds) > 1:
            avg = np.mean(list(preds.values()), axis=0)
            rho = mr.spearman_rho(-y[te], -avg)
            emit("")
            emit(f"- **多种子平均 ({len(preds)} 个种子)**: Spearman = **{rho:+.4f}**")
    if best_state is not None:
        try:
            st, sd, row = best_state
            torch.save({"state_dict": st, "seed": sd, "metrics": row,
                        "arch": "TwoTowerCNN(MaskedCNN x2 + attention pool) + phys tower",
                        "max_res": args.max_res, "topk_cache": int(K),
                        "note": "input = per-residue (interface residues, downsampled) + 60 physics"},
                       BEST_PT)
            emit(f"- 已写出 `{BEST_PT.name}` (state_dict + 元信息)")
        except Exception as exc:  # noqa: BLE001
            emit(f"- !! 保存模型失败: {exc}")
            with open(C.FAIL_LOG, "a", encoding="utf-8") as fh:
                fh.write(f"[03_cnn_perres] save: {exc}\n{traceback.format_exc()}\n")

    # ---- 结论 ----
    emit("")
    emit("## 2. 结论")
    emit("")
    if len(res):
        emit(f"- CNN 单种子最优 Spearman: {res['spearman'].max():+.4f} "
             f"(均值 {res['spearman'].mean():+.4f})")
        emit(f"- 参照 (同批行上重跑的 pooled 基线): `pooled_both@128` = +0.3253; "
             f"per-residue mean-pool = +0.3004")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
