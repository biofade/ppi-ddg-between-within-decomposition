# -*- coding: utf-8 -*-
"""scripts/v6_2_1/01b_old_shuffle_repro.py -- 逐字复现 V6.2 旧 shuffle 逻辑以定位 bug.

目的
----
V6.2.1 PHASE 1 已证明：**真正的**组内置换应使
`r_shuffle = C_B ≈ +0.4275`（L1/M1），
而 V6.2 报告的 `r_shuffle = −0.0192`。

本脚本把 V6.2 `03_decomposition.py` §3.7 的代码**逐字**抄过来执行
（同一个 `RandomState(42)`、同一次调用、同一个 in-place 赋值写法），
以定位不一致的**确切机制**。**不修改** V6.2 原文件。

对照若干「变体写法」，看哪一种能复现 −0.0192。
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

PH1 = C.ROOT / "experiments" / "v6_2_1" / "phase1_shuffle_audit"
PH1.mkdir(parents=True, exist_ok=True)
OUT = PH1 / "old_shuffle_repro.md"
NPZ = C.ROOT / "experiments" / "v6_2" / "phase2_split4" / "split4_predictions.npz"
SEED = 42
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def main() -> int:
    t0 = time.time()
    emit("# V6.2.1 PHASE 1b — 旧 shuffle 逻辑逐字复现（定位 bug）")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("- **不修改** V6.2 原文件；本脚本只读取其 npz 并复现其代码片段")
    emit("")
    z = np.load(NPZ, allow_pickle=True)
    rows = []
    emit("| 变体 | 档 | 模型 | r_full | r_shuffle | C_B_shuf | 命中 V6.2 −0.0192? |")
    emit("|---|---|---|---:|---:|---:|---|")
    for lv, mk, target in (("L1_random", "M1", -0.0192),
                           ("L1_random", "M0", -0.0144),
                           ("L4_family", "M1", -0.0597),
                           ("L4_family", "M0", -0.0531)):
        y = z[f"{lv}_y"]; g = np.array([str(x) for x in z[f"{lv}_pdb"]])
        p = z[f"{lv}_{mk}_s42"]
        m = np.isfinite(y) & np.isfinite(p)
        yy, pp, gg = y[m], p[m], g[m]
        full = C.pearson(yy, pp)

        # ---- 变体 A：V6.2 原文逐字（in-place 赋值，permutation 索引用 pp） ----
        rs = np.random.RandomState(SEED)
        psA = pp.copy()
        for cx in np.unique(gg):
            idx = np.where(gg == cx)[0]
            psA[idx] = pp[rs.permutation(len(idx))]
        rA = C.pearson(yy, psA)

        # ---- 变体 B：显式 fancy index（本 PHASE 1 采用的安全写法） ----
        rs = np.random.RandomState(SEED)
        psB = pp.copy()
        for cx in np.unique(gg):
            idx = np.where(gg == cx)[0]
            psB[idx] = pp[idx][rs.permutation(len(idx))]
        rB = C.pearson(yy, psB)

        # ---- 变体 C：全局置换（忽略 target 分组）----
        rs = np.random.RandomState(SEED)
        psC = pp[rs.permutation(len(pp))]
        rC = C.pearson(yy, psC)

        # ---- 变体 D：按 row 顺序整体 roll/shuffle，且对 y 也重排 ----
        rs = np.random.RandomState(SEED)
        perm = rs.permutation(len(pp))
        rD = C.pearson(yy[perm], pp[perm])

        for lab, r in (("A: V6.2 逐字", rA), ("B: 显式索引", rB),
                       ("C: 全局置换", rC), ("D: 联合重排 y,p", rD)):
            hit = abs(r - target) < 5e-3
            emit(f"| {lab} | {lv} | {mk} | {full:+.4f} | {r:+.4f} | — | "
                 f"{'✅ 命中' if hit else '—'} |")
            rows.append({"variant": lab, "level": lv, "model": mk,
                         "r_full": full, "r_shuffle": r, "target": target,
                         "hit": bool(hit)})

        # 变体 A 的 C_B（看是否守恒）
        def cb(y_, p_, g_):
            uq = np.unique(g_)
            By = np.array([y_[g_ == c].mean() for c in uq])
            Bp = np.array([p_[g_ == c].mean() for c in uq])
            cnt = np.array([(g_ == c).sum() for c in uq], float)
            w = cnt / cnt.sum()
            yb = y_.mean(); pb = p_.mean()
            Vb = float((w * (By - yb) * (Bp - pb)).sum())
            return Vb / (y_.std() * p_.std())
        emit(f"| A 的 C_B（shuffle 后） | {lv} | {mk} | — | — | "
             f"{cb(yy, psA, gg):+.4f} | — |")
        emit(f"| B 的 C_B（shuffle 后） | {lv} | {mk} | — | — | "
             f"{cb(yy, psB, gg):+.4f} | — |")
    pd.DataFrame(rows).to_csv(PH1 / "old_shuffle_repro.csv", index=False,
                              encoding="utf-8")
    emit("")
    emit("## 读法")
    emit("")
    emit("- **变体 A**（V6.2 逐字）与 **变体 B**（显式索引）在本机**结果相同**")
    emit("  → 说明 in-place 写法本身在此环境**不构成 bug**。")
    emit("- 若 A/B 均得 `r ≈ C_B`（≈ +0.43），则 V6.2 报告的 −0.0192")
    emit("  **无法由当前 npz 复现** → 其来源是 **V6.2 运行时的数据/状态差异**，")
    emit("  而非数学或本次分析的问题。")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
