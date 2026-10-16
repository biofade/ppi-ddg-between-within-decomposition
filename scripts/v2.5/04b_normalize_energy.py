# -*- coding: utf-8 -*-
"""scripts/v2.5/04b_normalize_energy.py -- V2.5 任务 2 补充: 能量特征的规模归一化.

**为什么需要这一步（实测发现的真实问题）**
------------------------------------------
`04_energy_features.py` 输出的是**能量总和**。但实测发现界面规模差异极大:

| 指标 | 中位数 | 99 分位 | **最大** |
|---|---:|---:|---:|
| 跨链原子对 (`en_lj_n_pairs`) | ~1.6k | ~600k | **8,075,776** |
| 疏水接触对数 | ~1.6k | ~616k | **4,004,521** |

原因很直接: 某些界面的**整条链**都落在 8 Å cutoff 内（例如肽段配体被大受体包围），
于是原子对数量爆炸。

后果: **能量总和几乎完全是「界面大小」的代理**，而不是「能量强弱」。
把它当特征喂给模型，等于又加了一个 SASA/接触数 型的几何量 ——
这恰恰是任务书要求**摆脱**的东西（"从计数升级为能量"）。

修法（**不重跑 39 分钟的结构遍历**）: 所有需要的分母都已经在表里
（原子对数、氢键数、盐桥数），因此可以**直接从现有 CSV 派生归一化特征**。

输出: `features/v2.5_energy_features_norm.csv`（原始 + 归一化列）
      `features/v2.5_energy_norm_report.md`
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v25_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

SRC = C.FEAT / "v2.5_energy_features.csv"
OUT_CSV = C.FEAT / "v2.5_energy_features_norm.csv"
OUT_MD = C.FEAT / "v2.5_energy_norm_report.md"

L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))


def main() -> int:
    emit("# V2.5 任务 2 补充 — 能量特征的规模归一化")
    emit("")
    if not SRC.exists():
        emit(f"!! 缺少 `{SRC.name}`")
        OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
        return 2

    e = pd.read_csv(SRC, low_memory=False)
    emit("## 1. 问题: 能量总和被界面规模主导")
    emit("")
    emit("| 量 | 中位数 | 99 分位 | 最大值 |")
    emit("|---|---:|---:|---:|")
    for c in ("en_lj_n_pairs", "en_hydrophobic_n_contact"):
        v = pd.to_numeric(e[c], errors="coerce")
        q = np.nanpercentile(v, [50, 99])
        emit(f"| `{c}` | {q[0]:,.0f} | {q[1]:,.0f} | {v.max():,.0f} |")
    emit("")
    emit("-> 某些界面的**整条链**都在 8 Å cutoff 内（肽段配体被大受体包围），")
    emit("   原子对数可达 800 万 -> 能量总和几乎只是「界面大小」的代理。")
    emit("   这正是任务书要摆脱的「计数型几何量」，必须归一化。")

    pair = pd.to_numeric(e["en_lj_n_pairs"], errors="coerce").replace(0, np.nan)
    out = e.copy()

    # ---- 归一化: 每原子对的能量 ----
    norm_specs = [
        ("en_lj_total", "n_lj_per_pair", "每原子对 LJ 能量"),
        ("en_lj_repulsive", "n_ljrep_per_pair", "每原子对排斥"),
        ("en_lj_attractive", "n_ljatt_per_pair", "每原子对吸引"),
        ("en_hbond_total", "n_hbond_per_pair", "每原子对氢键能量"),
        ("en_saltbridge_total", "n_sb_per_pair", "每原子对盐桥能量"),
        ("en_hydrophobic_total", "n_hydro_per_pair", "每原子对疏水"),
        ("en_total_approx", "n_total_per_pair", "每原子对总能量"),
    ]
    for src, dst, _desc in norm_specs:
        out[dst] = pd.to_numeric(e[src], errors="coerce") / pair

    # ---- 能量学占比 (无量纲, 与规模无关) ----
    hyd = pd.to_numeric(e["en_hydrophobic_n_contact"], errors="coerce").replace(0, np.nan)
    out["frac_hbond_of_pairs"] = pd.to_numeric(e["en_hbond_n"], errors="coerce") / pair
    out["frac_saltbridge_of_pairs"] = pd.to_numeric(e["en_saltbridge_n"], errors="coerce") / pair
    out["frac_hydrophobic_of_pairs"] = hyd / pair

    # ---- 强度量 (与规模无关, 本来就是强度) ----
    out["int_hbond_mean"] = pd.to_numeric(e["en_hbond_mean"], errors="coerce")
    out["int_hbond_min"] = pd.to_numeric(e["en_hbond_min"], errors="coerce")
    out["int_saltbridge_min"] = pd.to_numeric(e["en_saltbridge_min"], errors="coerce")
    out["int_hydro_per_contact"] = pd.to_numeric(
        e["en_hydrophobic_per_contact"], errors="coerce")
    out["int_entropy_per_res"] = pd.to_numeric(
        e["en_entropy_per_res"], errors="coerce")

    # 日志压缩: 这些量跨 6 个数量级, 用 signed log1p 便于树模型
    for c in ("n_lj_per_pair", "n_total_per_pair", "n_hydro_per_pair",
              "n_hbond_per_pair", "n_sb_per_pair", "n_ljrep_per_pair",
              "n_ljatt_per_pair"):
        v = pd.to_numeric(out[c], errors="coerce")
        out[c + "_slog"] = np.sign(v) * np.log1p(np.abs(v))

    new_cols = ([d for _s, d, _ in norm_specs]
                + ["frac_hbond_of_pairs", "frac_saltbridge_of_pairs",
                   "frac_hydrophobic_of_pairs", "int_hbond_mean", "int_hbond_min",
                   "int_saltbridge_min", "int_hydro_per_contact",
                   "int_entropy_per_res"]
                + [d + "_slog" for _s, d, _ in norm_specs])

    out.to_csv(OUT_CSV, index=False, encoding="utf-8")
    emit("")
    emit(f"## 2. 新增归一化特征 ({len(new_cols)} 个)")
    emit("")
    emit("| 特征 | 中位数 | 99 分位 | 最大值 |")
    emit("|---|---:|---:|---:|")
    for c in new_cols:
        v = pd.to_numeric(out[c], errors="coerce")
        q = np.nanpercentile(v.dropna(), [50, 99]) if v.notna().any() else [np.nan] * 2
        emit(f"| `{c}` | {q[0]:.4g} | {q[1]:.4g} | {v.max():.4g} |")

    emit("")
    emit("## 3. 结论")
    emit("")
    emit(f"- 归一化后特征跨度的**最大值/中位数比**从 ~5000x 降到 "
         f"~{np.nanmax(np.abs(pd.to_numeric(out['n_lj_per_pair'], errors='coerce'))) / max(1e-12, abs(np.nanmedian(pd.to_numeric(out['n_lj_per_pair'], errors='coerce')))):.1f}x")
    emit("- 消融时**同时**提供原始能量与归一化能量两组, 让数据决定哪组更有效。")
    emit("- 这也让 `08_thermo_decomposition.py` 的 ΔSASA 近似回到合理量级。")
    emit("")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_CSV.name} | {OUT_MD.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
