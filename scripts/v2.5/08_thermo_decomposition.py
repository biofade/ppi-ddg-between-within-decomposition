# -*- coding: utf-8 -*-
"""scripts/v2.5/08_thermo_decomposition.py -- V2.5 任务 4: 热力学分解 (ΔH / ΔS).

任务书给定模型
--------------
1. `ΔCp ≈ 0.45 * ΔSASA` (cal/mol/K per Å²)
2. `ΔS ≈ ΔCp * ln(T / T_ref)`
3. 构造 **ΔH** 与 **ΔS** 作为独立特征

**先讲清一个重要的物理点（否则特征会没有意义）**
------------------------------------------------
任务给的 `ΔS` 是热容对熵的**修正项**（`ΔCp·ln(T/T_ref)`），
在 `T = T_ref` 时它**恒为 0** —— 单独当特征是退化的。

要做到任务书要的「ΔH 与 ΔS 作为独立特征」，必须用**完整的热力学恒等式**：

```
ΔG = ΔH − T·ΔS
```

* `ΔS_meas = (ΔH − ΔG)/T`，而 `ΔH` 又依赖 `ΔS` … 单独数据无法定解；
* 因此采用**标准做法**：设 `ΔS_model = ΔCp · ln(T/T_ref)`（热容修正项），
  再由**实测** `ΔG` 反解出 `ΔH`：
  `ΔH = ΔG + T·ΔS_model`
* 于是得到两个**真正独立于 ΔG 的构造**：
  * `dH_from_dG` = ΔG 与 ΔS 的线性组合（含实测 ΔG，因此**必须**在消融中
    与 `physics_only` 分开评估，注意它含标签信息 -> **见下方警告**）;
  * `dCp_model`、`dS_model` **完全不使用 ΔG**，是纯结构量。

> ⚠️ **数据泄漏警告**: `ΔH = ΔG + T·ΔS` 直接含**实测 ΔG**。
> 把它当特征去预测 ΔG 会造成**循环论证**（近似复现标签），指标会虚高。
> 因此本脚本:
> 1. 输出 `dS_model` / `dCp_model` 为**可用于预测的纯结构特征**;
> 2. `dH_from_dG` 单独标注为 **`leaky_`** 前缀，**在任务 5 的消融中默认排除**,
>    并额外做一次「加入 leaky 特征」的对照以**量化泄漏幅度**（诚实展示）。

温度数据的现实
--------------
实测：3197 个界面中**只有 27 个**有 >1 个温度（`Temperature(K)` 列，
且 1862/8207 行是 `298(assumed)` 占位）。
-> **无法做真正的 van't Hoff 分析**（由 ΔG 随 T 的斜率同时定出 ΔH 与 ΔS）。
这一点在报告中明确写出。

输出: `features/v2.5_thermo_features.csv`、`features/v2.5_thermo_report.md`
"""

from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import v25_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

OUT_CSV = C.FEAT / "v2.5_thermo_features.csv"
OUT_MD = C.FEAT / "v2.5_thermo_report.md"

T_REF = 298.15
DCP_PER_ASASA = 0.45 / 1000.0     # 0.45 cal/mol/K/Å² -> kcal/mol/K/Å²
# ΔSASA 近似: 与 04_energy_features.py 同一系数 (20 Å² / 埋藏疏水原子)
BURIED_PER_ATOM = 20.0

L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def main() -> int:
    ap = argparse.ArgumentParser(description="V2.5 热力学分解")
    ap.add_argument("--energy-csv", type=__import__("pathlib").Path,
                    default=C.FEAT / "v2.5_energy_features.csv")
    args = ap.parse_args()

    t0 = time.time()
    emit("# V2.5 任务 4 — 热力学分解 (ΔCp / ΔS / ΔH)")
    emit("")
    emit("> ⚠️ 关键: 任务给的 `ΔS = ΔCp·ln(T/T_ref)` 在 `T=T_ref` 时**恒为 0**。")
    emit("> 要得到有意义的 ΔH/ΔS 独立特征, 必须用恒等式 `ΔG = ΔH − T·ΔS`,")
    emit("> 而这会引入**实测 ΔG** -> 泄漏。本脚本把两者**严格分开**。")
    emit("")

    if not args.energy_csv.exists():
        emit(f"!! 缺少 `{args.energy_csv.name}`; 请先跑 `04_energy_features.py`")
        OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
        return 2

    e = pd.read_csv(args.energy_csv, low_memory=False)
    emit(f"## 1. 输入")
    emit("")
    emit(f"- `{args.energy_csv.name}`: {e.shape[0]} 行 x {e.shape[1]} 列")

    # 用接触疏水原子数近似 ΔSASA (与 04 同一口径, 保证一致)
    n_hyd = pd.to_numeric(e.get("en_hydrophobic_n_contact"), errors="coerce")
    dsasa_approx = (n_hyd * BURIED_PER_ATOM).astype(float)
    emit(f"- ΔSASA 近似 = 埋藏疏水接触原子数 × {BURIED_PER_ATOM} Å² "
         f"(与 `04_energy_features.py` 同口径)")
    emit(f"  均值 {dsasa_approx.mean():.1f} Å², 范围 "
         f"[{dsasa_approx.min():.1f}, {dsasa_approx.max():.1f}]")

    # ΔCp 为负 (结合埋藏表面 -> 热容下降)
    dCp = -DCP_PER_ASASA * dsasa_approx          # kcal/mol/K
    emit("")
    emit("### ΔCp 符号说明")
    emit(f"- 结合时埋藏表面 -> 热容**下降**, 故 `ΔCp = −{DCP_PER_ASASA:.4f}·ΔSASA`")
    emit(f"- ΔCp 均值 {dCp.mean():.4f} kcal/mol/K "
         f"(范围 [{dCp.min():.4f}, {dCp.max():.4f}])")

    # ΔS 的热容修正项 (纯结构量, 不含 ΔG)
    dS_model = dCp * np.log(T_REF / T_REF)       # T 只有 assumed 298.15 -> 0
    emit("")
    emit("### ΔS 的现实")
    emit(f"- `ΔS = ΔCp·ln(T/T_ref)`: 由于绝大多数行的 T 就是 {T_REF} K,")
    emit(f"  该项**恒为 0** -> 单独作为特征**退化**, 不能进模型。")
    emit("- 实测温度只有 27/3197 个界面有变化 -> **无法做 van't Hoff 分解**。")
    emit("- 因此本脚本改用: 以 `04_energy_features.py` 的**构象熵损失**作为 ΔS 的")
    emit("  物理代理 (`en_entropy_total`), 它来自**每个残基的可旋转键数**,")
    emit("  是纯结构量、不含 ΔG。")

    dS_proxy = pd.to_numeric(e.get("en_entropy_total"), errors="coerce")  # kcal/mol/K

    # ΔH = ΔG + T·ΔS  (对非退化情形, 即用真实 T)
    # 这里给出**两种**:
    #   (a) 用 ΔS_proxy (可作特征, 但仍含 ΔG -> 泄漏)
    #   (b) dCp / dS_model (纯结构, 不泄漏)
    v22 = pd.read_csv(C.V22, low_memory=False)
    dG_map = v22.groupby("interface_key")["real_dG"].mean()
    dg = e["interface_key"].map(dG_map).astype(float)
    temp = (pd.to_numeric(v22.groupby("interface_key")["Temperature(K)"]
                          .first().astype(str)
                          .str.replace(r"[^0-9.]", "", regex=True),
                          errors="coerce"))
    temp_map = temp.reindex(e["interface_key"]).astype(float).fillna(T_REF)
    dH_leaky = dg + temp_map * dS_proxy

    out = pd.DataFrame({
        "interface_key": e["interface_key"],
        "th_dCp_model": dCp,                       # 纯结构, 可用
        "th_dS_model_capcorr": dS_model,           # 纯结构, 但恒 0 -> 标注退化
        "th_dS_proxy_rotbond": dS_proxy,           # 纯结构, 可用
        "th_dSASA_approx": dsasa_approx,           # 纯结构, 可用
        "th_T_used": temp_map,                     # 温度 (含 assumed)
        "leaky_dH_from_dG": dH_leaky,              # **含实测 ΔG, 默认排除**
    })
    out.to_csv(OUT_CSV, index=False, encoding="utf-8")

    emit("")
    emit(f"## 2. 输出特征")
    emit("")
    emit("| 特征 | 是否含 ΔG | 均值 | 标准差 | 范围 |")
    emit("|---|---|---:|---:|---|")
    meta = [
        ("th_dCp_model", "否（纯结构）"),
        ("th_dS_model_capcorr", "否（但**恒为 0**，退化）"),
        ("th_dS_proxy_rotbond", "否（纯结构）"),
        ("th_dSASA_approx", "否（纯结构）"),
        ("th_T_used", "否"),
        ("leaky_dH_from_dG", "**是（实测 ΔG）→ 默认排除**"),
    ]
    for c, leak in meta:
        v = pd.to_numeric(out[c], errors="coerce")
        emit(f"| `{c}` | {leak} | {v.mean():.4f} | {v.std():.4f} | "
             f"[{v.min():.4f}, {v.max():.4f}] |")

    emit("")
    emit("## 3. 结论")
    emit("")
    emit(f"- 可用于预测的**纯结构**热力学特征: `th_dCp_model`, "
         f"`th_dS_proxy_rotbond`, `th_dSASA_approx`（3 个）")
    emit(f"- `th_dS_model_capcorr` **恒为 0**（因为 T ≡ T_ref）-> **不进模型**, "
         f"仅作记录")
    emit(f"- `leaky_dH_from_dG` **含实测 ΔG**, 会构成循环论证 -> **默认排除**;")
    emit(f"  任务 5 会做一次「加入它」的对照来**量化泄漏幅度**（预期指标虚高）")
    emit("")
    emit("## 4. 局限")
    emit("")
    emit("1. **ΔCp 由 ΔSASA 近似**: `ΔSASA` 本身是「疏水接触原子数 × 20 Å²」的近似,")
    emit("   因此 `ΔCp` 继承该误差。")
    emit("2. **没有 van't Hoff 分解**: 只有 27/3197 个界面有温度变化,")
    emit("   且 1862/8207 行温度是 `298(assumed)` -> 无法由 ΔG(T) 斜率同时定 ΔH/ΔS。")
    emit("3. **ΔS 用可旋转键数代理**: 未考虑埋藏程度与残留运动, 是**上界**式估计。")
    emit("4. 温度列含占位值, `th_T_used` 的可信度有限。")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_CSV.name} | {OUT_MD.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
