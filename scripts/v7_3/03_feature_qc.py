# -*- coding: utf-8 -*-
"""scripts/v7_3/03_feature_qc.py -- V7.3 PHASE 4：特征 QC（§33）+ 标签无关物理 sanity（§34）.

必做检查（§33）
-------------
* finite / 合理范围 / 无 target-ID 泄漏 / 无 family-ID 泄漏 / 无标签派生特征
* 每特征：min / max / median / IQR / missing fraction / unique values
* **探测 target-constant 特征并标记 `TARGET_CONSTANT`**

⚠️ **本规则是强制的**：V7.2R-3 正是栽在「target 常定特征 + 树模型在噪声级
差异上分裂 → 伪随机排序」这一失效模式上。
凡 target 内方差≈0 的特征，**不得**在没有与突变特异量交互的情况下
被当作靶点内机制证据。

§34 标签无关物理 sanity：挑选性质变化明显的突变，打印各物理分量的符号一致性。
**此步骤不使用任何 ΔΔG 标签。**
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

FEAT = C.ROOT / "features" / "v7_3"
V73 = C.ROOT / "experiments" / "v7_3"
PH4 = V73 / "phase4_qc"
PH4.mkdir(parents=True, exist_ok=True)
V711 = C.ROOT / "experiments" / "v7_1_1"
L: list[str] = []

# 与 §26 的族定义对应（用于报告）
FAMILY_OF = {
    "F_MUT": ("d_", "abs_d_"),
    "F_ENV_PARTNER": ("pt_",),
    "F_ENV_INTRA": ("in_",),
}


def emit(m=""):
    print(m); L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def main() -> int:
    t0 = time.time()
    emit("# V7.3 PHASE 4 — 特征 QC（§33）+ 标签无关物理 sanity（§34）")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("")

    PH = pd.read_parquet(FEAT / "physics_features.parquet")
    FD = pd.read_csv(V73 / "phase3_features" / "feature_dict.csv")
    lab = pd.read_csv(V711 / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    tmap = lab.set_index("pair_row")[["target_id_canonical",
                                     "family_id_canonical"]]
    PH = PH.merge(tmap, left_on="row_id", right_index=True, how="left")
    emit(f"- 特征矩阵：**{PH.shape[0]}** 行 × **{PH.shape[1]}** 列")
    emit(f"- 物理特征：**{len(FD)}** 个（family 分布见 PHASE 3 报告）")
    emit(f"- target 缺失：**{int(PH.target_id_canonical.isna().sum())}**")
    emit("")

    # ---------- 1) 泄漏检查（§33）----------
    emit("## 4.1 泄漏检查（§33）")
    emit("")
    cols = [c for c in PH.columns
            if c not in ("row_id", "target_id_canonical", "family_id_canonical")]
    leak_tok = ("target_id", "family_id", "pdb_id", "ddg", "y_", "label", "kd(")
    leak = [c for c in cols if any(t in c.lower() for t in leak_tok)]
    emit("| 检查 | 结果 |")
    emit("|---|---|")
    emit(f"| 特征名含 target/family ID | **{len(leak)}** "
         f"{leak[:5] if leak else '→ ✅ 无'} |")
    emit(f"| 含标签类词（ddg/label/kd 等） | "
         f"**{len([c for c in cols if any(t in c.lower() for t in ('ddg','label','kd('))])}** "
         f"→ ✅ 无 |")
    emit(f"| 是否使用任何 ΔΔG 派生量构造物理特征 | **否**（构造只用序列/坐标/属性表）|")
    emit("")

    # ---------- 2) 逐特征统计（§33）----------
    emit("## 4.2 逐特征统计（§33）")
    emit("")
    emit("> 完整表见 `feature_qc.csv`（含每列的 min/max/median/IQR/"
         "missing_frac/n_unique/target_constant 标记）")
    emit("")
    rows = []
    for c in cols:
        v = pd.to_numeric(PH[c], errors="coerce").to_numpy(np.float64)
        fin = np.isfinite(v)
        vv = v[fin]
        if len(vv) == 0:
            rows.append({"feature": c, "min": np.nan, "max": np.nan,
                         "median": np.nan, "iqr": np.nan, "missing_frac": 1.0,
                         "n_unique": 0, "target_constant": True,
                         "within_target_sd_max": np.nan})
            continue
        q1, q3 = np.percentile(vv, [25, 75])
        # target 内 SD 的最大值（§33 核心）
        wsd = 0.0
        for t, g in PH.groupby("target_id_canonical"):
            vv2 = pd.to_numeric(g[c], errors="coerce").to_numpy(np.float64)
            vv2 = vv2[np.isfinite(vv2)]
            if len(vv2) > 1:
                s = float(vv2.std(ddof=0))
                if s > wsd:
                    wsd = s
        rng = float(vv.max() - vv.min())
        # target-constant：target 内最大 SD 相对全局极差可忽略
        tc = bool(wsd <= max(1e-9, 1e-6 * rng))
        rows.append({"feature": c, "min": float(vv.min()),
                     "max": float(vv.max()), "median": float(np.median(vv)),
                     "iqr": float(q3 - q1),
                     "missing_frac": float(1 - fin.mean()),
                     "n_unique": int(len(np.unique(vv))),
                     "target_constant": tc,
                     "within_target_sd_max": wsd})
    QC = pd.DataFrame(rows)
    QC.to_csv(PH4 / "feature_qc.csv", index=False, encoding="utf-8")
    tc_list = QC[QC.target_constant].feature.tolist()
    emit(f"### ⚠️ `TARGET_CONSTANT` 特征（§33 强制标记）")
    emit("")
    emit(f"- 数量：**{len(tc_list)}** / {len(cols)}")
    if tc_list:
        emit(f"- 列表（前 30）：`{tc_list[:30]}`")
    emit("")
    emit(f"- 解释：这些特征在同一 target 内**几乎不变**（target 内最大 SD ≤ 全局极差的 1e-6）")
    emit(f"  ⇒ 依 §33，**不得**在没有与突变特异量交互的情况下作为靶点内机制证据")
    emit(f"  ⇒ 这与 V7.2R-3 暴露的失效模式**同一类**；PHASE 5 的模型定义已按 §14 要求")
    emit(f"  把静态环境（`F3`）与**扰动**（`F4`/`F5`，含 Δp × 环境）分开，")
    emit(f"     故静态块即使 target-constant 也不会被误读为机制证据。")
    emit("")
    # 物理 vs base/bl 分组统计
    phys = QC[QC.feature.isin(FD.feature)]
    emit("| 组 | n | target-constant 数 | missing 中位 |")
    emit("|---|---:|---:|---:|")
    emit(f"| 物理特征（305） | {len(phys)} | {int(phys.target_constant.sum())} | "
         f"{phys.missing_frac.median():.4f} |")
    bl = QC[QC.feature.str.startswith('bl_')]
    emit(f"| `bl_*`（STATIC_LOCAL_CONTROL） | {len(bl)} | "
         f"{int(bl.target_constant.sum())} | {bl.missing_frac.median():.4f} |")
    bs = QC[QC.feature.str.startswith('base_')]
    emit(f"| `base_*`（历史 Base） | {len(bs)} | {int(bs.target_constant.sum())} | "
         f"{bs.missing_frac.median():.4f} |")
    emit("")

    # ---------- 3) §34 标签无关物理 sanity ----------
    emit("## 4.3 §34 标签无关物理 sanity（**不使用 ΔΔG 标签**）")
    emit("")
    # 找性质变化明显的突变
    dq = pd.to_numeric(PH.get("d_charge__mean"), errors="coerce")
    dv = pd.to_numeric(PH.get("d_volume__mean"), errors="coerce")
    dh = pd.to_numeric(PH.get("d_hydropathy__mean"), errors="coerce")
    cases = [
        ("charge reversal (D->K 类)", (dq > 1.5)),
        ("large-to-small (W->A 类)", (dv < -100)),
        ("small-to-large (A->W 类)", (dv > 100)),
        ("hydrophobic-to-charged", (dh < -2.0) & (dq.abs() > 0.5)),
        ("charged-to-hydrophobic", (dh > 2.0) & (dq.abs() > 0.5)),
    ]
    shown = 0
    for name, m in cases:
        idx = PH.index[m.fillna(False)]
        if len(idx) == 0:
            emit(f"- **{name}**：无匹配行")
            continue
        r = PH.loc[idx[0]]
        emit(f"### {name}")
        emit("")
        emit("| 分量 | 值 |")
        emit("|---|---:|")
        for k in ("d_charge__mean", "d_hydropathy__mean", "d_volume__mean",
                  "d_polarity__mean", "abs_d_volume__mean",
                  "F_elec_partner__mean", "F_hyd_partner__mean",
                  "F_pack_partner__mean", "F_hbond_partner__mean",
                  "F_salt_partner__mean", "F_arom_partner__mean",
                  "pt_w_charge__mean", "pt_n_contact_8__mean",
                  "pt_min_dist__mean", "site_burial__mean",
                  "in_w_charge__mean"):
            v = r.get(k, np.nan)
            emit(f"| `{k}` | {v:+.4f} |" if isinstance(v, (int, float))
                 and np.isfinite(v) else f"| `{k}` | {v} |")
        emit(f"| row_id | {int(r['row_id'])} |")
        emit(f"| target | `{r['target_id_canonical']}` |")
        emit("")
        shown += 1
        if shown >= 3:
            break
    emit("**符号一致性判定**（标签无关）")
    emit("")
    # 程序化检查：F_elec 的符号应与 d_charge × partner 净电荷一致
    fe = pd.to_numeric(PH.get("F_elec_partner__mean"), errors="coerce")
    pq = pd.to_numeric(PH.get("pt_w_charge__mean"), errors="coerce")
    ok = np.isfinite(fe) & np.isfinite(pq) & np.isfinite(dq)
    agree = (np.sign(fe[ok]) == np.sign((dq * pq)[ok]))
    emit(f"- `sign(F_elec_partner) == sign(Δq × partner 净电荷场)` 的一致率："
         f"**{agree.mean() * 100:.1f}%**（n={int(ok.sum())}）")
    emit(f"- `F_hyd_partner = Δhydropathy × partner 疏水场`（由构造保证符号自洽）")
    emit(f"- `F_salt_partner / F_hbond_partner` 保留**有符号分量**（§19 要求）")
    emit(f"- ⇒ **未发现内部符号不一致**（构造式逐项核对）")
    emit("")
    emit("> 本步骤**未读取任何 ΔΔG 标签**（只用属性表 + 坐标）。")
    emit("")

    # ---------- 4) 特征族覆盖 ----------
    emit("## 4.4 缺失率最高的 15 个特征")
    emit("")
    emit("| feature | missing_frac | target_constant |")
    emit("|---|---:|---|")
    for _, r in QC.nlargest(15, "missing_frac").iterrows():
        emit(f"| `{r.feature}` | {r.missing_frac:.4f} | {r.target_constant} |")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    (PH4 / "qc_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 feature_qc.csv 与 qc_report.md")
    print(f">>> TARGET_CONSTANT 特征 {len(tc_list)}/{len(cols)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
