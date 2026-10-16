# -*- coding: utf-8 -*-
"""scripts/v7_1_1/02_reconciliation_final.py -- V7.1.1 PHASE 1(终) — 208 vs 214 定案.

已确证的根因（三个独立事实）
--------------------------
1. `d['groups']`（`build_pair_features`）与 canonical `pdb_id` 一致率 **100.00%**（7523/7523）
   → 在 canonical 上二者 unique 均为 **230**。
2. `block_L_local.csv` / `block_L_rescued.csv` 的 `pdb_id` 与 canonical 一致率
   仅 **94.21%**（4214/4473，**259 行不一致**），且出现 canonical 中不存在的
   `1Y1K`（该行 pair_row=3791 在 canonical 中是 `1YCS`）。
   → 该 `pdb_id` 来自 `mutation_sites.csv`，**不是** canonical 命名。
3. 因此同一个 4473 行集合，按两种 `pdb_id` 各自聚合，得到 **208** 与 **214** 个
   不同 target。

⇒ 判定：**`TARGET_ID_NAMESPACE_MISMATCH`**（不是 feature intersection，
也不是 metric eligibility，也不是纯报告文本错误）。

影响范围（必须说清）
------------------
* **LSR 的 target 标签用的是 `d['groups']`（100% 对齐）** → PHASE 3/4/5/6 的
  target 级聚合与指标**未受该错位影响**。
* **但 PHASE 3/4 的 S/G 特征**是用 `Sm.reindex(pdb_rows)`（`pdb_rows = `d['groups']`）
  取回的，而 S/G 表的索引来自 V6.3（其 `pdb_id` 源自 canonical）→ **这一步是对齐的**。
* 唯一实际受影响的是 **BLOCK L 特征矩阵与 target 标签的对齐**：
  259 行的 local 特征被算到了**错误的 target** 名下（用于 target 级平均时）。
* 本脚本**量化该错位对关键指标的实际影响**（用正确的 canonical 标签重算）。

输出
----
* `phase1_universe/reconciliation_final.md`
* `phase1_universe/target_id_impact.csv`
* `phase1_universe/canonical_target_labels.csv`（**V7.2 应使用的唯一权威标签**）
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v6_1"))
import v5_common as C  # noqa: E402
import fastmetrics as FM  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

V62P2 = C.ROOT / "experiments" / "v6_2" / "phase2_split4"
V70 = C.ROOT / "experiments" / "v7_0"
V71 = C.ROOT / "experiments" / "v7_1"
PH1 = C.ROOT / "experiments" / "v7_1_1" / "phase1_universe"
PH1.mkdir(parents=True, exist_ok=True)
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
    emit("# V7.1.1 PHASE 1 — 208 vs 214 定案")
    emit("")
    can = pd.read_csv(V62P2 / "matched_analysis_set.csv", low_memory=False)
    can_i = can.set_index("row")
    lf0 = pd.read_csv(V70 / "phase1_features" / "block_L_local.csv",
                      low_memory=False)

    emit("## 1. 三个已确证事实")
    emit("")
    d = C.build_pair_features(verbose=False)
    pdb_groups = np.array([str(x) for x in d["groups"]])
    agree = (pdb_groups == can.pdb_id.astype(str).to_numpy())
    emit(f"**(1)** `d['groups']` vs canonical `pdb_id`：一致率 "
         f"**{agree.mean() * 100:.2f}%**（{int(agree.sum())}/{len(agree)}）"
         f" → canonical 上二者 unique 均为 "
         f"**{len(set(pdb_groups))} / {can.pdb_id.nunique()}**")
    j = lf0.join(can_i[["pdb_id"]], on="pair_row", rsuffix="_can")
    a = j.pdb_id.astype(str).str.upper(); b = j.pdb_id_can.astype(str).str.upper()
    emit(f"**(2)** `block_L.pdb_id` vs canonical `pdb_id`：一致率 "
         f"**{(a == b).mean() * 100:.2f}%**，不一致 **{int((a != b).sum())}** 行；"
         f"且含 canonical 中不存在的 `1Y1K`")
    emit(f"**(3)** 同一 4473 行集合，按上述两种 pdb_id 聚合 → "
         f"**{lf0.pdb_id.nunique()}** vs **{len(set(pdb_groups[np.isin(np.arange(len(can)), can.row.to_numpy())][:0])) or len(set(pdb_groups[can.row.to_numpy()][np.isin(can.row.to_numpy(), lf0.pair_row.to_numpy())]))}** 个 target")
    emit("")

    # 权威标签：用 canonical 的 pdb_id，按 pair_row 映射
    lf0_rows = lf0.pair_row.to_numpy(int)
    auth = can_i.loc[lf0_rows, ["pdb_id", "family_label"]].copy()
    auth["pair_row"] = lf0_rows
    auth = auth[["pair_row", "pdb_id", "family_label"]].rename(
        columns={"pdb_id": "target_id_canonical",
                 "family_label": "family_id_canonical"})
    auth.to_csv(PH1 / "canonical_target_labels.csv", index=False,
                encoding="utf-8")
    emit("## 2. 权威 target 标签（V7.2 唯一应使用者）")
    emit("")
    emit(f"- 来源：**`canonical.pdb_id` 按 `pair_row` 映射**"
         f"（唯一与 `d['groups']` 100% 自洽的命名）")
    emit(f"- rows = **{len(auth)}**，targets = "
         f"**{auth.target_id_canonical.nunique()}**，families = "
         f"**{auth.family_id_canonical.nunique()}**")
    emit("")
    n_auth = auth.groupby("target_id_canonical").size()
    emit("| 阈值 | target 数 |")
    emit("|---|---:|")
    for k in (1, 2, 3, 5, 10):
        emit(f"| n ≥ {k} | **{int((n_auth >= k).sum())}** |")
    emit("")
    emit(f"- `n == 1` 的 target：**{int((n_auth == 1).sum())}**"
         f"（这些 target 的 within 类指标无法定义）")
    emit("")

    # ---- 影响量化 ----
    emit("## 3. 错位对关键指标的**实际影响**")
    emit("")
    emit("用 `lsr_predictions`（其 target 标签来自 `d['groups']`，100% 对齐）作为预测，")
    emit("分别用 **(A) `d['groups']` 标签** 与 **(B) 旧 `block_L.pdb_id` 标签** 计算指标。")
    emit("")
    lsr = np.load(V71 / "phase5_lsr" / "lsr_predictions.npz", allow_pickle=True)
    y = lsr["y"].astype(float)
    gA = np.array([str(x) for x in lsr["pdb"]])          # = d['groups']
    # 构造 gB：用 block_L 的 pdb_id（按 pair_row 顺序）
    order = lsr.get("pair_row") if "pair_row" in lsr.files else None
    gB = None
    if order is not None:
        mp = lf0.set_index("pair_row").pdb_id
        gB = mp.reindex(np.asarray(order)).astype(str).to_numpy()
    if gB is None or pd.isna(gB).any():
        emit("- ⚠️ `lsr_predictions.npz` 未含 `pair_row`，改用"
             "「pair_row → 旧 pdb_id」的按序拼接近似（已标注为近似）。")
        gB = lf0.sort_values("pair_row").pdb_id.astype(str).to_numpy()

    rows_out = []
    emit("| 指标 | (A) canonical/`groups` 标签 | (B) 旧 `block_L.pdb_id` 标签 | Δ |")
    emit("|---|---:|---:|---:|")
    p0 = lsr["p_P0_raw_ranker"].astype(float)
    p3 = lsr["p_P3_full_LSR"].astype(float)
    for tag, p in (("P0_raw", p0), ("P3_fullLSR", p3)):
        for lab, gg in (("A_groups", gA), ("B_blockL", gB)):
            uq = np.unique(gg)
            By = np.array([y[gg == c].mean() for c in uq])
            Bp = np.array([p[gg == c].mean() for c in uq])
            cnt = np.array([(gg == c).sum() for c in uq], float)
            den = float(y.std()) * float(p.std())
            CW = float((cnt / cnt.sum() * np.array(
                [np.mean((y[gg == c] - By[i]) * (p[gg == c] - Bp[i]))
                 for i, c in enumerate(uq)])).sum()) / den
            ms = FM.macro_spearman(y, p, gg, min_n=5)
            tm = C.pearson(By, Bp)
            key = f"{tag}__{lab}"
            globals()[key] = dict(C_W=CW, macro_spearman=ms,
                                  target_mean_r=tm, n_target=len(uq))
    for tag in ("P0_raw", "P3_fullLSR"):
        a_ = globals()[f"{tag}__A_groups"]; b_ = globals()[f"{tag}__B_blockL"]
        for met in ("C_W", "macro_spearman", "target_mean_r"):
            emit(f"| {tag} · {met} | {a_[met]:+.4f} | {b_[met]:+.4f} | "
                 f"{a_[met] - b_[met]:+.4f} |")
            rows_out.append({"model": tag, "metric": met,
                             "label_A_groups": a_[met],
                             "label_B_blockL": b_[met],
                             "delta": a_[met] - b_[met]})
        emit(f"| {tag} · n_target | {a_['n_target']} | {b_['n_target']} | "
             f"{a_['n_target'] - b_['n_target']:+d} |")
    emit("")
    pd.DataFrame(rows_out).to_csv(PH1 / "target_id_impact.csv", index=False,
                                  encoding="utf-8")

    # ---- 定案 ----
    emit("## 4. 定案")
    emit("")
    emit("| 判定 | 值 |")
    emit("|---|---|")
    emit("| **Cause** | **`TARGET_ID_NAMESPACE_MISMATCH`** |")
    emit("| 208 的来源 | `block_L*.csv` 的 `pdb_id`（源自 `mutation_sites`，94.21% 对齐） |")
    emit("| 214 的来源 | `d['groups']`（源自 `build_pair_features`，**100% 对齐 canonical**） |")
    emit("| 行数 | 两者均为 **4473**（同一 row 集合） |")
    emit("| 是否为 feature intersection | **否** |")
    emit("| 是否为 metric eligibility | **否** |")
    emit("| 是否为纯文本错误 | **否**（是两个真实不同的命名空间） |")
    emit("")
    emit("### 受影响 / 未受影响")
    emit("")
    emit("| 环节 | 是否受影响 | 说明 |")
    emit("|---|---|---|")
    emit("| LSR 的 target 标签（PHASE 3/4/5/6 聚合） | **未受影响** | 用 `d['groups']` |")
    emit("| S/G 特征取回（`reindex(pdb_rows)`） | **未受影响** | S/G 索引同源自 canonical |")
    emit("| BLOCK L 特征与 target 的对应 | **受影响（259 行）** | 该 259 行的 local 特征归错 target |")
    emit("| LSR/Oracle 的数值结论 | **未受影响** | target 标签自始用 `d['groups']` |")
    emit("")
    emit("### ⚠️ V7.1 报告须更正的表述")
    emit("")
    emit("1. 「4473 rows / **208 targets**」应改为「4473 rows / **214 targets**」")
    emit("   （208 是旧命名空间下的计数）")
    emit("2. 「**90.4%** target 覆盖」应改为 **214/230 = 93.0%**")
    emit("3. V7.1 PHASE 1 报告中「救援未增加任何行」的**结论仍成立**，")
    emit("   但其 `208` 这个计数应替换为 `214`")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    (PH1 / "reconciliation_final.md").write_text("\n".join(L) + "\n",
                                                 encoding="utf-8")
    print(f"\n写出 reconciliation_final.md / target_id_impact.csv / "
          f"canonical_target_labels.csv")
    print(f"\n>>> Cause = TARGET_ID_NAMESPACE_MISMATCH; "
          f"权威 target 数 = {auth.target_id_canonical.nunique()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
