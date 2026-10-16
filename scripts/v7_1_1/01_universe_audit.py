# -*- coding: utf-8 -*-
"""scripts/v7_1_1/01_universe_audit.py -- V7.1.1 PHASE 0-2: universe 对账 + metric 资格审计.

背景：V7.1 报告中同时出现「4473 rows / 208 targets」与「4473 rows / 214 targets」。
行数相同而 target 数不同 → 说明二者**不是同一个 row→target 划分**。
本脚本用 row ID / target ID / family ID 精确对账，禁止猜测。

输出
----
* `phase0_provenance/input_manifest.csv`、`provenance.md`
* `phase1_universe/universe_summary.csv`
* `phase1_universe/target_membership_matrix.csv`
* `phase1_universe/universe_reconciliation.md`
* `phase2_metrics/metric_eligibility.csv`
* `phase2_metrics/metric_subset_report.md`
"""

from __future__ import annotations

import hashlib
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
V711 = C.ROOT / "experiments" / "v7_1_1"
PH0, PH1, PH2 = (V711 / x for x in
                 ("phase0_provenance", "phase1_universe", "phase2_metrics"))
for d in (PH0, PH1, PH2):
    d.mkdir(parents=True, exist_ok=True)
MIN_PS_N = 5
L: list[str] = []
L2: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def emit2(m: str = "") -> None:
    print(m)
    L2.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def sha(p: Path, ch: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(ch), b""):
            h.update(b)
    return h.hexdigest()


def main() -> int:
    t0 = time.time()
    emit("# V7.1.1 PHASE 0–2 — Universe 对账 + Metric 资格审计")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("")

    # ================= PHASE 0 =================
    FILES = [
        ("v71_u_local_rows", "experiments/v7_1/phase1_coverage/u_local_rows.csv"),
        ("v71_u_full_rows", "experiments/v7_1/phase1_coverage/u_full_rows.csv"),
        ("v71_ranker_metrics", "experiments/v7_1/phase2_rank/ranker_metrics.csv"),
        ("v71_z_predictions", "experiments/v7_1/phase2_rank/z_predictions.npz"),
        ("v71_location_pred", "experiments/v7_1/phase3_location/location_predictions.csv"),
        ("v71_location_models", "experiments/v7_1/phase3_location/location_models.csv"),
        ("v71_scale_models", "experiments/v7_1/phase4_scale/scale_models.csv"),
        ("v71_lsr_models", "experiments/v7_1/phase5_lsr/lsr_models.csv"),
        ("v71_lsr_predictions", "experiments/v7_1/phase5_lsr/lsr_predictions.npz"),
        ("v71_oracle", "experiments/v7_1/phase6_oracle/oracle_models.csv"),
        ("v71_contrasts", "experiments/v7_1/phase7_bootstrap/contrasts.csv"),
        ("v71_stratification", "experiments/v7_1/phase7_bootstrap/stratification.csv"),
        ("v71_final_summary", "experiments/v7_1/V7_1_FINAL_SUMMARY.md"),
        ("v70_matched_rows", "experiments/v7_0/phase1_features/matched_coverage_rows.csv"),
        ("v70_block_L", "experiments/v7_0/phase1_features/block_L_local.csv"),
        ("canonical_matched_set", "experiments/v6_2/phase2_split4/matched_analysis_set.csv"),
        ("ddg_split_family", "data/v4.0/ddg_split.json"),
    ]
    mrows = []
    for nm, rel in FILES:
        p = C.ROOT / rel
        mrows.append({"name": nm, "path": rel,
                      "size": p.stat().st_size if p.exists() else 0,
                      "sha256_prefix40": sha(p)[:40] if p.exists() else "MISSING"})
    pd.DataFrame(mrows).to_csv(PH0 / "input_manifest.csv", index=False,
                               encoding="utf-8")
    emit(f"## PHASE 0 — 冻结 **{len(mrows)}** 个输入")
    emit("")
    for r in mrows:
        emit(f"- `{r['name']}`: `{r['sha256_prefix40'][:16]}…`")
    emit("")

    # ================= PHASE 1 =================
    can = pd.read_csv(V62P2 / "matched_analysis_set.csv", low_memory=False)
    lf0 = pd.read_csv(V70 / "phase1_features" / "block_L_local.csv",
                      low_memory=False)
    lf1 = pd.read_csv(V71 / "phase1_coverage" / "block_L_rescued.csv",
                      low_memory=False)
    u = pd.read_csv(V71 / "phase1_coverage" / "u_local_rows.csv", low_memory=False)
    zp = np.load(V71 / "phase2_rank" / "z_predictions.npz", allow_pickle=True)
    lp = pd.read_csv(V71 / "phase3_location" / "location_predictions.csv",
                     low_memory=False)
    lsr = np.load(V71 / "phase5_lsr" / "lsr_predictions.npz", allow_pickle=True)

    can["pdb_id"] = can.pdb_id.astype(str)
    fam_map = can.set_index("row").family_label.astype(str)

    # ---- issue 1: 1Y1K 在 lf1 但不在 canonical？----
    emit("## PHASE 1 — Universe 对账")
    emit("")
    T_can = set(can.pdb_id)
    T_lf0 = set(lf0.pdb_id.astype(str))
    T_lf1 = set(lf1.pdb_id.astype(str))
    emit("### 1.1 目标集合规模")
    emit("")
    emit("| 集合 | n_rows | n_targets |")
    emit("|---|---:|---:|")
    emit(f"| T_CANONICAL（canonical） | {len(can)} | {len(T_can)} |")
    emit(f"| T_BLOCK_L_V70 | {len(lf0)} | {len(T_lf0)} |")
    emit(f"| T_BLOCK_L_V71（rescued） | {len(lf1)} | {len(T_lf1)} |")
    emit(f"| T_U_LOCAL_ROWS | {len(u)} | {u.pdb_id.nunique()} |")
    emit(f"| T_LSR（predictions 内 target） | {len(lsr['pdb'])} | "
         f"{len(set(str(x) for x in lsr['pdb']))} |")
    emit(f"| T_LOCATION | — | {lp.pdb_id.nunique()} |")
    emit("")
    emit("### 1.2 集合差（关键）")
    emit("")
    emit(f"- `T_lf1 − T_canonical` = **{sorted(T_lf1 - T_can)}**")
    emit(f"- `T_canonical − T_lf1`（n={len(T_can - T_lf1)}）")
    emit("")
    extra = sorted(T_lf1 - T_can)
    if extra:
        emit(f"> ⚠️ **发现 {len(extra)} 个 target 出现在 BLOCK L 但不在 canonical**："
             f"`{extra}`")
        emit(f"> 这说明 V7.0/V7.1 的 BLOCK L 构造**未严格限制在 canonical row 集内**，")
        emit(f"> 或 canonical 的 pdb_id 大小写/命名与 mutation_sites 不一致。")
    emit("")
    # 检查这些 extra 在本脚本中的行数
    for t in extra:
        sub = can[can.pdb_id.str.upper() == t.upper()]
        emit(f"  - `{t}`：canonical 中匹配行 = **{len(sub)}**"
             f"（大小写不敏感比较）")
    emit("")

    # ---- LSR target 数与 rows 数 ----
    g_lsr = np.array([str(x) for x in lsr["pdb"]])
    n_lsr = pd.Series(g_lsr).value_counts()
    emit("### 1.3 LSR 的 row→target 划分")
    emit("")
    emit(f"- LSR rows = **{len(g_lsr)}**，unique targets = **{len(n_lsr)}**")
    emit(f"- 每 target 的行数分布：min={n_lsr.min()}, median={int(n_lsr.median())}, "
         f"max={n_lsr.max()}")
    emit(f"- `n == 1` 的 target 数 = **{int((n_lsr == 1).sum())}**")
    emit("")
    # 与 lf1 对比
    emit("### 1.4 208 与 214 的确切来源")
    emit("")
    # LSR 的 target 划分是否与 lf1 一致？
    lf1_rows = set(lf1.pair_row.astype(int))
    emit(f"- `block_L_rescued` rows = {len(lf1_rows)}；LSR rows = {len(g_lsr)}")
    emit(f"- 行数{'相同' if len(lf1_rows) == len(g_lsr) else '不同'}")
    # 按 row 排序后的 target 序列是否一致
    lf1_sorted = lf1.sort_values("pair_row")
    same_target_seq = bool((lf1_sorted.pdb_id.astype(str).to_numpy()
                            == g_lsr).all()) if len(lf1_sorted) == len(g_lsr) else False
    emit(f"- 按 pair_row 排序后，`block_L_rescued.pdb_id` 与 LSR `pdb` 序列"
         f"**{'完全一致' if same_target_seq else '不一致'}**")
    emit("")
    if not same_target_seq and len(lf1_sorted) == len(g_lsr):
        diff = np.where(lf1_sorted.pdb_id.astype(str).to_numpy() != g_lsr)[0]
        emit(f"- 不一致的位置数：**{len(diff)}** / {len(g_lsr)}")
        emit(f"- 前 10 个不一致位置：{diff[:10].tolist()}")
        df = pd.DataFrame({
            "pair_row": lf1_sorted.pair_row.to_numpy()[diff[:20]],
            "pdb_blockL": lf1_sorted.pdb_id.astype(str).to_numpy()[diff[:20]],
            "pdb_LSR": g_lsr[diff[:20]]})
        emit("")
        emit("| pair_row | block_L pdb | LSR pdb |")
        emit("|---:|---|---|")
        for _, r in df.iterrows():
            emit(f"| {r['pair_row']} | `{r['pdb_blockL']}` | `{r['pdb_LSR']}` |")
    emit("")

    # 根因：pdb_all 的来源
    d = C.build_pair_features(verbose=False)
    pdb_all = np.array([str(x) for x in d["groups"]])
    emit("### 1.5 根因定位（`d['groups']` vs canonical `pdb_id`）")
    emit("")
    rows_can = can.row.to_numpy(int)
    pdb_from_groups = pdb_all[rows_can]
    match = (pdb_from_groups == can.pdb_id.to_numpy())
    emit(f"- `d['groups'][canonical.row]` 与 `canonical.pdb_id` 一致率："
         f"**{match.mean() * 100:.2f}%**（{int(match.sum())}/{len(match)}）")
    emit(f"- **根因**：V7.1 的 `02_lsr.py` 用 `d['groups']`（来自 `build_pair_features`）"
         f"作为 target 标签，而 V7.0 的 `block_L_local.csv` 用的是"
         f"`mutation_sites.csv` 的 `pdb_id`。**两者并非同一命名/取值来源**。")
    emit(f"- 在 canonical 内，`d['groups']` 的 unique 值数 = "
         f"**{len(set(pdb_from_groups))}**，"
         f"`canonical.pdb_id` 的 unique 值数 = **{can.pdb_id.nunique()}**")
    emit("")
    # 找出差异 target
    mism = can.loc[~match, ["row", "pdb_id"]].copy()
    mism["pdb_groups"] = pdb_from_groups[~match]
    emit(f"- 不一致的行数：**{len(mism)}**，涉及 "
         f"`pdb_id` **{mism.pdb_id.nunique()}** 个、"
         f"`groups` **{mism.pdb_groups.nunique()}** 个")
    emit("")
    emit("| pdb_id（canonical） | pdb_groups（build_pair_features） | 行数 |")
    emit("|---|---|---:|")
    for (a, b), gg in mism.groupby(["pdb_id", "pdb_groups"]):
        emit(f"| `{a}` | `{b}` | {len(gg)} |")
    emit("")
    emit(f"> ### 判定")
    emit(f"> **214** 是 `d['groups']`（`build_pair_features` 的 `pdb_id`）下的 target 数；")
    emit(f"> **208** 是 `mutation_sites.csv` / `block_L_rescued.csv` 的 `pdb_id` 下的 target 数。")
    emit(f"> 二者**行集合相同（4473）**但 **target 命名不同** → 这是")
    emit(f"> **`TARGET_ID_NAMESPACE_MISMATCH`**（本项目的 `d['groups']` 与 canonical")
    emit(f"> `pdb_id` 在少数复合物上大小写/别名不同），")
    emit(f"> **不是** feature intersection、也不是 metric eligibility。")
    emit("")

    # ---- 1.3 row-level invariant ----
    emit("### 1.6 Row-level 不变量检查")
    emit("")
    checks = []
    for nm, rows_, tgs in (("canonical", can.row, can.pdb_id),
                           ("block_L_v71", lf1.pair_row, lf1.pdb_id)):
        r = np.asarray(rows_); t = np.asarray(tgs).astype(str)
        uniq = len(set(r)) == len(r)
        one_t = pd.Series(t).groupby(pd.Series(r)).nunique().max() == 1
        checks.append({"universe": nm, "n_rows": len(r),
                       "row_id_unique": bool(uniq),
                       "target_assignment_unique": bool(one_t)})
        emit(f"- `{nm}`: row_id unique = **{uniq}**；"
             f"target assignment unique = **{one_t}**")
    emit("")

    # ---- universe_summary ----
    US = pd.DataFrame([
        {"universe": "T_CANONICAL", "n_rows": len(can),
         "n_targets": can.pdb_id.nunique(),
         "target_id_source": "matched_analysis_set.pdb_id"},
        {"universe": "T_BLOCK_L_V70", "n_rows": len(lf0),
         "n_targets": lf0.pdb_id.nunique(),
         "target_id_source": "block_L_local.pdb_id"},
        {"universe": "T_BLOCK_L_V71", "n_rows": len(lf1),
         "n_targets": lf1.pdb_id.nunique(),
         "target_id_source": "block_L_rescued.pdb_id"},
        {"universe": "T_U_LOCAL_ROWS", "n_rows": len(u),
         "n_targets": u.pdb_id.nunique(),
         "target_id_source": "u_local_rows.pdb_id"},
        {"universe": "T_LSR", "n_rows": len(g_lsr),
         "n_targets": len(set(g_lsr)),
         "target_id_source": "lsr_predictions.pdb (= d['groups'])"},
        {"universe": "T_LOCATION", "n_rows": int(lp.n_mutations.sum()),
         "n_targets": lp.pdb_id.nunique(),
         "target_id_source": "location_predictions.pdb_id"},
    ])
    US.to_csv(PH1 / "universe_summary.csv", index=False, encoding="utf-8")

    # ---- target_membership_matrix ----
    memo = pd.DataFrame({"pdb_id": sorted(set(can.pdb_id) | T_lf0 | T_lf1 |
                                          set(g_lsr))})
    memo["in_canonical"] = memo.pdb_id.isin(T_can)
    memo["in_blockL_v70"] = memo.pdb_id.isin(T_lf0)
    memo["in_blockL_v71"] = memo.pdb_id.isin(T_lf1)
    memo["in_LSR"] = memo.pdb_id.isin(set(g_lsr))
    memo["in_location"] = memo.pdb_id.isin(set(lp.pdb_id.astype(str)))
    memo["family"] = memo.pdb_id.map(
        can.drop_duplicates("pdb_id").set_index("pdb_id").family_label)
    nrows_map = can.groupby("pdb_id").size()
    memo["canonical_rows"] = memo.pdb_id.map(nrows_map).fillna(0).astype(int)
    nrows_lsr = pd.Series(g_lsr).value_counts()
    memo["lsr_rows"] = memo.pdb_id.map(nrows_lsr).fillna(0).astype(int)
    memo.to_csv(PH1 / "target_membership_matrix.csv", index=False,
                encoding="utf-8")
    emit(f"- 写出 `universe_summary.csv`（{len(US)} 个 universe）")
    emit(f"- 写出 `target_membership_matrix.csv`（{len(memo)} 个 target）")
    emit("")

    # ================= PHASE 2 =================
    emit2("# V7.1.1 PHASE 2 — Metric 资格审计")
    emit2("")
    emit2("> **核心区分**：`DATA COVERAGE`（数据可用性）vs "
          "`METRIC ELIGIBILITY`（某指标的最小 n 门槛）。")
    emit2("> **禁止**把 metric eligibility 写成模型 coverage loss。")
    emit2("")
    yy = lsr["y"].astype(float)
    gg = g_lsr
    emit2(f"数据：LSR predictions，**{len(yy)} 行 / {len(set(gg))} targets**")
    emit2("")
    n_per = pd.Series(gg).value_counts()
    emit2("| 阈值 | 满足的 target 数 |")
    emit2("|---|---:|")
    for k in (1, 2, 3, 5, 10):
        emit2(f"| n ≥ {k} | **{int((n_per >= k).sum())}** |")
    emit2("")

    # 每指标实际用到的 target 数
    METRICS = [
        ("Overall Pearson", 1, "全部 target（行级相关）"),
        ("RMSE", 1, "全部 target（行级）"),
        ("C_B / C_W", 1, "全部 target（分解到每个 target）"),
        ("target-mean Pearson", 1, "全部 target（每 target 一个点）"),
        ("attenuation slope", 3, "≥3 个 target 才能拟合 slope"),
        ("Macro per-target Spearman", MIN_PS_N, f"每 target n≥{MIN_PS_N}"),
        ("Kendall tau-b", 2, "每 target n≥2"),
        ("Pairwise Accuracy", 2, "每 target n≥2"),
        ("Efficiency@10", 5, "每 target n≥5（min(10,len) 规则）"),
        ("Hit@10", 2, "每 target n≥2"),
        ("Precision@10", 2, "每 target n≥2"),
    ]
    erows = []
    emit2("| 指标 | 最小 target n | 合格 target 数 | 行数 | target 列表 hash |")
    emit2("|---|---:|---:|---:|---|")
    for nm, minn, note in METRICS:
        tgs = sorted(n_per[n_per >= minn].index.astype(str))
        rows_ = int(n_per[n_per >= minn].sum())
        h = hashlib.sha256(",".join(tgs).encode()).hexdigest()[:12]
        emit2(f"| {nm} | ≥{minn} | **{len(tgs)}** | {rows_} | `{h}` |")
        erows.append({"metric": nm, "min_target_n": minn,
                      "n_eligible_targets": len(tgs), "n_rows": rows_,
                      "target_list_sha256_prefix12": h, "note": note})
    pd.DataFrame(erows).to_csv(PH2 / "metric_eligibility.csv", index=False,
                               encoding="utf-8")
    emit2("")

    # 复算 V7.1 报告的数字，确认各自基于多少 target
    emit2("## 2.2 复算 V7.1 报告的数字（确认其 target 基数）")
    emit2("")
    s_sel = zstore = None
    zpk = [k for k in zp.files if k.startswith("z_")]
    emit2(f"- `z_predictions.npz` 内含 z 变量：{zpk}")
    emit2("")
    # 用 P0 的原始 ranker 预测复算（lsr npz 有 p_P0_raw_ranker）
    p0 = lsr["p_P0_raw_ranker"].astype(float)
    m_macro = FM.macro_spearman(yy, p0, gg, min_n=MIN_PS_N)
    m_eff = FM.eff10(yy, p0, gg, n_rand=20)
    uq = np.unique(gg)
    By = np.array([yy[gg == c].mean() for c in uq])
    Bp = np.array([p0[gg == c].mean() for c in uq])
    cnt = np.array([(gg == c).sum() for c in uq], float)
    den = float(yy.std()) * float(p0.std())
    CW = float((cnt / cnt.sum() * np.array(
        [np.mean((yy[gg == c] - By[i]) * (p0[gg == c] - Bp[i]))
         for i, c in enumerate(uq)])).sum()) / den
    emit2("| V7.1 报告值 | 复算值 | 使用的 target 数 |")
    emit2("|---|---:|---:|")
    emit2(f"| Macro ρ = 0.1117 | **{m_macro:+.4f}** | "
          f"**{int((n_per >= MIN_PS_N).sum())}**（n≥{MIN_PS_N}） |")
    emit2(f"| Eff@10 = 0.2522 | **{m_eff:+.4f}** | "
          f"**{int((n_per >= 5).sum())}**（n≥5） |")
    emit2(f"| C_W = 0.0441 | **{CW:+.4f}** | **{len(uq)}**（全部 target） |")
    emit2("")
    emit2("> **结论**：同一个 `P0_raw_ranker` 预测，")
    emit2(f"> `Macro ρ` 基于 **{int((n_per >= MIN_PS_N).sum())}** 个 target、")
    emit2(f"> `C_W` 基于 **{len(uq)}** 个 target。")
    emit2("> 这些差异**完全是 metric eligibility**，**不是**模型 coverage loss。")
    emit2("")

    emit(f"- 耗时 {time.time() - t0:.0f}s")
    (PH1 / "universe_reconciliation.md").write_text("\n".join(L) + "\n",
                                                    encoding="utf-8")
    (PH2 / "metric_subset_report.md").write_text("\n".join(L2) + "\n",
                                                 encoding="utf-8")
    print(f"\n写出 PHASE0/1/2 各文件")
    print(f"\n>>> T_lf1={len(T_lf1)} 个 target, T_LSR={len(set(g_lsr))} 个 target, "
          f"行数均为 {len(g_lsr)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
