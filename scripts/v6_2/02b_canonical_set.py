# -*- coding: utf-8 -*-
"""scripts/v6_2/02b_canonical_set.py -- V6.2 PHASE 2.1: Canonical Row Universe.

目标
----
构建 `matched_analysis_set`：**四档（L1–L4）划分都能用**的行。
四档主比较全部基于该集合；另报 full-coverage sensitivity。

要求列（V6.2 §2.1）
-------------------
`pdb_id` / `mutation` / `ddG` / `family_label` / `seq40_cluster` / 特征

缺失任一 → 剔除并记录原因。

输出
----
* `experiments/v6_2/phase2_split4/matched_analysis_set.csv`
* `experiments/v6_2/phase2_split4/excluded_rows.csv`
* `experiments/v6_2/phase2_split4/coverage_report.md`
"""

from __future__ import annotations

import json
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

PH2 = C.ROOT / "experiments" / "v6_2" / "phase2_split4"
PH2.mkdir(parents=True, exist_ok=True)
OUT_SET = PH2 / "matched_analysis_set.csv"
OUT_EXC = PH2 / "excluded_rows.csv"
OUT_MD = PH2 / "coverage_report.md"
OUT_FAM = PH2 / "family_labels.csv"
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
    emit("# V6.2 PHASE 2.1 — Canonical Row Universe（matched_analysis_set）")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("")

    # ---------- 加载 V4/V5 已审计的特征装配（PHASE 2.0 已确认复用） ----------
    d = C.build_pair_features(verbose=True)
    df = d["df"].copy()
    y = d["y"]
    groups = np.array([str(x) for x in d["groups"]])
    n0 = len(y)
    emit("")
    emit(f"- `build_pair_features` 产出: **{n0}** 行（真 ΔΔG forward 配对）")
    emit(f"- 唯一复合物 (pdb_id): **{len(set(groups))}**")
    emit("")

    # ---------- 构建必需列 ----------
    M = pd.DataFrame({
        "row": np.arange(n0),
        "pdb_id": groups,
        "mutation": df["mutation"].astype(str).to_numpy()
        if "mutation" in df.columns else "",
        "ddG": y,
    })
    if "n_mut" in df.columns:
        M["n_mut"] = pd.to_numeric(df["n_mut"], errors="coerce").to_numpy()
    if "ligand_seq_changed" in df.columns:
        M["ligand_seq_changed"] = df["ligand_seq_changed"].astype(bool).to_numpy()
    if "receptor_seq_changed" in df.columns:
        M["receptor_seq_changed"] = df["receptor_seq_changed"].astype(bool).to_numpy()

    # ---------- state_id（V6.2 §2.4）：pdb|变异序列 哈希 ----------
    # 用 v2.2 的序列列构造（WT/MUT 各一侧）
    v22 = pd.read_csv(C.V22, low_memory=False)
    wt_i = d["fwd_v22_wt"]
    mut_i = d["fwd_v22_mut"]
    seq_cols = [c for c in ("Ligand Sequences", "Receptor Sequences")
                if c in v22.columns]

    def seq_of(idx):
        parts = []
        for c in seq_cols:
            parts.append(str(v22[c].iloc[int(idx)])[:400])
        return "|".join(parts)

    emit("## state_id 构造（V6.2 §2.4 要求的 state identity）")
    emit("")
    emit(f"- 序列列: {seq_cols}")
    wt_seq = [seq_of(i) for i in wt_i]
    mut_seq = [seq_of(i) for i in mut_i]
    import hashlib
    state_ids = [hashlib.sha256(
        f"{g}|{ws}|{ms}".encode("utf-8", errors="replace")).hexdigest()[:16]
        for g, ws, ms in zip(groups, wt_seq, mut_seq)]
    M["state_id"] = state_ids
    emit(f"- 唯一 `state_id`: **{len(set(state_ids))}** / {n0} 行")
    emit(f"  → 说明每个 `(pdb, WT序列, MUT序列)` 组合基本唯一（符合预期：每对独立）")
    emit("")

    # ---------- family_label（union-find，与 V4 一致：共享任一链簇即同族） ----------
    emit("## family_label（≤40% 序列同一性 + union-find）")
    emit("")
    fam_map = d.get("family_of_complex", {}) or {}
    emit(f"- V4 split 提供的 `family_of_complex`: **{len(fam_map)}** 个映射")
    if fam_map:
        M["family_label"] = [str(fam_map.get(g, f"UNK_{g}")) for g in groups]
    else:
        M["family_label"] = groups
    emit(f"- 唯一 family: **{M['family_label'].nunique()}**")
    emit("")
    # seq40_cluster：用 family 的父级（同族即同 cluster 的上界）
    M["seq40_cluster"] = M["family_label"]

    # ---------- forward/reverse twin 标记 ----------
    # 同一 (pdb, mutation) 若在 ddg_pairs 中同时有 forward/reverse，则互为 twin。
    # 本表已是 forward-only，故 twin 只可能来自「同一 pdb 的同一 mutation 重复记录」。
    key = M["pdb_id"].astype(str) + "|" + M["mutation"].astype(str)
    dup_mask = key.duplicated(keep=False)
    M["obs_key"] = key
    emit("## 重复观测 / twin")
    emit("")
    emit(f"- `(pdb, mutation)` 重复的行: **{int(dup_mask.sum())}**"
         f"（{key[dup_mask].nunique() if dup_mask.any() else 0} 个重复键）")
    emit(f"- 本表为 forward-only，故不含 forward/reverse twin；"
         f"重复键来自同一 (pdb, mutation) 的多条实验记录")
    emit("")

    # ---------- 过滤：四档都能用 ----------
    emit("## 过滤（四档都可用）")
    emit("")
    reasons = pd.Series("", index=M.index)
    ok = pd.Series(True, index=M.index)

    def drop(mask, why):
        nonlocal ok
        sel = ok & mask
        reasons.loc[sel] = why
        ok &= ~mask
        emit(f"- 剔除 `{why}`: **{int(sel.sum())}** 行")

    drop(~np.isfinite(M["ddG"].to_numpy(float)), "NaN_ddG")
    drop(M["mutation"].astype(str).str.strip() == "", "missing_mutation")
    drop(M["pdb_id"].astype(str).str.strip() == "", "missing_pdb_id")
    drop(M["family_label"].astype(str).str.startswith("UNK_"),
         "missing_family_label")
    # 特征完整性：phys/Z 块不得全为 NaN
    blk = d["X_blocks"]
    feat_ok = np.isfinite(blk["phys_diff"]).any(axis=1) & \
        np.isfinite(blk["Z_diff"]).any(axis=1)
    drop(~feat_ok, "missing_features")
    emit("")

    matched = M[ok].copy()
    excluded = M[~ok].copy()
    excluded["exclude_reason"] = reasons[~ok].to_numpy()
    emit(f"- **matched_analysis_set: {len(matched)} 行**（原始 {n0}）")
    emit(f"- excluded: **{len(excluded)} 行**")
    emit("")

    # 与 state_id 冲突检查（同 state_id 不得落在不同 fold）
    n_state_multi = int((matched.groupby("state_id").size() > 1).sum())
    emit(f"- 在 matched set 内，同一 `state_id` 出现多行的组数: **{n_state_multi}**"
         f"（L1 需保证它们同 fold）")
    emit("")

    # ---------- 写出 ----------
    keep = ["row", "pdb_id", "mutation", "ddG", "family_label", "seq40_cluster",
            "state_id", "obs_key"]
    for c in ("n_mut", "ligand_seq_changed", "receptor_seq_changed"):
        if c in matched.columns:
            keep.append(c)
    matched[keep].to_csv(OUT_SET, index=False, encoding="utf-8")
    if len(excluded):
        excluded[["row", "pdb_id", "mutation", "ddG", "exclude_reason"]].to_csv(
            OUT_EXC, index=False, encoding="utf-8")
    else:
        pd.DataFrame(columns=["row", "pdb_id", "mutation", "ddG",
                             "exclude_reason"]).to_csv(OUT_EXC, index=False,
                                                       encoding="utf-8")
    matched[["row", "pdb_id", "family_label", "seq40_cluster", "state_id",
             "obs_key"]].to_csv(OUT_FAM, index=False, encoding="utf-8")

    # ---------- coverage report ----------
    emit("## 覆盖统计")
    emit("")
    emit("| 项 | 值 |")
    emit("|---|---:|")
    emit(f"| 原始 forward 配对 | {n0} |")
    emit(f"| **matched_analysis_set** | **{len(matched)}** |")
    emit(f"| excluded | {len(excluded)} |")
    emit(f"| matched 唯一 pdb_id | {matched['pdb_id'].nunique()} |")
    emit(f"| matched 唯一 family | {matched['family_label'].nunique()} |")
    emit(f"| matched 唯一 state_id | {matched['state_id'].nunique()} |")
    emit(f"| ddG mean / std | {matched['ddG'].mean():+.4f} / "
         f"{matched['ddG'].std():.4f} |")
    emit(f"| ddG min / max | {matched['ddG'].min():+.4f} / "
         f"{matched['ddG'].max():+.4f} |")
    emit(f"| Base Rate (ΔΔG<0) | **{(matched['ddG'] < 0).mean():.4f}** |")
    emit("")
    emit("> **四档主比较全部基于 `matched_analysis_set`**；")
    emit("> 另报 full-coverage sensitivity（V6.2 §2.1）。")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.1f}s")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_SET}  ({len(matched)} 行)")
    print(f"写出 {OUT_EXC}  ({len(excluded)} 行)")
    print(f"写出 {OUT_FAM}")
    print(f"写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
