# -*- coding: utf-8 -*-
"""scripts/v8/03_model_registry.py -- V8 PHASE 3：模型注册表（§16/§17）。

盘点**内部**模型的每个预测键，判定：
    * 连接方式（row_id / 行对齐 canonical / 行对齐 matched_analysis_set）
    * 覆盖（rows / targets / families / n>=5）
    * 证据层级与溯源状态
输出 `experiments/v8/model_registry.csv` + `registry/keys_inventory.csv`。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
REG = R / "experiments" / "v8" / "registry"
REG.mkdir(parents=True, exist_ok=True)
OUTCSV = R / "experiments" / "v8" / "model_registry.csv"

lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)
mas = pd.read_csv(R / "experiments" / "v6_2" / "phase2_split4" /
                  "matched_analysis_set.csv", low_memory=False)
CANON = lab[["pair_row", "target_id_canonical", "family_id_canonical"]].rename(
    columns={"pair_row": "row_id", "target_id_canonical": "target_id",
             "family_id_canonical": "family_id"})
CANON_ROWS = CANON.row_id.to_numpy(int)
MAS_ROWS = mas.row.to_numpy(int)
MAS_TGT = mas.pdb_id.astype(str).to_numpy()
# canonical 的 target 序列（与 v7_0/v7_1 制品的 pdb 数组比对）
CANON_TGT_SEQ = CANON.target_id.to_numpy(str)

# artifact -> 连接方式
#   canonical_rowid : 制品内含 row_id / pair_row（逐位等于 canonical）
#   canonical_align : 无 row_id，但 pdb 序列与 canonical 100% 一致（行对齐）
#   mas_align       : 无 row_id，行对齐 matched_analysis_set(7523)
ARTIFACTS = [
    ("experiments/v7_1/phase2_rank/z_predictions.npz", "canonical_rowid"),
    ("experiments/v7_2r_1/phase2_full_eval/full_ranker_predictions.npz", "canonical_rowid"),
    ("experiments/v7_2r_2/phase5_eval/trackB_predictions.npz", "canonical_rowid"),
    ("experiments/v7_2r_3/phase5_eval/ppi_predictions.npz", "canonical_rowid"),
    ("experiments/v7_2r_3/phase8_noise_control/noise_control_predictions.npz", "canonical_rowid"),
    ("experiments/v7_3/phase5_rankers/family_predictions.npz", "canonical_rowid"),
    ("experiments/v7_2/phase5_nested_cv/ranker_predictions.npz", "canonical_rowid"),
    ("experiments/v7_0/phase2_ablation/ablation_predictions.npz", "canonical_align"),
    ("experiments/v7_0/phase3_losses/loss_predictions.npz", "canonical_align"),
    ("experiments/v7_1/phase5_lsr/lsr_predictions.npz", "canonical_align"),
    ("experiments/v6_2/phase2_split4/split4_predictions.npz", "mas_align"),
    ("experiments/v6_2_1/phase3_d_only/d_only_predictions.npz", "mas_align"),
]
# 非预测键
NONPRED = {"y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows", "z_raw",
           "y_cal", "family", "seed"}
# 家族/前缀 -> 模型族
def family_of(key: str) -> str:
    k = key
    for p, fam in (("F", "V7.3_physics"), ("R", "V7.2R1_representation"),
                   ("B", "V7.2R2_strong_repr"), ("A", "V7.2R2_small_plm"),
                   ("NC", "V7.2R3_noise_control"), ("L", "V6_2_split4"),
                   ("M", "V6_2_M"), ("p_", "V7.1_LSR"), ("z_", "V7.1_rank"),
                   ("s_", "V7.1_rank_scaled"), ("L4_family", "V7.0_ablation"),
                   ("L0_", "V7.0_losses")):
        if k.startswith(p):
            return fam
    return "other"


rows, keyrows = [], []
for rel, join in ARTIFACTS:
    p = R / rel
    if not p.exists():
        print("MISSING", rel); continue
    z = np.load(p, allow_pickle=True)
    keys = list(z.files)
    n = None
    # 基准行/目标
    if join == "canonical_rowid":
        rid = None
        for k in ("row_id", "pair_row", "rows", "row"):
            if k in keys:
                rid = np.asarray(z[k]).ravel().astype(int); break
        if rid is None and "pdb" in keys:
            rid = None
        base_rows = rid if rid is not None else CANON_ROWS[:len(np.asarray(z[keys[0]]))]
        tgt = np.asarray(z["pdb"]).astype(str) if "pdb" in keys else None
        fam = CANON.set_index("row_id").family_id.reindex(base_rows).to_numpy()
    elif join == "canonical_align":
        base_rows = CANON_ROWS
        tgt = np.asarray(z["pdb"]).astype(str) if "pdb" in keys else CANON_TGT_SEQ
        fam = CANON.family_id.to_numpy()
        # 校验对齐
        if "pdb" in keys:
            eq = (np.asarray(z["pdb"]).astype(str) == CANON_TGT_SEQ).mean()
            if eq < 0.999:
                print(f"  !! {rel} 行对齐校验失败 eq={eq:.4f}")
    else:  # mas_align
        base_rows = MAS_ROWS
        tgt = MAS_TGT
        fam = mas.family_label.astype(str).to_numpy()

    for k in keys:
        if k in NONPRED:
            continue
        a = np.asarray(z[k])
        if a.ndim != 1 or a.size != len(base_rows):
            continue
        # 只接受数值型预测（排除 per-fold 的字符串键，如 L1_random_pdb）
        try:
            af = a.astype(np.float64)
        except (ValueError, TypeError):
            continue
        fin = np.isfinite(af)
        if fin.sum() == 0:
            continue
        g = pd.Series(tgt).astype(str)
        sizes = g.value_counts()
        n5 = int((sizes >= 5).sum())
        keyrows.append({
            "artifact": rel, "key": k, "join": join,
            "n_rows": int(a.size), "n_finite": int(fin.sum()),
            "n_targets": int(g.nunique()), "n_families": int(pd.Series(fam).nunique()),
            "n_targets_ge5": n5,
            "model_family": family_of(k),
        })
        rows.append({
            "model_id": f"{Path(rel).stem}::{k}",
            "model_name": k, "model_family": family_of(k),
            "resource_tier": "TIER_1_INTERNAL",
            "provenance_status": "GREEN" if join.startswith("canonical") else "YELLOW",
            "prediction_file": rel, "prediction_key": k, "join_mode": join,
            "official_source": "internal historical artifact",
            "training_source": "internal",
            "supervised_on_skempi": "YES_SHARED_LABELS",
            "known_training_families": "see OUTER_FOLDS_V7_2",
            "split_provenance": "OUTER_FOLDS_V7_2 (family-OOD)",
            "prediction_semantics": "ddG (positive = weaker binding)",
            "ddg_sign_convention": "+ddG=weaker",
            "n_rows": int(a.size), "n_finite": int(fin.sum()),
            "n_targets": int(g.nunique()),
            "n_families": int(pd.Series(fam).nunique()),
            "n_targets_ge5": n5,
            "usable_primary": True, "usable_secondary": True,
            "reason_if_not": "",
        })

K = pd.DataFrame(keyrows)
Md = pd.DataFrame(rows)
K.to_csv(REG / "keys_inventory.csv", index=False, encoding="utf-8")
Md.to_csv(OUTCSV, index=False, encoding="utf-8")
print(f"预测键总数 = {len(K)}")
print(f"模型注册条目 = {len(Md)}")
print()
print("=== 按连接方式 ===")
print(K.groupby("join").agg(keys=("key", "size"),
                            models=("artifact", "nunique")).to_string())
print()
print("=== 按模型族 ===")
print(K.groupby("model_family").agg(keys=("key", "size"),
                                    max_rows=("n_rows", "max"),
                                    max_tgt=("n_targets", "max"),
                                    max_ge5=("n_targets_ge5", "max")).to_string())
print()
print("=== 最大宇宙的键（n>=5 target 最多）===")
print(K.nlargest(15, "n_targets_ge5")[["key", "join", "n_rows", "n_targets",
                                       "n_families", "n_targets_ge5"]].to_string(index=False))
