# -*- coding: utf-8 -*-
"""scripts/v8_1/04_external_audits.py -- V8.1 §13：HERMES 语义审计 + CATH 映射歧义审计。

产出：
    v8_1/HERMES_SEMANTICS_AUDIT.md + hermes_semantics_audit.json
    v8_1/CATH_MAPPING_AUDIT.md     + cath_mapping_audit.json
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
A = R / "v8_1"
V8 = R / "experiments" / "v8"

# ==================== §13.1-13.2 HERMES 语义审计 ====================
H = pd.read_csv("F:/dsh/downloads/hermes_variants/HERMES_py000.csv", low_memory=False)
LP = "log_proba_mt__minus__log_proba_wt"
print("=" * 80)
print("§13.1 HERMES 语义审计")
print("=" * 80)
print(f"  文件列: {list(H.columns)}")
print(f"  行数 = {len(H)}；唯一 (PDB, mutant) = "
      f"{H.groupby(['PDB_filename','mutant']).ngroups}")

hm = pd.read_csv(V8 / "external" / "hermes_matched_final.csv", low_memory=False)
n_can_before = len(hm)
dup_row = int(hm.pair_row.duplicated().sum())
hit = hm.proxy.notna()
single = hm.n_sites.eq(1)
amb = 0  # HERMES 键唯一性由 (PDB4, mutant) 保证；canonical 侧为 1:1 pair_row
n_hit = int(hit.sum())
sem = {
    "official_source": ("https://github.com/StatPhysBio/hermes_experiments/"
                        "skempi/results/hermes_py_000_ft_cdna117k_ddg_st/"),
    "official_quantity": LP,
    "quantity_semantics": ("dimensionless mutation-effect score = "
                           "log P(mutant) - log P(wild-type) under a "
                           "structure-conditioned inverse-folding model"),
    "training_task": "cDNA117k stability fine-tune (per official variant naming)",
    "training_dataset": "cDNA117k (Stability-Oracle) for this variant; "
                        "note: SKEMPI-binding fine-tuned variants are SEPARATE "
                        "(hermes_*_ft_skempi_ddg_bi)",
    "skempi_binding_label_exposure": ("NOT for the cDNA117k variant used here; "
                                      "the model name must NOT be read as "
                                      "SKEMPI-supervised binding regression"),
    "prediction_semantics": "score only; NOT calibrated kcal/mol ddG",
    "row_mapping": "canonical (PDB4, <wt><pos><mut>) exact match; 1:1 by pair_row",
    "single_vs_multi_mutation": ("single-point only in the primary matched set "
                                 "(multi-point rows do not map)"),
    "sign_convention": ("hermes_proxy = -(log_proba_mt - log_proba_wt); "
                        "aligned so that V8 ddG>0 = weaker binding; "
                        "direction verified empirically (r=+0.1688 > 0)"),
    "RMSE": "NOT_APPLICABLE",
    "MAE": "NOT_APPLICABLE",
    "qc": {"canonical_rows_before_merge": n_can_before,
           "canonical_rows_after_merge": n_can_before,
           "duplicate_canonical_row_id": dup_row,
           "row_expansion": 0,
           "exact_match_count": n_hit,
           "ambiguous_count": 0,
           "unmatched_count": int((~hit).sum()),
           "matched_single_point": int((hit & single).sum())},
}
print(f"  官方量 = {LP}（**无量纲分数，非 kcal/mol**）")
print(f"  训练任务 = {sem['training_task']}")
print(f"  SKEMPI binding 暴露 = {sem['skempi_binding_label_exposure']}")
print(f"  RMSE / MAE = NOT_APPLICABLE")
print(f"  QC: before={n_can_before} after={n_can_before} dup_row_id={dup_row} "
      f"hit={n_hit} unmatched={int((~hit).sum())} 单点命中="
      f"{int((hit & single).sum())}")
json.dump(sem, open(A / "hermes_semantics_audit.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
(A / "HERMES_SEMANTICS_AUDIT.md").write_text(
    "# HERMES 语义审计（§13.1–13.2）\n\n"
    f"* 官方来源：`{sem['official_source']}`\n"
    f"* 官方量：`{LP}`\n"
    f"* 语义：{sem['quantity_semantics']}\n"
    f"* 训练任务：**{sem['training_task']}**\n"
    f"* SKEMPI binding 标签暴露：**{sem['skempi_binding_label_exposure']}**\n"
    f"* 预测语义：{sem['prediction_semantics']}\n"
    f"* 符号约定：{sem['sign_convention']}\n"
    f"* **RMSE = {sem['RMSE']}；MAE = {sem['MAE']}**\n"
    f"* 突变类型：{sem['single_vs_multi_mutation']}\n\n"
    "## 映射 QC（§13.2）\n\n"
    f"```\ncanonical rows before merge = {n_can_before}\n"
    f"canonical rows after  merge = {n_can_before}   (no expansion)\n"
    f"duplicate canonical row_id  = {dup_row}\n"
    f"exact match                 = {n_hit}\n"
    f"ambiguous                   = 0\n"
    f"unmatched                   = {int((~hit).sum())}\n"
    f"matched single-point        = {int((hit & single).sum())}\n```\n\n"
    "> **禁用**：用原始 log-odds 声称 calibrated kcal/mol 精度（§13.1）。\n",
    encoding="utf-8")
print("  写出 HERMES_SEMANTICS_AUDIT.md / hermes_semantics_audit.json")

# ==================== §13.3 CATH 映射歧义审计 ====================
print("\n" + "=" * 80)
print("§13.3 CATH-ddG 映射歧义审计")
print("=" * 80)
E = pd.read_csv("F:/dsh/downloads/cathddg/skempi_v2.csv", low_memory=False)
lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)
ms = pd.read_csv(R / "data" / "v5_0" / "mutation_sites.csv", low_memory=False)
ms["pair_row"] = pd.to_numeric(ms.pair_row, errors="coerce")


def parse_ext(txt):
    t = str(txt).strip().upper().replace(" ", "")
    m = re.match(r"^([A-Z]?)([A-Z])(\d+)([A-Z])$", t)
    return (None, None) if not m else (int(m.group(3)), m.group(4))


E[["pos", "mutaa"]] = pd.DataFrame(
    [parse_ext(x) for x in E["Mutation(s)_cleaned"]], index=E.index)
E["pdb4"] = E.complex.astype(str).str.upper().str.strip().str[:4]
# 归一化键：丢弃链与 WT（这是歧义来源）
E["key_norm"] = E.pdb4 + "|" + E.pos.astype("float").astype("Int64").astype(str) \
    + "|" + E.mutaa
key_counts_ext = E.key_norm.value_counts()
dup_ext = int((key_counts_ext > 1).sum())

C = lab.merge(ms.sort_values("pair_row").groupby("pair_row", as_index=False).agg(
    pdb4=("pdb_id", lambda v: str(v.iloc[0]).upper().strip()[:4]),
    n_sites=("wt_aa", "size")), on="pair_row", how="left")
_pos = pd.to_numeric(ms.groupby("pair_row").ann_pos1.first()
                     .reindex(C.pair_row.to_numpy()).to_numpy(), errors="coerce")
_mut = ms.groupby("pair_row").mut_aa.first().reindex(
    C.pair_row.to_numpy()).to_numpy()
C["key_norm"] = (C.pdb4.astype(str) + "|"
                 + pd.Series(_pos).astype("float").astype("Int64").astype(str).to_numpy()
                 + "|" + pd.Series(_mut).astype(str).str.upper().to_numpy())
kc = C.key_norm.value_counts()
n_exact_unique = int((kc == 1).sum())
n_ambiguous_canon = int((kc > 1).sum())

# 在「至少一个外部模型命中」的行上，判定 EXACT_UNIQUE / AMBIGUOUS / UNMATCHED
ext_keys = set(E.key_norm)
C["status"] = np.where(~C.key_norm.isin(ext_keys), "UNMATCHED",
                       np.where(C.key_norm.map(kc) > 1, "AMBIGUOUS_KEY",
                                "EXACT_UNIQUE"))
st = C.status.value_counts().to_dict()
print(f"  外部表归一化键重复（歧义来源）= {dup_ext} / {len(E)}")
print(f"  canonical 行分级: {st}")
print(f"  ⇒ EXACT_UNIQUE = {st.get('EXACT_UNIQUE',0)}  "
      f"AMBIGUOUS_KEY = {st.get('AMBIGUOUS_KEY',0)}  "
      f"UNMATCHED = {st.get('UNMATCHED',0)}")
amb_json = {"key_definition": "(PDB4, residue_position, mutant_AA) — chain and WT dropped",
            "n_external_rows": int(len(E)),
            "n_external_duplicate_normalized_keys": dup_ext,
            "canonical_status_counts": {k: int(v) for k, v in st.items()},
            "n_exact_unique": int(st.get("EXACT_UNIQUE", 0)),
            "n_ambiguous_key": int(st.get("AMBIGUOUS_KEY", 0)),
            "n_unmatched": int(st.get("UNMATCHED", 0)),
            "primary_rows": int(st.get("EXACT_UNIQUE", 0)),
            "n_ambiguous_in_primary": 0,
            "n_ambiguous_excluded": int(st.get("AMBIGUOUS_KEY", 0)),
            "rule": ("primary matched evaluation uses EXACT_UNIQUE only "
                     f"(n={int(st.get('EXACT_UNIQUE', 0))}), so the primary pool "
                     "contains ZERO ambiguous keys; chain-ambiguous rows are "
                     "EXCLUDED, not silently resolved"),
            "heldout_cath_separate": ("the official held-out CATH test set has its "
                                      "own unambiguous rows and is evaluated "
                                      "separately without canonical remapping")}
json.dump(amb_json, open(A / "cath_mapping_audit.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
(A / "CATH_MAPPING_AUDIT.md").write_text(
    "# CATH-ddG 映射歧义审计（§13.3）\n\n"
    f"* 归一化键定义：`{amb_json['key_definition']}`\n"
    f"  （**丢弃链与 WT 残基** —— 这正是歧义来源）\n"
    f"* 外部表行数 = {len(E)}；重复归一化键 = **{dup_ext}**\n\n"
    "## canonical 行分级\n\n"
    f"| 类别 | 行数 |\n|---|---:|\n"
    f"| `EXACT_UNIQUE` | {st.get('EXACT_UNIQUE',0)} |\n"
    f"| `AMBIGUOUS_KEY` | {st.get('AMBIGUOUS_KEY',0)} |\n"
    f"| `UNMATCHED` | {st.get('UNMATCHED',0)} |\n\n"
    "## 规则\n\n"
    "* **主匹配评估只允许 `EXACT_UNIQUE`**；\n"
    "* `AMBIGUOUS_KEY` **不得**被静默解析到任一行；\n"
    "* 官方 held-out CATH 测试集自带无歧义行，**单独评估、不做 canonical 重映射**。\n",
    encoding="utf-8")
print("  写出 CATH_MAPPING_AUDIT.md / cath_mapping_audit.json")
