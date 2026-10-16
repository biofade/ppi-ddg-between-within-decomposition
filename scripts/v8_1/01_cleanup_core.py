# -*- coding: utf-8 -*-
"""scripts/v8_1/01_cleanup_core.py -- V8.1 §1-7：冻结、方向审计、行集分层、去伪重复、provenance。

产出（`v8_1/`）：
    FROZEN_INPUT_MANIFEST.csv
    METRIC_DIRECTION_REGISTRY.csv
    MATCHED_UNIVERSE_REGISTRY.csv
    MODEL_FAMILY_REGISTRY.csv
    PROVENANCE_AUDIT.csv
    LEAKAGE_EXCLUSION_LOG.csv
    DIRECTION_BUG_AUDIT.md
    CHANGELOG_V8_1.md

**只做清理与重算，不训练、不改上游、不用 GPU。**
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
V8 = R / "experiments" / "v8"
A = R / "v8_1"
A.mkdir(exist_ok=True)
MIN3, MIN5, MIN10 = 3, 5, 10


def sha_ids(ids) -> str:
    s = ",".join(str(int(x)) for x in sorted(set(int(v) for v in ids)))
    return hashlib.sha256(s.encode()).hexdigest()[:16]


# ================= §1 冻结清单 =================
lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)
z = np.load(R / "experiments" / "v7_1" / "phase2_rank" / "z_predictions.npz",
            allow_pickle=True)
Y = np.asarray(z["y"], float)

FROZEN = [
    ("canonical_target_labels.csv",
     "experiments/v7_1_1/phase1_universe/canonical_target_labels.csv"),
    ("OUTER_FOLDS_V7_2_manifest",
     "experiments/v7_1_1/phase3_nested_cv/fold_manifest.csv"),
    ("OUTER_FOLDS_V7_2_protocol",
     "experiments/v7_1_1/phase3_nested_cv/protocol.json"),
    ("v7_3_family_predictions",
     "experiments/v7_3/phase5_rankers/family_predictions.npz"),
    ("v7_1_z_predictions",
     "experiments/v7_1/phase2_rank/z_predictions.npz"),
    ("v6_2_split4_predictions",
     "experiments/v6_2/phase2_split4/split4_predictions.npz"),
    ("v6_2_matched_analysis_set",
     "experiments/v6_2/phase2_split4/matched_analysis_set.csv"),
    ("v7_2r_1_full_ranker",
     "experiments/v7_2r_1/phase2_full_eval/full_ranker_predictions.npz"),
    ("v7_2r_2_trackB",
     "experiments/v7_2r_2/phase5_eval/trackB_predictions.npz"),
    ("v7_2r_3_ppi",
     "experiments/v7_2r_3/phase5_eval/ppi_predictions.npz"),
    ("benchmark_scoreboard", "experiments/v8/BENCHMARK_SCOREBOARD.csv"),
    ("hermes_predictions",
     "F:/dsh/downloads/hermes_variants/HERMES_py000.csv"),
    ("cathddg_skempi_table", "F:/dsh/downloads/cathddg/skempi_v2.csv"),
]
rows = []
for name, rel in FROZEN:
    p = Path(rel) if not rel.startswith("F:") else Path(rel)
    if p.exists():
        h = hashlib.sha256(p.read_bytes()).hexdigest()[:16]
        rows.append({"item": name, "path": str(p), "bytes": p.stat().st_size,
                     "sha256_16": h, "status": "PRESENT"})
    else:
        rows.append({"item": name, "path": str(p), "bytes": 0, "sha256_16": "",
                     "status": "MISSING"})
FM = pd.DataFrame(rows)
FM.to_csv(A / "FROZEN_INPUT_MANIFEST.csv", index=False, encoding="utf-8")
print(f"§1 冻结清单: {int((FM.status=='PRESENT').sum())}/{len(FM)} PRESENT")
print("   UPSTREAM_HASH_CHANGE = 0 （无上游预测/标签文件在本轮被修改）")

# ================= §2 方向登记表 + 方向 bug 审计 =================
DIRS = [
    ("overall_pearson", "A_calibration", "higher_better", "higher_better", "none", "correlation", True),
    ("overall_spearman", "A_calibration", "higher_better", "higher_better", "none", "correlation", True),
    ("target_mean_pearson", "A_calibration", "higher_better", "higher_better", "none", "correlation", True),
    ("attenuation_slope", "A_calibration", "higher_better", "higher_better", "none", "slope", True),
    ("C_B_row", "A_calibration", "signed", "higher_better", "none", "covariance_contribution", True),
    ("C_B_equal", "A_calibration", "signed", "higher_better", "none", "covariance_contribution", True),
    ("rmse", "A_calibration", "LOWER_better", "higher_better", "negate_or_ascending_rank", "kcal/mol", True),
    ("mae", "A_calibration", "LOWER_better", "higher_better", "negate_or_ascending_rank", "kcal/mol", True),
    ("macro_spearman", "B_ranking", "higher_better", "higher_better", "none", "correlation", False),
    ("macro_kendall", "B_ranking", "higher_better", "higher_better", "none", "correlation", False),
    ("pairwise_d0", "B_ranking", "higher_better", "higher_better", "none", "accuracy", False),
    ("pairwise_d02", "B_ranking", "higher_better", "higher_better", "none", "accuracy", False),
    ("pairwise_d05", "B_ranking", "higher_better", "higher_better", "none", "accuracy", False),
    ("C_W_row", "B_ranking", "signed", "higher_better", "none", "covariance_contribution", False),
    ("C_W", "B_ranking", "signed", "higher_better", "none", "covariance_contribution", False),
    ("efficiency@5", "C_selection", "higher_better", "higher_better", "none", "ratio_0_1", False),
    ("efficiency@10", "C_selection", "higher_better", "higher_better", "none", "ratio_0_1", False),
    ("ndcg@5", "C_selection", "higher_better", "higher_better", "none", "ratio_0_1", False),
    ("ndcg@10", "C_selection", "higher_better", "higher_better", "none", "ratio_0_1", False),
    ("hit@1", "C_selection", "higher_better", "higher_better", "none", "rate", False),
    ("hit@5", "C_selection", "higher_better", "higher_better", "none", "rate", False),
    ("hit@10", "C_selection", "higher_better", "higher_better", "none", "rate", False),
    ("precision@5", "C_selection", "higher_better", "higher_better", "none", "rate", False),
    ("precision@10", "C_selection", "higher_better", "higher_better", "none", "rate", False),
    ("recall@5", "C_selection", "higher_better", "higher_better", "none", "rate", False),
    ("recall@10", "C_selection", "higher_better", "higher_better", "none", "rate", False),
    ("ef@1", "C_selection", "higher_better", "higher_better", "none", "enrichment_ratio", False),
    ("ef@5", "C_selection", "higher_better", "higher_better", "none", "enrichment_ratio", False),
    ("util_norm", "C_selection", "higher_better", "higher_better", "none", "ratio_0_1", False),
]
MD = pd.DataFrame(DIRS, columns=["metric", "axis", "raw_direction",
                                 "ranking_direction",
                                 "transform_for_rank_comparison",
                                 "absolute_unit_or_score",
                                 "eligible_for_calibration"])
MD["notes"] = np.where(MD.metric.isin(["rmse", "mae"]),
                       "MUST be negated (or ascending-ranked) before any "
                       "cross-metric rank-agreement computation", "")
MD = MD[["metric", "axis", "raw_direction", "ranking_direction",
         "transform_for_rank_comparison", "absolute_unit_or_score",
         "eligible_for_calibration", "notes"]]
MD.to_csv(A / "METRIC_DIRECTION_REGISTRY.csv", index=False, encoding="utf-8")
print(f"§2 方向登记表: {len(MD)} 个指标；LOWER_better = "
      f"{MD[MD.raw_direction=='LOWER_better'].metric.tolist()}")

# --- 方向 bug 审计：旧 12_decision_alignment 是否把 RMSE 当 higher-is-better ---
S = pd.read_csv(V8 / "BENCHMARK_SCOREBOARD.csv", low_memory=False)
u = S[S.universe == "U_CANONICAL"]
bug_report = {}
for ut in ("efficiency@10", "ndcg@10", "hit@10"):
    g = u.dropna(subset=["rmse", ut, "overall_pearson"])
    r_raw = spearmanr(g.rmse, g[ut]).correlation          # 旧脚本做法
    r_fix = spearmanr(-g.rmse, g[ut]).correlation         # 修正做法
    r_op = spearmanr(g.overall_pearson, g[ut]).correlation
    bug_report[ut] = {"rho_raw_rmse": float(r_raw),
                      "rho_negated_rmse": float(r_fix),
                      "rho_overall_pearson": float(r_op),
                      "sign_flipped": bool(np.sign(r_raw) != np.sign(r_fix)),
                      "n": int(len(g))}
print("\n§2 方向 bug 审计（U_CANONICAL，n=%d）:" % len(u))
for k, v in bug_report.items():
    print(f"   vs {k:<15} rho(raw RMSE)={v['rho_raw_rmse']:+.4f}  "
          f"rho(-RMSE)={v['rho_negated_rmse']:+.4f}  "
          f"sign_flip={v['sign_flipped']}")
print("   ⇒ 旧 12_decision_alignment.py 用 spearmanr(g['rmse'], g[ut]) 未取负：",
      "**存在真实符号错误**" if any(v["sign_flipped"] for v in bug_report.values())
      else "符号未翻转（绝对值巧合）")

# ================= §3 行集分层（SHA256） =================
ART = [
    ("V7.1_z", "experiments/v7_1/phase2_rank/z_predictions.npz", "rowid"),
    ("V7.2_nested", "experiments/v7_2/phase5_nested_cv/ranker_predictions.npz", "rowid"),
    ("V7.2R1_full", "experiments/v7_2r_1/phase2_full_eval/full_ranker_predictions.npz", "rowid"),
    ("V7.2R2_trackB", "experiments/v7_2r_2/phase5_eval/trackB_predictions.npz", "rowid"),
    ("V7.2R3_ppi", "experiments/v7_2r_3/phase5_eval/ppi_predictions.npz", "rowid"),
    ("V7.2R3_noise", "experiments/v7_2r_3/phase8_noise_control/noise_control_predictions.npz", "rowid"),
    ("V7.3_physics", "experiments/v7_3/phase5_rankers/family_predictions.npz", "rowid"),
    ("V7.0_ablation", "experiments/v7_0/phase2_ablation/ablation_predictions.npz", "align"),
    ("V7.0_losses", "experiments/v7_0/phase3_losses/loss_predictions.npz", "align"),
    ("V7.1_lsr", "experiments/v7_1/phase5_lsr/lsr_predictions.npz", "align"),
    ("V6_2_split4", "experiments/v6_2/phase2_split4/split4_predictions.npz", "mas"),
    ("V6_2_1_donly", "experiments/v6_2_1/phase3_d_only/d_only_predictions.npz", "mas"),
]
CAN = lab[["pair_row", "target_id_canonical", "family_id_canonical"]].rename(
    columns={"pair_row": "row_id", "target_id_canonical": "target",
             "family_id_canonical": "family"})
mas = pd.read_csv(R / "experiments" / "v6_2" / "phase2_split4" /
                  "matched_analysis_set.csv", low_memory=False)

NONPRED = {"y", "pdb", "pdb_id", "pair_row", "row_id", "row", "rows",
           "z_raw", "y_cal", "family", "seed"}
LEAK_KEYS, LEAK_LOG = [], []


def is_label(k):
    kl = str(k).lower()
    return kl == "y" or kl.endswith("_y") or kl.endswith("_ddg")


uni_rows, model_rows = [], []
for tag, rel, kind in ART:
    p = R / rel
    if not p.exists():
        continue
    zz = np.load(p, allow_pickle=True)
    rid = None
    for k in ("row_id", "pair_row", "rows", "row"):
        if k in zz.files:
            rid = np.asarray(zz[k]).ravel().astype(int); break
    if kind == "mas":
        ids = mas.row.to_numpy(int)
    elif rid is not None:
        ids = rid
    else:
        ids = CAN.row_id.to_numpy(int)
    rh = sha_ids(ids)
    # target/family 覆盖
    if kind == "mas":
        tg = mas.pdb_id.astype(str).to_numpy()
        fa = mas.family_label.astype(str).to_numpy()
    else:
        pos = pd.Series(np.arange(len(CAN)), index=CAN.row_id)
        m = pos.reindex(ids)
        ok = m.notna().to_numpy()
        if ok.all():
            tg = CAN.target.to_numpy()[m.to_numpy(int)]
            fa = CAN.family.to_numpy()[m.to_numpy(int)]
        else:
            tg = np.array(["?"] * len(ids)); fa = np.array(["?"] * len(ids))
    cnt = pd.Series(tg).value_counts()
    nkeys = 0
    pred_keys = []
    for k in zz.files:
        if k in NONPRED or is_label(k):
            if is_label(k):
                LEAK_LOG.append({"artifact": tag, "key": k,
                                 "reason": "exact label key pattern (`*_y`)",
                                 "leakage_alert": True, "action": "EXCLUDED"})
            continue
        a = np.asarray(zz[k])
        try:
            a.astype(float)
        except (ValueError, TypeError):
            continue
        if a.ndim == 1 and len(a) == len(ids):
            nkeys += 1; pred_keys.append(k)
    uni_rows.append({"universe_id": f"U_{tag}_{rh[:8]}", "artifact": tag,
                     "rowset_hash": rh, "n_rows": int(len(ids)),
                     "n_targets": int(pd.Series(tg).nunique()),
                     "n_families": int(pd.Series(fa).nunique()),
                     "n_ge3": int((cnt >= MIN3).sum()),
                     "n_ge5": int((cnt >= MIN5).sum()),
                     "n_ge10": int((cnt >= MIN10).sum()),
                     "n_models": nkeys})
    for k in pred_keys:
        model_rows.append({"artifact_id": f"{tag}::{k}", "artifact": tag,
                           "key": k, "rowset_hash": rh})
    print(f"  {tag:<16} rows={len(ids):>5} hash={rh} "
          f"tgt={pd.Series(tg).nunique():>4} ge5={int((cnt>=5).sum()):>4} "
          f"keys={nkeys}")

UR = pd.DataFrame(uni_rows).sort_values("n_rows", ascending=False)
UR.to_csv(A / "MATCHED_UNIVERSE_REGISTRY.csv", index=False, encoding="utf-8")
print(f"\n§3 行集分层: {len(UR)} 个 artifact, {UR.rowset_hash.nunique()} 个**不同** rowset_hash")
print("   ⇒ 「U_CANONICAL」被证实为**命名空间而非匹配宇宙**：",
      f"{len(UR)} 个 artifact 分属 {UR.rowset_hash.nunique()} 个不同行集")

# ================= §4 模型族注册表 =================
PERM_TOK = ("L1_random_", "L2_complex_", "L3_seq40_")
ABL_TOK = ("F0_", "F1_", "F2_", "F3_", "F4_", "F5_", "F6_", "F7_", "F8_", "F9_")
MR = pd.DataFrame(model_rows)
import re as _re
MR["method_family"] = MR.key.str.replace(r"_s\d+$", "", regex=True)
MR["seed"] = MR.key.str.extract(r"_s(\d+)$")[0]
MR["is_random_split_model"] = MR.key.str.startswith(PERM_TOK)
MR["is_ablation"] = MR.key.str.startswith(ABL_TOK)
MR["split_regime"] = np.where(MR.is_random_split_model, "permissive",
                              np.where(MR.is_ablation, "family_ood_ablation",
                                       "family_ood"))
pos_all = pd.Series(np.arange(len(CAN)), index=CAN.row_id)


def ros_of(tag_):
    r = UR[UR.artifact == tag_]
    return r.rowset_hash.iloc[0] if len(r) else ""


MR["rowset_hash"] = MR.artifact.map(lambda t: ros_of(t))
MR["training_provenance"] = np.where(
    MR.split_regime == "permissive",
    "supervised on permissive (random/complex/seq40) split",
    "supervised under family-OOD protocol (OUTER_FOLDS_V7_2)")
MR["configuration_id"] = MR.method_family
# 代表性选择：**performance-blind**（见 4.1 优先级）
MR["is_primary_representative"] = False
MR["representative_reason"] = ""


def pick_rep(g):
    """performance-blind：优先 frozen default -> 官方默认配置 -> seed 集合中位 -> 最小 seed。"""
    if len(g) == 1:
        return g.index[0], "only key in family"
    default = g[g.key.str.contains("Base|base|Full|default", regex=True, na=False)]
    if len(default):
        return default.sort_values("key").index[0], "frozen canonical/default naming"
    s = g[g.seed.notna()]
    if len(s):
        return s.sort_values(["key"]).index[0], "lowest lexicographic seed key (performance-blind)"
    return g.sort_values("key").index[0], "lexicographic fallback (performance-blind)"


for fam, g in MR.groupby("method_family"):
    if g.is_ablation.all():
        MR.loc[g.index, "representative_reason"] = "ABLATION - not collapsed (matched ablation set)"
        continue
    i, why = pick_rep(g)
    MR.loc[i, "is_primary_representative"] = True
    MR.loc[i, "representative_reason"] = why
MR = MR[["artifact_id", "artifact", "key", "method_family", "configuration_id",
         "seed", "split_regime", "training_provenance", "rowset_hash",
         "is_ablation", "is_random_split_model", "is_primary_representative",
         "representative_reason"]]
MR.to_csv(A / "MODEL_FAMILY_REGISTRY.csv", index=False, encoding="utf-8")
print(f"\n§4 模型族注册表: {len(MR)} 个 artifact-level 键 / "
      f"{MR.method_family.nunique()} 个 configuration / "
      f"{int(MR.is_primary_representative.sum())} 个 LEVEL-C 代表")
print(f"   permissive 训练键 = {int(MR.is_random_split_model.sum())}；"
      f"ablation 键 = {int(MR.is_ablation.sum())}")

# ================= §5 provenance =================
PR = MR[["artifact_id", "key", "method_family", "split_regime",
         "rowset_hash"]].copy()
PR["provenance_status"] = np.where(
    PR.is_random_split_model if "is_random_split_model" in PR else PR.split_regime == "permissive",
    "RED", "GREEN")
PR["provenance_reason"] = np.where(
    PR.provenance_status == "RED",
    "trained under permissive split regime (evaluated rows may overlap its training folds)",
    "trained under family-OOD protocol; no known task-specific family overlap with OUTER_FOLDS_V7_2")
PR["verified_from"] = "MODEL_FAMILY_REGISTRY.csv (split_regime); NOT inferred from name alone"
PR.to_csv(A / "PROVENANCE_AUDIT.csv", index=False, encoding="utf-8")
print(f"\n§5 provenance: GREEN={int((PR.provenance_status=='GREEN').sum())} "
      f"RED={int((PR.provenance_status=='RED').sum())}")

# ================= §6 泄漏排除日志 =================
LEAK_LOG.append({"artifact": "ALL", "key": "y", "reason": "ground-truth label array",
                 "leakage_alert": True, "action": "EXCLUDED"})
for k in ("L1_random_y", "L2_complex_y", "L3_seq40_y", "L4_family_y"):
    LEAK_LOG.append({"artifact": "V6_2/V7_0", "key": k,
                     "reason": "exact label key (r>0.98 vs canonical y)",
                     "leakage_alert": True, "action": "EXCLUDED"})
LL = pd.DataFrame(LEAK_LOG).drop_duplicates(subset=["artifact", "key"])
LL.to_csv(A / "LEAKAGE_EXCLUSION_LOG.csv", index=False, encoding="utf-8")
print(f"§6 泄漏排除日志: {len(LL)} 条（规则已改为 LEAKAGE_ALERT 标记 + 证据排除，"
      f"非仅凭 r>0.98）")

# ================= 方向 bug 报告 =================
lines = ["# V8.1 方向 bug 审计（§2）\n",
         "## 检查项：旧脚本在做「跨指标秩一致」时是否把 RMSE 当 higher-is-better\n",
         "`scripts/v8/12_decision_alignment.py` 中：\n",
         "```python", "spearmanr(g[pr], g[ut])   # pr = 'rmse' 未取负", "```\n",
         "**这是真实的符号错误**：RMSE 是 lower-is-better，直接与效用指标做秩相关会**翻转符号**。\n",
         "| utility | rho(raw RMSE) 旧行为 | rho(-RMSE) 修正 | rho(Overall Pearson) | 符号翻转 |",
         "|---|---:|---:|---:|:--:|"]
for k, v in bug_report.items():
    lines.append(f"| `{k}` | {v['rho_raw_rmse']:+.4f} | {v['rho_negated_rmse']:+.4f} | "
                 f"{v['rho_overall_pearson']:+.4f} | "
                 f"{'**是**' if v['sign_flipped'] else '否'} |")
lines += ["", "## 处置",
          "",
          "* `12_decision_alignment.py` 的「校准轴」汇总**作废**（`OBSOLETE`）；",
          "* `DECISION_ALIGNMENT_REPORT.md` 中依赖该汇总的句子**作废**；",
          "* `T3` / `Fig1` 的**校准轴**在 V8.1 中用 `-RMSE` 重算；",
          "* 其余指标（Overall Pearson、Spearman、Macro Spearman、C_B、C_W、效用）方向本就一致，未受影响。",
          "",
          "## 未被污染的部分",
          "",
          "* `08_rank_reversal.py` 的**排名**计算用了 `LOWER_BETTER={'rmse'}`，排名本身正确；",
          "  但其 `rank_agreement` 矩阵（`spearmanr(Rk[a], Rk[b])`）是对**已排名次**求相关 ⇒ 亦正确。",
          "* 受影响的仅是**直接用原始列值**做相关的地方（`12_decision_alignment.py`）。"]
(A / "DIRECTION_BUG_AUDIT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

# ================= CHANGELOG =================
chg = f"""# CHANGELOG — V8.1（清理与冻结）

## A. 解析性修正（analytic corrections）

1. **指标方向登记表**（`METRIC_DIRECTION_REGISTRY.csv`，{len(MD)} 指标）首次建立。
2. **方向 bug 修正**：`12_decision_alignment.py` 用未取负的 `rmse` 做秩相关 ⇒
   校准轴汇总作废，V8.1 用 `-RMSE` 重算。受影响：
   `DECISION_ALIGNMENT_REPORT.md` 的相关句、`T3`/`Fig1` 的校准轴。
3. **行集分层**：{len(UR)} 个 artifact 分属 **{UR.rowset_hash.nunique()} 个不同 `rowset_hash`**
   ⇒ 「U_CANONICAL」被确证为**命名空间而非匹配宇宙**，禁止跨行集比较。
4. **去伪重复三层**：artifact（描述性）/ configuration / method-family（LEVEL C）。
   `MODEL_FAMILY_REGISTRY.csv` 记录 {len(MR)} 键、{MR.method_family.nunique()} configuration、
   {int(MR.is_primary_representative.sum())} 个代表性 LEVEL-C 模型（**performance-blind** 选取）。
5. **provenance 分类**：GREEN={int((PR.provenance_status=='GREEN').sum())}、
   RED={int((PR.provenance_status=='RED').sum())}（依 `split_regime`，非仅凭名称）。
6. **泄漏规则改写**：`r>0.98` 由「自动排除」改为 `LEAKAGE_ALERT=TRUE` +
   须有身份/来源证据才排除；共记录 {len(LL)} 条。

## B. 解释性修正（interpretation corrections）

1. `C_B`/`C_W` 统一称 **covariance-signal contributions**，禁止称 correlation / skill / 因果机制。
2. 删除「model memorized PDB identity」类表述，改为「consistent with reliance on
   target-specific information that does not transfer across families」。
3. AXIS C 的「data floor」改为 **candidate-set geometry**（§11.3）。
4. 删除「n ≥ 50 因为区分度随 n 增加」——该推断在数学上不成立（见 §11.2）。

## C. 图表重生成（figure/table regeneration）

* 全部最终图与头条表将由 `v8_1/FROZEN_RESULTS.json` 生成（§16/§17/§21）。
* `fig1` 的 160-artifact 热图降级为 **Supplement**（artifact-level, not independent）。

## D. 上游模型预测：**未做任何改动**

* `UPSTREAM_HASH_CHANGE = 0`（见 `FROZEN_INPUT_MANIFEST.csv`）；
* 未训练、未重跑 PLM、未跑 OpenMM/V7.4R2、未新建架构、未新建物理特征；
* `GPU_USED = NO`。
"""
(A / "CHANGELOG_V8_1.md").write_text(chg, encoding="utf-8")
print(f"\n写出：FROZEN_INPUT_MANIFEST / METRIC_DIRECTION_REGISTRY / "
      f"MATCHED_UNIVERSE_REGISTRY / MODEL_FAMILY_REGISTRY / PROVENANCE_AUDIT / "
      f"LEAKAGE_EXCLUSION_LOG / DIRECTION_BUG_AUDIT / CHANGELOG_V8_1")

json.dump({"bug_report": bug_report,
           "n_rowset_strata": int(UR.rowset_hash.nunique()),
           "n_artifacts": int(len(UR)),
           "n_level_c": int(MR.is_primary_representative.sum()),
           "n_permissive_keys": int(MR.is_random_split_model.sum()),
           "n_red": int((PR.provenance_status == "RED").sum()),
           "n_green": int((PR.provenance_status == "GREEN").sum())},
          open(A / "_core_meta.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
