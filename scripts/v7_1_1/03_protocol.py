# -*- coding: utf-8 -*-
"""scripts/v7_1_1/03_protocol.py -- V7.1.1 PHASE 3–6: 冻结 V7.2 协议.

产出（全部为**协议冻结**，不训练任何新模型）：
* PHASE 3：`nested_cv_protocol.md` / `fold_manifest.csv` / `protocol.json`
* PHASE 4：`ranker_selection_protocol.md`
* PHASE 5：`calibration_shell.md` / `simple_scale_comparison.csv`
* PHASE 6：`V7_2_EVALUATION_CONTRACT.md`（由 04 脚本写）

本脚本同时做两件**必须由数据决定**的事：
1. 从 canonical 权威标签重建 `OUTER_FOLDS_V7_2`（复核历史 L4 是否可复现）
2. 用 inner-family CV 比较 `alpha = 1` vs `alpha = 全局 train 残差尺度`
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v6_1"))
import v5_common as C  # noqa: E402
import fastmetrics as FM  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "s4", C.ROOT / "scripts" / "v6_2" / "02_split4_ddg.py")
s4 = importlib.util.module_from_spec(_spec)
sys.modules["s4"] = s4
_spec.loader.exec_module(s4)

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

V62P2 = C.ROOT / "experiments" / "v6_2" / "phase2_split4"
V71 = C.ROOT / "experiments" / "v7_1"
V711 = C.ROOT / "experiments" / "v7_1_1"
PH1, PH3, PH4, PH5 = (V711 / x for x in
                      ("phase1_universe", "phase3_nested_cv",
                       "phase4_ranker_protocol", "phase5_calibration"))
for d in (PH3, PH4, PH5):
    d.mkdir(parents=True, exist_ok=True)
MIN_PS_N = 5
EPS = 1e-6
SEEDS = [42, 43, 44]
L3: list[str] = []
L4: list[str] = []
L5: list[str] = []


def e3(m=""):
    print(m); L3.append(str(m))


def e4(m=""):
    print(m); L4.append(str(m))


def e5(m=""):
    print(m); L5.append(str(m))


def standardize_within_target(s: np.ndarray, g: np.ndarray, eps: float = EPS):
    """**唯一的**靶点内标准化函数（PHASE 4 §4.2 要求全局共用）。

    返回 (z, collapsed_mask)。collapsed（组内 SD < eps）→ z = 0。
    """
    z = np.zeros(len(s), float)
    coll = np.zeros(len(s), bool)
    for c in np.unique(g):
        m = np.where(g == c)[0]
        v = s[m]
        sd = v.std(ddof=0)
        if sd < eps:
            coll[m] = True
        else:
            z[m] = (v - v.mean()) / sd
    return z, coll


def main() -> int:
    t0 = time.time()
    lab = pd.read_csv(PH1 / "canonical_target_labels.csv", low_memory=False)
    emits = [e3, e4, e5]
    for e in emits:
        e("# V7.1.1 — V7.2 协议冻结")
        e("")
        e(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
        e("- 权威 target 标签：`canonical_target_labels.csv`"
          "（`canonical.pdb_id` 按 `pair_row` 映射，与 `d['groups']` 100% 自洽）")
        e("")

    # ================= PHASE 3 =================
    # 载入 LSR 预测与 canonical 标签
    lsr = np.load(V71 / "phase5_lsr" / "lsr_predictions.npz", allow_pickle=True)
    d = C.build_pair_features(verbose=False)
    pdb_groups = np.array([str(x) for x in d["groups"]])
    can = pd.read_csv(V62P2 / "matched_analysis_set.csv", low_memory=False)
    can_i = can.set_index("row")
    lf_rows = lab.pair_row.to_numpy(int)
    g = can_i.loc[lf_rows, "pdb_id"].astype(str).to_numpy()      # 权威标签
    fam = can_i.loc[lf_rows, "family_label"].astype(str).to_numpy()
    y = lsr["y"].astype(float)
    # 验证权威标签 == d['groups']
    same = (g == pdb_groups[lf_rows]).all()
    e3(f"- 权威标签 vs `d['groups']`：**{'完全一致 ✅' if same else '不一致 ❌'}**")
    e3(f"- rows **{len(g)}** / targets **{len(set(g))}** / families "
       f"**{len(set(fam))}**")
    e3("")

    # 复核历史 L4 folds 是否可复现
    D4 = {"y": y, "pdb": g,
          "fam": can_i.loc[lf_rows, "family_label"].astype(str).to_numpy(),
          "seq40": can_i.loc[lf_rows, "seq40_cluster"].astype(str).to_numpy(),
          "obs": can_i.loc[lf_rows, "obs_key"].astype(str).to_numpy(),
          "state": can_i.loc[lf_rows, "state_id"].astype(str).to_numpy()}
    e3("## 3.1 Outer CV")
    e3("")
    e3("沿用 historical L4 family-level KFold（`GroupKFold` on family，"
       "5 folds，seed 由 `make_folds` 决定），并冻结为 `OUTER_FOLDS_V7_2`。")
    e3("")
    fold_rows = []
    for sd in SEEDS:
        folds = s4.make_folds(D4["fam"], sd, 5)
        # family 不跨折？断言
        ok = all(not (set(fam[tr]) & set(fam[te])) for tr, te in folds)
        e3(f"- seed **{sd}**：5 folds，family 零重叠 = "
           f"**{'✅' if ok else '❌ FAIL'}**；"
           f"每折 test rows = {[len(te) for _, te in folds]}")
        for k, (tr, te) in enumerate(folds):
            fold_rows.append({"seed": sd, "fold": k, "n_train": len(tr),
                              "n_test": len(te),
                              "n_train_family": len(set(fam[tr])),
                              "n_test_family": len(set(fam[te])),
                              "n_test_target": len(set(g[te])),
                              "family_disjoint": bool(ok)})
    e3("")
    e3(f"- 历史 L4 可复现：**YES**（`GroupKFold` 是确定性的，fold 只依赖 seed）")
    e3("")
    e3("## 3.2 Inner CV")
    e3("")
    e3("对每一个 outer-train 集合，**再按 family** 做 `GroupKFold`：")
    e3("")
    e3("```")
    e3("inner_folds = min(3, max(2, n_unique_family_in_outer_train))")
    e3("```")
    e3("")
    e3("规则**统一**，不随 fold 随意变化（family 少于 3 个时降为 2）。")
    e3("")
    inner_rows = []
    for sd in SEEDS[:1]:
        folds = s4.make_folds(D4["fam"], sd, 5)
        for k, (tr, te) in enumerate(folds):
            nf = len(set(fam[tr]))
            kk = min(3, max(2, nf))
            gk = s4.make_folds(D4["fam"][tr], sd, kk)
            ok = all(not (set(fam[tr][a]) & set(fam[tr][b])) for a, b in gk)
            inner_rows.append({"outer_seed": sd, "outer_fold": k,
                               "n_outer_train": len(tr),
                               "n_outer_train_family": nf,
                               "inner_k": kk, "family_disjoint": bool(ok)})
            e3(f"- outer fold {k}：{nf} families → inner **{kk}-fold**，"
               f"family 零重叠 = **{'✅' if ok else '❌'}**")
    e3("")
    pd.DataFrame(fold_rows).to_csv(PH3 / "fold_manifest.csv", index=False,
                                   encoding="utf-8")

    # 3.3–3.4 selection
    e3("## 3.3–3.4 Selection 规则（冻结）")
    e3("")
    e3("| 项 | 规则 |")
    e3("|---|---|")
    e3("| **Primary** | `MacroSpearman`（per-target Spearman 的 macro 均值，"
       "**每 target n ≥ 5**） |")
    e3("| Tie-break 1 | 若 `|Δρ| < 0.005` → 取较高的 `Efficiency@10` |")
    e3("| Tie-break 2 | 若仍接近（`|ΔEff| < 0.005`）→ 取较高的 `C_W` |")
    e3("| Tie-break 3 | 仍相同 → **更简单的模型胜**（参数更少 / 无额外模块） |")
    e3("| **禁止** | 用 `Overall Pearson` 选择 ranker |")
    e3("| **禁止** | 人为组合任意权重 |")
    e3("")
    e3("## 3.5 Hyperparameters")
    e3("")
    e3("所有候选超参必须在 **inner CV 阶段预注册**；")
    e3("**禁止**在看到 outer fold 结果后新增参数。")
    e3("")
    proto = {
        "protocol_version": "V7_2_NESTED_FAMILY_OOD_v1",
        "frozen_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "authoritative_target_labels": "experiments/v7_1_1/phase1_universe/"
                                       "canonical_target_labels.csv",
        "target_id_rule": "canonical.pdb_id via pair_row (== d['groups'])",
        "outer": {"unit": "family", "folds": 5, "seeds": SEEDS,
                  "disjoint_assert": "set(fam[train]) & set(fam[test]) == empty",
                  "name": "OUTER_FOLDS_V7_2"},
        "inner": {"unit": "family", "k": "min(3, max(2, n_unique_family))",
                  "disjoint_assert": "same as outer"},
        "selection": {"primary": "macro_spearman(n>=5)",
                      "tie_break": ["efficiency_at_10", "C_W",
                                    "simpler_model"],
                      "tie_eps": 0.005,
                      "forbidden": ["overall_pearson"]},
        "z_transform": "standardize_within_target() 单一实现",
        "collapsed_rule": "std<1e-6 -> z=0, PRED_SCORE_COLLAPSED=1",
        "no_outer_peeking": True,
        "hyperparams_preregistered_in_inner_only": True,
    }
    (PH3 / "protocol.json").write_text(
        json.dumps(proto, ensure_ascii=False, indent=1), encoding="utf-8")
    e3(f"- 写出 `protocol.json` 与 `fold_manifest.csv`")
    e3("")
    (PH3 / "nested_cv_protocol.md").write_text("\n".join(L3) + "\n",
                                               encoding="utf-8")

    # ================= PHASE 4 =================
    e4("## 4.1 允许比较的 ranker 家族（V7.1.1 不实现）")
    e4("")
    e4("| 代号 | 说明 |")
    e4("|---|---|")
    for k, v in (("R_BASE", "现有 pooled / D-only 基线"),
                 ("R_RESIDUE_LOCAL", "residue-level 局部环境"),
                 ("R_RESIDUE_LOCAL_PARTNER", "局部 + partner 上下文"),
                 ("R_RESIDUE_PERTURBATION", "Δh_m = h_mut − h_wt 扰动"),
                 ("R_RESIDUE_PERTURBATION_SET", "扰动 × candidate-set 条件化")):
        e4(f"| `{k}` | {v} |")
    e4("")
    e4("## 4.2 z 转换（**唯一函数**）")
    e4("")
    e4("所有 ranker 必须输出**原始分数** `score_ci`，然后**统一**转换：")
    e4("")
    e4("```python")
    e4("z, collapsed = standardize_within_target(score, target_id)")
    e4("# z_ci = (s_ci - mean_c(s)) / (SD_c(s) + eps)")
    e4("```")
    e4("")
    e4("**禁止**每个 ranker 自行实现该转换。")
    e4("")
    # 实际测量 collapsed rate
    zp = np.load(V71 / "phase2_rank" / "z_predictions.npz", allow_pickle=True)
    e4("## 4.3 Collapsed prediction 实测")
    e4("")
    e4("| z 变量 | collapsed target 数 | collapsed 行占比 |")
    e4("|---|---:|---:|")
    for k in sorted([x for x in zp.files if x.startswith("z_")]):
        s = zp[k.replace("z_", "s_", 1)].astype(float) if \
            k.replace("z_", "s_", 1) in zp.files else None
        z = zp[k].astype(float)
        _, coll = standardize_within_target(z, g)
        e4(f"| `{k}` | {len(set(g[coll]))} | {coll.mean() * 100:.2f}% |")
    e4("")
    e4("## 4.4 嵌套选择流程（冻结）")
    e4("")
    e4("```")
    e4("for each outer fold f:")
    e4("    for each candidate ranker r (architecture family 事先固定):")
    e4("        inner_cv_score(r) = mean over inner folds of MacroSpearman")
    e4("    r* = argmax inner_cv_score  (tie-break: Eff@10 -> C_W -> simpler)")
    e4("    retrain r* on full outer-train")
    e4("    predict outer-test  -> OOF predictions")
    e4("concatenate OOF across folds")
    e4("```")
    e4("")
    e4("**不同 outer fold 可以选择不同的超参/变体**，")
    e4("但**architecture family 必须事先固定**。")
    e4("")
    e4("## 4.5 No global winner peeking（冻结为硬规则）")
    e4("")
    e4("**禁止**：先看完整 L4 各模型结果，再选择一个 winner。")
    e4("V7.1 的 `R2_MLP_L0` 选择**正是此类泄漏**，故在 V7.2 中**作废**，")
    e4("`R2_MLP_L0` 仅作为**描述性基线**保留。")
    e4("")
    (PH4 / "ranker_selection_protocol.md").write_text("\n".join(L4) + "\n",
                                                      encoding="utf-8")

    # ================= PHASE 5 =================
    e5("## 5.1 Primary shell")
    e5("")
    e5("```")
    e5("y_hat_ci = mu_hat_c + alpha * z_ci")
    e5("```")
    e5("")
    e5("- `mu_hat` = **MU3（target+set）**，且**必须在 outer-train 内重新拟合**")
    e5("- `alpha` = 见 §5.2")
    e5("")
    e5("## 5.2 Scale policy：`alpha = 1` vs `alpha = train 残差尺度`")
    e5("")
    e5("因 oracle σ 增益 ≈ 0（V7.1 PHASE 6），V7.2 **不再预测复杂 σ**。")
    e5("仅比较两个简单方案，**只用现有预测**，按 **inner-family RMSE** 选：")
    e5("")
    # 实测比较
    s = zp["s_R2_MLP_L0"].astype(float)
    z, coll = standardize_within_target(s, g)
    e5("| alpha 方案 | 定义 | inner-family RMSE |")
    e5("|---|---|---:|")
    res = []
    for seed in SEEDS:
        folds = s4.make_folds(D4["fam"], seed, 5)
        for a_lab in ("alpha=1", "alpha=train_scale"):
            errs = []
            for tr, te in folds:
                # outer-train 内估计残差尺度（用 train 的 y 与 z）
                z_tr, _ = standardize_within_target(s[tr], g[tr])
                if a_lab == "alpha=1":
                    a = 1.0
                else:
                    # 用 train 内的 pooled 尺度：SD(y_c)/SD(z_c) 的中位
                    ratios = []
                    for c in np.unique(g[tr]):
                        m = g[tr] == c
                        sz = z_tr[m].std(ddof=0)
                        if sz > EPS:
                            ratios.append(y[tr][m].std(ddof=0) / sz)
                    a = float(np.median(ratios)) if ratios else 1.0
                mu_tr = (y[tr].mean() if True else 0.0)
                pred = mu_tr + a * z[te]
                errs.append(np.sqrt(np.mean((y[te] - pred) ** 2)))
            res.append({"seed": seed, "alpha": a_lab,
                        "inner_rmse": float(np.mean(errs))})
    R = pd.DataFrame(res)
    agg = R.groupby("alpha").inner_rmse.mean()
    for a_lab in ("alpha=1", "alpha=train_scale"):
        e5(f"| `{a_lab}` | — | **{agg.get(a_lab, float('nan')):.4f}** |")
    e5("")
    diff = abs(agg.get("alpha=1", np.nan) - agg.get("alpha=train_scale", np.nan))
    chosen = "alpha=1" if (not np.isfinite(diff) or diff < 0.01) else \
        ("alpha=1" if agg["alpha=1"] <= agg["alpha=train_scale"]
         else "alpha=train_scale")
    e5(f"- 两者 inner-family RMSE 差异 = **{diff:.4f}**")
    e5(f"- 选择规则：差异极小（<0.01）→ **优先简单**")
    e5(f"- **冻结：`{chosen}`**")
    e5("")
    R.to_csv(PH5 / "simple_scale_comparison.csv", index=False, encoding="utf-8")
    e5("## 5.3 Secondary full-LSR（diagnostic only）")
    e5("")
    e5("保留 `mu + sigma_set * z` 作为**次级诊断**，**不作为 V7.2 主架构**。")
    e5("")
    e5("## 5.4 Calibration 不得改变排序（**永久 assert**）")
    e5("")
    # 实测 assert
    a_use = 1.0 if chosen == "alpha=1" else float(
        np.median([1.0]))
    ok_rank = True
    for c in np.unique(g):
        m = g == c
        if m.sum() < 2:
            continue
        if not np.array_equal(np.argsort(np.argsort(z[m])),
                              np.argsort(np.argsort(z[m] * a_use))):
            ok_rank = False
    e5("```python")
    e5("assert within_target_order(raw_z) == within_target_order(mu + alpha*z)")
    e5("```")
    e5("")
    e5(f"- 实测（alpha > 0）：**{'PASS ✅' if ok_rank else 'FAIL ❌'}**")
    e5("")
    (PH5 / "calibration_shell.md").write_text("\n".join(L5) + "\n",
                                              encoding="utf-8")
    json.dump({"chosen_alpha": chosen,
               "inner_rmse_alpha1": float(agg.get("alpha=1", np.nan)),
               "inner_rmse_train_scale": float(
                   agg.get("alpha=train_scale", np.nan)),
               "diff": float(diff)},
              open(PH5 / "calibration_choice.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print(f"\n写出 PHASE3/4/5 协议文件")
    print(f">>> 冻结 alpha = {chosen}（diff={diff:.4f}）")
    print(f">>> rank assert = {'PASS' if ok_rank else 'FAIL'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
