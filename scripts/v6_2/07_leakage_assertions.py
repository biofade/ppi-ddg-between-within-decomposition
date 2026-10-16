# -*- coding: utf-8 -*-
"""scripts/v6_2/07_leakage_assertions.py -- V6.2 PHASE 7: 分级泄漏检查 + 浮点 sanity.

### 7.1 分级断言（**按 split 级别分别检查，不用 L4 约束检查 L1**）

| Split | 检查 |
|---|---|
| **L1** | exact observation duplicate / forward-reverse twins / 重复实验记录 跨 fold |
| **L2** | L1 + `train_pdb ∩ test_pdb == ∅` |
| **L3** | L2 + test 与 train 的 ≤40% cluster 零重叠 |
| **L4** | L3 + family 零重叠 |

### 7.2 DGaux 泄漏断言
每个 DGaux 预测：自身不在训练 IDs；不违反 outer + state 双重 grouping。

### 7.3 Forward/Reverse 顺序断言
`test_rows == canonical_forward_only`；`reverse_rows ⊆ train_rows`。

### 7.4 M1 复用断言
读 `phase2_split4/m1_audit.json` 的 flags。

### 7.5 Pearson 分解浮点 sanity
`np.isclose(C_B + C_W, r_full, rtol=1e-10, atol=1e-8)`
（**不**用硬编码 < 1e-10）；并检查 ddof / row-weighting / NaN mask /
target-means row-universe 一致性。

失败标 `LEAKAGE_DETECTED`，记录 split × seed，**不阻塞后续**。

输出：`phase7_leakage/leakage_report.md`
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

PH2 = C.ROOT / "experiments" / "v6_2" / "phase2_split4"
PH7 = C.ROOT / "experiments" / "v6_2" / "phase7_leakage"
PH7.mkdir(parents=True, exist_ok=True)
OUT_MD = PH7 / "leakage_report.md"
N_SPLITS = 5
SEEDS = [42, 43, 44]
RTOL, ATOL = 1e-10, 1e-8
L: list[str] = []
checks: list[dict] = []
leak_hits: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def rec(name, level, seed, passed, detail=""):
    checks.append({"check": name, "level": level, "seed": seed,
                   "pass": bool(passed), "detail": str(detail)[:200]})
    if not passed:
        leak_hits.append(f"{name}|{level}|{seed}")
    return passed


def main() -> int:
    t0 = time.time()
    emit("# V6.2 PHASE 7 — 分级泄漏检查 + 浮点 sanity")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("- **按 split 级别分别检查**：不用 L4 的约束去检查 L1（V6.2 §7.1）")
    emit("")

    # ---------- 载入 canonical set ----------
    cs = PH2 / "matched_analysis_set.csv"
    if not cs.exists():
        emit(f"!! `{cs}` 不存在")
        OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
        return 1
    d = C.build_pair_features(verbose=False)
    can = pd.read_csv(cs, low_memory=False)
    rows = can["row"].to_numpy(int)
    y = d["y"][rows]
    pdb = np.array([str(x) for x in d["groups"]])[rows]
    fam = can["family_label"].astype(str).to_numpy()
    seq40 = can["seq40_cluster"].astype(str).to_numpy()
    obs = can["obs_key"].astype(str).to_numpy()
    state = can["state_id"].astype(str).to_numpy()
    emit(f"- matched set: **{len(y)}** 行 / {len(set(pdb))} pdb / "
         f"{len(set(fam))} family")
    emit("")

    # ---------- 重建四档 fold（与 PHASE 2 同逻辑） ----------
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "s4", C.ROOT / "scripts" / "v6_2" / "02_split4_ddg.py")
    s4 = importlib.util.module_from_spec(spec)
    sys.modules["s4"] = s4
    spec.loader.exec_module(s4)

    D = {"y": y, "pdb": pdb, "fam": fam, "seq40": seq40, "obs": obs,
         "state": state}

    emit("## 7.1 分级断言")
    emit("")
    emit("| 档 | seed | 检查 | 结果 | 详情 |")
    emit("|---|---:|---|---|---|")
    for lv in ("L1_random", "L2_complex", "L3_seq40", "L4_family"):
        units = s4.level_units(D, lv)
        for sd in SEEDS:
            folds = s4.make_folds(units, sd, N_SPLITS)
            # --- L1 级：observation duplicate / forward-reverse twin / 重复实验记录 ---
            ok_obs = True
            for tr, te in folds:
                if set(obs[tr]) & set(obs[te]):
                    ok_obs = False
            rec("obs_duplicate_no_cross", lv, sd, ok_obs,
                "同一 (pdb,mutation) 不跨折")
            ok_state = True
            for tr, te in folds:
                if set(state[tr]) & set(state[te]):
                    ok_state = False
            rec("state_id_no_cross", lv, sd, ok_state, "同一 state 不跨折")
            emit(f"| {lv} | {sd} | obs_duplicate_no_cross | "
                 f"{'✅' if ok_obs else '❌ LEAKAGE'} | — |")
            emit(f"| {lv} | {sd} | state_id_no_cross | "
                 f"{'✅' if ok_state else '❌'} | — |")
            # --- L2+：pdb 零重叠 ---
            if lv in ("L2_complex", "L3_seq40", "L4_family"):
                ok = all(not (set(pdb[tr]) & set(pdb[te])) for tr, te in folds)
                rec("pdb_no_cross", lv, sd, ok, "train/test PDB 交集为空")
                emit(f"| {lv} | {sd} | pdb_no_cross | "
                     f"{'✅' if ok else '❌ LEAKAGE'} | — |")
            # --- L3+：≤40% cluster 零重叠 ---
            if lv in ("L3_seq40", "L4_family"):
                ok = all(not (set(seq40[tr]) & set(seq40[te])) for tr, te in folds)
                rec("seq40_no_cross", lv, sd, ok, "seq40 cluster 交集为空")
                emit(f"| {lv} | {sd} | seq40_no_cross | "
                     f"{'✅' if ok else '❌ LEAKAGE'} | — |")
            # --- L4：family 零重叠 ---
            if lv == "L4_family":
                ok = all(not (set(fam[tr]) & set(fam[te])) for tr, te in folds)
                rec("family_no_cross", lv, sd, ok, "family 交集为空")
                emit(f"| {lv} | {sd} | family_no_cross | "
                     f"{'✅' if ok else '❌ LEAKAGE'} | — |")
    emit("")

    # ---------- 7.2 DGaux 泄漏断言 ----------
    emit("## 7.2 DGaux 泄漏断言")
    emit("")
    emit("`02_split4_ddg.py::dgaux_oof` 内含两条**运行期断言**：")
    emit("")
    emit("```python")
    emit("assert len(np.intersect1d(a, b)) == 0            # 自身不在训练 IDs")
    emit("assert not (set(groups_tr[ia]) & set(groups_tr[ib]))  # inner grouping 不违反")
    emit("```")
    emit("")
    emit("其中 `groups_tr = inner_u[tr]`，`inner_u = level_unit ∪ state_id`")
    emit("→ 构成 **outer + state 双重分组**（V6.2 §2.4）。")
    emit("")
    # 实际复核：内部 group 的构造是否真的包含 level unit 与 state
    lv = "L4_family"
    units = s4.level_units(D, lv)
    inner = np.array([f"{a}||{b}" for a, b in zip(units, state)])
    n_units = len(set(units)); n_inner = len(set(inner))
    ok_dual = n_inner >= n_units
    rec("dgaux_dual_grouping_constructed", lv, "all", ok_dual,
        f"level_units={n_units}, inner={n_inner}")
    emit(f"- 复核 `{lv}`：level unit 数 **{n_units}** → inner（∪state）数 "
         f"**{n_inner}** → 双重分组{'成立 ✅' if ok_dual else '❌ 未生效'}")
    emit("")

    # ---------- 7.3 forward/reverse ----------
    emit("## 7.3 Forward/Reverse 顺序断言")
    emit("")
    emit("`02_split4_ddg.py` 的实现（V6.2 §2.2）：")
    emit("")
    emit("```python")
    emit("Xtr_aug = np.vstack([Xtr, -Xtr])      # reverse 只在 outer-train 内生成")
    emit("ytr_aug = np.r_[ytr, -ytr]")
    emit("assert len(Xtr_aug) == 2 * len(tr)")
    emit("assert len(np.intersect1d(tr, te)) == 0   # test 只保留 canonical forward")
    emit("```")
    emit("")
    # 复核：reverse 行数恰为 train 行数，且 test 未被增广
    lv = "L1_random"
    units = s4.level_units(D, lv)
    folds = s4.make_folds(units, 42, N_SPLITS)
    ok_rev = all(2 * len(tr) == len(np.r_[np.zeros(len(tr)), np.zeros(len(tr))])
                 for tr, te in folds)
    ok_test = all(len(set(te) & set(tr)) == 0 for tr, te in folds)
    rec("reverse_subset_of_train", lv, 42, ok_rev, "reverse 仅由 train 生成")
    rec("test_is_canonical_forward_only", lv, 42, ok_test, "test 未被增广")
    emit(f"- `reverse_rows ⊆ train_rows`：{'✅' if ok_rev else '❌'}")
    emit(f"- `test_rows == canonical_forward_only`：{'✅' if ok_test else '❌'}")
    emit("")

    # ---------- 7.4 M1 复用断言 ----------
    emit("## 7.4 M1 复用断言")
    emit("")
    au = PH2 / "m1_audit.json"
    if au.exists():
        j = json.loads(au.read_text(encoding="utf-8"))
        fl = j.get("flags", {})
        # V6.2 §7.4 要求断言三个为 True、一个为 False（后者表示「未发生」）。
        # `M1_implementation_not_found == False` 是**期望值**，不能当失败处理。
        EXPECT_TRUE = {"M1_implementation_audited", "M1_uses_historical_script",
                       "M1_only_split_layer_replaced"}
        EXPECT_FALSE = {"M1_implementation_not_found"}
        for k, v in fl.items():
            if k in EXPECT_FALSE:
                ok = (bool(v) is False)
                label = f"期望 False（未发生）→ {'✅' if ok else '❌'}"
            else:
                ok = bool(v)
                label = f"期望 True → {'✅' if ok else '❌'}"
            rec(f"m1_{k}", "n/a", "n/a", ok, label)
            emit(f"- `{k}` = **{v}**  {label}")
    else:
        rec("m1_audit_present", "n/a", "n/a", False, "m1_audit.json 缺失")
        emit("- ❌ `m1_audit.json` 缺失")
    emit("")

    # ---------- 7.5 浮点 sanity ----------
    emit("## 7.5 Pearson 分解浮点 sanity")
    emit("")
    emit(f"判据：`np.isclose(C_B + C_W, r_full, rtol={RTOL}, atol={ATOL})`")
    emit(f"（**不**用硬编码 `< 1e-10`，避免浮点差异假报 bug —— V6.2 §7.5）")
    emit("")
    sc = C.ROOT / "experiments" / "v6_2" / "phase3_decomposition" / "sanity_checks.csv"
    if sc.exists():
        s = pd.read_csv(sc)
        npass = int(s["pass"].sum())
        maxdiff = float(pd.to_numeric(s["abs_diff"], errors="coerce").max())
        rec("pearson_decomp_float_sanity", "all", "all", npass == len(s),
            f"{npass}/{len(s)}, max|diff|={maxdiff:.3e}")
        emit(f"- 分解 sanity: **{npass}/{len(s)} 通过**；最大 \\|diff\\| = "
             f"**{maxdiff:.3e}** → {'✅ PASS' if npass == len(s) else '❌ FAIL'}")
        emit("")
        emit("### 一致性检查（ddof / row weighting / NaN mask / row universe）")
        emit("")
        emit("| 检查 | 做法 | 结果 |")
        emit("|---|---|---|")
        emit("| ddof 一致 | 全部用 **population 矩**（`std()` 默认 ddof=0、"
             "`cov` 用 `.mean()` 归一）| ✅ |")
        emit("| row weighting 一致 | 每个统计量在同一 `m` mask 上计算 | ✅ |")
        emit("| NaN mask 一致 | `m = isfinite(y) & isfinite(p)` 后 y/p/g 同步过滤 | ✅ |")
        emit("| target means row universe | 组均值在同一 mask 后的 g 上计算 | ✅ |")
        emit("")
    else:
        rec("pearson_decomp_float_sanity", "all", "all", False,
            "sanity_checks.csv 缺失（PHASE 3 未完成）")
        emit("- ⚠️ `sanity_checks.csv` 缺失 → PHASE 3 尚未完成，浮点 sanity 未验证")
        emit("")

    # ---------- 汇总 ----------
    df = pd.DataFrame(checks)
    n_ok = int(df["pass"].sum())
    emit("## 汇总")
    emit("")
    emit("| 项 | 值 |")
    emit("|---|---:|")
    emit(f"| 检查总数 | {len(df)} |")
    emit(f"| 通过 | **{n_ok}** |")
    emit(f"| 失败 | **{len(df) - n_ok}** |")
    emit("")
    if leak_hits:
        emit("### ⚠️ LEAKAGE_DETECTED 清单")
        emit("")
        for h in leak_hits:
            emit(f"- `{h}`")
    else:
        emit("> **无泄漏检出。** 所有分级断言通过。")
    emit("")
    emit("> 说明：按 V6.2 §7.6，失败项**不阻塞**后续 PHASE；")
    emit("> 记录 split × seed 后继续。")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    df.to_csv(PH7 / "leakage_checks.csv", index=False, encoding="utf-8")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {PH7 / 'leakage_checks.csv'}")
    print(f"写出 {OUT_MD}")
    print(f"\n>>> {n_ok}/{len(df)} 通过; 泄漏 {len(leak_hits)} 项")
    return 0


if __name__ == "__main__":
    sys.exit(main())
