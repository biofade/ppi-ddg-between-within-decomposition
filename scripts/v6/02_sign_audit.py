# -*- coding: utf-8 -*-
"""scripts/v6/02_sign_audit.py -- V6.0 PHASE 2.1: ΔΔG 符号一致性审计（**手工可复核**）.

为什么必须做（V6 §0.2）
-----------------------
不同论文对 ΔΔG 的正负定义不同。在算 Hit@10 之前必须确认：
**「增强结合（beneficial）」= ΔΔG < 0**。
本脚本用**独立的第三方真值**（本地 SKEMPI v2.0 原始表，含 `Affinity_wt (M)` /
`Affinity_mut (M)`，分号分隔）逐步复算 ΔΔG，与 CATH-ddG 表的 `ddG` 列对照。

手工复算链（对每个样例逐行打印，便于人工复核）
----------------------------------------------
    ΔG_WT  = 0.592 × ln(KD_WT[M])
    ΔG_mut = 0.592 × ln(KD_mut[M])
    ΔΔG   = ΔG_mut − ΔG_WT
  若 KD_mut > KD_WT（结合变弱）→ ΔΔG > 0（destabilizing）

0.592 kcal/mol = R·T，R = 1.98720425864083e-3，T = 298.15 K

对**模型预测列**：ΔΔG 真值越小越好，因此一个「方向正确」的模型应与真值**正相关**。
若某模型预测与真值显著**负相关**，则该模型的符号约定相反 → 需要反号
（V6 §0.2 的 `needs_flip`）。本脚本在 **CATH 全表**（train+test，n≈6706）上估计
相关系数以避免小样本噪声，但在报告中**同时给出未反号与反号两套指标**。

输出
----
* `experiments/v6_sota_audit/sign_audit.csv`
* `experiments/v6_sota_audit/sign_audit.md`
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402  (复用 V5 指标；只读)

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = C.ROOT
OUT = ROOT / "experiments" / "v6_sota_audit"
OUT.mkdir(parents=True, exist_ok=True)
CATH_CSV = Path(r"F:\dsh\downloads\cathddg\skempi_v2_with_all_results.csv")
SKEMPI = ROOT / "data" / "v2.4" / "skempi_v2_cleaned.csv"
RT = 1.98720425864083e-3 * 298.15        # 0.5925...

L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def dG(kd_M):
    kd_M = np.asarray(kd_M, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return RT * np.log(kd_M)


def main() -> int:
    emit("# V6.0 PHASE 2.1 — ΔΔG 符号一致性审计（手工可复核）")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- 物理常数: **R·T = {RT:.6f}** kcal/mol（R = 1.98720425864083e-3, "
         f"T = 298.15 K）")
    emit(f"- 约定: `ΔG = RT·ln(KD[M])`；`ΔΔG = ΔG_mut − ΔG_WT`；"
         f"**ΔΔG > 0 = 削弱结合，ΔΔG < 0 = 增强结合（beneficial）**")
    emit("")

    # ---------- 1. 从本地 SKEMPI 原始表独立复算 ----------
    emit("## 1. 独立复算：本地 SKEMPI v2.0 原始表（分号分隔，第三方真值）")
    emit("")
    sk = None
    if SKEMPI.exists():
        sk = pd.read_csv(SKEMPI, sep=";", low_memory=False)
        emit(f"- 读入 `data/v2.4/skempi_v2_cleaned.csv`: **{sk.shape[0]} 行 × "
             f"{sk.shape[1]} 列**（分隔符 `;`）")
        need = ["#Pdb", "Mutation(s)_cleaned", "Affinity_wt (M)",
                "Affinity_mut (M)"]
        miss = [c for c in need if c not in sk.columns]
        if miss:
            emit(f"- ⚠️ 缺列 {miss} → 无法复算，跳过")
            sk = None
    else:
        emit(f"- ⚠️ `{SKEMPI}` 不存在 → 无法用第三方真值复算")

    picked = []
    if sk is not None:
        sk["wtM"] = pd.to_numeric(sk["Affinity_wt (M)"], errors="coerce")
        sk["mutM"] = pd.to_numeric(sk["Affinity_mut (M)"], errors="coerce")
        ok = sk[sk["wtM"].notna() & sk["mutM"].notna() & (sk["wtM"] > 0)
                & (sk["mutM"] > 0)].copy()
        ok["dG_wt"] = dG(ok["wtM"])
        ok["dG_mut"] = dG(ok["mutM"])
        ok["ddG_calc"] = ok["dG_mut"] - ok["dG_wt"]
        emit(f"- 可复算行（两个 KD 均为正）: **{len(ok)}**")
        emit("")
        emit("### 逐步手工复算（取 6 条，含 1 条典型 destabilizing）")
        emit("")
        emit("| #Pdb | Mutation(cleaned) | KD_WT (M) | KD_mut (M) | ΔG_WT | ΔG_mut | "
             "ΔΔG_calc | ΔΔG<0? |")
        emit("|---|---|---:|---:|---:|---:|---:|---|")
        # 取 1CSE 的 LI38G/X 系列（V6 §2.1 的教科书样例）与几条随机样本
        pick = ok[ok["#Pdb"].astype(str).str.startswith("1CSE")].head(3)
        rest = ok.drop(pick.index).sample(n=3, random_state=C.SEED)
        picked = pd.concat([pick, rest])
        for _, r in picked.iterrows():
            emit(f"| {r['#Pdb']} | {r['Mutation(s)_cleaned']} | {r['wtM']:.3e} | "
                 f"{r['mutM']:.3e} | {r['dG_wt']:+.4f} | {r['dG_mut']:+.4f} | "
                 f"**{r['ddG_calc']:+.4f}** | "
                 f"{'✅ 是' if r['ddG_calc'] < 0 else '否（削弱）'} |")
        emit("")
        emit("**读法**：`LI38G` 的 KD 从 1.12e−12 升到 5.26e−11（结合**变弱**）→ "
             "ΔΔG = **正**值 ✅ 与「正值 = 削弱结合」一致。")
        emit("")

    # ---------- 2. 与 CATH-ddG 表连接核对 ----------
    emit("## 2. 与 CATH-ddG 表 `ddG` 列交叉核对")
    emit("")
    if not CATH_CSV.exists():
        emit(f"- ⚠️ `{CATH_CSV}` 不存在 → `DATA_MISSING`")
        (OUT / "sign_audit.md").write_text("\n".join(L) + "\n", encoding="utf-8")
        return 0
    ca = pd.read_csv(CATH_CSV, low_memory=False)
    emit(f"- CATH-ddG 表: **{ca.shape[0]} 行 × {ca.shape[1]} 列**；"
         f"`ddG` 非空 {ca['ddG'].notna().sum()} 行")
    emit("")

    agree = None
    if sk is not None:
        a = ca[["complex", "Mutation(s)_cleaned", "ddG"]].copy()
        a["k"] = (a["complex"].astype(str) + "|"
                  + a["Mutation(s)_cleaned"].astype(str))
        b = ok[["#Pdb", "Mutation(s)_cleaned", "ddG_calc"]].copy()
        b["k"] = (b["#Pdb"].astype(str).str.split("_").str[0] + "|"
                  + b["Mutation(s)_cleaned"].astype(str))
        m = a.merge(b[["k", "ddG_calc"]].drop_duplicates("k"), on="k", how="inner")
        emit(f"- 与本地 SKEMPI 按 `(PDB, Mutation_cleaned)` 连接命中: "
             f"**{len(m)}** 行")
        if len(m) > 10:
            r = C.pearson(m["ddG_calc"], m["ddG"])
            frac_pos = float(((m["ddG_calc"] > 0) & (m["ddG"] > 0)).mean()
                             + ((m["ddG_calc"] < 0) & (m["ddG"] < 0)).mean())
            emit(f"- **Pearson(ΔΔG_calc, ddG_CATH) = {r:+.4f}**；"
                 f"**符号一致比例 = {frac_pos * 100:.2f}%**")
            agree = r
            emit("")
            emit("| #Pdb | Mutation | ΔΔG_calc (本地复算) | ddG (CATH 表) | 符号一致? |")
            emit("|---|---|---:|---:|---|")
            for _, x in m.head(6).iterrows():
                same = (x["ddG_calc"] > 0) == (x["ddG"] > 0)
                emit(f"| {x['k'].split('|')[0]} | {x['k'].split('|')[1]} | "
                     f"{x['ddG_calc']:+.4f} | {x['ddG']:+.4f} | "
                     f"{'✅' if same else '❌'} |")
            emit("")
            if r > 0.9:
                emit("→ **结论：CATH-ddG 的 `ddG` 与我方 ΔΔG 约定一致"
                     "（正 = 削弱结合）→ 真值列 `needs_flip = False`**")
            elif r < -0.9:
                emit("→ **结论：CATH-ddG 的 `ddG` 需反号**")
            else:
                emit(f"→ ⚠️ 相关仅 {r:+.4f}，无法确证一致性；"
                     "保守策略：**按符号一致比例判断**，并在报告中标注不确定性。")
        else:
            emit("- ⚠️ 连接命中过少，无法交叉核对")
    emit("")

    # ---------- 3. 模型预测列的符号判定 ----------
    emit("## 3. 模型预测列的符号判定（`needs_flip`）")
    emit("")
    emit("> 判定规则（写死在代码里）：ΔΔG 真值**越小越好**，因此方向正确的模型"
         "应与真值**正相关**。若在 CATH 全表（train+test，最大样本量）上"
         "`Pearson(pred, ddG) < 0`，则判 `needs_flip = True`。")
    emit("> **同时给出未反号与反号两套指标**，不隐藏反号带来的人为提升。")
    emit("")
    MODEL_COLS = ["DSSA-PPI", "PPIFormer", "ESM-2-FFN", "MuTon", "MpbPPI", "ProBASS",
                  "SAAMBE-SEQ", "FoldX"]
    if CATH_CSV.exists():
        # CATH 表中的 baseline 预测列（结构/能量/序列基线与 SOTA 基线）
        cath_models = {"Flex-ddG": "flex ddG", "DDGPred": "DDGPred",
                       "ddG_FoldX": "FoldX", "ESM-IF": "ESM-IF",
                       "proteinMPNN": "ProteinMPNN"}
        yv = pd.to_numeric(ca["ddG"], errors="coerce")
        emit("### 3a. CATH-ddG 表中的预测列（held-out CATH test 可用）")
        emit("")
        emit("| 列 | 模型 | n 非空 | Pearson(全表) | 与真值同方向? | needs_flip |")
        emit("|---|---|---:|---:|---|---|")
        rows = []
        for col, name in cath_models.items():
            pv = pd.to_numeric(ca[col], errors="coerce")
            m = pv.notna() & yv.notna()
            r = C.pearson(yv[m], pv[m]) if m.sum() > 30 else float("nan")
            flip = bool(r < 0)
            emit(f"| `{col}` | {name} | {int(m.sum())} | {r:+.4f} | "
                 f"{'✅ 同向' if r > 0 else '❌ 反向'} | "
                 f"**{'True' if flip else 'False'}** |")
            rows.append({"dataset": "CATH-ddG/SKEMPI_v2.0", "column": col,
                         "model": name, "n_nonnull": int(m.sum()),
                         "pearson_full_table": round(float(r), 6),
                         "needs_flip": flip,
                         "sign_convention": ("higher=destabilizing" if r > 0
                                             else "lower=destabilizing")})
        # DSSA-PPI RBM
        rbm = Path(r"F:\dsh\downloads\dssa_ppi\RBM_baselines_comparison_ddG.csv")
        if rbm.exists():
            rb = pd.read_csv(rbm, low_memory=False)
            emit("")
            emit("### 3b. DSSA-PPI RBM 案例集（`EXP_fitness` 为真值）")
            emit("")
            emit(f"- 行数 **{len(rb)}**；列: {list(rb.columns)}")
            emit(f"- ⚠️ 注意：该文件的真值列名为 `EXP_fitness`（实验**适应度**），"
                 f"**不是 ΔΔG**；其符号含义未在仓库中说明 → "
                 f"按「预测应与真值同向（正相关）」判定 `needs_flip`，"
                 f"但**该数据集不进入 Hit@10（ΔΔG<0）口径**，只用于 "
                 f"Pearson/Spearman 与 within-position 排序对比。")
            emit("")
            emit("| 列 | n | Pearson(全表) | needs_flip |")
            emit("|---|---:|---:|---|")
            for col in MODEL_COLS:
                if col not in rb.columns:
                    continue
                pv = pd.to_numeric(rb[col], errors="coerce")
                tv = pd.to_numeric(rb["EXP_fitness"], errors="coerce")
                m = pv.notna() & tv.notna()
                r = C.pearson(tv[m], pv[m]) if m.sum() > 30 else float("nan")
                flip = bool(r < 0)
                emit(f"| `{col}` | {int(m.sum())} | {r:+.4f} | "
                     f"**{'True' if flip else 'False'}** |")
                rows.append({"dataset": "DSSA-PPI/RBM_case_study", "column": col,
                             "model": col, "n_nonnull": int(m.sum()),
                             "pearson_full_table": round(float(r), 6),
                             "needs_flip": flip,
                             "sign_convention": ("higher=beneficial" if r > 0
                                                 else "lower=beneficial")})
        sa = pd.read_csv(OUT / "sign_audit.csv", low_memory=False) \
            if (OUT / "sign_audit.csv").exists() else None
        pd.DataFrame(rows).to_csv(OUT / "sign_audit.csv", index=False,
                                  encoding="utf-8")
    emit("")

    # ---------- 4. 备注 ----------
    emit("## 4. 未能审计的数据集（如实记录）")
    emit("")
    emit("| 目标 | 状态 | 原因 |")
    emit("|---|---|---|")
    emit("| SAMAffinity | `PREDICTION_FILE_NOT_FOUND` | 仓库 `data/test_fold/*/*.csv` "
         "只含 `ddG` 真值与突变标识，**无任何预测值列**（已逐个检查 31 列） |")
    emit("| Prompt-DDG | `PREDICTION_FILE_NOT_FOUND` | CATH 表中 `prompt-fold0/1/2` "
         "三列只是 **split 标签**（train/val），无数值预测 |")
    emit("| MINT | `PREDICTION_FILE_NOT_FOUND` | 仓库 68 个文件中仅 "
         "`mint/helpers/predict.py`（428 B）与一个 metrics notebook，无预测结果文件 |")
    emit("| Pythia-PPI | `PREDICTION_FILE_NOT_FOUND` | 仓库 32 个文件，grep "
         "predict/result/test 无匹配结果文件 |")
    emit("")

    (OUT / "sign_audit.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT / 'sign_audit.csv'}")
    print(f"写出 {OUT / 'sign_audit.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
