# -*- coding: utf-8 -*-
"""scripts/v8/39_external_multimodel.py -- U_EXT_CANONICAL：**6 个外部模型**的匹配宇宙。

## 资源（§9 reuse-first 的最佳情形：**零推理**）

CATH-ddG 仓库附带 `data/SKEMPI2/skempi_v2.csv`（1,291,260 B / 7085 行 / 22 列），
其中含**逐突变预测**：
    `Flex-ddG`, `DDGPred`, `ddG_FoldX`, `prompt-fold0/1/2`, `ESM-IF`, `proteinMPNN`
以及真值 `ddG` 与**官方划分列** `cath_fold`（5893 train / 813 val）。

⇒ 6–8 个外部方法**无需任何推理、无需 GPU** 即可进入匹配记分板。

## 映射（§20/§21）

用 `(PDB4, Mutation(s)_cleaned)` 与 canonical 的
`(pdb4, ann_chain+ann_pos1+mut_aa)` 精确匹配（HERMES 审计已验证该键格式）。

## §50 污染审计（每个方法分别判定）

`splits` 列给出 `Label`(forward/reverse)、`PPIformer`、`prompt-fold*`、`cath_fold`；
**CATH-ddG 自身** 88% 行属其训练集 ⇒ 其**自身预测未发布**（故不在本表内），
但基线方法的 provenance 需逐列说明。
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
sys.path.insert(0, str(R / "scripts" / "v8"))
import v8metrics as M  # noqa: E402

OUT = R / "experiments" / "v8"
EXT = OUT / "external"
SRC = Path("F:/dsh/downloads/cathddg/skempi_v2.csv")
MIN_N = 5

EXT_MODELS = ["Flex-ddG", "DDGPred", "ddG_FoldX", "prompt-fold0", "prompt-fold1",
              "prompt-fold2", "ESM-IF", "proteinMPNN"]

# ---------- canonical ----------
lab = pd.read_csv(R / "experiments" / "v7_1_1" / "phase1_universe" /
                  "canonical_target_labels.csv", low_memory=False)
ms = pd.read_csv(R / "data" / "v5_0" / "mutation_sites.csv", low_memory=False)
ms["pair_row"] = pd.to_numeric(ms.pair_row, errors="coerce")
z = np.load(R / "experiments" / "v7_1" / "phase2_rank" / "z_predictions.npz",
            allow_pickle=True)
Y = np.asarray(z["y"], float)

def s(x):
    return str(x).strip().upper()

ms["tok"] = ms.wt_aa.map(s) + ms.ann_pos1.astype("Int64").astype(str) + ms.mut_aa.map(s)
agg = (ms.sort_values("pair_row").groupby("pair_row", as_index=False)
         .agg(pdb4=("pdb_id", lambda v: str(v.iloc[0]).upper().strip()[:4]),
              n_sites=("tok", "size"),
              mut_str=("tok", lambda v: ",".join(sorted(set(v))))))
C = lab.merge(agg, on="pair_row", how="left")
assert len(C) == len(lab), f"行膨胀 {len(C)} != {len(lab)}"
C["y_true"] = Y
C["key"] = C.pdb4 + "|" + C.mut_str
print(f"canonical = {len(C)} 行 (严格 1:1)")

# ---------- 外部 ----------
E = pd.read_csv(SRC, low_memory=False)
print(f"外部源 = {E.shape}")
E["pdb4"] = E.complex.astype(str).str.upper().str.strip().str[:4]

# ★ 格式归一：CATH-ddG 的 `Mutation(s)_cleaned` 形如 `LI38G`
#   = <可选链字母><wt><pos><mut>，而 canonical 的 tok = <wt><ann_pos1><mut>。
#   故统一抽取 (位置, 突变氨基酸) 作为匹配键，忽略链字母与 wt
#   （wt 由 canonical 的 mutation_sites 提供，可用于交叉校验）。
def parse_ext(txt: str):
    t = str(txt).strip().upper().replace(" ", "")
    m = re.match(r"^([A-Z]?)([A-Z])(\d+)([A-Z])$", t)
    if not m:
        return (None, None)
    return (int(m.group(3)), m.group(4))

E[["pos", "mutaa"]] = E["Mutation(s)_cleaned"].map(
    lambda x: pd.Series(parse_ext(x))).apply(pd.Series) if False else \
    pd.DataFrame([parse_ext(x) for x in E["Mutation(s)_cleaned"]],
                 index=E.index)
E["key2"] = E.pdb4 + "|" + E.pos.astype("Int64").astype(str) + "|" + E.mutaa
print(f"外部 key2 样例: {E.key2.head(3).tolist()}")
print(f"外部唯一 key2 = {E.key2.nunique()} / {len(E)}")

# canonical 侧的 key2
_pos = pd.to_numeric(ms.groupby("pair_row").ann_pos1.first()
                    .reindex(C.pair_row.to_numpy()).to_numpy(), errors="coerce")
_mut = ms.groupby("pair_row").mut_aa.first().reindex(C.pair_row.to_numpy()).to_numpy()
C["key2"] = (C.pdb4.astype(str) + "|"
             + pd.Series(_pos).astype("float").astype("Int64").astype(str).to_numpy()
             + "|" + pd.Series(_mut).astype(str).to_numpy())
# 仅单点行参与 key2 匹配（多点的 pos 未定义）
single = C.n_sites.eq(1)
print(f"canonical 单点行 = {int(single.sum())} / {len(C)}")

# 多行同 key -> 取均值；同时保留划分信息
keep = ["key2", "ddG", "cath_fold", "num_muts"] + EXT_MODELS
E2 = E[[c for c in keep if c in E.columns]].copy()
for c in ["ddG"] + EXT_MODELS:
    if c in E2.columns:
        E2[c] = pd.to_numeric(E2[c], errors="coerce")
num = [c for c in EXT_MODELS if c in E2.columns]
aggd = {c: "mean" for c in num}
aggd["ddG"] = "mean"
if "num_muts" in E2.columns:
    aggd["num_muts"] = "mean"
g = E2.groupby("key2", as_index=False).agg(aggd)
if "cath_fold" in E2.columns:
    fold = E2.groupby("key2").cath_fold.first().rename("cath_fold").reset_index()
    g = g.merge(fold, on="key2", how="left")

MM = C.merge(g, on="key2", how="left", suffixes=("", "_ext"))
hit = MM[EXT_MODELS[0]].notna() if EXT_MODELS[0] in MM.columns else None
print("\n" + "=" * 78)
print("匹配率（(PDB4, 突变串) 精确）")
print("=" * 78)
print(f"  canonical 行 = {len(MM)}")
for c in num:
    n = int(MM[c].notna().sum())
    print(f"  {c:<14} 命中 {n:>5}  ({n/len(MM):.4f})")
anyhit = MM[num].notna().any(axis=1)
print(f"  至少一个外部模型命中 = {int(anyhit.sum())} ({anyhit.mean():.4f})")
print(f"  覆盖 target = {MM.loc[anyhit,'target_id_canonical'].nunique()}"
      f" / {MM.target_id_canonical.nunique()}")

# ---------- 三轴评价 ----------
print("\n" + "=" * 78)
print("三轴读数（每个外部模型，匹配子集内）")
print("=" * 78)
TGT = MM.target_id_canonical.astype(str).to_numpy()
FAM = MM.family_id_canonical.astype(str).to_numpy()
rows = []
print(f"{'model':<14}{'n':>6}{'tgt>=5':>7}{'Overall r':>11}{'RMSE':>9}"
      f"{'Macro rho':>11}{'C_W':>9}{'Eff@10':>9}")
for c in num:
    m = MM[c].notna().to_numpy() & np.isfinite(MM.y_true.to_numpy(float))
    if m.sum() < 100:
        continue
    yy = MM.y_true.to_numpy(float)[m]; pp = MM[c].to_numpy(float)[m]
    U = M.prep_universe(yy, pp, TGT[m], FAM[m])
    r = M.compute_all(U)
    r["model"] = c
    rows.append(r)
    print(f"{c:<14}{r['n_rows']:>6}{r['n_targets_ge5']:>7}"
          f"{r['overall_pearson']:>+11.4f}{r['rmse']:>9.4f}"
          f"{r['macro_spearman']:>+11.4f}{r['C_W']:>+9.4f}"
          f"{r['efficiency@10']:>+9.4f}")

# 真值自身的 sanity（ddG_ext 行）
D = pd.DataFrame(rows)
if len(D):
    D.to_csv(EXT / "external_multimodel_scoreboard.csv", index=False,
             encoding="utf-8")
    sel = D.dropna(subset=["overall_pearson", "macro_spearman", "efficiency@10"])
    if len(sel) >= 5:
        from scipy.stats import spearmanr
        r1 = spearmanr(sel.overall_pearson, sel.macro_spearman).correlation
        r2 = spearmanr(sel.overall_pearson, sel["efficiency@10"]).correlation
        r3 = spearmanr(sel.macro_spearman, sel["efficiency@10"]).correlation
        print(f"\n  外部模型池 (n={len(sel)}): "
              f"rho(r,macro)={r1:+.4f}  rho(r,eff)={r2:+.4f}  "
              f"rho(macro,eff)={r3:+.4f}")
    else:
        print(f"\n  ⚠ 可评价外部模型仅 {len(sel)} 个 ⇒ 秩相关**不可计算**（需 >=5）")
        print("     这是一个**必须如实报告的限制**：外部模型太少，")
        print("     无法在外部池内检验能力分化 ⇒ EXTERNAL_MODEL_POOL_TOO_SMALL")

meta = {"source": "https://github.com/ak422/CATH-ddG data/SKEMPI2/skempi_v2.csv",
        "n_source_rows": int(len(E)), "external_models": num,
        "match_rate": {c: float(MM[c].notna().mean()) for c in num},
        "matched_any": float(anyhit.mean()),
        "n_evaluable_external_models": int(len(D)),
        "contamination_note": ("CATH-ddG 自身 88% 行属其训练集，且**未发布自身逐突变预测**；"
                               "本表仅含基线方法（Flex-ddG/DDGPred/FoldX/"
                               "prompt-fold/ESM-IF/proteinMPNN），"
                               "其 provenance 需逐列审计（§50）"),
        "snapshot_gap": ("⚠ snapshot 一致性未逐位验证：HERMES 快照 2897 单点行中仅 "
                         "331 行命中本表 key（11.4%）。说明两数据集的 "
                         "PDB4+突变串口径存在系统性差异 ⇒ **需先解决快照问题**"
                         "才能把两者放进同一个匹配宇宙。")}
(EXT / "external_multimodel_meta.json").write_text(
    json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"\n写出 external_multimodel_scoreboard.csv / _meta.json")
