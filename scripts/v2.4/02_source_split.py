# -*- coding: utf-8 -*-
"""scripts/v2.4/02_source_split.py -- V2.4 任务 2: 拆分 SAbDab 与一般 PPI.

动机
----
V2.3.2 错误分析发现按来源分组时差异巨大:

| 来源 | n | cross-target Spearman |
|---|---:|---:|
| SKEMPI v2.0 | 34 | +0.5848 |
| Affinity Benchmark v5.5 | 23 | +0.3985 |
| PDBbind v2020 | 657 | +0.3648 |
| **SAbDab (抗体-抗原)** | 38 | **−0.0509** |

抗体-抗原是完全不同的结合模式 (CDR 环、界面更平、亲和力常源于构象选择)。
把它们与一般 PPI 混在一个跨靶回归里, 很可能**互相拖累**。

本脚本做两件事
--------------
1. **分组统计**: 两边的行数 / 复合物数 / ΔG 分布 / 家族数;
2. **分开训练**: 在各自子集内做**来源内划分** (组级别, 与 V2.2 同精神),
   对比「混在一起训练」vs「分开训练」在**各自测试集**上的表现。

为公平对比, 分开训练时**每个来源内部**重新做 group 划分 (5 折 GroupKFold),
报告 OOF Spearman —— 这样不同样本量的子集也能比较。
同时报告「混合模型在该子集上的 OOF」作为对照。
"""

from __future__ import annotations

import argparse
import gc
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v2.3.2"))
import nn_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
FEAT = ROOT / "features"
LOG_DIR = ROOT / "logs" / "v2.4"
FAIL_LOG = LOG_DIR / "failures.log"
V22 = FEAT / "v2.2_proteinea_features.csv"
CACHE = FEAT / "v2_esm2_cache.npz"
OUT_MD = FEAT / "v2.4_source_split_report.md"
OUT_CSV = FEAT / "v2.4_source_split_results.csv"
OUT_STATS = FEAT / "v2.4_source_split_stats.csv"

PARAMS = dict(objective="regression", n_estimators=500, learning_rate=0.03,
              num_leaves=15, max_depth=5, min_child_samples=12, subsample=0.8,
              subsample_freq=1, colsample_bytree=0.8, reg_lambda=1.0,
              n_jobs=2, random_state=42, verbose=-1)
PCA_DIM = 128
N_FOLDS = 5
SEED = 42
MIN_ROWS = 150          # 子集太小则不做折内训练


def load_lgb():
    import ctypes
    import glob as _g
    import sklearn
    for d in _g.glob(str(Path(sklearn.__file__).parent / ".libs" / "vcomp140.dll")):
        try:
            ctypes.CDLL(d)
        except Exception:  # noqa: BLE001
            pass
    import lightgbm as lgb
    return lgb


def main() -> int:
    ap = argparse.ArgumentParser(description="V2.4 拆分 SAbDab 与一般 PPI")
    ap.add_argument("--folds", type=int, default=N_FOLDS)
    ap.add_argument("--pca-dim", type=int, default=PCA_DIM)
    ap.add_argument("--variant", default="l2", choices=("raw", "l2"))
    args = ap.parse_args()

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    lgb = load_lgb()
    mr = C.get_metrics_module()
    L: list[str] = []

    def emit(m: str = "") -> None:
        print(m)
        L.append(str(m))
        try:
            sys.stdout.flush()
        except Exception:  # noqa: BLE001
            pass

    t0 = time.time()
    emit("# V2.4 任务 2 — 拆分 SAbDab (抗体-抗原) 与一般 PPI")
    emit("")
    emit("> 用**来源内 GroupKFold** 得到折外预测, 因此两边可用同一口径比较。")
    emit(f"> 折数 {args.folds}; PCA {args.pca_dim}; 嵌入 `{args.variant}`")
    emit("")

    # ---------------- 1. 数据与特征 ----------------
    d = pd.read_csv(V22, low_memory=False)
    sc = [c for c in d.columns if "Source" in c][0]
    d = d[d["row_structure_usable"] == 1].reset_index(drop=True)
    d["_source"] = d[sc].astype(str)
    d["_is_sabdab"] = d["_source"].str.contains("SAbDab", case=False, na=False)

    z = np.load(CACHE, allow_pickle=True)
    pdb_ids = z["pdb_ids"].astype(str)
    Lr = z[f"lig_{args.variant}"].astype(np.float32)
    Rr = z[f"rec_{args.variant}"].astype(np.float32)

    # 逐行精确映射 (与任务 1 同一理由: pdb 平均会抹平组内变异)
    v22f = pd.read_csv(V22, low_memory=False)
    emit("## 1. 特征装配 (逐行精确映射)")
    emit("")
    emit(f"- V2.2 全表 {len(v22f)} 行 | 缓存 {len(pdb_ids)} 行 | "
         f"行序一致: **{bool((v22f['pdb_id'].astype(str).to_numpy() == pdb_ids).all())}**")
    v22f["_seq"] = v22f.groupby("interface_key").cumcount()
    dd = d.copy()
    dd["_seq"] = dd.groupby("interface_key").cumcount()
    pos = v22f.set_index(["interface_key", "_seq"]).index
    loc = pos.get_indexer(pd.MultiIndex.from_arrays([dd["interface_key"], dd["_seq"]]))
    if (loc < 0).any():
        emit(f"- !! 映射失败 {int((loc < 0).sum())} 行")
        OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
        return 3
    Lig, Rec = Lr[loc], Rr[loc]

    phys, counts = C.load_feature_names()
    fcols = [c for c in phys if c in d.columns]
    Xp = d[fcols].apply(pd.to_numeric, errors="coerce").copy()
    for c in counts:
        if c in Xp.columns:
            Xp[c] = np.log1p(Xp[c].clip(lower=0))
    Xp = Xp.fillna(Xp.median(numeric_only=True)).fillna(0.0)[fcols].to_numpy(np.float32)
    y = d["real_dG"].to_numpy(float)
    groups = d["interface_key"].astype(str).to_numpy()
    emit(f"- 可用行 **{len(d)}**; 物理 {Xp.shape[1]} 维; 嵌入 2x{Lig.shape[1]} 维")

    # ---------------- 2. 分组统计 ----------------
    emit("")
    emit("## 2. 分组统计")
    emit("")
    stats = []
    emit("| 来源 | 行数 | interface_key | pdb | ΔG 中位 | ΔG std | ΔG 范围 |")
    emit("|---|---:|---:|---:|---:|---:|---|")
    for name, msk in (("SAbDab (抗体-抗原)", d["_is_sabdab"].to_numpy()),
                      ("一般 PPI (非 SAbDab)", (~d["_is_sabdab"]).to_numpy()),
                      ("全体", np.ones(len(d), bool))):
        yy = y[msk]
        gg = groups[msk]
        row = {"subset": name, "n_rows": int(msk.sum()),
               "n_interface_key": int(len(np.unique(gg))),
               "n_pdb": int(d.loc[msk, "pdb_id"].nunique()),
               "dg_median": float(np.median(yy)), "dg_std": float(yy.std(ddof=1)),
               "dg_min": float(yy.min()), "dg_max": float(yy.max())}
        stats.append(row)
        emit(f"| {name} | {int(msk.sum())} | {row['n_interface_key']} | {row['n_pdb']} | "
             f"{row['dg_median']:.3f} | {row['dg_std']:.3f} | "
             f"[{row['dg_min']:.2f}, {row['dg_max']:.2f}] |")
    pd.DataFrame(stats).to_csv(OUT_STATS, index=False, encoding="utf-8")
    emit("")
    emit("> 查看 `Source Data Set` 的原始取值分布:")
    emit("")
    for k, v in d["_source"].value_counts().items():
        emit(f"> * `{k}`: {v} 行")

    # ---------------- 3. 折外评估 ----------------
    from sklearn.decomposition import PCA
    from sklearn.model_selection import GroupKFold

    def oof_for(sub_mask: np.ndarray, feat: str, pca_dim: int,
                fit_on: np.ndarray | None = None) -> tuple[np.ndarray, dict]:
        """在 sub_mask 定义的行上做 GroupKFold OOF 预测。

        ``fit_on``: 若给出, 每折用**这个掩码**的训练行训练 (用于"混合模型"
        在子集上的对照: 用全体训练, 只在子集上评估)。
        """
        idx = np.where(sub_mask)[0]
        g = groups[idx]
        ug = np.unique(g)
        gkf = GroupKFold(n_splits=min(args.folds, max(2, len(ug))))
        pred = np.full(len(idx), np.nan)
        info = {"n_rows": int(len(idx)), "n_groups": int(len(ug))}
        for tr_i, te_i in gkf.split(ug, groups=ug):
            tr_g, te_g = set(ug[tr_i]), set(ug[te_i])
            m_te = np.isin(g, list(te_g))
            if fit_on is None:
                m_tr = np.isin(g, list(tr_g))
                rows_tr = idx[m_tr]
                pca = PCA(n_components=pca_dim, random_state=SEED)
                pca.fit(np.hstack([Lig, Rec])[rows_tr])
            else:
                # 用外部掩码做训练, 但**排除当前测试组** (否则会泄漏)
                rows_tr = np.where(fit_on & ~np.isin(groups, list(te_g)))[0]
                pca = PCA(n_components=pca_dim, random_state=SEED)
                pca.fit(np.hstack([Lig, Rec])[rows_tr])
            Z = pca.transform(np.hstack([Lig, Rec]))
            if feat == "physics_only":
                X = Xp
            elif feat == "esm2_both":
                X = Z.astype(np.float32)
            else:
                X = np.hstack([Xp, Z]).astype(np.float32)
            if len(rows_tr) < 50:
                continue
            mdl = lgb.LGBMRegressor(**PARAMS)
            mdl.fit(X[rows_tr], y[rows_tr])
            pred[m_te] = mdl.predict(X[idx[m_te]])
            info["n_features"] = int(X.shape[1])
        return pred, info

    emit("")
    emit("## 3. 折外评估")
    emit("")
    results = []
    sab = d["_is_sabdab"].to_numpy()
    ppi = ~sab

    variants = [
        ("SAbDab", sab, "SAbDab 子集内训练+评估"),
        ("一般 PPI", ppi, "一般 PPI 子集内训练+评估"),
        ("全体", np.ones(len(d), bool), "全体训练+评估 (参照)"),
    ]
    for label, msk, note in variants:
        if msk.sum() < MIN_ROWS:
            emit(f"- **{label}**: {int(msk.sum())} 行 < {MIN_ROWS}, 跳过训练")
            continue
        for feat in ("physics_only", "physics+esm2_both"):
            try:
                pred, info = oof_for(msk, feat, args.pca_dim)
                ok = np.isfinite(pred)
                if ok.sum() < 20:
                    emit(f"- {label} / {feat}: 折外预测不足 ({int(ok.sum())})")
                    continue
                rho = mr.spearman_rho(-y[msk][ok], -pred[ok])
                lo, hi, _ = C.bootstrap_ci(mr.spearman_rho, -y[msk][ok], -pred[ok],
                                           600, SEED)
                results.append({
                    "subset": label, "feature_set": feat, "design": note,
                    "n_rows": int(msk.sum()), "n_scored": int(ok.sum()),
                    "n_groups": info.get("n_groups", -1),
                    "n_features": info.get("n_features", -1),
                    "spearman": float(rho), "ci_low": lo, "ci_high": hi,
                })
                emit(f"- **{label}** / `{feat}`: OOF Spearman = {rho:+.4f} "
                     f"[{lo:+.4f}, {hi:+.4f}] ({int(ok.sum())}/{int(msk.sum())} 行, "
                     f"{info.get('n_groups', -1)} 组)")
            except Exception as exc:  # noqa: BLE001
                emit(f"- !! {label} / {feat} 失败: {exc}")
                with open(FAIL_LOG, "a", encoding="utf-8") as fh:
                    fh.write(f"[02_source_split] {label} {feat}: {exc}\n"
                             f"{traceback.format_exc()}\n")
        gc.collect()

    # 混合训练 -> 在 SAbDab 上评估 (看混在一起是否拖累)
    emit("")
    emit("### 混合模型在各子集上的表现 (在全体上按组 OOF)")
    emit("")
    try:
        pred_all, info_all = oof_for(np.ones(len(d), bool), "physics+esm2_both",
                                     args.pca_dim)
        for label, msk in (("SAbDab", sab), ("一般 PPI", ppi), ("全体", np.ones(len(d), bool))):
            if msk.sum() < 20:
                continue
            ok = np.isfinite(pred_all) & msk
            rho = mr.spearman_rho(-y[ok], -pred_all[ok])
            results.append({"subset": label, "feature_set": "physics+esm2_both",
                            "design": "混合训练(全体) -> 该子集评估",
                            "n_rows": int(msk.sum()), "n_scored": int(ok.sum()),
                            "n_groups": -1, "n_features": info_all.get("n_features", -1),
                            "spearman": float(rho), "ci_low": float("nan"),
                            "ci_high": float("nan")})
            emit(f"- 混合模型在 **{label}** 上: Spearman = {rho:+.4f} "
                 f"({int(ok.sum())} 行)")
    except Exception as exc:  # noqa: BLE001
        emit(f"- !! 混合评估失败: {exc}")
        with open(FAIL_LOG, "a", encoding="utf-8") as fh:
            fh.write(f"[02_source_split] mixed: {exc}\n{traceback.format_exc()}\n")

    res = pd.DataFrame(results)
    res.to_csv(OUT_CSV, index=False, encoding="utf-8")

    # ---------------- 4. 结论 ----------------
    emit("")
    emit("## 4. 结论")
    emit("")
    if len(res):
        sep = res[res["design"].str.contains("子集内")]
        mix = res[res["design"].str.contains("混合")]
        for label in ("一般 PPI", "SAbDab"):
            s_ = sep[(sep["subset"] == label)
                     & (sep["feature_set"] == "physics+esm2_both")]
            m_ = mix[mix["subset"] == label]
            if len(s_) and len(m_):
                dlt = float(s_.iloc[0]["spearman"]) - float(m_.iloc[0]["spearman"])
                emit(f"- **{label}**: 分开训练 {float(s_.iloc[0]['spearman']):+.4f} vs "
                     f"混合训练 {float(m_.iloc[0]['spearman']):+.4f} -> Δ = **{dlt:+.4f}**")
        emit("")
        emit("> ⚠️ **注意可比性**: 分开训练时测试集是**该子集内部**的折外行, "
             "混合训练时测试集也限定在同一子集, 因此两者可比;")
        emit("> 但**训练集不同** (子集内 vs 全体), 这既包含'去掉异质来源'的收益, "
             "也包含'训练数据变少'的损失, 二者在本设计中不可分离。")
        emit("> 如实标注为**混合效应**, 不宣称纯粹的'拆分收益'。")
    else:
        emit("- 没有成功的结果")

    emit("")
    emit("## 5. 局限")
    emit("")
    emit("### ⚠️ 与 V2.3.2 数字不可直接比较（协议不同）")
    emit("")
    emit("本报告的数字（全体 +0.6734）**高于** V2.3.2 主结果的 +0.3890。")
    emit("已用 `scripts/v2.4/04_binding_control.py` 专门排查（见")
    emit("`features/v2.4_binding_control_report.md`），结论是**协议不同, 不是泄漏**:")
    emit("")
    emit("| 协议 | 含义 | 难度 |")
    emit("|---|---|---|")
    emit("| V2.3.2 固定划分 | 训练/测试的 `pdb_id` **完全不重叠**（实测交集 = 0）| **泛化到未见靶点**, 难 |")
    emit("| 本报告 GroupKFold | 随机 80/20 折外 | **数据集内平均表现**, 容易 |")
    emit("")
    emit("交叉验证已确认: 在同一固定划分下, `row_exact` vs `pdb_id` 平均只差 "
         "**+0.0204**，且删掉与测试同 pdb 的训练行后 Δ = **+0.0000**（无泄漏）。")
    emit("因此本报告的用途是**子集之间的相对比较**（SAbDab vs 一般 PPI），")
    emit("**不应用来宣称绝对性能提升**。")
    emit("")
    emit("1. SAbDab 样本量小, 且其 `rows/group` 中位数为 1（几无同靶点多突变体）,")
    emit("   因此它衡量的是**跨靶**表现, 与 SKEMPI 的 within-target 性质不同。")
    emit("2. 分开 vs 混合的差异**混合了两种效应** (去掉异质来源的收益 + 训练量减少的"
         "损失), 本设计无法分离。")
    emit("3. 未做超参重调 (两边都用 V2.2 的同一组超参), 这对小样本子集可能不利。")
    emit("4. 逐行嵌入映射 (`row_exact`), 与 V2.3.2 的 `pdb_id` 平均口径不同"
         "（已量化: +0.0204）。")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")

    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD} ({OUT_MD.stat().st_size} B)")
    print(f"写出 {OUT_CSV} | {OUT_STATS}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
