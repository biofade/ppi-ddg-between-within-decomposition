# -*- coding: utf-8 -*-
"""scripts/v2.4/01_skempi_within_target.py -- V2.4 任务 1: SKEMPI within-target 基准.

为什么用「复合物级 GroupKFold + 折外预测」
------------------------------------------
V2.2 的 test 划分里 SKEMPI 只有 **32 行** (见 `v2.4_within_target_diagnosis.md`),
无法支撑 within-target 基准。但 SKEMPI 子集总共有 **5109 行 / 288 个 interface_key**。
用复合物级 GroupKFold 可以得到**每一行的折外 (out-of-fold) 预测**,
再在复合物内部排序算 per-complex Spearman。这样:

* 用上全部 5109 行而不是 32 行;
* `group = interface_key` -> **同一复合物的突变体绝不跨折**, 因此折外预测对该复合物
  是"没见过的", **与 cross-target 划分精神一致, 不泄漏**;
* 每个复合物内部排序 = 真正的 within-target 任务形态。

评估协议
--------
1. 仅保留突变体数 >= `--min-conditions` 的复合物 (否则组内无排序可言);
2. 每个复合物内: `spearman(pred, -ΔG)`  —— **ΔG 越负越强, 我们的分数越大越强**,
   因此对标签取负后再算 Spearman (与全项目口径一致);
3. 汇总: 中位数 / IQR / 均值 / 按突变体数加权均值 / 有效复合物数;
4. 与 V2.3.2 的 cross-target Spearman (0.3890) 直接对比。

模型
----
* LightGBM (V2.2 超参, physics_only / esm2_both_raw / physics+esm2_both)
* NN 双塔 MLP (V2.3.2 规格, 同一 OOF 框架)
* rank 平均集成

**不伪造**: 任何模型若在某复合物内预测全并列 (Spearman 未定义) 则记 NaN 并排除,
并在报告中明确给出被排除的复合物数。
"""

from __future__ import annotations

import argparse
import gc
import json
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
OUT_MD = FEAT / "v2.4_skempi_within_target_report.md"
OUT_CSV = FEAT / "v2.4_skempi_within_target_results.csv"
OUT_PER = FEAT / "v2.4_skempi_per_complex.csv"
OUT_DATA = FEAT / "v2.4_skempi_within_target.csv"

# V2.2 超参 (逐字沿用)
PARAMS = dict(objective="regression", n_estimators=500, learning_rate=0.03,
              num_leaves=15, max_depth=5, min_child_samples=12, subsample=0.8,
              subsample_freq=1, colsample_bytree=0.8, reg_lambda=1.0,
              n_jobs=2, random_state=42, verbose=-1)
PCA_DIM = 128
N_FOLDS = 5
SEED = 42
CROSS_TARGET_RHO = 0.3890     # V2.3.2 的 physics+esm2_both@l2@128 结果, 用于对比
RT_LITERAL = 0.592


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
    ap = argparse.ArgumentParser(description="V2.4 SKEMPI within-target 基准")
    ap.add_argument("--min-conditions", type=int, default=3,
                    help="复合物内最少突变体数 (默认 3, 组内才可算 Spearman)")
    ap.add_argument("--folds", type=int, default=N_FOLDS)
    ap.add_argument("--pca-dim", type=int, default=PCA_DIM)
    ap.add_argument("--variant", default="l2", choices=("raw", "l2"))
    ap.add_argument("--binding", default="row_exact",
                    choices=("row_exact", "pdb_avg"),
                    help="嵌入绑定方式: row_exact=逐行精确映射 (within-target 唯一可行); "
                         "pdb_avg=按 pdb_id 平均 (与 V2.3.2 同口径, 但组内会全并列)")
    ap.add_argument("--nn", action="store_true", help="同时跑 NN 双塔")
    ap.add_argument("--nn-seeds", type=int, nargs="*", default=[42, 43])
    ap.add_argument("--n-perm", type=int, default=200,
                    help="置换对照次数 (建立 within-target 的 null 分布)")
    ap.add_argument("--limit-folds", type=int, default=0,
                    help=">0 时只用前 N 折 (调试/时限用)")
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
    emit("# V2.4 任务 1 — SKEMPI 2.0 within-target ranking 基准")
    emit("")
    emit("> 协议: **复合物级 GroupKFold** -> 每行折外预测 -> 复合物内部 Spearman。")
    emit(f"> 组内最少突变体数: **{args.min_conditions}**; 折数: **{args.folds}**; "
         f"PCA 维度: {args.pca_dim}; 嵌入版本: `{args.variant}`")
    emit("")

    # ---------------- 1. 数据 ----------------
    d = pd.read_csv(V22, low_memory=False)
    sc = [c for c in d.columns if "Source" in c][0]
    s = d[d[sc].astype(str).str.contains("SKEMPI", case=False, na=False)].copy()
    s = s[s["row_structure_usable"] == 1].reset_index(drop=True)
    emit(f"## 1. 数据")
    emit("")
    emit(f"- SKEMPI 子集 (结构可用): **{len(s)} 行**, "
         f"{s['interface_key'].nunique()} 个 interface_key, {s['pdb_id'].nunique()} 个 pdb")

    # 复核标签 (从 KD 重算, 与 V2.2 的 real_dG 对齐)
    k = pd.to_numeric(s["KD(M)"], errors="coerce")
    rec = RT_LITERAL * np.log(k)
    ok = k.notna() & s["real_dG"].notna()
    emit(f"- 标签复核: 由 `KD(M)` 重算 ΔG 与 `real_dG` 最大差 "
         f"**{float((rec[ok] - s['real_dG'][ok]).abs().max()):.2e}** (符号约定一致)")

    s["_n_cond"] = s.groupby("interface_key")["interface_key"].transform("size")
    keep = s[s["_n_cond"] >= args.min_conditions].reset_index(drop=True)
    emit(f"- 保留 (突变体 >= {args.min_conditions}): **{len(keep)} 行**, "
         f"**{keep['interface_key'].nunique()} 个复合物**")
    emit(f"- 被排除: {len(s) - len(keep)} 行 "
         f"({s['interface_key'].nunique() - keep['interface_key'].nunique()} 个复合物突变体过少)")

    # 保存基准表
    keep.to_csv(OUT_DATA, index=False, encoding="utf-8")
    emit(f"- 基准表: `{OUT_DATA.name}`")

    # ---------------- 2. 特征装配 ----------------
    emit("")
    emit("## 2. 特征装配")
    emit("")
    z = np.load(CACHE, allow_pickle=True)
    pdb_ids = z["pdb_ids"].astype(str)
    lmd5, rmd5 = z["lig_md5"].astype(str), z["rec_md5"].astype(str)
    E = {
        "raw": (z["lig_raw"].astype(np.float32), z["rec_raw"].astype(np.float32)),
        "l2": (z["lig_l2"].astype(np.float32), z["rec_l2"].astype(np.float32)),
    }
    Lr_all, Rr_all = E[args.variant]

    # ---- 嵌入绑定方式 (决定性选择, 见诊断文档) ----
    # * "pdb_avg": 按 pdb_id 取平均 —— 与 V2.3.2 主结果同口径。
    #   但实测会把 4987 行压成 **197 个唯一嵌入**, 同一复合物的所有突变体特征完全相同
    #   -> 组内预测必然全并列 -> within-target Spearman 恒为 NaN。
    # * "row_exact": 缓存与 filtered.csv/V2.2 表**行序一致**, 做逐行 1:1 映射
    #   (2071 个唯一嵌入)。这是 within-target 唯一可行的绑定方式。
    if args.binding == "pdb_avg":
        emb = pd.DataFrame({"pdb_id": pdb_ids, "lig": list(Lr_all),
                            "rec": list(Rr_all)})
        agg = emb.groupby("pdb_id", sort=False).agg(
            lig=("lig", lambda c: np.mean(np.stack(c.to_numpy()), axis=0)),
            rec=("rec", lambda c: np.mean(np.stack(c.to_numpy()), axis=0)),
        ).reset_index()
        m = keep.merge(agg, on="pdb_id", how="left", validate="many_to_one")
        hit = m["lig"].notna().to_numpy()
        emit(f"- 绑定方式: **`pdb_avg`** (按 `pdb_id` 取平均, 与 V2.3.2 同口径)")
        emit(f"- ESM-2 join 命中率: **{hit.mean() * 100:.2f}%** ({int(hit.sum())}/{len(m)})")
        if hit.mean() < 0.95:
            emit("- !! 命中率 < 95%, 停止")
            OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
            return 3
        m = m[hit].reset_index(drop=True)
        Lig = np.stack(m["lig"].to_numpy()).astype(np.float32)
        Rec = np.stack(m["rec"].to_numpy()).astype(np.float32)
    else:
        # 逐行精确映射: 以 (interface_key, 组内出现次序) 为键, 在**未过滤的 V2.2 全表**
        # (8207 行, 与 ESM-2 缓存逐行同序) 上定位全局行号 -> 直接索引缓存的对应行。
        # 注意: 不能用 `row_structure_usable==1` 过滤后的表做键 —— 那是 8199 行,
        # 行号会整体偏移 (实测踩过: 8199 vs 8207 对不上)。
        v22 = pd.read_csv(V22, low_memory=False)
        emit(f"- 绑定方式: **`row_exact`** (逐行 1:1 精确映射)")
        emit(f"- V2.2 全表: {len(v22)} 行 | ESM-2 缓存: {len(pdb_ids)} 行")
        if len(v22) != len(pdb_ids):
            emit(f"- !! 行数不等 ({len(v22)} vs {len(pdb_ids)}), 无法逐行映射")
            OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
            return 4
        emit(f"- 行序一致 (pdb_id 逐位置相同): "
             f"**{bool((v22['pdb_id'].astype(str).to_numpy() == pdb_ids).all())}**")
        tmp = keep.copy()
        tmp["_seq"] = tmp.groupby("interface_key").cumcount()
        v22["_seq"] = v22.groupby("interface_key").cumcount()
        pos = v22.set_index(["interface_key", "_seq"]).index
        idx = pd.MultiIndex.from_arrays([tmp["interface_key"], tmp["_seq"]])
        loc = pos.get_indexer(idx)
        emit(f"- 映射命中: **{int((loc >= 0).sum())} / {len(tmp)}**")
        if (loc < 0).any():
            emit("- !! 有行无法定位, 停止")
            OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
            return 5
        m = keep.reset_index(drop=True)
        Lig = Lr_all[loc]
        Rec = Rr_all[loc]
        u = np.unique(Lig, axis=0).shape[0]
        emit(f"- 映射后唯一配体嵌入行数: **{u} / {len(Lig)}** "
             f"(对比 `pdb_avg` 会把 4987 行压成 197 行)")
        # 逐复合物核对: 组内是否有可区分的嵌入
        lg = pd.DataFrame({"g": m["interface_key"].to_numpy(),
                           "r": [hash(x.tobytes()) for x in Lig],
                           "c": [hash(x.tobytes()) for x in Rec]})
        per_g = lg.groupby("g").agg(nl=("r", "nunique"), nc=("c", "nunique"))
        emit(f"- 组内**配体**嵌入有变化的复合物: "
             f"**{int((per_g['nl'] > 1).sum())} / {len(per_g)}**")
        emit(f"- 组内**受体**嵌入有变化的复合物: "
             f"**{int((per_g['nc'] > 1).sum())} / {len(per_g)}**")
        emit(f"- 组内**两边都不变** (任何模型都无法区分) 的复合物: "
             f"**{int(((per_g['nl'] == 1) & (per_g['nc'] == 1)).sum())}**")

    phys, counts = C.load_feature_names()
    fcols = [c for c in phys if c in m.columns]
    Xp = m[fcols].apply(pd.to_numeric, errors="coerce").copy()
    for c in counts:
        if c in Xp.columns:
            Xp[c] = np.log1p(Xp[c].clip(lower=0))
    # 中位数兜底: 用**全体**中位数 (折内再算会因小样本不稳; 此处无标签泄漏)
    Xp = Xp.fillna(Xp.median(numeric_only=True)).fillna(0.0)[fcols]
    Xp = Xp.to_numpy(np.float32)
    y = m["real_dG"].to_numpy(float)
    groups = m["interface_key"].astype(str).to_numpy()
    emit(f"- 物理特征 {Xp.shape[1]} 维; 配体嵌入 {Lig.shape[1]} 维; 受体嵌入 {Rec.shape[1]} 维")
    # 物理特征的组内变异 (关键诊断)
    _cmp = pd.DataFrame({"g": groups, "v": Xp[:, 0]})
    n_vary = 0
    Xdf = pd.DataFrame(Xp)
    for c in range(Xp.shape[1]):
        if Xdf.groupby(groups)[c].nunique().gt(1).any():
            n_vary += 1
    emit(f"- **组内会变化的物理特征数: {n_vary} / {Xp.shape[1]}** "
         + ("(-> physics_only 在原理上无法区分同一复合物的突变体)"
            if n_vary == 0 else ""))

    # PCA (拟合在**训练折**内 —— 由下面的循环对每折单独 fit, 避免任何泄漏)
    from sklearn.decomposition import PCA

    # ---------------- 3. 折外预测 ----------------
    from sklearn.model_selection import GroupKFold
    ug = np.unique(groups)
    gkf = GroupKFold(n_splits=args.folds)
    folds = list(gkf.split(ug, groups=ug))
    if args.limit_folds:
        folds = folds[:args.limit_folds]
    emit("")
    emit(f"## 3. 折外预测 ({len(folds)} 折, group = interface_key)")
    emit("")
    for i, (tr_i, te_i) in enumerate(folds):
        emit(f"- 折 {i}: 训练 {len(tr_i)} 复合物 / 测试 {len(te_i)} 复合物")

    FEATURE_SETS = ["physics_only", "esm2_both", "physics+esm2_both"]

    def build(fs: str, tr_groups: set, pca: PCA | None) -> np.ndarray:
        if fs == "physics_only":
            return Xp
        if pca is None:
            raise ValueError("需要 PCA")
        zb = pca.transform(np.hstack([Lig, Rec]))
        if fs == "esm2_both":
            return zb.astype(np.float32)
        return np.hstack([Xp, zb]).astype(np.float32)

    oof: dict[str, np.ndarray] = {}
    oof_pred_store: dict[str, np.ndarray] = {fs: np.full(len(m), np.nan)
                                             for fs in FEATURE_SETS}
    fit_rows: list[dict] = []
    for fi, (tr_i, te_i) in enumerate(folds):
        tr_mask = np.isin(groups, ug[tr_i])
        te_mask = np.isin(groups, ug[te_i])
        # PCA 只在**该折的训练部分 + 该复合物的训练行**上 fit
        pca = PCA(n_components=args.pca_dim, random_state=SEED)
        pca.fit(np.hstack([Lig, Rec])[tr_mask])
        for fs in FEATURE_SETS:
            try:
                X = build(fs, set(ug[tr_i]), pca)
                mdl = lgb.LGBMRegressor(**PARAMS)
                mdl.fit(X[tr_mask], y[tr_mask])
                oof_pred_store[fs][te_mask] = mdl.predict(X[te_mask])
                fit_rows.append({"fold": fi, "feature_set": fs,
                                 "n_train": int(tr_mask.sum()),
                                 "n_test": int(te_mask.sum()),
                                 "n_features": int(X.shape[1])})
            except Exception as exc:  # noqa: BLE001
                emit(f"  !! 折 {fi} {fs} 失败: {exc}")
                with open(FAIL_LOG, "a", encoding="utf-8") as fh:
                    fh.write(f"[01_skempi_within_target] fold{fi} {fs}: {exc}\n"
                             f"{traceback.format_exc()}\n")
        gc.collect()
        emit(f"  折 {fi} 完成 ({time.time() - t0:.0f}s)")

    oof = oof_pred_store

    # ---------------- 4. 可选: NN 双塔 ----------------
    if args.nn:
        emit("")
        emit("## 4. NN 双塔 (同一 OOF 框架, CPU only)")
        emit("")
        try:
            import torch
            import torch.nn as nn

            d_emb = args.pca_dim
            class TwoTower(nn.Module):
                def __init__(self, dl, dr, dp):
                    super().__init__()
                    self.lig = nn.Sequential(nn.Linear(dl, 256), nn.ReLU(),
                                             nn.Dropout(0.3), nn.Linear(256, 64))
                    self.rec = nn.Sequential(nn.Linear(dr, 256), nn.ReLU(),
                                             nn.Dropout(0.3), nn.Linear(256, 64))
                    self.phys = nn.Sequential(nn.Linear(dp, 64), nn.ReLU())
                    self.head = nn.Sequential(nn.Linear(192, 64), nn.ReLU(),
                                              nn.Linear(64, 1))
                def forward(self, b):
                    h = torch.cat([self.lig(b["lig"]), self.rec(b["rec"]),
                                   self.phys(b["phys"])], dim=1)
                    return self.head(h)

            nn_oof = {f"nn_s{sd}": np.full(len(m), np.nan) for sd in args.nn_seeds}
            import torch as _t
            import torch.nn as nn

            for fi, (tr_i, te_i) in enumerate(folds):
                # 训练折内再切出一个**按复合物分组**的验证集做早停监测。
                # 注意: 早停监测用训练分布 (与 V2.3.2 一致), 评估仍只看 test 折;
                # 验证集与训练集在**复合物层面**不重叠, 避免早停利用测试信息。
                rs = np.random.RandomState(SEED + fi)
                perm = rs.permutation(len(tr_i))
                n_va = max(1, int(round(len(tr_i) * 0.15)))
                va_g = set(np.asarray(ug)[tr_i][perm[:n_va]])
                tr_g = set(np.asarray(ug)[tr_i][perm[n_va:]])
                m_tr = np.isin(groups, list(tr_g))
                m_va = np.isin(groups, list(va_g))
                m_te = np.isin(groups, ug[te_i])
                pca = PCA(n_components=args.pca_dim, random_state=SEED)
                pca.fit(np.hstack([Lig, Rec])[m_tr])
                Z = pca.transform(np.hstack([Lig, Rec])).astype(np.float32)
                Zl, Zr = Z[:, :args.pca_dim], Z[:, args.pca_dim:]
                for sd in args.nn_seeds:
                    try:
                        r = C.train_nn(
                            lambda n: TwoTower(n["dl"], n["dr"], n["dp"]),
                            {},
                            {"dl": Zl.shape[1], "dr": Zr.shape[1], "dp": Xp.shape[1]},
                            y[m_tr], y[m_va],
                            {"lig": Zl[m_tr], "rec": Zr[m_tr], "phys": Xp[m_tr]},
                            {"lig": Zl[m_va], "rec": Zr[m_va], "phys": Xp[m_va]},
                            {"lig": Zl[m_te], "rec": Zr[m_te], "phys": Xp[m_te]},
                            seed=sd)
                        nn_oof[f"nn_s{sd}"][m_te] = r["test_pred"]
                    except Exception as exc:  # noqa: BLE001
                        emit(f"  !! 折 {fi} NN seed{sd} 失败: {exc}")
                        with open(FAIL_LOG, "a", encoding="utf-8") as fh:
                            fh.write(f"[01_skempi_within_target] NN fold{fi} s{sd}: {exc}\n"
                                     f"{traceback.format_exc()}\n")
                emit(f"  折 {fi} NN 完成 ({time.time() - t0:.0f}s)")
            ok_nn = [k for k, v in nn_oof.items() if np.isfinite(v).all()]
            if ok_nn:
                oof["nn_two_tower"] = np.mean([nn_oof[k] for k in ok_nn], axis=0)
                emit(f"- NN 折外预测完整, 用 {len(ok_nn)} 个种子平均")
            else:
                emit("- !! NN 折外预测不完整, 跳过")
        except Exception as exc:  # noqa: BLE001
            emit(f"- !! NN 分支失败: {exc}")
            with open(FAIL_LOG, "a", encoding="utf-8") as fh:
                fh.write(f"[01_skempi_within_target] NN branch: {exc}\n"
                         f"{traceback.format_exc()}\n")

    # 集成 (rank 平均)
    if "nn_two_tower" in oof and np.isfinite(oof["physics+esm2_both"]).all():
        from scipy.stats import rankdata
        oof["ensemble_rankmean"] = (rankdata(oof["physics+esm2_both"])
                                    + rankdata(oof["nn_two_tower"])) / 2.0

    # ---------------- 5. per-complex Spearman ----------------
    emit("")
    emit("## 5. Within-target 评估 (复合物内部排序)")
    emit("")
    per_rows: list[dict] = []
    for model, pred in oof.items():
        if not np.isfinite(pred).all():
            emit(f"- `{model}`: 折外预测不完整, **跳过** "
                 f"(缺失 {int((~np.isfinite(pred)).sum())} 行)")
            continue
        for g in np.unique(groups):
            msk = groups == g
            if msk.sum() < args.min_conditions:
                continue
            yv, pv = y[msk], pred[msk]
            if len(np.unique(yv)) < 2:
                per_rows.append({"model": model, "interface_key": g, "n": int(msk.sum()),
                                 "spearman": np.nan, "excluded": "标签全并列"})
                continue
            if len(np.unique(pv)) < 2:
                per_rows.append({"model": model, "interface_key": g, "n": int(msk.sum()),
                                 "spearman": np.nan, "excluded": "预测全并列"})
                continue
            rho = mr.spearman_rho(-yv, -pv)
            per_rows.append({"model": model, "interface_key": g, "n": int(msk.sum()),
                             "spearman": float(rho) if rho == rho else np.nan,
                             "excluded": ""})

    per = pd.DataFrame(per_rows, columns=["model", "interface_key", "n",
                                          "spearman", "excluded"])
    per.to_csv(OUT_PER, index=False, encoding="utf-8")

    summary_rows: list[dict] = []
    emit("| 模型 | 复合物数 | 有效 | 排除 | **中位 Spearman** | IQR | 均值 | 加权均值 |")
    emit("|---|---:|---:|---:|---:|---|---:|---:|")
    for model in oof:
        if not len(per):
            break
        sub = per[per["model"] == model]
        if not len(sub):
            continue
        val = sub["spearman"].dropna()
        exc = int(sub["spearman"].isna().sum())
        if not len(val):
            emit(f"| {model} | {len(sub)} | 0 | {exc} | — | — | — | — |")
            continue
        w = sub.loc[val.index, "n"].to_numpy(float)
        wmean = float(np.average(val.to_numpy(float), weights=w))
        q1, q3 = np.percentile(val, [25, 75])
        summary_rows.append({
            "model": model, "n_complex": int(len(sub)), "n_valid": int(len(val)),
            "n_excluded": exc,
            "spearman_median": float(np.median(val)), "iqr_low": float(q1),
            "iqr_high": float(q3), "spearman_mean": float(val.mean()),
            "spearman_weighted_mean": wmean,
            "spearman_std": float(val.std(ddof=1)) if len(val) > 1 else float("nan"),
            "median_n_conditions": float(sub["n"].median()),
        })
        emit(f"| {model} | {len(sub)} | {len(val)} | {exc} | "
             f"**{np.median(val):+.4f}** | [{q1:+.3f}, {q3:+.3f}] | "
             f"{val.mean():+.4f} | {wmean:+.4f} |")

    sdf = pd.DataFrame(summary_rows).sort_values("spearman_median", ascending=False)
    sdf.to_csv(OUT_CSV, index=False, encoding="utf-8")

    # ---------------- 5b. 置换对照: 建立 null ----------------
    # within-target 的组内 Spearman 在"随机预测"下期望值**不是 0** 的近邻吗?
    # 实际上随机预测的组内 Spearman 期望是 0, 但**噪声很大** (组内只有 3-10 个点)。
    # 因此必须给出 null 分布, 才能判断 +0.10 是否有意义。
    emit("")
    emit("## 5b. 置换对照 (null 分布)")
    emit("")
    rs = np.random.RandomState(SEED)
    n_perm = args.n_perm
    null_med: list[float] = []
    for _ in range(n_perm):
        yperm = rs.permutation(y)
        vals = []
        for g in np.unique(groups):
            msk = groups == g
            if msk.sum() < args.min_conditions:
                continue
            yv, pv = yperm[msk], oof["physics+esm2_both"][msk]
            if len(np.unique(yv)) < 2 or len(np.unique(pv)) < 2:
                continue
            rr = mr.spearman_rho(-yv, -pv)
            if rr == rr:
                vals.append(rr)
        if vals:
            null_med.append(float(np.median(vals)))
    if null_med:
        na = np.asarray(null_med)
        obs = float(sdf.iloc[0]["spearman_median"]) if len(sdf) else float("nan")
        p_emp = float((np.sum(na >= obs) + 1) / (len(na) + 1))
        emit(f"- 固定特征、**打乱标签** {n_perm} 次, 每次重算全体复合物的中位 Spearman")
        emit(f"- null 分布: 中位 {np.median(na):+.4f}, "
             f"95% 区间 [{np.percentile(na, 2.5):+.4f}, {np.percentile(na, 97.5):+.4f}]")
        emit(f"- 观测值 {obs:+.4f} -> 经验单侧 p = **{p_emp:.4f}**")
        if p_emp < 0.05:
            emit(f"- -> **观测显著高于随机** (但幅度仍然很小, 见下)")
        else:
            emit(f"- -> **观测与随机不可区分** —— within-target 上模型基本没有学到东西")
        perm_res = {"n_perm": n_perm, "null_median": float(np.median(na)),
                    "null_lo": float(np.percentile(na, 2.5)),
                    "null_hi": float(np.percentile(na, 97.5)),
                    "observed": obs, "p_empirical": p_emp}
    else:
        perm_res = {}
        emit("- 置换对照无法计算 (没有有效复合物)")

    # ---------------- 6. 与 cross-target 对比 + 显著性 ----------------
    emit("")
    emit("## 6. 与 cross-target 对比")
    emit("")
    emit(f"- V2.3.2 cross-target (Proteinea test, n=752) 最优: "
         f"**{CROSS_TARGET_RHO:+.4f}** (physics+esm2_both@l2@128)")
    if len(sdf):
        best = sdf.iloc[0]
        emit(f"- 本基准 within-target 中位 Spearman: **{best['spearman_median']:+.4f}** "
             f"(`{best['model']}`, {int(best['n_valid'])} 个复合物)")
        emit(f"- 差值: **{best['spearman_median'] - CROSS_TARGET_RHO:+.4f}**")
        emit("")
        emit("> ⚠️ **这两个数字不可直接做统计检验**: 一个是对 752 个**样本**的单一"
             "排序相关, 另一个是对数百个**复合物**各自相关的分布中位数。")
        emit("> 样本单位不同 (行 vs 复合物), 因此此处只做**描述性对比**, ")
        emit("> 不编造 p 值。要检验需要设计配对方案 (例如把 cross-target 也按复合物"
             "分组重算), 已列入下一步。")

        # 关键对照: physics_only 是否真的"瞎"
        po = sdf[sdf["model"] == "physics_only"]
        if len(po):
            emit("")
            emit(f"- **关键对照 — `physics_only` 的中位 Spearman = "
                 f"{float(po['spearman_median'].iloc[0]):+.4f}** "
                 f"(排除 {int(po['n_excluded'].iloc[0])} 个复合物)")
            emit("  物理特征在同一界面内**完全相同** (诊断: 0/63 个特征有变化),")
            emit("  因此它在原理上**无法区分同一复合物的突变体**。")
            emit("  若该值接近 0 或为 NaN 主导, 就**直接证明**了「序列信息是必要的」。")

    # ---------------- 7. 结论 ----------------
    emit("")
    emit("## 7. 结论")
    emit("")
    if len(sdf):
        best = sdf.iloc[0]
        med = float(best["spearman_median"])
        bm = best["model"]
        n_rows_covered = int(per[(per["model"] == bm)
                                 & per["spearman"].notna()]["n"].sum())
        emit(f"* 最优模型: `{bm}` -> 中位 within-target Spearman "
             f"**{med:+.4f}** (IQR [{best['iqr_low']:+.3f}, {best['iqr_high']:+.3f}]), "
             f"覆盖 {int(best['n_valid'])} 个复合物 / {n_rows_covered} 行")
        if med >= 0.55:
            emit(f"* **达到 V2.4 任务 5.1 的「分支 A」阈值 (>= 0.55)** "
                 f"-> within-target 是更合适的任务形态。")
        elif med >= 0.30:
            emit(f"* 中等水平 (0.30–0.55): 优于随机但未达分支 A 阈值。")
        else:
            emit(f"* 偏低 (< 0.30): within-target 并未像 V2.3.2 的 n=34 结果暗示的那么好。")
            emit(f"  -> 说明 V2.3.2 的 ρ=+0.5848 **很可能是 n=34 的偶然波动**, ")
            emit(f"     本基准用 {int(best['n_valid'])} 个复合物给出了更可靠的估计。")
    else:
        emit("* 没有任何模型产出可评估的折外预测。")

    emit("")
    emit("## 8. 局限")
    emit("")
    emit("1. **PCA 与中位数兜底**: PCA 按折内训练部分 fit (无泄漏); 物理特征缺失值用"
         "全表中位数兜底 (无标签泄漏, 但严格说应折内计算)。")
    emit(f"2. **嵌入绑定用逐行精确映射 (`row_exact`)**, 不是 V2.3.2 主结果的 "
         f"`pdb_id` 平均。原因是 `pdb_id` 平均会把 4987 行压成 **197 个唯一嵌入**,")
    emit("   同一复合物的突变体特征完全相同 -> 组内 Spearman 恒为 NaN (实测确认)。")
    emit("   因此本基准**必须**用逐行映射; 这也意味着本基准的嵌入粒度比 V2.3.2 更细,")
    emit("   两者数字不可混用。")
    emit(f"3. **per-complex Spearman 在 n=3 时极不稳定** (只有 3 个点), 这也是为什么"
         f"同时报告加权均值与 IQR, 并给出置换 null。")
    emit("4. 与 cross-target 的对比是**描述性**的, 未做统计检验 (样本单位不同: "
         "行 vs 复合物)。")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")

    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD} ({OUT_MD.stat().st_size} B)")
    print(f"写出 {OUT_CSV} | {OUT_PER} | {OUT_DATA}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
