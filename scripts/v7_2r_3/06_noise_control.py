# -*- coding: utf-8 -*-
"""scripts/v7_2r_3/06_noise_control.py -- V7.2R-3 机制判别（决定性收尾实验）.

科学动机
--------
V7.2R-3 发现 `B3a`（base + Δr + **静态 target 级 PPI 上下文** `G_tgt`）
相对 `B2`（base + Δr）显著提升 Macro ρ（+0.0710 [+0.0147, +0.1345]）。

但 `G_tgt` 对同一 target 内所有突变**数值上几乎恒定**
（实测 target 内最大逐维标准差 **2.34e-05**，float32 精度量级）。
理论上 target-常定特征在**加性模型**中不可能改变靶点内排序；
而实测 B3a 与 B2 的靶点内预测 Spearman 中位仅 0.774 →
说明 LightGBM 在**噪声级差异**上找到了分裂点。

判别设计
--------
构造一个**同等方差、零信息**的对照特征 `G_noise`：
逐 target、逐维地用真实 `G_tgt` 的**该 target 内标准差**作为尺度，
生成独立高斯噪声（seed=42 固定）。然后跑同一冻结嵌套 family-CV。

判据（§14 的机制问题是否解决）
---------------------------
* **若 `NC_noise` ≈ `B3a`（差异不显著）** → B3a 的增益是
  **噪声级伪随机排序**，不是 PPI 信息
  ⇒ `CONTINUE_REPRESENTATION_ROUTE = NO`
* **若 `B3a` 显著 > `NC_noise`** → 增益含有真实 target 级 PPI 信息（仍需续跑确认覆盖率）
"""

from __future__ import annotations

import gc
import importlib.util
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

V711 = C.ROOT / "experiments" / "v7_1_1"
V723 = C.ROOT / "experiments" / "v7_2r_3"
PH8 = V723 / "phase8_noise_control"
PH8.mkdir(parents=True, exist_ok=True)
IMP_B = C.ROOT / "features" / "v7_2r_2" / "trackB_imported"
IMP_P = C.ROOT / "features" / "v7_2r_3" / "ppi_imported"
TREES = 300
PCA_NEW = 32
MIN_PS_N = 5
SEEDS = [42, 43, 44]
N_BOOT = 1500
L: list[str] = []

MODELS = {
    "B2_dr":            ["base", "B_dr"],
    "B3a_real_PPI":     ["base", "B_dr", "G_real"],
    "NC_noise_PPI":     ["base", "B_dr", "G_noise"],
    "B_shuf_target":    ["base", "B_dr", "G_shuf"],
}


def emit(m=""):
    print(m); L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def mset(y, p, g):
    uq = np.unique(g)
    By = np.array([y[g == c].mean() for c in uq])
    Bp = np.array([p[g == c].mean() for c in uq])
    cnt = np.array([(g == c).sum() for c in uq], float)
    den = float(y.std()) * float(p.std())
    if den <= 0:
        den = np.nan
    return {"macro_spearman": FM.macro_spearman(y, p, g, min_n=MIN_PS_N),
            "eff10": FM.eff10(y, p, g, n_rand=20),
            "C_W": float((cnt / cnt.sum() * np.array(
                [np.mean((y[g == c] - By[i]) * (p[g == c] - Bp[i]))
                 for i, c in enumerate(uq)])).sum()) / den,
            "overall_r": float(C.pearson(y, p)),
            "rmse": float(np.sqrt(np.mean((y - p) ** 2)))}


def main() -> int:
    t0 = time.time()
    emit("# V7.2R-3 PHASE 8 — §14 机制判别：同方差噪声对照")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    pred_p = V723 / "phase5_eval" / "ppi_predictions.npz"
    if not pred_p.exists():
        emit("## 状态：**WAITING_FOR_V7_2R_3_EVAL**")
        (PH8 / "noise_control.md").write_text("\n".join(L) + "\n",
                                              encoding="utf-8")
        return 2
    z = np.load(pred_p, allow_pickle=True)
    rows = z["row_id"].astype(np.int64)
    y = z["y"].astype(float)
    g = np.array([str(x) for x in z["pdb"]])
    fam = np.array([str(x) for x in z["family"]])

    d = C.build_pair_features(verbose=False)
    can = pd.read_csv(C.ROOT / "experiments" / "v6_2" / "phase2_split4" /
                      "matched_analysis_set.csv", low_memory=False)
    can_i = can.set_index("row")
    blk = d["X_blocks"]
    base = np.hstack([blk["phys_diff"], blk["Z_diff"]])[rows].astype(np.float32)

    def load(dirp, keys):
        out = {}
        for p in sorted(dirp.glob("*.npz")):
            try:
                q = np.load(p, allow_pickle=True)
            except Exception:  # noqa: BLE001
                continue
            if not all(k in q.files for k in keys):
                continue
            for i, r in enumerate(q["row_id"].astype(np.int64)):
                out[int(r)] = {k: q[k][i].astype(np.float32) for k in keys}
        return out
    B = load(IMP_B, ["delta_r"])
    P = load(IMP_P, ["r_wt"])
    Bdr = np.stack([B[int(r)]["delta_r"] for r in rows]).astype(np.float32)
    Pwt = np.stack([P[int(r)]["r_wt"] for r in rows]).astype(np.float32)

    # target 级真实 g_c
    Greal = np.zeros_like(Pwt)
    tgt_idx = {}
    for i, t in enumerate(g):
        tgt_idx.setdefault(t, []).append(i)
    for t, ii in tgt_idx.items():
        Greal[ii] = Pwt[ii].mean(0)

    # ---- 同方差噪声对照 ----
    # 用「每个 target 内 G_tgt 的**逐维标准差」作为噪声尺度（逐目标、逐维校准）。
    # 补一个最小地板，使噪声在 float16 存储精度之上仍然可被树模型利用：
    #   G_tgt 本身以 float16 存盘，其可分辨步长 ≈ 2^-11 × |值| → 地板取 1e-4
    # 这对 MINT 是**更保守**的对照（噪声更有机会被利用 → 更严格地检验 PPI 特异性）。
    NOISE_FLOOR = 1e-4
    keys = sorted(tgt_idx)
    rs = np.random.RandomState(42)
    Gnoise = np.zeros_like(Greal)
    per_tgt_sd = []
    for t, ii in tgt_idx.items():
        sd_t = Greal[ii].std(0).astype(np.float32) if len(ii) > 1 \
            else np.zeros(Greal.shape[1], np.float32)
        sd_use = np.maximum(sd_t, NOISE_FLOOR)
        per_tgt_sd.append(float(sd_use.max()))
        Gnoise[ii] = rs.normal(0.0, 1.0, size=(len(ii), Greal.shape[1])) \
            * sd_use
    emit("## 8.1 对照尺度校准")
    emit("")
    emit("| 量 | 值 |")
    emit("|---|---:|")
    emit(f"| `G_real` target 内最大逐维 SD（全 target max） | "
         f"**{Greal[tgt_idx[keys[0]]].std(0).max():.3e}** … 全局 max "
         f"**{max(float(np.max(Greal[ii].std(0))) for ii in tgt_idx.values() if len(ii) > 1):.3e}** |")
    emit(f"| 噪声地板 `NOISE_FLOOR` | **{NOISE_FLOOR:.1e}** |")
    emit(f"| `G_noise` 逐维 SD 中位（target 内） | "
         f"**{np.median([float(np.median(Gnoise[ii].std(0))) for ii in tgt_idx.values() if len(ii) > 1]):.3e}** |")
    emit(f"| `G_real` 逐维 SD 中位（target 内） | "
         f"**{np.median([float(np.median(Greal[ii].std(0))) for ii in tgt_idx.values() if len(ii) > 1]):.3e}** |")
    emit("")
    emit("> `G_noise` 是**零信息**对照：逐 target、逐维匹配 `G_real` 的可用方差尺度")
    emit("> （并施加 1e-4 地板，使其**对 MINT 更有利** → 更严格）。")
    emit("> 差别仅在于 **real 的扰动与复合物身份相关，noise 的与任何东西无关**。")
    emit("")

    # ---- target-shuffle 对照 ----
    perm = np.random.RandomState(42).permutation(len(keys))
    t2s = {keys[j]: keys[perm[j]] for j in range(len(keys))}
    Gshuf = np.zeros_like(Greal)
    for t, ii in tgt_idx.items():
        Gshuf[ii] = Greal[tgt_idx[t2s[t]][0]]     # 取该 target 的常量向量

    RAW = {"base": base, "B_dr": Bdr, "G_real": Greal, "G_noise": Gnoise,
           "G_shuf": Gshuf}

    # ---- 冻结嵌套 family-CV ----
    from sklearn.decomposition import PCA
    D4 = {"y": y, "pdb": g, "fam": fam,
          "seq40": can_i.loc[rows, "seq40_cluster"].astype(str).to_numpy(),
          "obs": can_i.loc[rows, "obs_key"].astype(str).to_numpy(),
          "state": can_i.loc[rows, "state_id"].astype(str).to_numpy()}
    emit("## 8.2 冻结嵌套 family-CV")
    emit("")
    emit("| Model | Macro ρ | Eff@10 | C_W | Overall r |")
    emit("|---|---:|---:|---:|---:|")
    res, preds = [], {}
    for sd in SEEDS:
        outer = s4.make_folds(D4["fam"], sd, 5)
        for mk, blocks in MODELS.items():
            yhat = np.full(len(y), np.nan)
            for tr, te in outer:
                assert not (set(fam[tr]) & set(fam[te])), "outer 泄漏"
                mtr, mte = [], []
                for b in blocks:
                    X = RAW[b]
                    if b == "base":
                        mtr.append(X[tr]); mte.append(X[te]); continue
                    mu = np.nanmean(X[tr], axis=0)
                    sv = np.nanstd(X[tr], axis=0); sv[sv == 0] = 1
                    Ztr = (X[tr] - mu) / sv; Zte = (X[te] - mu) / sv
                    kk = min(PCA_NEW, Ztr.shape[1], max(2, len(tr) - 1))
                    pc = PCA(n_components=kk, random_state=sd).fit(Ztr)
                    mtr.append(pc.transform(Ztr).astype(np.float32))
                    mte.append(pc.transform(Zte).astype(np.float32))
                Xtr = np.hstack(mtr).astype(np.float32)
                Xte = np.hstack(mte).astype(np.float32)
                lgb = C.load_lgb()
                mdl = lgb.LGBMRegressor(**{**C.LGB_PARAMS,
                                           "n_estimators": TREES,
                                           "random_state": sd})
                mdl.fit(np.vstack([Xtr, -Xtr]), np.r_[y[tr], -y[tr]])
                yhat[te] = mdl.predict(Xte)
                del mdl; gc.collect()
            ok = np.isfinite(yhat)
            m = mset(y[ok], yhat[ok], g[ok])
            res.append({"seed": sd, "model": mk, "n_rows": int(ok.sum()),
                        "n_target": len(set(g[ok])), **m})
            preds[f"{mk}_s{sd}"] = yhat.astype(np.float32)
        emit(f"| （seed {sd} 完成） | | | | |")
    R = pd.DataFrame(res)
    R.to_csv(PH8 / "noise_control_results.csv", index=False, encoding="utf-8")
    preds["y"] = y.astype(np.float32); preds["pdb"] = g.astype(str)
    preds["family"] = fam.astype(str); preds["row_id"] = rows.astype(np.int64)
    np.savez(PH8 / "noise_control_predictions.npz", **preds)

    emit("")
    emit("### seed-mean 与 ensemble 两口径")
    emit("")
    sm = R.groupby("model")[["macro_spearman", "eff10", "C_W",
                             "overall_r"]].mean()
    ens = {}
    for mk in MODELS:
        ks = [f"{mk}_s{s}" for s in SEEDS if f"{mk}_s{s}" in preds]
        Pm = np.mean([preds[k].astype(float) for k in ks], axis=0)
        ok = np.isfinite(Pm)
        ens[mk] = mset(y[ok], Pm[ok], g[ok])
    emit("| Model | seed-mean Macro ρ | ensemble Macro ρ | Eff@10 | C_W |")
    emit("|---|---:|---:|---:|---:|")
    for mk in MODELS:
        emit(f"| {mk} | {sm.loc[mk,'macro_spearman']:+.4f} | "
             f"{ens[mk]['macro_spearman']:+.4f} | {ens[mk]['eff10']:+.4f} | "
             f"{ens[mk]['C_W']:+.4f} |")
    emit("")

    # ---- 关键对比的 bootstrap ----
    emit("## 8.3 关键对比（paired cluster bootstrap，1500）")
    emit("")
    uq = np.unique(g)
    cidx = [np.where(g == c)[0] for c in uq]
    CONTRASTS = [("B3a_real_PPI", "NC_noise_PPI", "B3a_real − NC_noise（判别核心）"),
                 ("B3a_real_PPI", "B_shuf_target", "B3a_real − B_shuf_target"),
                 ("NC_noise_PPI", "B2_dr", "NC_noise − B2_dr（噪声是否也「有效」）")]
    out = []
    for a, b, lab in CONTRASTS:
        PA = np.mean([preds[f"{a}_s{s}"].astype(float) for s in SEEDS
                      if f"{a}_s{s}" in preds], axis=0)
        PB = np.mean([preds[f"{b}_s{s}"].astype(float) for s in SEEDS
                      if f"{b}_s{s}" in preds], axis=0)
        m = np.isfinite(PA) & np.isfinite(PB)
        yy, pa, pb, gg = y[m], PA[m], PB[m], g[m]
        A0, B0 = mset(yy, pa, gg), mset(yy, pb, gg)
        names = ["macro_spearman", "eff10", "C_W", "overall_r", "rmse"]
        pt = {n: A0[n] - B0[n] for n in names}
        acc = {n: [] for n in names}
        rs2 = np.random.RandomState(42)
        for _ in range(N_BOOT):
            pick = rs2.randint(0, len(uq), len(uq))
            r_ = np.concatenate([cidx[j] for j in pick])
            rep = np.concatenate([np.full(len(cidx[j]), i)
                                  for i, j in enumerate(pick)])
            try:
                A = mset(y[r_], pa[r_], rep); Bx = mset(y[r_], pb[r_], rep)
            except Exception:  # noqa: BLE001
                continue
            for n in names:
                v = A[n] - Bx[n]
                if v == v:
                    acc[n].append(v)
        emit(f"### {lab}")
        emit("")
        emit("| 指标 | point | 95% CI | 排除 0? |")
        emit("|---|---:|---|---|")
        for n in names:
            arr = np.asarray(acc[n], float)
            if len(arr) < 20:
                lo = hi = np.nan; ex = False
            else:
                lo, hi = np.percentile(arr, [2.5, 97.5])
                ex = bool((lo > 0) or (hi < 0))
            emit(f"| Δ{n} | {pt[n]:+.4f} | [{lo:+.4f}, {hi:+.4f}] | "
                 f"{'✅ 排除 0' if ex else '❌ 含 0'} |")
            out.append({"contrast": lab, "metric": f"d_{n}", "point": pt[n],
                        "ci_low": float(lo), "ci_high": float(hi),
                        "n_valid": len(arr), "ci_excludes_zero": ex})
        emit("")
        pd.DataFrame(out).to_csv(PH8 / "noise_control_contrasts.csv",
                                 index=False, encoding="utf-8")

    # ---- 判定 ----
    d1 = [o for o in out if o["contrast"].startswith("B3a_real − NC_noise")
          and o["metric"] == "d_macro_spearman"]
    dn = [o for o in out if o["contrast"].startswith("NC_noise − B2_dr")
          and o["metric"] == "d_macro_spearman"]
    emit("## 8.4 机制判定（决定性）")
    emit("")
    if d1:
        r1 = d1[0]
        emit(f"- `B3a_real − NC_noise` 的 ΔMacro ρ = **{r1['point']:+.4f}** "
             f"[{r1['ci_low']:+.4f}, {r1['ci_high']:+.4f}]")
        if r1["ci_excludes_zero"] and r1["point"] > 0:
            emit("")
            emit("```")
            emit("NOISE_CONTROL_SEPARATED = YES")
            emit("CONTINUE_REPRESENTATION_ROUTE = CONDITIONAL_YES (需续跑至 >=80% 覆盖)")
            emit("```")
            emit("")
            emit("⇒ 真实 PPI 上下文**显著优于**同方差噪声 ⇒ 增益含真实 target 级信息。")
        else:
            emit("")
            emit("```")
            emit("NOISE_CONTROL_SEPARATED = NO")
            emit("CONTINUE_REPRESENTATION_ROUTE = NO")
            emit("```")
            emit("")
            emit("⇒ 真实 PPI 上下文**未能显著优于**同方差噪声 ⇒")
            emit("**B3a 的增益不能排除「噪声级伪随机排序」**（§14 的机制问题**未解决**）。")
            emit("⇒ 依 §25：**关闭表示路线**，进入 V7.3 physics / structure-informed ranking。")
    if dn:
        rn = dn[0]
        emit("")
        emit(f"- `NC_noise − B2_dr` 的 ΔMacro ρ = **{rn['point']:+.4f}** "
             f"[{rn['ci_low']:+.4f}, {rn['ci_high']:+.4f}] → "
             f"{'显著' if rn['ci_excludes_zero'] else '不显著'}")
        if rn["ci_excludes_zero"] and rn["point"] > 0:
            emit("  ⇒ **纯噪声特征本身也能显著提升 Macro ρ** ——")
            emit("  这直接证明该提升**不是 PPI 信息**，而是 LightGBM 在")
            emit("  target 内噪声上做伪随机排序的结果。")
    emit("")
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (PH8 / "noise_control.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {PH8 / 'noise_control_results.csv'} / "
          f"noise_control_contrasts.csv / noise_control.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
