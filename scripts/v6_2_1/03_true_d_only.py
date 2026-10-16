# -*- coding: utf-8 -*-
"""scripts/v6_2_1/03_true_d_only.py -- V6.2.1 PHASE 3: 真正的 M2 = D-only.

问题
----
V6.2 中 `M2` 与 `M0` 完全相同（`yhat["M2"][te] = p0`），
因此 `Gain_W(M2 − M0) = 0` **没有科研含义**。

本脚本实现**真正的** D-only 模型。

历史 D-only 定义（**复用，非重写**）
-----------------------------------
来源：`scripts/v5/02_baseline_repro.py`（V5 Gate A 复现，SHA256 已冻结）

```python
Dx = np.hstack([blk["phys_diff"], blk["Z_diff"]])      # 反对称块（D = x_mut − x_wt）
Dn, _ = C.standardize_train(Dx, train_mask)            # 只用 train 行拟合标准化
# forward/reverse 增广（V4/V5 一致）：
Xc = np.vstack([Daux[mtr], Daux_rev[mtr]])             # Daux_rev 的 D 部分取负
y  = np.r_[y[mtr], -y[mtr]]
```

**M2（D-only）** = 只用 `Dn`（**不含** DGaux 辅助列、**不含** `S` 对称块），
并保留同样的 **forward/reverse 增广** 语义。

| 模型 | 特征 | 增广 |
|---|---|---|
| **M0** pooled baseline | `[phys_diff, Z_diff]`（原始，未标准化） | fwd+rev |
| **M2** D-only（本脚本） | `standardize_train(D)` = `[phys_diff, Z_diff]` 标准化 | fwd+rev |
| **M1** DGaux | `[D, DGaux_OOF]` | fwd+rev（D 取负、aux 不取负） |

> M0 与 M2 的**唯一**差别是是否用 train 统计量标准化（`standardize_train`）。
> 这是历史 D-only 的确切定义（Gate A 的 `Dn`），因此 M2 ≠ M0 是**实质性**的。

只跑 L1 与 L4（任务书 §3.3 指定优先档），3 seeds。
"""

from __future__ import annotations

import gc
import hashlib
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v6_1"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import v5_common as C  # noqa: E402
import fastmetrics as FM  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

V62P2 = C.ROOT / "experiments" / "v6_2" / "phase2_split4"
PH3 = C.ROOT / "experiments" / "v6_2_1" / "phase3_d_only"
PH3.mkdir(parents=True, exist_ok=True)
OUT_NPZ = PH3 / "d_only_predictions.npz"
OUT_MET = PH3 / "d_only_metrics.csv"
OUT_GAIN = PH3 / "gain_w_bootstrap.csv"
OUT_MD = PH3 / "d_only_report.md"
FAIL = C.ROOT / "logs" / "v6_2_1" / "failures.log"
FAIL.parent.mkdir(parents=True, exist_ok=True)

TREES = 300          # 与 V6.2 完全一致（PHASE 2 已审计）
SEEDS = [42, 43, 44]
LEVELS = ["L1_random", "L4_family"]
N_BOOT = 5000
MIN_PS_N = 5
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def fail(where: str, exc: BaseException) -> None:
    import traceback
    with open(FAIL, "a", encoding="utf-8") as fh:
        fh.write(f"[{time.strftime('%H:%M:%S')}] [03_d_only] {where}: "
                 f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}\n")


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ---------------------------------------------------------------- 分解
def decompose(y, p, g):
    y = np.asarray(y, float); p = np.asarray(p, float); g = np.asarray(g)
    order = np.argsort(g, kind="mergesort")
    gs = g[order]
    uniq, _ = np.unique(gs, return_index=True)
    codes = np.searchsorted(uniq, gs)
    cnt = np.bincount(codes, minlength=len(uniq))
    sy = np.bincount(codes, weights=y[order], minlength=len(uniq))
    sp = np.bincount(codes, weights=p[order], minlength=len(uniq))
    By = np.empty(len(y)); Bp = np.empty(len(p))
    By[order] = sy[codes] / cnt[codes]
    Bp[order] = sp[codes] / cnt[codes]
    Wy, Wp = y - By, p - Bp

    def cov(a, b):
        return float(((a - a.mean()) * (b - b.mean())).mean())
    den = float(y.std()) * float(p.std())
    if den <= 0:
        return dict(r_full=np.nan, C_B=np.nan, C_W=np.nan)
    return dict(r_full=cov(y, p) / den, C_B=cov(By, Bp) / den,
                C_W=cov(Wy, Wp) / den)


def boot_paired(fn, arrays, g, n_boot=N_BOOT, seed=42):
    """配对 cluster bootstrap（抽 complex，同一索引同时算所有量）。"""
    uniq = np.unique(g)
    cidx = {u: np.where(g == u)[0] for u in uniq}
    rs = np.random.RandomState(seed)
    pt = fn(*arrays, g)
    vals = []
    for _ in range(n_boot):
        pick = rs.randint(0, len(uniq), len(uniq))
        rows = np.concatenate([cidx[uniq[j]] for j in pick])
        rep = np.concatenate([np.full(len(cidx[uniq[j]]), i)
                              for i, j in enumerate(pick)])
        v = fn(*[a[rows] for a in arrays], rep)
        if v == v:
            vals.append(v)
    a = np.asarray(vals, float)
    if len(a) < 20:
        return {"point": pt, "ci_low": np.nan, "ci_high": np.nan,
                "n_valid": len(a)}
    lo, hi = np.percentile(a, [2.5, 97.5])
    return {"point": float(pt), "ci_low": float(lo), "ci_high": float(hi),
            "n_valid": len(a), "p_gt0": float((a > 0).mean())}


def f_CW(y, p, g):
    return decompose(y, p, g)["C_W"]


def main() -> int:
    t0 = time.time()
    emit("# V6.2.1 PHASE 3 — 真正的 M2 = D-only")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- **历史 D-only 定义复用**：`scripts/v5/02_baseline_repro.py` 的 `Dn`")
    emit(f"  - SHA256: `{sha(C.ROOT / 'scripts' / 'v5' / '02_baseline_repro.py')[:32]}…`")
    emit(f"- 树数 **{TREES}**（与 V6.2 一致，PHASE 2 已审计）")
    emit(f"- 档: {LEVELS}；seeds: {SEEDS}")
    emit("")

    # ---------- 复用 V6.2 的 fold 生成（绝不重写 split 逻辑） ----------
    spec = importlib.util.spec_from_file_location(
        "s4", C.ROOT / "scripts" / "v6_2" / "02_split4_ddg.py")
    s4 = importlib.util.module_from_spec(spec)
    sys.modules["s4"] = s4
    spec.loader.exec_module(s4)
    emit(f"- 复用 V6.2 的 `level_units` / `make_folds`："
         f"`{sha(C.ROOT / 'scripts' / 'v6_2' / '02_split4_ddg.py')[:32]}…`")
    emit("")

    d = C.build_pair_features(verbose=False)
    can = pd.read_csv(V62P2 / "matched_analysis_set.csv", low_memory=False)
    rows = can["row"].to_numpy(int)
    blk = d["X_blocks"]
    Dtra = {
        "rows": rows,
        "y": d["y"][rows],
        "pdb": np.array([str(x) for x in d["groups"]])[rows],
        "fam": can["family_label"].astype(str).to_numpy(),
        "seq40": can["seq40_cluster"].astype(str).to_numpy(),
        "obs": can["obs_key"].astype(str).to_numpy(),
        "state": can["state_id"].astype(str).to_numpy(),
        "mut": can["mutation"].astype(str).to_numpy(),
        "all_rows": d,
    }
    Dx = np.hstack([blk["phys_diff"], blk["Z_diff"]])[rows].astype(np.float32)
    emit(f"- matched set: **{len(Dtra['y'])}** 行 / "
         f"{len(set(Dtra['pdb']))} 复合物")
    emit(f"- D（反对称块）维度: **{Dx.shape[1]}**"
         f"（phys_diff 60 + Z_diff 128）")
    emit("")

    # ---------- 训练 M2（D-only） ----------
    lgb = C.load_lgb()
    preds: dict[str, np.ndarray] = {}
    met_rows: list[dict] = []
    emit("## 3.4–3.5 训练与指标")
    emit("")
    emit("| 档 | seed | 模型 | r_full | C_B | C_W | Macro ρ | Kendall τ | Eff@10 |")
    emit("|---|---:|---|---:|---:|---:|---:|---:|---:|")
    z62 = np.load(V62P2 / "split4_predictions.npz", allow_pickle=True)
    for lv in LEVELS:
        units = s4.level_units(Dtra, lv)
        n_u = len(np.unique(units))
        for sd in SEEDS:
            try:
                folds = s4.make_folds(units, sd, 5)
                n = len(Dtra["y"])
                yhat = np.full(n, np.nan)
                for tr, te in folds:
                    # ---- 关键：standardize_train 只用 **本折 train 行** 拟合 ----
                    m_tr = np.zeros(n, dtype=bool)
                    m_tr[tr] = True
                    Dn, _sc = C.standardize_train(Dx, m_tr)
                    # forward/reverse 增广（历史语义：D 取负、y 取负）
                    Xtr = np.vstack([Dn[tr], -Dn[tr]])
                    ytr = np.r_[Dtra["y"][tr], -Dtra["y"][tr]]
                    m = lgb.LGBMRegressor(**{**C.LGB_PARAMS,
                                             "n_estimators": TREES,
                                             "random_state": sd})
                    m.fit(Xtr, ytr)
                    yhat[te] = m.predict(Dn[te])
                    assert len(np.intersect1d(tr, te)) == 0, "fold 泄漏"
                    del m
                    gc.collect()
                ok = np.isfinite(yhat) & np.isfinite(Dtra["y"])
                yy = Dtra["y"][ok]; gg = Dtra["pdb"][ok]; pp = yhat[ok]
                dd = decompose(yy, pp, gg)
                mr = FM.macro_spearman(yy, pp, gg, min_n=MIN_PS_N)
                kt = _macro_kendall(yy, pp, gg)
                e10 = FM.eff10(yy, pp, gg, n_rand=20)
                emit(f"| {lv} | {sd} | M2(D-only) | {dd['r_full']:+.4f} | "
                     f"{dd['C_B']:+.4f} | {dd['C_W']:+.4f} | {mr:+.4f} | "
                     f"{kt:+.4f} | {e10:+.4f} |")
                met_rows.append({"level": lv, "seed": sd, "model": "M2_D_only",
                                 "r_full": dd["r_full"], "C_B": dd["C_B"],
                                 "C_W": dd["C_W"], "macro_spearman": mr,
                                 "macro_kendall": kt, "eff10": e10,
                                 "n_rows": int(ok.sum()),
                                 "n_target": int(len(set(gg))),
                                 "status": "OK"})
                preds[f"{lv}_M2_s{sd}"] = yhat.astype(np.float32)
                preds[f"{lv}_y"] = Dtra["y"].astype(np.float32)
                preds[f"{lv}_pdb"] = gg.astype(str)
                pd.DataFrame(met_rows).to_csv(OUT_MET, index=False,
                                              encoding="utf-8")
            except MemoryError as exc:
                fail(f"{lv}/seed{sd} OOM", exc)
                emit(f"| {lv} | {sd} | M2 | OOM（跳过） | | | | | |")
                met_rows.append({"level": lv, "seed": sd, "model": "M2_D_only",
                                 "status": "OOM_SKIPPED"})
                gc.collect()
            except Exception as exc:  # noqa: BLE001
                fail(f"{lv}/seed{sd}", exc)
                emit(f"| {lv} | {sd} | M2 | FAILED {type(exc).__name__} | | | | | |")
                met_rows.append({"level": lv, "seed": sd, "model": "M2_D_only",
                                 "status": f"FAILED_{type(exc).__name__}"})
                gc.collect()
    emit("")

    # ---------- M0 / M1 指标（直接读 V6.2 保存的预测，不重训） ----------
    emit("## M0 / M1（复用 V6.2 已保存预测，**不重训**）")
    emit("")
    emit("| 档 | seed | 模型 | r_full | C_B | C_W | Macro ρ | Eff@10 |")
    emit("|---|---:|---|---:|---:|---:|---:|---:|")
    for lv in LEVELS:
        y = z62[f"{lv}_y"]
        g = np.array([str(x) for x in z62[f"{lv}_pdb"]])
        for mk in ("M0", "M1"):
            for sd in SEEDS:
                k = f"{lv}_{mk}_s{sd}"
                if k not in z62.files:
                    continue
                p = z62[k]
                m = np.isfinite(y) & np.isfinite(p)
                yy = y[m].astype(float); gg = g[m]; pp = p[m].astype(float)
                dd = decompose(yy, pp, gg)
                mr = FM.macro_spearman(yy, pp, gg, min_n=MIN_PS_N)
                e10 = FM.eff10(yy, pp, gg, n_rand=20)
                emit(f"| {lv} | {sd} | {mk} | {dd['r_full']:+.4f} | "
                     f"{dd['C_B']:+.4f} | {dd['C_W']:+.4f} | {mr:+.4f} | "
                     f"{e10:+.4f} |")
                met_rows.append({"level": lv, "seed": sd, "model": mk,
                                 "r_full": dd["r_full"], "C_B": dd["C_B"],
                                 "C_W": dd["C_W"], "macro_spearman": mr,
                                 "eff10": e10,
                                 "n_rows": int(m.sum()),
                                 "n_target": int(len(set(gg))),
                                 "status": "OK(REUSED)"})
    pd.DataFrame(met_rows).to_csv(OUT_MET, index=False, encoding="utf-8")
    emit("")

    # ---------- Gain_W bootstrap（配对） ----------
    emit("## 3.5 Gain_W（配对 cluster bootstrap，5000）")
    emit("")
    emit("| 档 | 对比 | point | 95% CI | 有效 | 排除 0? |")
    emit("|---|---|---:|---|---:|---|")
    gain_rows = []
    for lv in LEVELS:
        y = z62[f"{lv}_y"]
        g = np.array([str(x) for x in z62[f"{lv}_pdb"]])
        P = {}
        for mk in ("M0", "M1"):
            ks = [f"{lv}_{mk}_s{s}" for s in SEEDS if f"{lv}_{mk}_s{s}" in z62.files]
            if ks:
                P[mk] = np.mean([z62[k] for k in ks], axis=0)
        ks2 = [f"{lv}_M2_s{s}" for s in SEEDS if f"{lv}_M2_s{s}" in preds]
        if ks2:
            P["M2"] = np.mean([preds[k] for k in ks2], axis=0)
        msk = np.isfinite(y)
        for v in P.values():
            msk &= np.isfinite(v)
        yy = y[msk].astype(float); gg = g[msk]
        for lab, a, b in (("Gain_W(M2−M0) D-only", "M2", "M0"),
                          ("Gain_W(M1−M2) DGaux", "M1", "M2"),
                          ("Gain_W(M1−M0) total", "M1", "M0")):
            if a not in P or b not in P:
                continue
            pa, pb = P[a][msk].astype(float), P[b][msk].astype(float)
            r = boot_paired(
                lambda yy_, x_, z_, gg_: f_CW(yy_, x_, gg_) - f_CW(yy_, z_, gg_),
                [yy, pa, pb], gg)
            excl = (r["ci_low"] > 0) or (r["ci_high"] < 0)
            emit(f"| {lv} | {lab} | {r['point']:+.4f} | "
                 f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}] | {r['n_valid']} | "
                 f"{'✅ 排除 0' if excl else '❌ 含 0'} |")
            gain_rows.append({"level": lv, "contrast": lab, "point": r["point"],
                              "ci_low": r["ci_low"], "ci_high": r["ci_high"],
                              "n_valid": r["n_valid"],
                              "ci_excludes_zero": bool(excl),
                              "unit": "complex", "n_boot": N_BOOT})
    pd.DataFrame(gain_rows).to_csv(OUT_GAIN, index=False, encoding="utf-8")
    emit("")

    np.savez(OUT_NPZ, **preds)
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_NPZ}")
    print(f"写出 {OUT_MET}")
    print(f"写出 {OUT_GAIN}")
    print(f"写出 {OUT_MD}")
    return 0


def _macro_kendall(y, p, g, min_n=2):
    """Macro Kendall τ_b（含 tie 校正）。"""
    vals = []
    for c in np.unique(g):
        m = np.where(g == c)[0]
        if len(m) < min_n:
            continue
        a, b = y[m], p[m]
        n = len(a)
        conc = disc = ta = tb = 0
        for i in range(n):
            for j in range(i + 1, n):
                da = a[i] - a[j]; db = b[i] - b[j]
                if da == 0 and db == 0:
                    ta += 1; tb += 1
                elif da == 0:
                    ta += 1
                elif db == 0:
                    tb += 1
                elif (da > 0) == (db > 0):
                    conc += 1
                else:
                    disc += 1
        n0 = conc + disc
        den = np.sqrt((n0 + ta) * (n0 + tb))
        if den > 0:
            vals.append((conc - disc) / den)
    return float(np.mean(vals)) if vals else float("nan")


if __name__ == "__main__":
    sys.exit(main())
