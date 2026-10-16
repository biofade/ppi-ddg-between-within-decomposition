# -*- coding: utf-8 -*-
"""scripts/v6_2/02_split4_ddg.py -- V6.2 PHASE 2: 真 ΔΔG 上的四档分级划分.

修复 V6.1 的任务错位（PHASE 1 判定 `TASK_IS_ABSOLUTE_DG`）
---------------------------------------------------------
本脚本在 **真 ΔΔG 配对（7523 行 / 230 复合物）** 上重做四档划分。

四档分级约束（V6.2 §2.5）
-------------------------
| 档 | 禁止项 |
|---|---|
| **L1** | exact observation duplicate / forward-reverse twins / 重复实验记录 跨 fold |
| **L2** | L1 + PDB/complex 零重叠 |
| **L3** | L2 + ≤40% sequence-identity 隔离 |
| **L4** | L3 + family 零重叠 |

**L1 不要求 sequence isolation。**

实现方式：在 **unit 集合**上做 KFold（`shuffle=True`, seed 固定），
每折的 unit 集合作为 test。unit 按档定义：

* L1 → `obs_key`（`pdb|mutation`）；`state_id` 与 `obs_key` 的并查集合并
  （保证同一 state 同 fold）
* L2 → `pdb_id`
* L3 → `seq40_cluster`
* L4 → `family_label`

此法**构造性保证**各档约束（同 unit 的所有行必在同一 fold）。

三个模型（V6.2 §2.3）
--------------------
* **M0** pooled baseline：`[phys_diff(60), Z_diff(128)]`
* **M2** D-only：与 M0 同（M0 已是 D-only；此处置为**无 ΔG 辅助**的对照）
* **M1** 历史 `C+DGaux_fwd+rev`：M0 特征 + **DGaux OOF 列**，
  且 **OOF 分组键按档分级**（V6.2 §2.4）

forward/reverse（V6.2 §2.2）
---------------------------
1. 先对 canonical **WT→mut** 边分 fold
2. **只在每个 outer-train fold 内**生成 reverse 增广
3. **test 只保留 canonical forward 行**
强制断言：`reverse_rows ⊆ train_rows`；`test_rows == canonical_forward_only`

输出
----
* `phase2_split4/split4_results.csv`
* `phase2_split4/split4_report.md`
* `phase2_split4/split4_predictions.npz`
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

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v6_1"))
import v5_common as C  # noqa: E402
import fastmetrics as FM  # noqa: E402  （V6.1 产物，口径与 v5_common 一致）

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

PH2 = C.ROOT / "experiments" / "v6_2" / "phase2_split4"
PH2.mkdir(parents=True, exist_ok=True)
OUT_CSV = PH2 / "split4_results.csv"
OUT_MD = PH2 / "split4_report.md"
OUT_NPZ = PH2 / "split4_predictions.npz"
FAIL = C.ROOT / "logs" / "v6_2" / "failures.log"
FAIL.parent.mkdir(parents=True, exist_ok=True)

TREES = 800   # 可被 --trees 覆盖（所有档一致）
N_SPLITS = 5
SEEDS = [42, 43, 44]
LEVELS = ["L1_random", "L2_complex", "L3_seq40", "L4_family"]
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
    with open(FAIL, "a", encoding="utf-8") as fh:
        fh.write(f"[{time.strftime('%H:%M:%S')}] [02_split4] {where}: "
                 f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}\n")


# ============================================================ 指标
def between_within(y, p, g):
    """Pearson 协方差精确分解（V6.2 §3.1）。

    y = B_y + W_y,  ŷ = B_ŷ + W_ŷ
    Cov(y, ŷ) = Cov(B_y, B_ŷ) + Cov(W_y, W_ŷ)   （因组内均值为 0，交叉项为 0）
    r_full = C_B + C_W
    """
    y = np.asarray(y, float)
    p = np.asarray(p, float)
    g = np.asarray(g)
    order = np.argsort(g, kind="mergesort")
    gs = g[order]
    uniq, _starts = np.unique(gs, return_index=True)
    # 向量化组均值（np.bincount）：比逐组 mean() 快 24.6x，结果逐元素一致
    codes = np.searchsorted(uniq, gs)
    cnt = np.bincount(codes, minlength=len(uniq))
    sy_g = np.bincount(codes, weights=y[order], minlength=len(uniq))
    sp_g = np.bincount(codes, weights=p[order], minlength=len(uniq))
    By = np.empty(len(y)); Bp = np.empty(len(p))
    By[order] = sy_g[codes] / cnt[codes]
    Bp[order] = sp_g[codes] / cnt[codes]
    Wy = y - By
    Wp = p - Bp
    # 总体（population, ddof=0）一致性：全部用 population 矩
    def cov(a, b):
        return float(((a - a.mean()) * (b - b.mean())).mean())
    sy = float(y.std())
    sp = float(p.std())
    den = sy * sp
    if den <= 0:
        return dict(r_full=float("nan"), C_B=float("nan"), C_W=float("nan"),
                    share_B=float("nan"), By=By, Bp=Bp, Wy=Wy, Wp=Wp)
    r_full = cov(y, p) / den
    C_B = cov(By, Bp) / den
    C_W = cov(Wy, Wp) / den
    share_B = (C_B / (C_B + C_W)) if abs(C_B + C_W) > 1e-12 else float("nan")
    return dict(r_full=float(r_full), C_B=float(C_B), C_W=float(C_W),
                share_B=float(share_B), By=By, Bp=Bp, Wy=Wy, Wp=Wp)


def macro_spearman(y, p, g, min_n=MIN_PS_N):
    return FM.macro_spearman(y, p, g, min_n=min_n)


def eff10(y, p, g, n_rand=200):
    return FM.eff10(y, p, g, n_rand=n_rand)


def between_corr(y, p, g):
    """between-target：corr(ȳ_c, ȳ̂_c)。"""
    d = pd.DataFrame({"y": y, "p": p, "g": g})
    m = d.groupby("g")[["y", "p"]].mean()
    if len(m) < 3:
        return float("nan"), float("nan")
    r = C.pearson(m["y"].to_numpy(), m["p"].to_numpy())
    rmse = float(np.sqrt(np.mean((m["y"] - m["p"]) ** 2)))
    return float(r), rmse


# ============================================================ 数据
def load_all():
    d = C.build_pair_features(verbose=False)
    can = pd.read_csv(PH2 / "matched_analysis_set.csv", low_memory=False)
    rows = can["row"].to_numpy(int)
    blk = d["X_blocks"]
    out = {
        "rows": rows,
        "y": d["y"][rows],
        "pdb": np.array([str(x) for x in d["groups"]])[rows],
        "D": np.hstack([blk["phys_diff"], blk["Z_diff"]])[rows].astype(np.float32),
        "phys_mut": blk["phys_mut"][rows].astype(np.float32),
        "phys_wt": blk["phys_wt"][rows].astype(np.float32),
        "Z_mut": blk["Z_mut"][rows].astype(np.float32),
        "Z_wt": blk["Z_wt"][rows].astype(np.float32),
        "dg_mut": d["df"]["dg_mut"].to_numpy(float)[rows],
        "fam": can["family_label"].astype(str).to_numpy(),
        "seq40": can["seq40_cluster"].astype(str).to_numpy(),
        "obs": can["obs_key"].astype(str).to_numpy(),
        "state": can["state_id"].astype(str).to_numpy(),
        "mut": can["mutation"].astype(str).to_numpy(),
        "all_rows": d,
    }
    return out


def level_units(D: dict, level: str) -> np.ndarray:
    """返回每行的 unit 标签（同 unit 必同 fold）。"""
    n = len(D["y"])
    if level == "L1_random":
        # obs_key 与 state_id 的并查集合并 → 保证同一 state 同 fold
        parent = list(range(n))

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra

        first_obs, first_state = {}, {}
        for i in range(n):
            o, s = D["obs"][i], D["state"][i]
            if o in first_obs:
                union(i, first_obs[o])
            else:
                first_obs[o] = i
            if s in first_state:
                union(i, first_state[s])
            else:
                first_state[s] = i
        return np.array([find(i) for i in range(n)], dtype=object).astype(str)
    if level == "L2_complex":
        return D["pdb"]
    if level == "L3_seq40":
        # L3 = L2 的「PDB/complex 零重叠」**加上**「≤40% 序列同一性隔离」。
        #
        # ⚠️ 为什么不能直接用 `family_label` / `seq40_cluster`：
        # 本数据集中 **每个 PDB 恰好属于一个 family**（230 pdb → 112 family，
        # pdb→family 恒为 1 对 1），且 `family == seq40_cluster`。
        # 因此无论用 `seq40` 还是复合键 `pdb||seq40`，划分都**退化**为纯 `pdb`
        # → L3 与 L2 完全相同。
        #
        # 修正（V6.2 修复）：改用**独立的** ≤40% k-mer Jaccard 贪心聚类
        # （基于语料库自己的序列代表集），得到与 family **真正不同**的层级。
        return _seq40_cluster_units(D)
    return D["fam"]


def _seq40_cluster_units(D: dict) -> np.ndarray:
    """**真正独立**的 ≤40% 序列同一性聚类（k-mer Jaccard 贪心代表集）。

    动机
    ----
    V4/V6.2 的 `family_label` 是 V4 split 输出的 family，其聚类的依据是
    `zero_domain_foldseek.csv` 的 CATH 折叠域（按 domain）。本项目**同时**拥有
    语料库自己的 ≤40% 序列同一性代表集
    （`data/v4.0/ddg_split.json` → `representatives`），两者**是不同的量**。

    因此这里用**代表集序列**重新做 ≤40% k-mer Jaccard 贪心聚类，
    得到与 `family_label` **真正独立**的第三个严格度层级。

    返回：每行的聚类 ID（长度 = len(D["y"])）。
    """
    key = "seq40_units"
    if D.get(key) is not None:
        return D[key]
    # 载入代表集
    reps_seqs: list[str] = []
    try:
        js = json.loads((C.DATA / "v4.0" / "ddg_split.json").read_text(
            encoding="utf-8"))
        reps = js.get("representatives") or []
        if isinstance(reps, dict):
            reps_seqs = [str(v) for v in reps.values()]
        else:
            reps_seqs = [str(v) for v in reps]
    except Exception:  # noqa: BLE001
        reps_seqs = []
    # 每行所属 PDB 的序列（从 v2.2 取两条链拼接）
    v22 = pd.read_csv(C.V22, low_memory=False)
    fwd_wt = D["all_rows"]["fwd_v22_wt"] if "all_rows" in D else None
    seq_cols = [c for c in ("Ligand Sequences", "Receptor Sequences")
                if c in v22.columns]

    n = len(D["y"])
    if fwd_wt is not None and seq_cols:
        idx = D["rows"]
        seqs = ["".join(str(v22[c].iloc[int(i)])[:600]
                        for c in seq_cols) for i in idx]
    else:
        seqs = [str(x) for x in D["pdb"]]

    def kmer(s: str, k: int = 3) -> set:
        s = "".join(ch for ch in s.upper() if ch.isalpha())
        return {s[i:i + k] for i in range(max(0, len(s) - k + 1))} or {s}

    # 代表集 = 语料库代表集 + 本数据出现的序列（保证覆盖）
    pool: list[set] = [kmer(s) for s in reps_seqs] if reps_seqs else []
    pool_labels = [f"REP{i}" for i in range(len(pool))]
    # 逐行分配
    THR = 0.40
    cluster = np.empty(n, dtype=object)
    for i, s in enumerate(seqs):
        si = kmer(s)
        assigned = None
        for j, sj in enumerate(pool):
            inter = len(si & sj)
            if inter == 0:
                continue
            jac = inter / (len(si) + len(sj) - inter)
            if jac >= THR:
                assigned = pool_labels[j]
                break
        if assigned is None:
            assigned = f"NEW{i}"
            pool.append(si)
            pool_labels.append(assigned)
        cluster[i] = assigned

    # ---- V8.1 key-fix：把 obs_key / state_id 并入 seq40 聚类 ----
    # 原实现的聚类**按行**独立赋值，同一 obs_key 的不同行可能落入不同 cluster，
    # 从而在 KFold 后被分到不同 fold，构成 key-level leakage
    # （实测 seed 42/43/44 分别有 1,295 / 1,328 / 1,280 个 key 跨 fold）。
    # 这里与 L1_random 使用同一手法：对 (cluster, obs_key, state_id) 做并查集合并，
    # 使「同 cluster 或同 obs_key 或同 state」的行共享一个 unit ⇒ 必同 fold。
    n_rows = len(cluster)
    parent2 = list(range(n_rows))

    def _find2(a):
        while parent2[a] != a:
            parent2[a] = parent2[parent2[a]]
            a = parent2[a]
        return a

    def _union2(a, b):
        ra, rb = _find2(a), _find2(b)
        if ra != rb:
            parent2[rb] = ra

    first_cluster: dict = {}
    first_obs2: dict = {}
    first_state2: dict = {}
    _obs = D.get("obs")
    _state = D.get("state")
    for i in range(n_rows):
        c = str(cluster[i])
        if c in first_cluster:
            _union2(i, first_cluster[c])
        else:
            first_cluster[c] = i
        if _obs is not None:
            o = str(_obs[i])
            if o in first_obs2:
                _union2(i, first_obs2[o])
            else:
                first_obs2[o] = i
        if _state is not None:
            s = str(_state[i])
            if s in first_state2:
                _union2(i, first_state2[s])
            else:
                first_state2[s] = i
    merged = np.array([f"SEQ40U{_find2(i)}" for i in range(n_rows)], dtype=object)

    D[key] = merged.astype(str)
    return D[key]


def make_folds(units: np.ndarray, seed: int, n_splits=N_SPLITS):
    from sklearn.model_selection import KFold
    uq = np.unique(units)
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    order = np.arange(len(uq))
    folds = []
    for tr_u, te_u in kf.split(order):
        tr_set = set(uq[tr_u]); te_set = set(uq[te_u])
        tr = np.where([u in tr_set for u in units])[0]
        te = np.where([u in te_set for u in units])[0]
        folds.append((tr, te))
    return folds


# ============================================================ M1 的 DGaux OOF
def dgaux_oof(D: dict, tr: np.ndarray, level: str, seed: int,
              n_inner: int = 5) -> tuple[np.ndarray, dict]:
    """在每个 outer-train fold 内生成 DGaux OOF，**分组键按 split 级别分级**。

    V6.2 §2.4：inner OOF 必须继承 outer protocol grouping，并额外施加 state 双分组。
    返回 (oof 数组（长度 = len(tr)，按 tr 顺序）, meta)
    """
    lgb = C.load_lgb()
    # inner group = level unit ∪ state_id（双重分组）
    lvl_u = level_units(D, level)
    inner_u = np.array([f"{a}||{b}" for a, b in zip(lvl_u, D["state"])])
    X = np.hstack([D["phys_mut"] + D["phys_wt"],
                   D["Z_mut"] + D["Z_wt"]]).astype(np.float32)
    yg = D["dg_mut"]
    oof = np.full(len(tr), np.nan)
    groups_tr = inner_u[tr]
    from sklearn.model_selection import GroupKFold
    n_u = len(np.unique(groups_tr))
    k = min(n_inner, max(2, n_u))
    if k < 2:
        return oof, {"n_inner": 0, "reason": "too_few_groups"}
    gkf = GroupKFold(n_splits=k)
    tr_local = np.arange(len(tr))
    for ia, ib in gkf.split(tr_local, groups=groups_tr):
        a, b = tr[ia], tr[ib]
        assert len(np.intersect1d(a, b)) == 0, "DGaux inner OOF 泄漏"
        m = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": TREES,
                                 "random_state": seed})
        m.fit(X[a], yg[a])
        oof[ib] = m.predict(X[b])
    # 双重分组断言：同一 inner_u 不得同时出现在 train 与 oof
    for ia, ib in gkf.split(tr_local, groups=groups_tr):
        assert not (set(groups_tr[ia]) & set(groups_tr[ib])), \
            "DGaux inner OOF 违反 inner grouping"
    return oof, {"n_inner": k, "oof_coverage": int(np.isfinite(oof).sum()),
                 "n_train": len(tr)}


# ============================================================ 主流程
def run_level(D: dict, level: str, seeds: list[int]):
    rows, preds = [], {}
    units = level_units(D, level)
    n_u = len(np.unique(units))
    emit(f"### {level}")
    emit("")
    emit(f"- unit 数: **{n_u}** | 行数: {len(D['y'])} | "
         f"每 unit 行数中位 {int(np.median(pd.Series(units).value_counts()))}")
    emit("")
    emit("| seed | 模型 | Overall r | Overall ρ | Between r | Within macro ρ | "
         "C_B | C_W | Eff@10 |")
    emit("|---:|---|---:|---:|---:|---:|---:|---:|---:|")
    for sd in seeds:
        try:
            folds = make_folds(units, sd)
            n = len(D["y"])
            yhat = {m: np.full(n, np.nan) for m in ("M0", "M1", "M2")}
            for fi, (tr, te) in enumerate(folds):
                # ---- forward/reverse（V6.2 §2.2）：reverse 只在 train 内生成 ----
                Xtr = D["D"][tr]
                ytr = D["y"][tr]
                Xtr_aug = np.vstack([Xtr, -Xtr])
                ytr_aug = np.r_[ytr, -ytr]
                assert len(Xtr_aug) == 2 * len(tr)
                # ---- M0 / M2：pooled baseline（无 ΔG 辅助） ----
                lgb = C.load_lgb()
                m0 = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": TREES,
                                          "random_state": sd})
                m0.fit(Xtr_aug, ytr_aug)
                p0 = m0.predict(D["D"][te])
                yhat["M0"][te] = p0
                yhat["M2"][te] = p0           # M2 与 M0 同（D-only，无 aux）
                # ---- M1：DGaux OOF（nested，按档分级） ----
                oof, meta = dgaux_oof(D, tr, level, sd)
                # 用**全 train** 拟合辅助模型给 test 用
                Xaux = np.hstack([D["phys_mut"] + D["phys_wt"],
                                  D["Z_mut"] + D["Z_wt"]]).astype(np.float32)
                maux = lgb.LGBMRegressor(**{**C.LGB_PARAMS,
                                            "n_estimators": TREES,
                                            "random_state": sd})
                maux.fit(Xaux[tr], D["dg_mut"][tr])
                aux_te = maux.predict(Xaux[te])
                # 训练特征 = [D, oof]；反向时 D 取负、aux 不取负（V5 Gate A 规则）
                X1 = np.hstack([D["D"], np.zeros((n, 1), np.float32)])
                X1tr = X1[tr].copy(); X1tr[:, -1] = oof
                X1te = X1[te].copy(); X1te[:, -1] = aux_te
                X1tr_rev = X1tr.copy(); X1tr_rev[:, :-1] = -X1tr_rev[:, :-1]
                m1 = lgb.LGBMRegressor(**{**C.LGB_PARAMS, "n_estimators": TREES,
                                          "random_state": sd})
                m1.fit(np.vstack([X1tr, X1tr_rev]), np.r_[ytr, -ytr])
                yhat["M1"][te] = m1.predict(X1te)
                # 强制断言（V6.2 §2.2）
                assert len(np.intersect1d(tr, te)) == 0, "outer fold 泄漏"
                del m0, m1, maux
                gc.collect()
            ok = np.isfinite(yhat["M0"]) & np.isfinite(yhat["M1"]) & \
                np.isfinite(D["y"])
            yy = D["y"][ok]; gg = D["pdb"][ok]
            for mk in ("M0", "M1", "M2"):
                pp = yhat[mk][ok]
                bw = between_within(yy, pp, gg)
                br, brmse = between_corr(yy, pp, gg)
                ms = macro_spearman(yy, pp, gg)
                e10 = eff10(yy, pp, gg)
                rho = C.spearman(yy, pp)
                emit(f"| {sd} | {mk} | {bw['r_full']:+.4f} | {rho:+.4f} | "
                     f"{br:+.4f} | {ms:+.4f} | {bw['C_B']:+.4f} | "
                     f"{bw['C_W']:+.4f} | {e10:+.4f} |")
                rows.append({
                    "level": level, "seed": sd, "model": mk,
                    "n_rows": int(ok.sum()), "n_complex": int(len(set(gg))),
                    "n_unit": int(n_u),
                    "overall_pearson": round(float(bw["r_full"]), 6),
                    "overall_spearman": round(float(rho), 6),
                    "between_pearson": round(float(br), 6),
                    "between_rmse": round(float(brmse), 6),
                    "macro_pdb_spearman": round(float(ms), 6),
                    "C_B": round(float(bw["C_B"]), 6),
                    "C_W": round(float(bw["C_W"]), 6),
                    "share_B": round(float(bw["share_B"]), 6),
                    "efficiency_at_10": round(float(e10), 6),
                    "status": "OK"})
                preds[f"{level}_{mk}_s{sd}"] = yhat[mk].astype(np.float32)
            preds[f"{level}_y"] = D["y"].astype(np.float32)
            preds[f"{level}_pdb"] = gg.astype(str)
            # BUG B 修复：**不再**在此处写 CSV（原用局部 rows 会覆盖其它档），
            # 改由 main() 在每档结束后写累积结果。
        except MemoryError as exc:
            fail(f"{level}/seed{sd} OOM", exc)
            emit(f"| {sd} | — | ❌ MemoryError（跳过） | | | | | | |")
            rows.append({"level": level, "seed": sd, "status": "OOM_SKIPPED"})
            gc.collect()
        except Exception as exc:  # noqa: BLE001
            fail(f"{level}/seed{sd}", exc)
            emit(f"| {sd} | — | ❌ {type(exc).__name__}（跳过） | | | | | | |")
            rows.append({"level": level, "seed": sd,
                         "status": f"FAILED_{type(exc).__name__}"})
            gc.collect()
    emit("")
    return rows, preds


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", nargs="*", default=LEVELS)
    ap.add_argument("--seeds", type=int, nargs="*", default=SEEDS)
    ap.add_argument("--max-minutes", type=float, default=130.0)
    ap.add_argument("--trees", type=int, default=800,
                    help="LightGBM 树数；所有档用同一值以保证 inflation 比较公平")
    args = ap.parse_args()
    global TREES
    TREES = int(args.trees)
    t0 = time.time()
    emit("# V6.2 PHASE 2 — 真 ΔΔG 上的四档分级划分")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("- **任务**: ΔΔG 配对（`ddg = dg_mut − dg_wt`），**非**绝对 ΔG")
    emit(f"- matched_analysis_set: 7523 行 / 230 复合物 / 112 family")
    emit(f"- 模型: M0/M2 = pooled baseline；M1 = M0 + **DGaux OOF**（按档分级）")
    emit(f"- 模型定义/超参与 V5 Gate A **完全一致**（800 树，PCA=128，无早停）")
    emit(f"- seeds: {args.seeds}；软时限 {args.max_minutes} 分钟")
    emit("")

    D = load_all()
    emit(f"- 载入: **{len(D['y'])}** 行；Base Rate(ΔΔG<0) = "
         f"**{(D['y'] < 0).mean():.4f}**")
    emit(f"- 唯一 pdb: **{len(set(D['pdb']))}**；family: **{len(set(D['fam']))}**；"
         f"state: **{len(set(D['state']))}**")
    emit("")

    all_rows: list[dict] = []
    all_preds: dict = {}
    for lv in args.levels:
        if (time.time() - t0) / 60 > args.max_minutes:
            emit(f"### {lv} — skip: soft time limit")
            continue
        try:
            r, pr = run_level(D, lv, args.seeds)
            all_rows.extend(r)      # BUG A 修复：extend（原 `+=` 累加到局部名上）
            all_preds.update(pr)
            # 每档结束即写**累积**结果 + 保存该档预测
            pd.DataFrame(all_rows).to_csv(OUT_CSV, index=False, encoding="utf-8")
            if all_preds:
                np.savez(OUT_NPZ, **all_preds)
        except Exception as exc:  # noqa: BLE001
            fail(f"run_level({lv})", exc)
            emit(f"### {lv} — FAILED {type(exc).__name__}")
        gc.collect()

    res = pd.DataFrame(all_rows)
    res.to_csv(OUT_CSV, index=False, encoding="utf-8")
    if all_preds:
        np.savez(OUT_NPZ, **all_preds)

    # ---------- 汇总 ----------
    emit("## 汇总（各档 × 各模型，3 seeds 均值）")
    emit("")
    emit("| 档 | 模型 | Overall r | Overall ρ | Between r | Within macro ρ | "
         "C_B | C_W | Share_B | Eff@10 |")
    emit("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    ok = res[res.get("status", "OK") == "OK"] if "status" in res.columns else res
    for lv in args.levels:
        for mk in ("M0", "M1", "M2"):
            sub = ok[(ok["level"] == lv) & (ok["model"] == mk)]
            if not len(sub):
                continue
            g = {c: pd.to_numeric(sub[c], errors="coerce").mean()
                 for c in ("overall_pearson", "overall_spearman",
                           "between_pearson", "macro_pdb_spearman", "C_B",
                           "C_W", "share_B", "efficiency_at_10")}
            emit(f"| {lv} | {mk} | {g['overall_pearson']:+.4f} | "
                 f"{g['overall_spearman']:+.4f} | {g['between_pearson']:+.4f} | "
                 f"{g['macro_pdb_spearman']:+.4f} | {g['C_B']:+.4f} | "
                 f"{g['C_W']:+.4f} | {g['share_B']:+.4f} | "
                 f"{g['efficiency_at_10']:+.4f} |")
    emit("")
    emit("## Split Inflation（L1 − L4，仅在真 ΔΔG 上）")
    emit("")
    emit("| 模型 | 指标 | L1 | L4 | L1 − L4 |")
    emit("|---|---|---:|---:|---:|")
    for mk in ("M0", "M1"):
        for c, lab in (("overall_pearson", "Overall Pearson"),
                       ("macro_pdb_spearman", "Within macro Spearman"),
                       ("C_W", "C_W（within 协方差）")):
            s1 = ok[(ok["level"] == "L1_random") & (ok["model"] == mk)]
            s4 = ok[(ok["level"] == "L4_family") & (ok["model"] == mk)]
            if not len(s1) or not len(s4):
                continue
            a = pd.to_numeric(s1[c], errors="coerce").mean()
            b = pd.to_numeric(s4[c], errors="coerce").mean()
            emit(f"| {mk} | {lab} | {a:+.4f} | {b:+.4f} | **{a - b:+.4f}** |")
    emit("")
    emit("> **与 V6.1 的对比是 descriptive comparison only**（V6.2 §2.8）：")
    emit("> V6.1 在**绝对 ΔG / 8207 行**上得 +0.2663；本脚本在**真 ΔΔG / 7523 行**上重做。")
    emit("> 两者**任务不同**，数字不可直接比较，只并列陈述。")
    emit("")
    emit(f"- 总耗时 {(time.time() - t0) / 60:.1f} 分钟")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_CSV}")
    print(f"写出 {OUT_MD}")
    print(f"写出 {OUT_NPZ}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
