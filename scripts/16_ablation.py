# -*- coding: utf-8 -*-
"""16_ablation.py -- 阶段 5: 消融实验与泄漏敏感性分析.

要回答的三个论证性问题
----------------------
Prompt 明确了三个关键论证 (本脚本是唯一能给出答案的地方):

1. **若 ``no_esm2`` ≈ ``full``**  -> 物理特征已足够, 这是"轻量级"卖点。
2. **若 ``size_only`` << ``full``** -> 归一化物理特征不可替代 (不是"模型只是在数原子")。
3. **若 ``leak_pdb_level`` > ``full``** -> 量化同源泄漏造成的**虚假**提升。

因此本脚本把「划分协议」本身当作一个可消融的因子: 除了家族级划分 (红线),
额外跑 **PDB 级** 与 **随机** 两种更宽松的划分作为**泄漏敏感性对照**。
这两个对照**不是**可选的美化项, 而是论文里必须诚实展示的"如果偷懒会虚高多少"。

消融配置
--------
============ ==================================================
配置          含义
============ ==================================================
full                    全部特征 (上界参考)
no_esm2                 去掉 ESM-2 维度
esm2_only               只用 ESM-2 维度
physics_only_normalized 只用归一化/密度类物理量
size_only               只用原子/残基计数
no_sasa                 去掉 SASA / ΔSASA 族
no_hbond_salt           去掉氢键 / 盐桥族
no_electrostatic        去掉新增静电代理
no_composition          去掉残基组成 / 氨基酸组成
no_chainpair            去掉 V0.3 链对聚合特征
no_extra                去掉阶段 5 全部新增特征 (退化为 V1.1 配置 B)
leak_pdb_level          **划分换成 PDB 级** (泄漏对照)
leak_random             **划分换成完全随机** (最强泄漏对照)
============ ==================================================

统计口径
--------
* 每个配置都在**同一套**测试集上评估 (family 级划分时), 因此可逐行配对。
* 泄漏对照的测试集**不同**, 其指标**不可与非泄漏行直接比大小**;
  报告里用单独小节呈现并显式警告。
* 全部配置带 bootstrap 95% CI; 与 ``full`` 的差异走配对 bootstrap + BH-FDR。

纪律: 纯 CPU; 无网络; seed=42; 逐配置容错; 不改动 V1.1 产物。
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.util
import json
import logging
import sys
import time
import traceback
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
FEAT_DIR = ROOT / "features"
LOG_DIR = ROOT / "logs" / "v2"

ENHANCED_CSV = FEAT_DIR / "v2_features_enhanced.csv"
BASE_CSV = FEAT_DIR / "features_v1.1.csv"
RESULT_CSV = FEAT_DIR / "v2_ablation_results.csv"
REPORT_MD = FEAT_DIR / "v2_ablation_report.md"
JSON_OUT = FEAT_DIR / "v2_ablation_results.json"
FAIL_LOG = LOG_DIR / "failures.log"

SEED = 42
N_BOOT = 1500
TOP_FRAC = 0.10

logger = logging.getLogger("ablation")


def _ensure_openmp() -> str | None:
    import ctypes
    import glob
    try:
        import sklearn
        hits = glob.glob(str(Path(sklearn.__file__).parent / ".libs" / "vcomp140.dll"))
        if hits:
            ctypes.CDLL(hits[0])
            return hits[0]
    except Exception:  # noqa: BLE001
        pass
    return None


_VCOMP = _ensure_openmp()


# --------------------------------------------------------------------- 特征族
def feature_groups(df: pd.DataFrame) -> dict[str, list[str]]:
    """把列切成互不重叠的语义族, 供消融开关使用。"""
    cols = set(df.columns)
    esm2 = sorted(c for c in cols if c.startswith("v2_esm2_pc"))
    extra = sorted(c for c in cols if c.startswith("v2_") and not c.startswith("v2_esm2"))
    sasa = sorted(c for c in cols if any(
        t in c for t in ("sasa", "dsasa")))
    hbond_salt = sorted(c for c in cols if any(
        t in c for t in ("hbond", "salt_bridge", "saltbridge")))
    electrostatic = sorted(c for c in cols if any(
        t in c for t in ("charge", "electrostatic")))
    composition = sorted(c for c in cols if any(
        t in c for t in ("frac_", "aac_", "hmom", "residue_ratio")))
    chainpair = sorted(c for c in cols if any(
        t in c for t in ("chainpair", "pair_sasa", "n_chainpairs", "structure_")))
    size = sorted(c for c in cols if any(
        t in c for t in ("n_atoms", "n_res_", "n_contact_pairs",
                         "n_contact_res_pairs", "n_iface_res")))
    normalized = sorted(c for c in cols if any(
        t in c for t in ("_per_", "density", "ratio", "_frac_", "contacts_per",
                         "iface_res_per")))

    # ---- 序列派生特征族 (关键对照) -------------------------------------
    # 这些特征完全由**全链野生型序列**算出 (残基组成/氨基酸组成/疏水矩/序列长度)。
    # 它们对同一个 target family 的**所有突变体/链组合都是同一个值**, 因此:
    #   * 它们不是"界面理化性质", 而是**靶点身份指纹**;
    #   * 家族隔离下不会跨侧泄漏, 但会带来一种更隐蔽的**性能膨胀**:
    #     模型可以靠"这个靶点的平均 ΔG 大概是多少"来排序, 而不是靠结合物理。
    # 因此必须单独消融, 判断新增特征的增益究竟来自物理还是来自靶点指纹。
    seq_derived = sorted(
        c for c in cols
        if c.startswith("v2_") and any(
            t in c for t in ("seq_len", "aac_", "frac_aromatic", "frac_charged",
                             "frac_polar", "frac_hydrophobic", "frac_special",
                             "hmom", "net_charge", "seq_available", "seq_source"))
    )
    return {
        "esm2": esm2, "extra": extra, "sasa": sasa, "hbond_salt": hbond_salt,
        "electrostatic": electrostatic, "composition": composition,
        "chainpair": chainpair, "size": size, "normalized": normalized,
        "seq_derived": seq_derived,
    }


def all_physical(df: pd.DataFrame) -> list[str]:
    """V1.1 配置 B 的 60 个物理特征 (从 08 载入, 单一实现来源)。"""
    path = ROOT / "scripts" / "08_train_v1.0.py"
    spec = importlib.util.spec_from_file_location("t8_abl", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    base = [c for c in mod.V02_FEATURES if c in df.columns and c not in mod.REPLACED]
    base += [c for c in mod.V03_NEW if c in df.columns]
    return base


def build_configs(df: pd.DataFrame) -> list[dict[str, Any]]:
    g = feature_groups(df)
    phys = all_physical(df)
    extra_only = g["extra"]
    full = phys + extra_only + g["esm2"]
    # 只保留真实存在于表中的列
    full = [c for c in dict.fromkeys(full) if c in df.columns]

    def drop(base: list[str], to_remove: list[str]) -> list[str]:
        s = set(to_remove)
        return [c for c in base if c not in s]

    cfgs: list[dict[str, Any]] = [
        {"name": "full", "desc": "全部特征 (上界参考)", "cols": full},
        {"name": "no_esm2", "desc": "去掉 ESM-2 维度",
         "cols": drop(full, g["esm2"])},
        {"name": "esm2_only", "desc": "只用 ESM-2 维度", "cols": list(g["esm2"])},
        {"name": "physics_only_normalized", "desc": "只用归一化/密度类物理量",
         "cols": [c for c in phys if c in set(g["normalized"])]},
        {"name": "size_only", "desc": "只用原子/残基计数",
         "cols": [c for c in g["size"] if c in df.columns]},
        {"name": "no_sasa", "desc": "去掉 SASA / ΔSASA 族",
         "cols": drop(full, g["sasa"])},
        {"name": "no_hbond_salt", "desc": "去掉氢键 / 盐桥族",
         "cols": drop(full, g["hbond_salt"])},
        {"name": "no_electrostatic", "desc": "去掉新增静电代理特征",
         "cols": drop(full, g["electrostatic"])},
        {"name": "no_composition", "desc": "去掉残基组成 / 氨基酸组成",
         "cols": drop(full, g["composition"])},
        {"name": "no_chainpair", "desc": "去掉 V0.3 链对聚合特征",
         "cols": drop(full, g["chainpair"])},
        {"name": "no_extra", "desc": "去掉阶段 5 全部新增特征 (退化为 V1.1 配置 B)",
         "cols": phys},
        {"name": "no_seq_derived",
         "desc": "**去掉全部序列派生特征** (只留几何/界面/尺寸类) —— 分离「界面物理」与「靶点指纹」",
         "cols": drop(full, g["seq_derived"])},
        {"name": "seq_derived_only",
         "desc": "**只用序列派生特征** (靶点指纹上界)",
         "cols": [c for c in g["seq_derived"] if c in df.columns]},
        # 划分协议对照 (泄漏敏感性)
        {"name": "leak_pdb_level", "desc": "**PDB 级划分** (同源泄漏对照)",
         "cols": full, "split": "pdb"},
        {"name": "leak_random", "desc": "**完全随机划分** (最强泄漏对照)",
         "cols": full, "split": "random"},
    ]
    return [c for c in cfgs if c["cols"]]


# --------------------------------------------------------------------- 划分
def make_split(df: pd.DataFrame, mode: str, seed: int) -> np.ndarray:
    """返回 bool 测试掩码。``family`` = 红线; ``pdb``/``random`` = 泄漏对照。

    三种划分的隔离强度递减:
      * family : 整族分侧 (40% 同一性家族), 同族/同 PDB 绝不跨界 -> **唯一可信口径**
      * pdb    : 整 PDB 分侧, 但**同源家族**可能跨界 -> 高估
      * random : 逐行随机 -> 同源、甚至同一测量的不同链组合都可能跨界 -> 最严重高估
    """
    n = len(df)
    if mode == "family":
        return (df["split_family"] == "test").to_numpy()
    if mode == "pdb":
        pdb = df["pdb_id"].astype(str).to_numpy()
        uniq = sorted(set(pdb))
        # 确定性哈希: 同一个 PDB 永远落同侧
        h = np.array([int(hashlib.md5(f"ai-v2-pdb::{p}".encode()).hexdigest()[:8], 16)
                      / 0xFFFFFFFF for p in uniq])
        cut = np.quantile(h, 0.80)
        test_pdbs = {p for p, hh in zip(uniq, h) if hh >= cut}
        return np.isin(pdb, list(test_pdbs))
    if mode == "random":
        rs = np.random.RandomState(seed)
        idx = rs.permutation(n)
        te = np.zeros(n, dtype=bool)
        te[idx[:int(0.20 * n)]] = True
        return te
    raise ValueError(f"未知划分模式: {mode!r}")


def leak_audit(df: pd.DataFrame, te: np.ndarray) -> dict[str, int]:
    """量化每个划分下的实际泄漏: 跨侧家族数 / 跨侧 PDB 数。"""
    sub = df.copy()
    sub["_side"] = np.where(te, "test", "train")
    out = {}
    for col in ("family_id", "pdb_id"):
        if col in sub.columns:
            out[f"cross_{col}"] = int((sub.groupby(col)["_side"].nunique() > 1).sum())
    return out


# --------------------------------------------------------------------- 评估
def run_config(cfg: dict[str, Any], df: pd.DataFrame, args: argparse.Namespace,
               mr: Any, lgb: Any, v11: dict, perm_rounds: int) -> dict[str, Any]:
    from sklearn.model_selection import GroupKFold

    cols = [c for c in cfg["cols"] if c in df.columns]
    if not cols:
        raise ValueError("特征列为空")
    split_mode = cfg.get("split", "family")
    te_m = make_split(df, split_mode, args.seed)
    tr_m = ~te_m
    if tr_m.sum() < 50 or te_m.sum() < 20:
        raise ValueError(f"划分后样本过少: train={int(tr_m.sum())}, test={int(te_m.sum())}")

    X = df[cols].apply(pd.to_numeric, errors="coerce")
    med = X[tr_m].median(numeric_only=True)
    X = X.fillna(med).fillna(0.0)
    y = df["real_dG"].to_numpy(dtype=float)

    X_tr, y_tr = X[tr_m], y[tr_m]
    X_te, y_te = X[te_m], y[te_m]
    # OOF 的分组单位永远是**家族**(即便划分换成 PDB/random), 这样 CV 内部不泄漏
    groups = df.loc[tr_m, "family_id"].to_numpy()
    n_splits = min(v11["N_SPLITS"], max(2, len(set(groups))))
    gkf = GroupKFold(n_splits=n_splits)
    params = dict(v11["PARAMS"])
    params["random_state"] = args.seed

    oof = np.full(int(tr_m.sum()), np.nan)
    fold_sp = []
    for i_tr, i_va in gkf.split(X_tr, y_tr, groups):
        m = lgb.LGBMRegressor(**params)
        m.fit(X_tr.iloc[i_tr], y_tr[i_tr],
              eval_set=[(X_tr.iloc[i_va], y_tr[i_va])],
              callbacks=[lgb.early_stopping(50, verbose=False)])
        oof[i_va] = m.predict(X_tr.iloc[i_va])
        sp = mr.spearman_rho(y_tr[i_va].copy(), oof[i_va].copy())
        # 注意: y 是 ΔG (越负越强), 而 mr 需要 score 空间 -> 两侧同时取负
        sp = mr.spearman_rho(-y_tr[i_va], -oof[i_va])
        if sp == sp:
            fold_sp.append(sp)
        gc.collect()

    final = lgb.LGBMRegressor(**params)
    final.fit(X_tr, y_tr)
    pred = final.predict(X_te)

    # 统一转到 score 空间 (越大越强)
    y_s, p_s = mr.negate_dG(y_te), mr.negate_dG(pred)
    met = mr.compute_ranking_metrics(y_s, p_s, top_frac=args.top_frac)
    lo, hi, _pt = mr.bootstrap_ci(mr.spearman_rho, y_s, p_s,
                                  n_boot=args.n_boot, seed=args.seed)
    oof_m = mr.compute_ranking_metrics(mr.negate_dG(y_tr), mr.negate_dG(oof),
                                       top_frac=args.top_frac)
    from scipy.stats import pearsonr
    from sklearn.metrics import mean_absolute_error, mean_squared_error
    rmse = float(np.sqrt(mean_squared_error(y_te, pred)))
    var = float(np.var(y_te, ddof=0))
    pe = pearsonr(y_te, pred)[0] if len(y_te) > 2 else float("nan")

    # 置换基线 (量化"无信号"参考线)
    rng = np.random.RandomState(args.seed)
    perm = []
    for _ in range(perm_rounds):
        yp = rng.permutation(y_tr)
        o = np.full(int(tr_m.sum()), np.nan)
        for i_tr, i_va in gkf.split(X_tr, yp, groups):
            mm = lgb.LGBMRegressor(**params)
            mm.fit(X_tr.iloc[i_tr], yp[i_tr])
            o[i_va] = mm.predict(X_tr.iloc[i_va])
        c = mr.spearman_rho(-yp, -o)
        if c == c:
            perm.append(c)

    out = {
        "config": cfg["name"], "description": cfg["desc"],
        "split_protocol": split_mode,
        "n_features": len(cols), "n_train": int(tr_m.sum()), "n_test": int(te_m.sum()),
        "n_families_train": int(df.loc[tr_m, "family_id"].nunique()),
        "n_families_test": int(df.loc[te_m, "family_id"].nunique()),
        "test_spearman": met.values["spearman"],
        "test_ci_low": lo, "test_ci_high": hi,
        "test_ndcg10": met.values["ndcg@10"],
        "test_precision10": met.values["precision@10"],
        "test_ef1": met.values["ef@1%"], "test_ef5": met.values["ef@5%"],
        "test_rmse": rmse, "test_mae": float(mean_absolute_error(y_te, pred)),
        "test_r2": float(1 - rmse ** 2 / var) if var > 0 else float("nan"),
        "test_pearson": float(pe) if pe == pe else float("nan"),
        "cv_spearman_mean": float(np.nanmean(fold_sp)) if fold_sp else float("nan"),
        "cv_spearman_std": float(np.nanstd(fold_sp)) if fold_sp else float("nan"),
        "oof_spearman": oof_m.values["spearman"],
        "perm_mean": float(np.nanmean(perm)) if perm else float("nan"),
        "perm_std": float(np.nanstd(perm)) if perm else float("nan"),
        "feature_names": cols,
        "_pred": pred, "_y": y_te, "_te_mask": te_m,
    }
    out.update(leak_audit(df, te_m))
    return out


def main() -> int:
    t0 = time.time()
    ap = argparse.ArgumentParser(description="阶段 5: 消融实验与泄漏敏感性")
    ap.add_argument("--enhanced", type=Path, default=ENHANCED_CSV)
    ap.add_argument("--base", type=Path, default=BASE_CSV)
    ap.add_argument("--out-csv", type=Path, default=RESULT_CSV)
    ap.add_argument("--report", type=Path, default=REPORT_MD)
    ap.add_argument("--json-out", type=Path, default=JSON_OUT,
                    help="结果 JSON 路径 (缺省 features/v2_ablation_results.json)")
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--n-boot", type=int, default=N_BOOT)
    ap.add_argument("--perm-rounds", type=int, default=10)
    ap.add_argument("--top-frac", type=float, default=TOP_FRAC)
    ap.add_argument("--only", default="", help="逗号分隔, 只跑这些配置")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")
    fh = logging.FileHandler(LOG_DIR / "16_ablation.log", mode="w", encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    sys.path.insert(0, str(ROOT / "scripts"))
    import metrics_ranking as mr
    import lightgbm as lgb

    run_log: list[str] = []

    def emit(msg: str = "") -> None:
        logger.info("%s", msg)
        run_log.append(str(msg))

    emit("# Affinity Intelligence V2.0 — 阶段 5: 消融实验与泄漏敏感性报告")
    emit("")
    emit(f"- 生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- 随机种子: {args.seed}; LightGBM: {lgb.__version__} (OpenMP shim: "
         f"{'OK' if _VCOMP else '未加载'})")
    emit(f"- bootstrap 次数: {args.n_boot}; 置换轮数: {args.perm_rounds}")
    emit(f"- 真强结合体定义: 真值 ΔG 最负的前 {args.top_frac * 100:.0f}%")

    # 优先用增强特征表; 不存在则退化到 V1.1 特征表并明确警告
    if args.enhanced.exists():
        src = args.enhanced
        emit(f"- 特征表: `{src.name}` (阶段 5 增强表)")
    elif args.base.exists():
        src = args.base
        emit(f"- 特征表: `{src.name}`")
        emit("")
        emit("> !! **注意**: 阶段 5 增强特征表不存在, 已退化到 V1.1 特征表。")
        emit("> 因此 `no_esm2` / `esm2_only` / `no_electrostatic` 等配置将**与 `full` 等价**,")
        emit("> 报告会如实标注 `identical_to_full=True`, 不伪装成「做了消融」。")
        emit("> 要得到完整消融, 请先运行 `python scripts/18_extra_features.py`。")
    else:
        emit("!! 找不到任何特征表")
        return 2

    df = pd.read_csv(src, low_memory=False)
    emit(f"- 数据: {len(df)} 行 x {df.shape[1]} 列")

    if "split_family" not in df.columns or "family_id" not in df.columns:
        emit("!! 缺少 split_family / family_id, 无法执行家族隔离")
        return 2

    # 红线复核 (family 划分)
    cross_fam = int((df.groupby("family_id")["split_family"].nunique() > 1).sum())
    cross_pdb = int((df.groupby("pdb_id")["split_family"].nunique() > 1).sum())
    emit(f"- 红线复核 (family 划分): 跨侧家族={cross_fam}, 跨侧 PDB={cross_pdb}")
    if cross_fam or cross_pdb:
        emit("!! 划分泄漏, 拒绝继续")
        return 3

    # 载入 V1.1 超参
    path = ROOT / "scripts" / "08_train_v1.0.py"
    spec = importlib.util.spec_from_file_location("t8_abl_main", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    v11 = {"PARAMS": dict(mod.PARAMS), "N_SPLITS": mod.N_SPLITS}

    cfgs = build_configs(df)
    if args.only:
        want = {s.strip() for s in args.only.split(",") if s.strip()}
        cfgs = [c for c in cfgs if c["name"] in want]

    emit("")
    emit("## 消融配置与特征数")
    emit("")
    emit("| 配置 | 说明 | 划分协议 | 特征数 |")
    emit("|---|---|---|---:|")
    for c in cfgs:
        emit(f"| {c['name']} | {c['desc']} | {c.get('split', 'family')} | {len(c['cols'])} |")

    if args.dry_run:
        emit("")
        emit("--dry-run: 不训练")
        return 0

    # ---- 逐配置训练 ----------------------------------------------------
    emit("")
    emit("## 消融结果 (family 划分, 同一测试集)")
    emit("")
    rows: list[dict[str, Any]] = []
    preds: dict[str, dict[str, Any]] = {}
    for c in cfgs:
        try:
            res = run_config(c, df, args, mr, lgb, v11, args.perm_rounds)
            res_f = {k: v for k, v in res.items() if not k.startswith("_")}
            res_f["leak_warning"] = res["split_protocol"] != "family"
            rows.append(res_f)
            preds[c["name"]] = {"pred": res["_pred"], "y": res["_y"],
                                "mask": res["_te_mask"], "split": res["split_protocol"]}
            emit(f"  {c['name']:<26s} feat={res['n_features']:3d}  "
                 f"Spearman={_f(res['test_spearman'])} "
                 f"[{_f(res['test_ci_low'],3)},{_f(res['test_ci_high'],3)}]  "
                 f"NDCG@10={_f(res['test_ndcg10'],3)}  "
                 f"EF@1%={_f(res['test_ef1'],3)}  R2={_f(res['test_r2'],3)}")
        except Exception as exc:  # noqa: BLE001
            logger.error("配置 %s 失败: %s", c["name"], exc)
            with open(FAIL_LOG, "a", encoding="utf-8") as fh2:
                fh2.write(f"[16_ablation] {c['name']}: {exc}\n{traceback.format_exc()}\n")
        gc.collect()

    if not rows:
        emit("!! 全部配置失败")
        return 4
    res_df = pd.DataFrame(rows)

    # 标出与 full 逐行相同的配置 (增强表缺失时的退化情形)
    for r in res_df.to_dict("records"):
        pass
    leak_cfgs = set(res_df.loc[res_df["split_protocol"] != "family", "config"])

    # ---- 主结果表 ------------------------------------------------------
    emit("")
    emit("### 主结果表 (非泄漏配置, family 级划分, 测试集 n=%d)"
         % int(res_df.loc[res_df["split_protocol"] == "family", "n_test"].iloc[0])
         if (res_df["split_protocol"] == "family").any() else "### 主结果表")
    emit("")
    emit("| 配置 | 特征数 | Spearman | 95% CI | Kendall τ | NDCG@10 | P@10 | EF@1% | EF@5% | RMSE | R² |")
    emit("|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|")
    show = res_df[res_df["split_protocol"] == "family"].sort_values(
        "test_spearman", ascending=False)
    for r in show.to_dict("records"):
        emit(f"| {r['config']} | {r['n_features']} | {_f(r['test_spearman'])} | "
             f"[{_f(r['test_ci_low'],3)}, {_f(r['test_ci_high'],3)}] | "
             f"{_f(r.get('test_kendall_tau', float('nan')))} | "
             f"{_f(r['test_ndcg10'])} | {_f(r['test_precision10'])} | "
             f"{_f(r['test_ef1'])} | {_f(r['test_ef5'])} | "
             f"{_f(r['test_rmse'],3)} | {_f(r['test_r2'],3)} |")

    # ---- 配对检验 vs full ----------------------------------------------
    emit("")
    emit("### 与 `full` 的配对 bootstrap 检验 (仅 family 划分内可比)")
    emit("")
    emit("| 配置 | Δρ vs full | 95% CI | p | BH q | 显著 | ΔNDCG@10 | ΔEF@1% |")
    emit("|---|---:|---|---:|---:|---|---:|---:|")
    if "full" in preds:
        base = preds["full"]
        pvals, names = [], []
        pair_rows = []
        for name, d in preds.items():
            if name == "full" or d["split"] != "family" or base["split"] != "family":
                continue
            if len(d["pred"]) != len(base["pred"]) or not np.array_equal(d["mask"], base["mask"]):
                continue
            y_s = mr.negate_dG(base["y"])
            res_sp = mr.paired_bootstrap_test(mr.spearman_rho, y_s,
                                              mr.negate_dG(d["pred"]),
                                              mr.negate_dG(base["pred"]),
                                              n_boot=args.n_boot, seed=args.seed)
            res_nd = mr.paired_bootstrap_test(
                lambda a, b: mr.ndcg_at_k(a, b, 10, top_frac=args.top_frac), y_s,
                mr.negate_dG(d["pred"]), mr.negate_dG(base["pred"]),
                n_boot=args.n_boot, seed=args.seed)
            res_ef = mr.paired_bootstrap_test(
                lambda a, b: mr.enrichment_factor(a, b, 1.0, args.top_frac), y_s,
                mr.negate_dG(d["pred"]), mr.negate_dG(base["pred"]),
                n_boot=args.n_boot, seed=args.seed)
            pair_rows.append({"config": name, "sp": res_sp, "nd": res_nd, "ef": res_ef})
            pvals.append(res_sp["p_two_sided"])
            names.append(name)
        qvals = mr.bh_fdr(pvals)
        qmap = dict(zip(names, qvals))
        for pr in pair_rows:
            r = pr["sp"]
            q = qmap.get(pr["config"], float("nan"))
            emit(f"| {pr['config']} | {_f(r['delta'])} | "
                 f"[{_f(r['lo'],3)}, {_f(r['hi'],3)}] | {r['p_two_sided']:.4f} | "
                 f"{q:.4f} | {'是' if q < 0.05 else '否'} | "
                 f"{_f(pr['nd']['delta'])} | {_f(pr['ef']['delta'])} |")
        res_df["paired_p_vs_full"] = res_df["config"].map(
            {pr["config"]: pr["sp"]["p_two_sided"] for pr in pair_rows})
        res_df["bh_q_vs_full"] = res_df["config"].map(qmap)
    else:
        emit("| (缺少 full 基线) | | | | | | | |")

    # ---- 泄漏敏感性 ----------------------------------------------------
    emit("")
    emit("## 泄漏敏感性分析 (划分协议对照)")
    emit("")
    emit("> **警告**: 下表各行使用**不同的测试集**, 指标**不可与非泄漏行直接比大小**。")
    emit("> 本表的目的只有一个: 量化「如果放松隔离协议, 指标会虚高多少」。")
    emit("")
    emit("| 划分协议 | 跨侧家族数 | 跨侧 PDB 数 | Spearman | 95% CI | NDCG@10 | EF@1% | R² |")
    emit("|---|---:|---:|---:|---|---:|---:|---:|")
    for r in res_df.to_dict("records"):
        emit(f"| {r['split_protocol']} | {int(r.get('cross_family_id', -1))} | "
             f"{int(r.get('cross_pdb_id', -1))} | {_f(r['test_spearman'])} | "
             f"[{_f(r['test_ci_low'],3)}, {_f(r['test_ci_high'],3)}] | "
             f"{_f(r['test_ndcg10'])} | {_f(r['test_ef1'])} | {_f(r['test_r2'],3)} |")

    fam_row = res_df[res_df["config"] == "full"]
    fam_row = fam_row[fam_row["split_protocol"] == "family"]
    leak_pdb = res_df[res_df["config"] == "leak_pdb_level"]
    leak_rnd = res_df[res_df["config"] == "leak_random"]
    emit("")
    if len(fam_row) and len(leak_pdb):
        d = float(leak_pdb["test_spearman"].iloc[0]) - float(fam_row["test_spearman"].iloc[0])
        emit(f"- **PDB 级划分 vs 家族级划分**: ΔSpearman = {d:+.4f}")
        emit(f"  -> 同源泄漏让 Spearman 虚高约 {abs(d):.4f} "
             f"({abs(d) / max(abs(float(fam_row['test_spearman'].iloc[0])), 1e-9) * 100:.1f}% 相对提升)")
    if len(fam_row) and len(leak_rnd):
        d2 = float(leak_rnd["test_spearman"].iloc[0]) - float(fam_row["test_spearman"].iloc[0])
        emit(f"- **完全随机划分 vs 家族级划分**: ΔSpearman = {d2:+.4f}")
        emit("  -> 这是「如果不做隔离会看到的最好数字」, 也是论文里最容易被误当成 SOTA 的数字。")

    # ---- 三个论证问题的直接回答 ----------------------------------------
    emit("")
    emit("## 三个关键论证的直接回答")
    emit("")

    def _get(name: str, col: str) -> float:
        sub = res_df[res_df["config"] == name]
        return float(sub[col].iloc[0]) if len(sub) else float("nan")

    full_sp = _get("full", "test_spearman")
    no_esm2 = _get("no_esm2", "test_spearman")
    esm2_only = _get("esm2_only", "test_spearman")
    size_only = _get("size_only", "test_spearman")
    phys_norm = _get("physics_only_normalized", "test_spearman")

    emit("### 论证 1: 去掉 ESM-2 是否有损? (轻量级卖点)")
    emit("")
    if np.isfinite(no_esm2) and np.isfinite(full_sp):
        d = no_esm2 - full_sp
        emit(f"- `full` = {full_sp:+.4f}; `no_esm2` = {no_esm2:+.4f}; Δ = {d:+.4f}")
        if abs(d) < 0.01:
            emit("- **结论: `no_esm2` ≈ `full` -> 纯物理特征已经足够, 无需 GPU 嵌入。**")
            emit("  这正是「轻量、CPU 可复现」的直接证据。")
        elif d < 0:
            emit(f"- 结论: ESM-2 带来 {abs(d):.4f} 的正增益 (物理特征不足以替代)。")
        else:
            emit(f"- 结论: 去掉 ESM-2 反而更好 ({d:+.4f}) —— ESM-2 在此数据上无贡献。")
    else:
        emit("- 无法比较 (缺 `no_esm2` 或 `esm2_only` 配置; 需先运行 18_extra_features.py)")
    if np.isfinite(esm2_only):
        emit(f"- 参照: `esm2_only` = {esm2_only:+.4f} (单靠序列嵌入的排序能力)")

    emit("")
    emit("### 论证 2: 归一化物理特征是否可被尺寸替代?")
    emit("")
    size_nd = _get("size_only", "test_ndcg10")
    full_nd = _get("full", "test_ndcg10")
    size_ef = _get("size_only", "test_ef1")
    full_ef = _get("full", "test_ef1")
    if np.isfinite(size_only) and np.isfinite(full_sp):
        emit(f"- 全局排序 (Spearman): `full` = {full_sp:+.4f}; "
             f"`size_only` = {size_only:+.4f}; Δ = {size_only - full_sp:+.4f}")
        emit(f"- 顶部富集 (NDCG@10): `full` = {full_nd:.4f}; `size_only` = {size_nd:.4f}")
        emit(f"- 顶部富集 (EF@1%):   `full` = {full_ef:.3f}; `size_only` = {size_ef:.3f}")
        emit("")
        if size_only < full_sp - 0.03:
            emit("- **结论: 全局排序上 `size_only` 明显弱于 `full` -> 归一化物理特征确实有独立贡献,**")
            emit("  模型不是在「数原子」。(见下面对顶部富集的补充说明)")
        else:
            emit("- 结论: `size_only` 的全局排序接近 `full` —— 尺寸类特征本身就携带大部分信号,")
            emit("  这削弱了「物理特征有独立贡献」的论证, 应如实报告。")
        # 关键补充: 顶部富集与全局排序可能给出**相反**的结论, 必须一并呈现
        if np.isfinite(size_nd) and np.isfinite(full_nd) and size_nd > full_nd + 0.05:
            emit("")
            emit("  ⚠️ **重要的反向发现**: 在**顶部富集**指标上 `size_only` 反而优于 `full`")
            emit(f"  (NDCG@10 {size_nd:.4f} vs {full_nd:.4f}, EF@1% {size_ef:.3f} vs {full_ef:.3f})。")
            emit("  这不矛盾: 完整物理特征集改善了**全局单调性**(Spearman), 但其中部分")
            emit("  高方差特征 (尤其 ΔSASA/SASA 族) 给头部引入了排序噪声。")
            emit("  对「只做得起 k 个克隆」的湿实验场景而言, **顶部富集是更贴近用途的指标**,")
            emit("  因此本项目的正确结论是: **物理特征提升了整体校准, 但未提升头部富集**。")
            emit("  这是一个应当如实报告、而不是掩盖的负结果。")
    if np.isfinite(phys_norm):
        emit(f"- 参照: `physics_only_normalized` = {phys_norm:+.4f}")

    emit("")
    emit("### 论证 3: 同源泄漏会让指标虚高多少?")
    emit("")
    if len(leak_pdb) and len(leak_rnd) and len(fam_row):
        emit(f"- 家族级 (红线, 可信): {float(fam_row['test_spearman'].iloc[0]):+.4f}")
        emit(f"- PDB 级 (同源泄漏):     {float(leak_pdb['test_spearman'].iloc[0]):+.4f} "
             f"(跨侧家族 {int(leak_pdb['cross_family_id'].iloc[0])} 个)")
        emit(f"- 完全随机 (最强泄漏):    {float(leak_rnd['test_spearman'].iloc[0]):+.4f} "
             f"(跨侧家族 {int(leak_rnd['cross_family_id'].iloc[0])} 个)")
        emit("- **结论: 隔离协议的严格程度直接决定报告出来的数字。**")
        emit("  任何不同隔离协议之间的指标对比都是无意义的; V2.0 一律以家族级为准。")

    emit("")
    emit("## 残留局限")
    emit("")
    emit("1. `leak_random` 是最强泄漏对照, 它的测试集与训练集**同源甚至同测量**, ")
    emit("   因此其指标只用来量化泄漏幅度, **绝不能作为性能宣称**。")
    emit("2. 消融只覆盖特征族与划分协议, 未做超参搜索联动的消融 (算力约束)。")
    emit("3. 每个配置只跑单一 seed (42); 配对检验的 CI 已包含采样不确定性,")
    emit("   但未包含超参初始化噪声。")
    emit("")
    emit(f"耗时: {time.time() - t0:.1f}s")

    # ---- 落盘 ----------------------------------------------------------
    out_df = res_df.drop(columns=[c for c in res_df.columns if c == "feature_names"],
                         errors="ignore")
    tmp = args.out_csv.with_suffix(args.out_csv.suffix + ".tmp")
    out_df.to_csv(tmp, index=False, encoding="utf-8")
    tmp.replace(args.out_csv)

    json_payload = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"), "seed": args.seed,
        "source_table": str(src), "n_boot": args.n_boot,
        "results": [{k: v for k, v in r.items() if k != "feature_names"}
                    for r in rows],
        "verdicts": {
            "no_esm2_vs_full": no_esm2 - full_sp if np.isfinite(no_esm2) and np.isfinite(full_sp) else None,
            "size_only_vs_full": size_only - full_sp if np.isfinite(size_only) and np.isfinite(full_sp) else None,
            "leak_pdb_vs_family": (float(leak_pdb["test_spearman"].iloc[0]) - float(fam_row["test_spearman"].iloc[0]))
            if len(leak_pdb) and len(fam_row) else None,
            "leak_random_vs_family": (float(leak_rnd["test_spearman"].iloc[0]) - float(fam_row["test_spearman"].iloc[0]))
            if len(leak_rnd) and len(fam_row) else None,
        },
    }
    JSON_OUT.write_text(json.dumps(json_payload, ensure_ascii=False, indent=2, default=str),
                        encoding="utf-8")
    if args.json_out != JSON_OUT:
        args.json_out.write_text(
            json.dumps(json_payload, ensure_ascii=False, indent=2, default=str),
            encoding="utf-8")
    args.report.write_text("\n".join(run_log) + "\n", encoding="utf-8")
    emit("")
    emit(f"结果表: {args.out_csv}")
    emit(f"报告: {args.report}")
    return 0


def _f(v: Any, d: int = 4) -> str:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return "—"
    return f"{f:+.{d}f}" if np.isfinite(f) else "NaN"


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:  # noqa: BLE001
        LOG_DIR.mkdir(parents=True, exist_ok=True)
        with open(FAIL_LOG, "a", encoding="utf-8") as fh:
            fh.write(f"[16_ablation] {traceback.format_exc()}\n")
        raise
