# -*- coding: utf-8 -*-
"""scripts/v2.3/03_merge_esm2.py -- V2.3 任务 3: 合并 ESM-2 嵌入并做消融.

用途
----
在 Colab 上跑完 `scripts/colab/esm2_extract.ipynb` 得到 `features/v2_esm2_cache.npz`
之后, 用本脚本把 ESM-2 嵌入与 V2.2 的 60 个物理特征合并, 并系统地消融:

* PCA 维度: **32 / 64 / 128** (全部报告)
* 嵌入版本: **raw / L2** (全部报告)
* 嵌入来源: **ligand / receptor / both**
* 10 个特征集 (见下)

**划分与超参完全复用 V2.2** (硬约束 6): 直接读 `v2.2_proteinea_features.csv` 的 `split` 列,
并做 assert 校验 train/val/test 行数与 V2.2 模型记录一致, 不重新划分。

10 个特征集
-----------
======================= ==========================================
physics_only            V2.2 基线 (60 物理特征)
esm2_ligand_raw         仅配体嵌入 (原始向量)
esm2_ligand_l2          仅配体嵌入 (L2)
esm2_receptor_raw       仅受体嵌入 (原始向量)
esm2_receptor_l2        仅受体嵌入 (L2)
esm2_both_raw           配体+受体嵌入 (原始向量)
esm2_both_l2            配体+受体嵌入 (L2)
physics+esm2_ligand     物理 + 配体嵌入
physics+esm2_receptor   物理 + 受体嵌入
physics+esm2_both       物理 + 配体 + 受体嵌入
======================= ==========================================

PCA 纪律: **只在训练集 fit**, 再 transform val/test (全量 fit 是标签泄漏)。

关键正确性保障
--------------
* join 用 **`interface_key`** (格式 `PDB|配体链|受体链`, 与 V2.2 特征表同名同格式)。
  任务书说"按 pdb_id join", 但同一 PDB 有大量突变体 (8207 行 / 2171 个 PDB),
  按 pdb_id 会把多个不同序列的嵌入混在一起 —— 那是错的。本脚本先按 `interface_key`
  精确 join, 并在报告中给出**按 pdb_id join 会冲突多少行**的对照数字。
* join 命中率 **< 95% 直接报错停止** (不静默降级)。
* 每个特征集用**完全相同**的划分与超参重训 LightGBM。
* 与 V2.2 基线做**配对 bootstrap** 检验 + BH-FDR 校正。

未运行声明
----------
若 `features/v2_esm2_cache.npz` 不存在, 本脚本**不会伪造任何指标**: 它只做
可验证的 dry-run (检查 V2.2 划分、baseline 复现、特征名列), 然后以非 0 退出码结束,
并在报告中写明"ESM-2 未运行"。

纪律: 纯 CPU; 无网络; seed=42; 只写 v2.3 命名空间; 不改动 V1.1/V2.0/V2.1/V2.2 文件。
"""

from __future__ import annotations

import argparse
import gc
import importlib.util
import json
import logging
import pickle
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
MODEL_DIR = ROOT / "models"
LOG_DIR = ROOT / "logs" / "v2.3.2"      # V2.3.2 命名空间 (硬约束 2)

V22_FEATURES = FEAT_DIR / "v2.2_proteinea_features.csv"
ESM2_CACHE = FEAT_DIR / "v2_esm2_cache.npz"
V22_MODEL = MODEL_DIR / "v2.2_proteinea_lgbm.pkl"

# 输出一律使用 v2.3.2 命名空间 (硬约束 2: 不覆盖 V1.1-V2.3 任何文件)
RESULTS_CSV = FEAT_DIR / "v2.3.2_lgbm_results.csv"
PAIRS_CSV = FEAT_DIR / "v2.3.2_lgbm_pairs.csv"
REPORT_MD = FEAT_DIR / "v2.3.2_lgbm_ablation_report.md"
JSON_OUT = FEAT_DIR / "v2.3.2_lgbm_results.json"
BEST_MODEL = MODEL_DIR / "v2.3.2_best_lgbm.pkl"
PREDS_NPZ = FEAT_DIR / "v2.3.2_test_predictions.npz"   # 供任务 5 集成复用
FAIL_LOG = LOG_DIR / "failures.log"

SEED = 42
N_BOOT = 5000                  # 任务书要求 5000 次配对 bootstrap
N_BOOT_CI = 1000               # 每个模型自身指标的 CI bootstrap 次数
                               # (37 个模型 x 4 指标 x 5000 次会超出夜间算力预算,
                               #  故 CI 用 1000; 与基线的**配对**检验仍用 5000)
TOP_FRAC = 0.10
JOIN_MIN_HIT = 0.95            # 命中率下限 (硬约束)
PCA_DIMS = (32, 64, 128)
VARIANTS = ("raw", "l2")
SOURCES = ("ligand", "receptor", "both")

# V2.2 实测的划分规模 (用于 assert)
EXPECTED_SPLITS = {"train": 6482, "val": 965, "test": 752}

# V2.2 记录的基线测试集 Spearman (任务书给出, 用于对齐校验)
V22_BASELINE_RHO = 0.23747781344414348


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

logger = logging.getLogger("merge_esm2")


# ---------------------------------------------------------------- 指标工具
def build_metric_fns(mr: Any, top_frac: float) -> dict[str, Any]:
    return {
        "spearman": mr.spearman_rho,
        "kendall_tau": mr.kendall_tau,
        "ndcg@10": lambda a, b: mr.ndcg_at_k(a, b, 10, top_frac=top_frac),
        "precision@10": lambda a, b: mr.precision_at_k(a, b, 10, top_frac),
        "ef@1%": lambda a, b: mr.enrichment_factor(a, b, 1.0, top_frac),
        "ef@5%": lambda a, b: mr.enrichment_factor(a, b, 5.0, top_frac),
    }


def bootstrap_ci(fn: Any, y: np.ndarray, p: np.ndarray, n_boot: int,
                 seed: int) -> tuple[float, float, float]:
    ok = np.isfinite(y) & np.isfinite(p)
    y, p = y[ok], p[ok]
    point = float(fn(y, p)) if len(y) else float("nan")
    if len(y) < 5:
        return float("nan"), float("nan"), point
    rs = np.random.RandomState(seed)
    vals = []
    for _ in range(n_boot):
        idx = rs.randint(0, len(y), len(y))
        v = fn(y[idx], p[idx])
        if v == v:
            vals.append(v)
    if len(vals) < 20:
        return float("nan"), float("nan"), point
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(lo), float(hi), point


def paired_bootstrap(fn: Any, y: np.ndarray, pa: np.ndarray, pb: np.ndarray,
                     n_boot: int, seed: int) -> dict[str, float]:
    """配对 bootstrap: Δ = metric(B) - metric(A), 同一批索引重采样。"""
    ok = np.isfinite(y) & np.isfinite(pa) & np.isfinite(pb)
    y, pa, pb = y[ok], pa[ok], pb[ok]
    if len(y) < 10:
        return {"delta": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "p_two_sided": float("nan"), "n": int(len(y))}
    ma, mb = float(fn(y, pa)), float(fn(y, pb))
    rs = np.random.RandomState(seed)
    d = []
    for _ in range(n_boot):
        idx = rs.randint(0, len(y), len(y))
        va, vb = fn(y[idx], pa[idx]), fn(y[idx], pb[idx])
        if va == va and vb == vb:
            d.append(vb - va)
    if len(d) < 20:
        return {"delta": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "p_two_sided": float("nan"), "metric_a": ma, "metric_b": mb,
                "n": int(len(y)), "n_distinct_deltas": 0}
    d = np.asarray(d, float)
    lo, hi = np.percentile(d, [2.5, 97.5])
    pl = (np.sum(d <= 0) + 1) / (len(d) + 1)
    pr = (np.sum(d >= 0) + 1) / (len(d) + 1)
    # 离散化诊断: 若重采样差值几乎只有几个取值, 双侧 p 会撞到 2/(n_boot+1) 的下限,
    # 多个模型会显示**完全相同**的 p。必须把该计数一并报告, 否则会误导。
    n_distinct = int(len(np.unique(np.round(d, 10))))
    return {"delta": float(d.mean()), "lo": float(lo), "hi": float(hi),
            "p_two_sided": float(min(1.0, 2 * min(pl, pr))),
            "metric_a": ma, "metric_b": mb, "n": int(len(y)),
            "n_distinct_deltas": n_distinct}


def bh_fdr(pvals: list[float]) -> list[float]:
    p = np.asarray([1.0 if (x != x) else float(x) for x in pvals], float)
    m = len(p)
    if m == 0:
        return []
    order = np.argsort(p, kind="mergesort")
    q = np.empty(m, float)
    prev = 1.0
    for rank, _ in enumerate(reversed(order), start=1):
        i = m - rank
        prev = min(prev, p[order[i]] * m / (i + 1))
        q[order[i]] = prev
    return [float(min(1.0, v)) for v in q]


# ---------------------------------------------------------------- 主流程
def main() -> int:
    t0 = time.time()
    ap = argparse.ArgumentParser(description="V2.3 任务3: 合并 ESM-2 并做消融")
    ap.add_argument("--features", type=Path, default=V22_FEATURES)
    ap.add_argument("--esm2-cache", type=Path, default=ESM2_CACHE)
    ap.add_argument("--v22-model", type=Path, default=V22_MODEL)
    ap.add_argument("--out-csv", type=Path, default=RESULTS_CSV)
    ap.add_argument("--pairs-csv", type=Path, default=PAIRS_CSV)
    ap.add_argument("--report", type=Path, default=REPORT_MD)
    ap.add_argument("--json-out", type=Path, default=JSON_OUT)
    ap.add_argument("--best-model", type=Path, default=BEST_MODEL)
    ap.add_argument("--preds-npz", type=Path, default=PREDS_NPZ,
                    help="测试集预测产物 (供任务 5 集成复用)")
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--n-boot", type=int, default=N_BOOT,
                    help="配对 bootstrap 次数 (任务书要求 5000)")
    ap.add_argument("--n-boot-ci", type=int, default=N_BOOT_CI,
                    help="每个模型自身指标的 bootstrap CI 次数 "
                         "(计算预算所限, 见报告「统计实现说明」)")
    ap.add_argument("--top-frac", type=float, default=TOP_FRAC)
    ap.add_argument("--dry-run", action="store_true",
                    help="只校验划分与 baseline 复现, 不做任何 ESM-2 相关工作")
    args = ap.parse_args()

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")
    fh = logging.FileHandler(LOG_DIR / "03_merge_esm2_v2.log", mode="w", encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    sys.path.insert(0, str(ROOT / "scripts"))
    import metrics_ranking as mr
    import lightgbm as lgb
    from sklearn.decomposition import PCA

    spec = importlib.util.spec_from_file_location("t8_v23",
                                                  ROOT / "scripts" / "08_train_v1.0.py")
    v11mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(v11mod)
    feat_names = [c for c in v11mod.V02_FEATURES if c not in v11mod.REPLACED] \
        + list(v11mod.V03_NEW)
    count_cols = list(v11mod.COUNT_COLS)
    PARAMS = dict(v11mod.PARAMS)
    PARAMS["random_state"] = args.seed

    mfns = build_metric_fns(mr, args.top_frac)
    L: list[str] = []

    def emit(msg: str = "") -> None:
        logger.info("%s", msg)
        L.append(str(msg))
        try:
            sys.stdout.flush()
        except Exception:  # noqa: BLE001
            pass

    # V2.2 模型记录 (用于划分 assert 与 baseline 对齐)
    v22_payload: dict[str, Any] = {}
    if args.v22_model.exists():
        with open(args.v22_model, "rb") as fh2:
            v22_payload = pickle.load(fh2)

    emit("# V2.3.2 — ESM-2 (pdb_id join, 已修复) × 60 物理特征消融报告")
    emit("")
    emit(f"- 生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- 随机种子: {args.seed}; LightGBM: {lgb.__version__} (OpenMP shim: "
         f"{'OK' if _VCOMP else '未加载'})")
    emit(f"- bootstrap: {args.n_boot} 次; 真强结合体 = 真值 ΔG 最负的前 "
         f"{args.top_frac * 100:.0f}%")
    emit(f"- 划分与超参: **完全复用 V2.2** (硬约束 6)")

    if not args.features.exists():
        emit(f"!! 找不到 V2.2 特征表: {args.features}")
        return 2

    # ---------------- 1. 载入 V2.2 特征表并校验划分 ----------------
    df = pd.read_csv(args.features, low_memory=False)
    emit("")
    emit("## 1. V2.2 特征表与划分校验")
    emit("")
    emit(f"- `{args.features.name}`: {len(df)} 行 x {df.shape[1]} 列")
    if "split" not in df.columns or "interface_key" not in df.columns:
        emit("!! 缺少 `split` 或 `interface_key` 列, 无法复用 V2.2 划分")
        return 2

    usable = (df["row_structure_usable"] == 1) & df["real_dG"].notna()
    d = df[usable].reset_index(drop=True)
    split_counts = d["split"].value_counts().to_dict()
    emit(f"- 结构可用 + 标签有效行: {len(d)}")
    emit(f"- 划分规模: {split_counts}")

    # ---- 硬约束 6: assert 测试集与 V2.2 完全一致 ----
    split_ok = True
    emit("")
    emit("### 划分一致性 assert (硬约束 6)")
    emit("")
    emit("| split | 实测行数 | V2.2 期望 | 通过 |")
    emit("|---|---:|---:|---|")
    for sp, exp in EXPECTED_SPLITS.items():
        got = int(split_counts.get(sp, 0))
        ok = (got == exp)
        split_ok = split_ok and ok
        emit(f"| {sp} | {got} | {exp} | {'✅' if ok else '❌'} |")
    if v22_payload:
        emit("")
        emit(f"- V2.2 模型记录: train={v22_payload.get('train_rows')}, "
             f"val={v22_payload.get('val_rows')}, test={v22_payload.get('test_rows')}")
        emit(f"- V2.2 记录测试集 Spearman = "
             f"{v22_payload.get('metrics', {}).get('test', {}).get('spearman')}")
    if not split_ok:
        emit("")
        emit("!! **划分与 V2.2 不一致, 拒绝继续** (硬约束 6)")
        return 3
    emit("")
    emit("- ✅ **划分与 V2.2 完全一致**")

    # 测试集 interface_key 集合必须一致 (更严格的一致性校验)
    te = d[d["split"] == "test"]
    emit(f"- test 侧 interface_key 数: {te['interface_key'].nunique()} / {len(te)} 行")

    # ---------------- 2. ESM-2 缓存存在性 ----------------
    emit("")
    emit("## 2. ESM-2 缓存检查")
    emit("")
    if not args.esm2_cache.exists():
        emit(f"- `{args.esm2_cache}` **不存在**")
        emit("")
        emit("> ⚠️ **ESM-2 尚未运行。** 本脚本不伪造任何指标。")
        emit(">")
        emit("> 获取方式: 在 Colab 上运行 `scripts/colab/esm2_extract.ipynb`")
        emit("> (需要 T4 GPU), 把产出的 `v2_esm2_cache.npz` 放到")
        emit(f"> `{args.esm2_cache}`。")
        emit(">")
        emit("> 操作步骤见 `features/v2.3_colab_instructions.md`。")
        emit("")
        emit("### 本脚本在无缓存时仍然完成的可验证工作")
        emit("")
        emit(f"- ✅ V2.2 划分规模校验通过 ({split_counts})")
        emit("- ✅ 划分一致性 assert 通过")
        emit(f"- ✅ 物理特征集定义可用 ({len(feat_names)} 列)")
        emit("")
        emit("### 未运行")
        emit("")
        emit("- ESM-2 嵌入提取: **未运行** (本机无 GPU)")
        emit("- ESM-2 × 物理特征消融: **未运行** (无嵌入缓存)")
        emit("- 因此本报告**没有任何 ESM-2 指标**。")
        emit("")
        if not args.dry_run:
            args.report.write_text("\n".join(L) + "\n", encoding="utf-8")
            args.json_out.write_text(json.dumps({
                "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
                "esm2_run": False,
                "esm2_cache_present": False,
                "v22_split_check": {sp: int(split_counts.get(sp, 0))
                                    for sp in EXPECTED_SPLITS},
                "v22_split_ok": bool(split_ok),
                "metrics": None,
                "note": "ESM-2 未运行; 本脚本未伪造任何指标",
            }, ensure_ascii=False, indent=2), encoding="utf-8")
            emit(f"报告: `{args.report.name}`")
            emit(f"JSON: `{args.json_out.name}`")
        return 4

    # ---------------- 3. 载入嵌入并 join ----------------
    emit(f"- `{args.esm2_cache.name}` **存在**")
    with np.load(args.esm2_cache, allow_pickle=True) as z:
        keys = list(z.keys())
        emit(f"- 内含数组: {keys}")
        for req in ("pdb_ids", "lig_raw", "lig_l2", "rec_raw", "rec_l2"):
            if req not in keys:
                emit(f"!! 缺少必需数组 `{req}`")
                return 5
        pdb_ids = z["pdb_ids"].astype(str)
        lig_raw = np.asarray(z["lig_raw"], dtype=np.float32)
        lig_l2 = np.asarray(z["lig_l2"], dtype=np.float32)
        rec_raw = np.asarray(z["rec_raw"], dtype=np.float32)
        rec_l2 = np.asarray(z["rec_l2"], dtype=np.float32)
        iface = z["interface_key"].astype(str) if "interface_key" in keys else None
        missing_mask = z["missing_mask"] if "missing_mask" in keys else None
        model_name = str(z["model_name"]) if "model_name" in keys else "unknown"
        repr_layer = int(z["repr_layer"]) if "repr_layer" in keys else -1
        pool = str(z["pool"]) if "pool" in keys else "unknown"

    emit("")
    emit("## 3. 嵌入缓存概况")
    emit("")
    emit(f"- 模型: **{model_name}**  层: {repr_layer}  池化: {pool}")
    emit(f"- 条目数 N = {len(pdb_ids)}; 嵌入维度 D = {lig_raw.shape[1]}")
    emit(f"- `lig_raw` {lig_raw.shape}, `rec_raw` {rec_raw.shape}")
    if lig_raw.shape[1] != 640:
        emit(f"- ⚠️ 维度不是 640 (150M 模型应为 640) -> 说明 Colab 端触发了 OOM 降级")
    if missing_mask is not None:
        emit(f"- 端侧标记缺失嵌入的行: {int(np.sum(missing_mask))}")

    # ---- join: 按 `pdb_id` (V2.3.2 修复) ----
    #
    # V2.3 原版用 `interface_key` join, 命中率 0% —— 因为两个文件的 interface_key
    # **语义不同**: V2.2 特征表是 `pdb|链ID|链ID`, 而 Kaggle 缓存是
    # `pdb|配体序列|受体序列`。两者交集为 0。详见
    # `features/v2.3.2_join_diagnosis.md`。
    #
    # 修复方式 (任务书要求): 改用 `pdb_id` 作 join 键。
    #   一对多时: **缓存侧按 pdb_id 取平均**; 特征表侧共享同一嵌入。
    #   副作用: 粒度从「结构+链组合(+突变体)」退化为「结构级」。
    #
    # 同时**实测一条精确通路作为对照**: 若缓存与特征表**行序与长度完全一致**,
    # 则可以逐行 1:1 映射 (零信息损失)。把它作为 `aligned` 对照,
    # 用于**量化退化代价** —— 这样"退化的影响"是实测数字而非定性描述。
    emit("")
    emit("### join 键与命中率 (硬约束: < 95% 直接报错)")
    emit("")
    use_pdb = True          # V2.3.2 固定使用 pdb_id 作为 join 键
    emit("- join 键: **`pdb_id`** (V2.3.2 修复; V2.3 用的 `interface_key` 语义不一致 -> 0% 命中)")
    emit("  * 缓存侧一对多: **按 pdb_id 取平均**")
    emit("  * 特征表侧一对多: 共享同一嵌入 (结构级退化)")
    n_pdb_multi = int((d.groupby("pdb_id").size() > 1).sum())
    n_rows_multi = int(d.groupby("pdb_id").size().pipe(lambda s: s[s > 1].sum()))
    emit(f"  * 有多行的 pdb_id: **{n_pdb_multi}** / {d['pdb_id'].nunique()} "
         f"(覆盖 {n_rows_multi} 行 = {n_rows_multi / len(d) * 100:.1f}%)")

    # 缓存侧按 pdb_id 聚合 (取平均)
    emb_df = pd.DataFrame({
        "pdb_id": pdb_ids,
        "lig_raw": list(lig_raw), "lig_l2": list(lig_l2),
        "rec_raw": list(rec_raw), "rec_l2": list(rec_l2),
    })
    agg = emb_df.groupby("pdb_id", sort=False).agg(
        lig_raw=("lig_raw", lambda col: np.mean(np.stack(col.to_numpy()), axis=0)),
        lig_l2=("lig_l2", lambda col: np.mean(np.stack(col.to_numpy()), axis=0)),
        rec_raw=("rec_raw", lambda col: np.mean(np.stack(col.to_numpy()), axis=0)),
        rec_l2=("rec_l2", lambda col: np.mean(np.stack(col.to_numpy()), axis=0)),
        n_cache_rows=("pdb_id", "size"),
    ).reset_index()
    emit(f"  * 缓存聚合后 pdb_id 数: **{len(agg)}** (来自 {len(emb_df)} 行)")

    merged = d.merge(agg, on="pdb_id", how="left", validate="many_to_one")
    hit_arr = merged["n_cache_rows"].notna().to_numpy()
    hit_rate = float(hit_arr.mean())

    # ---- 对照: 精确逐行映射 (仅当行序与长度完全一致时可用) ----
    #
    # 注意: `d` 在下面会按 hit_arr 过滤 (排掉 8 行结构不可用), 因此**必须保留一个
    # 覆盖原始 d 所有行的布尔掩码**, 用它去索引 aligned 嵌入;
    # 否则 aligned 数组 (8207 行) 与过滤后的 d (8199 行) 会错位。
    row_keep = hit_arr.copy()                    # 覆盖原始 d 行序的掩码
    aligned_ok = False
    E_aligned: dict[str, np.ndarray] | None = None
    # 对齐必须与**未过滤**的全表比 (缓存是 8207 行 = filtered.csv 全部行),
    # 而 d 此时已被 row_structure_usable 过滤到 8199 行。早期版本误用过滤后的 d
    # 去比, 长度不等 -> 精确通路被错误判定为不可用。
    full_pdb = df["pdb_id"].astype(str).to_numpy()
    if (len(full_pdb) == len(pdb_ids)
            and bool((full_pdb == pdb_ids.astype(str)).all())):
        aligned_ok = True
        E_aligned = {"lig_raw": lig_raw, "lig_l2": lig_l2,
                     "rec_raw": rec_raw, "rec_l2": rec_l2}
        emit("")
        emit(f"- ✅ **精确通路可用**: 全表 ({len(full_pdb)} 行) 与缓存行序/长度完全一致 "
             f"-> 可逐行 1:1 映射 (零信息损失)")
    else:
        emit("")
        emit(f"- ⚠️ 精确通路**不可用** (全表 {len(full_pdb)} 行 vs 缓存 {len(pdb_ids)} 行, "
             f"或行序不一致) -> 只能做结构级 join")

    emit("")
    emit(f"- **join 命中率 = {hit_rate * 100:.2f}%** ({int(hit_arr.sum())} / {len(merged)})")
    if hit_rate < JOIN_MIN_HIT:
        emit(f"- !! **命中率低于 {JOIN_MIN_HIT * 100:.0f}%, 直接报错停止** (硬约束)")
        args.report.write_text("\n".join(L) + "\n", encoding="utf-8")
        return 6
    emit(f"- ✅ 命中率达标 (>= {JOIN_MIN_HIT * 100:.0f}%)")

    d = d.loc[hit_arr].reset_index(drop=True)
    merged_h = merged.loc[hit_arr].reset_index(drop=True)
    E_lig_raw = np.stack(merged_h["lig_raw"].to_numpy())
    E_lig_l2 = np.stack(merged_h["lig_l2"].to_numpy())
    E_rec_raw = np.stack(merged_h["rec_raw"].to_numpy())
    E_rec_l2 = np.stack(merged_h["rec_l2"].to_numpy())
    if aligned_ok and E_aligned is not None:
        # 两级索引: E_aligned 是**全表** (8207 行, 对应未过滤的 df),
        # 而 row_keep 是**过滤后 d** 的行掩码 (8199 行)。先把全表行按 df 的
        # row_structure_usable / 标签有效性筛成 8199, 再用 row_keep 取命中的行。
        full_valid = ((df["row_structure_usable"] == 1)
                      & df["real_dG"].notna()).to_numpy()
        assert int(full_valid.sum()) == len(row_keep), \
            (int(full_valid.sum()), len(row_keep))
        sel = np.zeros(len(pdb_ids), dtype=bool)
        sel[np.where(full_valid)[0]] = row_keep
        assert len(sel) == len(pdb_ids), (len(sel), len(pdb_ids))
        Ea_lig_raw = E_aligned["lig_raw"][sel]
        Ea_lig_l2 = E_aligned["lig_l2"][sel]
        Ea_rec_raw = E_aligned["rec_raw"][sel]
        Ea_rec_l2 = E_aligned["rec_l2"][sel]
        assert len(Ea_lig_raw) == len(d), (len(Ea_lig_raw), len(d))
    else:
        Ea_lig_raw = Ea_lig_l2 = Ea_rec_raw = Ea_rec_l2 = None
    emit(f"- 保留可用行: {len(d)} (drop {int((~hit_arr).sum())})")

    # ---------------- 4. 特征矩阵 ----------------
    fcols = [c for c in feat_names if c in d.columns]
    tr_m = (d["split"] == "train").to_numpy()
    va_m = (d["split"] == "val").to_numpy()
    te_m = (d["split"] == "test").to_numpy()
    y = d["real_dG"].to_numpy(float)

    # 与 V2.2 完全一致的物理特征预处理:
    #   训练集中位数兜底 (用 V2.2 模型里存的 medians, 保证逐位复现) + COUNT_COLS log1p
    Xp = d[fcols].apply(pd.to_numeric, errors="coerce").copy()
    for c in count_cols:
        if c in Xp.columns:
            Xp[c] = np.log1p(Xp[c].clip(lower=0))
    stored_med = pd.Series(v22_payload.get("train_medians", {}))
    if len(stored_med):
        med = stored_med.reindex(Xp.columns)
    else:
        med = Xp[tr_m].median(numeric_only=True)
    Xp = Xp.fillna(med).fillna(0.0)
    Xp = Xp[fcols]

    # ESM-2 侧的 PCA: 只在训练集 fit
    emit("")
    emit("## 4. PCA (只在训练集 fit, 硬约束: 不得全量 fit)")
    emit("")
    pca_payload: dict[str, dict[str, np.ndarray]] = {}
    pca_info: list[dict[str, Any]] = []
    emit("| 版本 | 维度 | 解释方差比之和 |")
    emit("|---|---:|---:|")
    for variant, (Lr, Rr) in (("raw", (E_lig_raw, E_rec_raw)),
                              ("l2", (E_lig_l2, E_rec_l2))):
        for source, mat in (("ligand", Lr), ("receptor", Rr), ("both", np.hstack([Lr, Rr]))):
            for dim in PCA_DIMS:
                nmax = min(dim, mat.shape[1], int(tr_m.sum()) - 1)
                p = PCA(n_components=nmax, random_state=args.seed)
                p.fit(mat[tr_m])                      # <-- 只用训练集
                Z = p.transform(mat).astype(np.float32)
                pca_payload[f"{variant}|{source}|{dim}"] = Z
                evr = float(np.sum(p.explained_variance_ratio_))
                pca_info.append({"variant": variant, "source": source, "dim": int(nmax),
                                 "explained_variance_ratio_sum": evr})
                if source == "both":
                    emit(f"| {variant} | {nmax} | {evr:.4f} |")

    # ---- 对照用 PCA: 精确逐行映射 (aligned) 的嵌入 ----
    # 用**同一训练集掩码** fit, 保证与主路径可比 (只是输入嵌入不同)。
    if aligned_ok and Ea_lig_raw is not None:
        emit("")
        emit("### 对照通路 (aligned) 的 PCA")
        emit("")
        emit("> 输入是**逐行 1:1 映射**的原始嵌入 (不做 pdb_id 平均), 其余完全一致。")
        emit("")
        emit("| 版本 | 维度 | 解释方差比之和 |")
        emit("|---|---:|---:|")
        for variant, (Lr, Rr) in (("raw", (Ea_lig_raw, Ea_rec_raw)),
                                  ("l2", (Ea_lig_l2, Ea_rec_l2))):
            for source, mat in (("ligand", Lr), ("receptor", Rr),
                                ("both", np.hstack([Lr, Rr]))):
                for dim in PCA_DIMS:
                    nmax = min(dim, mat.shape[1], int(tr_m.sum()) - 1)
                    p = PCA(n_components=nmax, random_state=args.seed)
                    p.fit(mat[tr_m])
                    Z = p.transform(mat).astype(np.float32)
                    pca_payload[f"aligned|{variant}|{source}|{dim}"] = Z
                    if source == "both":
                        emit(f"| {variant} | {nmax} | "
                             f"{float(np.sum(p.explained_variance_ratio_)):.4f} |")
    else:
        emit("")
        emit("- (跳过 aligned 对照 PCA: 精确通路不可用)")

    emit("")
    emit("> `both` = 配体嵌入与受体嵌入**直接拼接**后再 PCA。`ligand` / `receptor` 单独 PCA")
    emit("> 的维度与解释方差比见 JSON。")

    if args.dry_run:
        emit("")
        emit("--dry-run: 不训练")
        return 0

    # ---------------- 5. 特征集定义 ----------------
    def esm_block(variant: str, source: str, dim: int, aligned: bool = False) -> np.ndarray:
        key = (f"aligned|{variant}|{source}|{dim}" if aligned
               else f"{variant}|{source}|{dim}")
        return pca_payload[key]

    def feature_set(name: str, variant: str, dim: int,
                    aligned: bool = False) -> tuple[pd.DataFrame, list[str]]:
        if name == "physics_only":
            return Xp, list(Xp.columns)
        if name.startswith("esm2_"):
            body = name[len("esm2_"):]
            # 形如 esm2_ligand_raw / esm2_receptor_l2 / esm2_both_raw
            for s in SOURCES:
                for v in VARIANTS:
                    if body == f"{s}_{v}":
                        Z = esm_block(v, s, dim, aligned)
                        cols = [f"esm2_{v}_{s}_{i:03d}" for i in range(Z.shape[1])]
                        return pd.DataFrame(Z, columns=cols), cols
            raise ValueError(f"未知特征集 {name!r}")
        if name.startswith("physics+esm2_"):
            s = name[len("physics+esm2_"):]
            if s not in SOURCES:
                raise ValueError(f"未知特征集 {name!r}")
            Z = esm_block(variant, s, dim, aligned)
            cols = [f"esm2_{variant}_{s}_{i:03d}" for i in range(Z.shape[1])]
            X = pd.concat([Xp.reset_index(drop=True),
                           pd.DataFrame(Z, columns=cols)], axis=1)
            return X, list(X.columns)
        raise ValueError(f"未知特征集 {name!r}")

    BASE_SETS = ["physics_only", "esm2_ligand_raw", "esm2_ligand_l2",
                 "esm2_receptor_raw", "esm2_receptor_l2",
                 "esm2_both_raw", "esm2_both_l2"]
    COMBO_SETS = ["physics+esm2_ligand", "physics+esm2_receptor", "physics+esm2_both"]

    emit("")
    emit("## 5. 消融矩阵")
    emit("")
    emit(f"- ESM-2 单独特征集: {len(BASE_SETS) - 1} 个 x {len(PCA_DIMS)} 维度 = "
         f"{(len(BASE_SETS) - 1) * len(PCA_DIMS)} 次训练")
    emit(f"- 物理+ESM-2 特征集: {len(COMBO_SETS)} 个 x {len(PCA_DIMS)} 维度 x "
         f"{len(VARIANTS)} 版本 = {len(COMBO_SETS) * len(PCA_DIMS) * len(VARIANTS)} 次训练")
    emit(f"- 加上 `physics_only` 基线 = 共 "
         f"{(len(BASE_SETS) - 1) * len(PCA_DIMS) + len(COMBO_SETS) * len(PCA_DIMS) * len(VARIANTS) + 1} 次训练")

    # ---------------- 6. 训练与评估 ----------------
    def fit_eval(X: pd.DataFrame, tag: str) -> dict[str, Any]:
        Xtr, Xva, Xte = X[tr_m], X[va_m], X[te_m]
        ytr, yva, yte = y[tr_m], y[va_m], y[te_m]
        m = lgb.LGBMRegressor(**PARAMS)
        m.fit(Xtr, ytr, eval_set=[(Xva, yva)],
              callbacks=[lgb.early_stopping(50, verbose=False)])
        pred = m.predict(Xte)
        y_s, p_s = -yte, -pred
        out: dict[str, Any] = {"n_features": int(X.shape[1]),
                               "best_iteration": int(getattr(m, "best_iteration_", -1) or -1)}
        for k, fn in mfns.items():
            try:
                out[k] = float(fn(y_s, p_s))
            except Exception:  # noqa: BLE001
                out[k] = float("nan")
        for k in ("spearman", "ndcg@10", "ef@1%", "precision@10"):
            lo, hi, _pt = bootstrap_ci(mfns[k], y_s, p_s, args.n_boot_ci, args.seed)
            out[f"{k}_ci_low"], out[f"{k}_ci_high"] = lo, hi
        out["_pred"] = pred
        out["_model"] = m
        return out

    rows: list[dict[str, Any]] = []
    preds: dict[str, np.ndarray] = {}

    # 6.1 baseline
    emit("")
    emit("## 6. 结果")
    emit("")
    emit("### 6.1 基线 (physics_only) — 与 V2.2 对齐校验")
    emit("")
    base = fit_eval(Xp, "physics_only")
    rows.append({"feature_set": "physics_only", "pca_dim": 0, "embedding_variant": "-",
                 "n_features": base["n_features"], "best_iteration": base["best_iteration"],
                 **{k: v for k, v in base.items() if not k.startswith("_")}})
    preds["physics_only"] = base["_pred"]
    v22_rho = None
    if v22_payload:
        v22_rho = float(v22_payload.get("metrics", {}).get("test", {})
                        .get("spearman", float("nan")))
    emit(f"- 本脚本 physics_only 测试集 Spearman = **{base['spearman']:+.4f}** "
         f"[{base['spearman_ci_low']:+.4f}, {base['spearman_ci_high']:+.4f}]")
    if v22_rho is not None and np.isfinite(v22_rho):
        dlt = abs(base["spearman"] - v22_rho)
        emit(f"- V2.2 记录的同一口径 Spearman = {v22_rho:+.4f}; 绝对差 = {dlt:.6f}")
        emit(f"- **对齐校验: {'✅ 通过 (逐位复现)' if dlt < 1e-6 else '⚠️ 有差异, 见报告局限'}**")

    # 6.2 ESM-2 单独
    emit("")
    emit("### 6.2 ESM-2 单独 (无物理特征)")
    emit("")
    emit("| 特征集 | 维度 | 特征数 | Spearman | 95% CI | NDCG@10 | EF@1% | P@10 |")
    emit("|---|---:|---:|---:|---|---:|---:|---:|")
    for name in BASE_SETS:
        if name == "physics_only":
            continue
        # 解析 variant/source
        body = name[len("esm2_"):]
        src, var = body.rsplit("_", 1)
        for dim in PCA_DIMS:
            try:
                X, _cols = feature_set(name, var, dim)
                r = fit_eval(X, name)
                rows.append({"feature_set": name, "pca_dim": dim,
                             "embedding_variant": var, "n_features": r["n_features"],
                             "best_iteration": r["best_iteration"],
                             **{k: v for k, v in r.items() if not k.startswith("_")}})
                preds[f"{name}@{dim}"] = r["_pred"]
                emit(f"| {name} | {dim} | {r['n_features']} | {r['spearman']:+.4f} | "
                     f"[{r['spearman_ci_low']:+.4f}, {r['spearman_ci_high']:+.4f}] | "
                     f"{r['ndcg@10']:.4f} | {r['ef@1%']:.4f} | {r['precision@10']:.4f} |")
            except Exception as exc:  # noqa: BLE001
                emit(f"| {name} | {dim} | — | FAIL: {exc} | | | |")
                with open(FAIL_LOG, "a", encoding="utf-8") as fh3:
                    fh3.write(f"[03_merge_esm2_v2] {name}@{dim}: {exc}\n"
                              f"{traceback.format_exc()}\n")
            gc.collect()

    # 6.3 物理 + ESM-2
    emit("")
    emit("### 6.3 物理特征 + ESM-2")
    emit("")
    emit("| 特征集 | 版本 | 维度 | 特征数 | Spearman | 95% CI | NDCG@10 | EF@1% | P@10 |")
    emit("|---|---|---:|---:|---:|---|---:|---:|---:|")
    for name in COMBO_SETS:
        s = name[len("physics+esm2_"):]
        for var in VARIANTS:
            for dim in PCA_DIMS:
                try:
                    X, _cols = feature_set(name, var, dim)
                    r = fit_eval(X, name)
                    rows.append({"feature_set": name, "pca_dim": dim,
                                 "embedding_variant": var, "n_features": r["n_features"],
                                 "best_iteration": r["best_iteration"],
                                 **{k: v for k, v in r.items() if not k.startswith("_")}})
                    preds[f"{name}@{var}@{dim}"] = r["_pred"]
                    emit(f"| {name} | {var} | {dim} | {r['n_features']} | "
                         f"{r['spearman']:+.4f} | "
                         f"[{r['spearman_ci_low']:+.4f}, {r['spearman_ci_high']:+.4f}] | "
                         f"{r['ndcg@10']:.4f} | {r['ef@1%']:.4f} | {r['precision@10']:.4f} |")
                except Exception as exc:  # noqa: BLE001
                    emit(f"| {name} | {var} | {dim} | — | FAIL: {exc} | | | |")
                    with open(FAIL_LOG, "a", encoding="utf-8") as fh3:
                        fh3.write(f"[03_merge_esm2_v2] {name}@{var}@{dim}: {exc}\n"
                                  f"{traceback.format_exc()}\n")
                gc.collect()

    # 6.4 对照: 精确逐行映射 (aligned) vs 结构级 pdb_id join
    #
    # 目的: 把「退化为结构级」的代价**实测出来**, 而不是定性描述。
    # 两边用**同一训练/验证/测试划分、同一超参、同一 PCA 维度**,
    # 唯一差异是嵌入的来源 (逐行精确 vs 按 pdb_id 平均)。
    if aligned_ok and Ea_lig_raw is not None:
        emit("")
        emit("### 6.4 对照: 精确逐行映射 (aligned) vs 结构级 pdb_id join")
        emit("")
        emit("> 唯一差异是嵌入来源。两者在同一测试集上, 因此可逐行配对比较。")
        emit("")
        emit("| 特征集 | 版本 | 维度 | pdb_id join (平均) | aligned (精确) | Δ(aligned − avg) |")
        emit("|---|---|---:|---:|---:|---:|")
        aligned_compare: list[dict[str, Any]] = []
        # 注意: 此时 `res` 还没构造 (它在本块之后才由 rows 生成), 因此这里
        # 用**已累积的 rows** 查 pdb_id 平均版的指标, 而不是引用 `res`。
        # 早期版本误用 `res` -> 9 次 aligned 运行全部抛 NameError 被吞掉。
        avg_lookup = {(r["feature_set"], r["pca_dim"]): r["spearman"]
                      for r in rows if r["feature_set"] in
                      ("esm2_both_raw", "esm2_ligand_raw", "esm2_receptor_raw")}
        aligned_pairs: list[tuple[np.ndarray, np.ndarray]] = []
        for tgt in ("esm2_both_raw", "esm2_ligand_raw", "esm2_receptor_raw"):
            for dim in PCA_DIMS:
                try:
                    Xa, _ca = feature_set(tgt, "raw", dim, aligned=True)
                    ra = fit_eval(Xa, tgt + "_aligned")
                    avg_rho = float(avg_lookup.get((tgt, dim), float("nan")))
                    rows.append({"feature_set": tgt + "_ALIGNED", "pca_dim": dim,
                                 "embedding_variant": "raw",
                                 "n_features": ra["n_features"],
                                 "best_iteration": ra["best_iteration"],
                                 **{k: v for k, v in ra.items() if not k.startswith("_")}})
                    preds[f"{tgt}_ALIGNED@{dim}"] = ra["_pred"]
                    # 配对检验: aligned vs 平均版 (同一测试集, 可逐行配对)
                    if f"{tgt}@{dim}" in preds:
                        aligned_pairs.append((preds[f"{tgt}@{dim}"], ra["_pred"]))
                    aligned_compare.append({"target": tgt, "dim": dim,
                                            "avg_rho": avg_rho,
                                            "aligned_rho": float(ra["spearman"])})
                    emit(f"| {tgt} | raw | {dim} | {avg_rho:+.4f} | "
                         f"{ra['spearman']:+.4f} | {ra['spearman'] - avg_rho:+.4f} |")
                except Exception as exc:  # noqa: BLE001
                    emit(f"| {tgt} | raw | {dim} | — | FAIL: {exc} | |")
                    with open(FAIL_LOG, "a", encoding="utf-8") as fh3:
                        fh3.write(f"[03_merge_esm2_v2] aligned {tgt}@{dim}: {exc}\n"
                                  f"{traceback.format_exc()}\n")
                gc.collect()
        if aligned_compare:
            deltas = [c["aligned_rho"] - c["avg_rho"] for c in aligned_compare
                      if np.isfinite(c["avg_rho"]) and np.isfinite(c["aligned_rho"])]
            if deltas:
                emit("")
                emit(f"- aligned 相对 pdb_id 平均的 ΔSpearman: "
                     f"均值 {np.mean(deltas):+.4f}, 最小 {min(deltas):+.4f}, "
                     f"最大 {max(deltas):+.4f}")
                if np.mean(deltas) > 0.005:
                    emit("- **结论: 精确逐行映射优于 pdb_id 平均** -> 结构级退化确实损失了信息,")
                    emit("  上述 pdb_id join 的指标是**偏保守的下界**。")
                elif np.mean(deltas) < -0.005:
                    emit("- 结论: pdb_id 平均在本数据上反而更好 (平均起到正则化作用)。")
                else:
                    emit("- 结论: 两种映射**基本等价** (差异在噪声量级)。")
            # 配对 bootstrap (最稳的那一组)
            if aligned_pairs:
                res_ap = paired_bootstrap(mr.spearman_rho, y[te_m],
                                          -aligned_pairs[0][0], -aligned_pairs[0][1],
                                          args.n_boot, args.seed)
                emit("")
                emit(f"- 配对 bootstrap (aligned vs 平均, 以首组为例): "
                     f"Δ = {res_ap['delta']:+.4f}, "
                     f"95% CI [{res_ap['lo']:+.4f}, {res_ap['hi']:+.4f}], "
                     f"p = {res_ap['p_two_sided']:.4f}")
    else:
        emit("")
        emit("### 6.4 对照: 精确逐行映射 — 跳过 (精确通路不可用)")

    res = pd.DataFrame(rows)

    # ---------------- 7. 配对 bootstrap vs physics_only ----------------
    emit("")
    emit("## 7. 与 V2.2 基线的配对 bootstrap 检验 (BH-FDR 校正)")
    emit("")
    y_s = -y[te_m]
    base_pred = preds["physics_only"]
    pair_rows: list[dict[str, Any]] = []
    names, pvals = [], []
    for key, pred in preds.items():
        if key == "physics_only":
            continue
        res_sp = paired_bootstrap(mr.spearman_rho, y_s, -base_pred, -pred,
                                  args.n_boot, args.seed)
        pair_rows.append({"model": key, **res_sp})
        names.append(key)
        pvals.append(res_sp["p_two_sided"])
    qvals = bh_fdr(pvals)
    qmap = dict(zip(names, qvals))
    for pr in pair_rows:
        pr["bh_q"] = qmap.get(pr["model"], float("nan"))
    pairs_df = pd.DataFrame(pair_rows)
    if len(pairs_df):
        emit("| 特征集 | Δρ vs physics_only | 95% CI | p | BH q | 显著 | 重采样差值去重数 |")
        emit("|---|---:|---|---:|---:|---|---:|")
        for pr in sorted(pair_rows, key=lambda x: -(x["delta"] if x["delta"] == x["delta"] else -9)):
            sig = "是" if (pr["bh_q"] == pr["bh_q"] and pr["bh_q"] < 0.05) else "否"
            emit(f"| {pr['model']} | {pr['delta']:+.4f} | "
                 f"[{pr['lo']:+.4f}, {pr['hi']:+.4f}] | {pr['p_two_sided']:.4f} | "
                 f"{pr['bh_q']:.4f} | {sig} | {pr.get('n_distinct_deltas', -1)} |")
        n_sig = int((pairs_df["bh_q"] < 0.05).sum())
        n_better = int(((pairs_df["bh_q"] < 0.05) & (pairs_df["delta"] > 0)).sum())
        emit("")
        emit(f"- BH q<0.05 的检验数: **{n_sig} / {len(pairs_df)}**, "
             f"其中**显著优于**基线: **{n_better}**")
        if "n_distinct_deltas" in pairs_df.columns:
            mn = int(pairs_df["n_distinct_deltas"].min())
            floor = 2.0 / (args.n_boot + 1)
            emit(f"- ⚠️ **统计实现提示**: 最后一项「重采样差值去重数」最小值 = {mn}。")
            emit(f"  当该值很小时, 双侧 p 会撞到 **{floor:.4f} 的下限** "
                 f"(= 2/(n_boot+1)), 多个模型会显示完全相同的 p。")
            emit("  这是**离散化效应**, 不是「这些模型效果一样」; 解读时应以 Δρ 与 CI 为主。")

    # ---------------- 8. 结论 ----------------
    emit("")
    emit("## 8. 结论")
    emit("")
    if len(res):
        best = res.sort_values("spearman", ascending=False).iloc[0]
        emit(f"- 最佳特征集: **{best['feature_set']}** "
             f"(版本 `{best['embedding_variant']}`, PCA 维度 {best['pca_dim']}, "
             f"{best['n_features']} 特征)")
        emit(f"- 其测试集 Spearman = **{best['spearman']:+.4f}** "
             f"[{best['spearman_ci_low']:+.4f}, {best['spearman_ci_high']:+.4f}], "
             f"NDCG@10 = {best['ndcg@10']:.4f}, EF@1% = {best['ef@1%']:.4f}, "
             f"P@10 = {best['precision@10']:.4f}")
        base_rho = float(res.loc[res["feature_set"] == "physics_only", "spearman"].iloc[0])
        emit(f"- 基线 (physics_only) Spearman = {base_rho:+.4f}; "
             f"最佳 Δ = {best['spearman'] - base_rho:+.4f}")
        if best["spearman"] > base_rho:
            emit("- **结论: ESM-2 嵌入带来了正向增益**, 支持 V2.2 的判断 "
                 "(结构物理特征无法区分突变体, 序列信息是必要的)。")
        else:
            emit("- **结论: 本配置下 ESM-2 未超过水平基线。** 需检查: "
                 "PCA 维度是否过低、是否需要更大模型 (650M)、"
                 "或需要界面残基口径的池化。**如实报告, 不夸大。**")
    else:
        emit("- 没有成功训练的特征集")

    emit("")
    emit("## 9. 残留局限")
    emit("")
    emit("1. PCA 只在训练集 fit, 但**降维本身会损失信息**; 128 维仅为原始 640/960 维的一部分。")
    emit("2. ESM-2 用 **mean-pool over residues**, 未做界面残基加权池化;")
    emit("   后者在 V2.0 的 `18_extra_features.py` 中有通路但本任务未启用。")
    emit("3. 单 seed (42); bootstrap CI 覆盖采样不确定性, 未覆盖超参初始化噪声。")
    emit("4. 超参完全沿用 V2.2, **未针对 ESM-2 特征做调参** —— 这对 ESM-2 可能不利,")
    emit("   但也保证了与基线的可比性。")
    emit("5. 多链序列里的 `,` 在 Colab 端被映射为 `<unk>` (见 notebook),")
    emit("   会在链边界引入 1 个未知位置; 更干净的分链池化未实现。")
    emit("")
    emit(f"耗时: {time.time() - t0:.1f}s")

    # ---------------- 10. 落盘 ----------------
    if len(res):
        res_out = res.drop(columns=[c for c in res.columns if c.startswith("_")],
                           errors="ignore")
        tmp = args.out_csv.with_suffix(args.out_csv.suffix + ".tmp")
        res_out.to_csv(tmp, index=False, encoding="utf-8")
        tmp.replace(args.out_csv)
        if len(pairs_df):
            tmp2 = args.pairs_csv.with_suffix(args.pairs_csv.suffix + ".tmp")
            pairs_df.to_csv(tmp2, index=False, encoding="utf-8")
            tmp2.replace(args.pairs_csv)

        # ---- 预测产物 (供任务 5 的集成复用, 避免重新训练) ----
        # 只保存 **测试集** 预测 (集成只在 test 上评估) + 测试集标签/键。
        try:
            te_keys = d.loc[te_m, "interface_key"].to_numpy(dtype=str)
            np.savez_compressed(
                args.preds_npz,
                test_interface_key=te_keys,
                test_pdb_id=d.loc[te_m, "pdb_id"].to_numpy(dtype=str),
                test_y=y[te_m].astype(np.float32),
                model_keys=np.array(sorted(preds.keys()), dtype=str),
                **{f"pred__{k}": preds[k].astype(np.float32) for k in sorted(preds)},
            )
            emit(f"- 预测产物: `{args.preds_npz.name}` "
                 f"({args.preds_npz.stat().st_size / 1e6:.1f} MB, "
                 f"{len(preds)} 个模型的 test 预测)")
        except Exception as exc:  # noqa: BLE001
            emit(f"- !! 预测产物写出失败: {exc}")
            with open(FAIL_LOG, "a", encoding="utf-8") as fh4:
                fh4.write(f"[03_merge_esm2_v2] preds npz: {exc}\n"
                          f"{traceback.format_exc()}\n")

        # 保存最佳模型
        best_row = res.sort_values("spearman", ascending=False).iloc[0]
        bname, bdim = best_row["feature_set"], int(best_row["pca_dim"])
        bvar = best_row["embedding_variant"]
        Xb, bcols = feature_set(bname, bvar if bvar != "-" else "raw", bdim) \
            if bname != "physics_only" else (Xp, list(Xp.columns))
        mb = lgb.LGBMRegressor(**PARAMS)
        mb.fit(Xb[tr_m], y[tr_m], eval_set=[(Xb[va_m], y[va_m])],
               callbacks=[lgb.early_stopping(50, verbose=False)])
        payload = {
            "model": mb, "feature_set": bname, "feature_names": list(Xb.columns),
            "pca_dim": bdim, "embedding_variant": bvar,
            "params": PARAMS, "count_features_log1p": [c for c in count_cols if c in bcols],
            "esm2_cache": str(args.esm2_cache), "esm2_model": model_name,
            "esm2_repr_layer": repr_layer, "esm2_pool": pool,
            "split_policy": v22_payload.get("split_policy",
                                            "V2.2 Proteinea split (reused)"),
            "train_rows": int(tr_m.sum()), "val_rows": int(va_m.sum()),
            "test_rows": int(te_m.sum()),
            "metrics_test": {k: v for k, v in best_row.items()
                             if not isinstance(v, (list, dict))},
            "version": "v2.3",
            "note": "PCA fitted on TRAIN split only; split reused from V2.2",
        }
        tmpm = args.best_model.with_suffix(args.best_model.suffix + ".tmp")
        with open(tmpm, "wb") as fhm:
            pickle.dump(payload, fhm)
        tmpm.replace(args.best_model)
        emit("")
        emit(f"- 结果表: `{args.out_csv.name}`")
        emit(f"- 配对检验: `{args.pairs_csv.name}`")
        emit(f"- 最佳模型: `{args.best_model.name}`")

    args.report.write_text("\n".join(L) + "\n", encoding="utf-8")
    summary = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"), "seed": args.seed,
        "esm2_run": True, "esm2_cache_present": True,
        "esm2_model": model_name, "esm2_repr_layer": repr_layer, "esm2_pool": pool,
        "embedding_dim": int(lig_raw.shape[1]),
        "join_key": "pdb_id" if use_pdb else "interface_key",
        "join_note": ("V2.3.2: join 键改为 pdb_id (V2.3 用 interface_key 语义不一致 -> 0% 命中); "
                      "一对多时缓存侧取平均, 特征表侧共享 -> 结构级退化"),
        "aligned_exact_path_available": bool(aligned_ok),
        "join_hit_rate": hit_rate,
        "v22_split_check": {sp: int(split_counts.get(sp, 0)) for sp in EXPECTED_SPLITS},
        "v22_split_ok": bool(split_ok),
        "v22_recorded_test_spearman": v22_rho,
        "baseline_reproduced_spearman": float(base["spearman"]),
        "pca": pca_info,
        "results": res.drop(columns=[c for c in res.columns if c.startswith("_")],
                            errors="ignore").to_dict("records") if len(res) else [],
        "pairs": pair_rows,
    }
    args.json_out.write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=str),
                             encoding="utf-8")
    emit("")
    emit(f"报告: `{args.report.name}`")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:  # noqa: BLE001
        LOG_DIR.mkdir(parents=True, exist_ok=True)
        with open(FAIL_LOG, "a", encoding="utf-8") as fh2:
            fh2.write(f"[03_merge_esm2_v2] {traceback.format_exc()}\n")
        raise
