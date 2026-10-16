# -*- coding: utf-8 -*-
"""scripts/v2.5/09_enzyme_data_prep.py -- V2.5 任务 6.2: 酶数据接入骨架（**不下载**）.

设计原则
--------
1. **只写框架，不实际下载**（任务书要求）。所有下载函数带 `dry_run=True` 默认值，
   调用即打印计划并返回空 DataFrame。
2. **统一 schema**：任何来源都必须映射到同一套列，否则无法合并。
3. **三层隔离**：EC 号 > CATH > 40% 序列同一性（与 V2.0 `13_enzyme_split.py` 同精神）。
4. **Km 绝不进入主回归**：用断言强制，而不是靠注释提醒。

运行方式
--------
```powershell
# 只做自检与 schema 演示 (不联网, 不下载)
python scripts/v2.5/09_enzyme_data_prep.py --self-test
```

主流程 `main()` 里所有数据获取都被包在 `if not args.dry_run` 之外 ——
**默认 dry_run**，因此直接运行也不会产生任何网络请求。
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
FEAT = ROOT / "features"
LOG_DIR = ROOT / "logs" / "v2.5"
LOG_DIR.mkdir(parents=True, exist_ok=True)

RT_LITERAL = 0.592          # kcal/mol @ 298.15 K (与 V2.2/V2.3.2 一致)

# ----------------------------------------------------------------- schema
UNIFIED_SCHEMA = [
    "source",           # 数据来源标识 (sabio_rk / catpred / skempi / v22_enzyme / ...)
    "pdb_id",           # 结构 ID (可空)
    "chain",            # 链 ID (可空)
    "sequence",         # 酶序列 (单字母)
    "ec_number",        # EC 号 (形如 3.4.21.4; 可空)
    "uniprot_id",       # UniProt AC (可空)
    "ligand_id",        # 底物/抑制剂标识 (可空)
    "label_type",       # **Kd | Ki | Km | kcat | kcat_Km | ddG** —— 决定能否进主回归
    "value",            # 原始数值 (单位见 unit)
    "unit",             # M / s^-1 / M^-1 s^-1 / kcal/mol
    "temperature_K",    # 实验温度 (可空 -> 用 298.15)
    "ph",               # 实验 pH (可空)
    "assay_method",     # 测定方法 (可空)
    "real_dG",          # **换算后的 ΔG (kcal/mol), 仅对 Kd/Ki/ddG 有定义**
    "group_id",         # 组内排序的单元 (同一酶+底物)
    "family_id",        # 家族 (用于隔离: EC 或 CATH 或序列簇)
]

# 允许进入**亲和力主回归**的标签类型
AFFINITY_LABELS = {"Kd", "Ki", "ddG"}
# 明确**禁止**进入主回归的标签类型 (动力学)
KINETIC_LABELS = {"Km", "kcat", "kcat_Km"}


# ----------------------------------------------------------- 换算与断言
def kd_to_dG(k: float, temperature_K: float | None = None,
             rt_literal: float = RT_LITERAL) -> float:
    """把 Kd/Ki (M) 换算为 ΔG (kcal/mol)。

    约定与 V2.2 完全一致: 用**固定** RT = 0.592 kcal/mol (298.15 K),
    不用每个样本的实际温度 —— 这是为了与既有标签体系可比。

    >>> round(kd_to_dG(1e-9), 4)
    -8.1747
    """
    k = float(k)
    if not np.isfinite(k) or k <= 0:
        return float("nan")
    return rt_literal * float(np.log(k))


def dG_to_kd(dg: float, rt_literal: float = RT_LITERAL) -> float:
    """ΔG (kcal/mol) -> Kd (M)。与 `kd_to_dG` 互为逆运算。"""
    return float(np.exp(float(dg) / rt_literal))


def assert_affinity_only(df: pd.DataFrame, label_col: str = "label_type") -> None:
    """**硬断言**: 主回归数据里不得出现动力学标签。

    这是本项目的红线之一（V2.2 已因 KD/Ki/Km 混用做过标签纯化）。
    任何把 Km/kcat 混入亲和力回归的行为都会在这里直接抛错。
    """
    bad = set(df[label_col].astype(str).unique()) & KINETIC_LABELS
    if bad:
        raise ValueError(
            f"亲和力回归数据里出现动力学标签 {sorted(bad)} —— 禁止混用。"
            f"动力学数据必须走独立的第二回归。"
        )


def build_group_family_ids(df: pd.DataFrame) -> pd.DataFrame:
    """生成 `group_id` 与 `family_id`。

    * `group_id`: 同一「酶 + 底物」下可相互排序的单元
      （优先 uniprot_id+ligand_id；退化到 sequence 哈希 + ligand_id）
    * `family_id`: 隔离用（优先 EC 号；缺失则用序列前缀占位并**标记**）
    """
    out = df.copy()
    uni = out.get("uniprot_id", pd.Series([""] * len(out))).fillna("").astype(str)
    lig = out.get("ligand_id", pd.Series([""] * len(out))).fillna("").astype(str)
    seq = out.get("sequence", pd.Series([""] * len(out))).fillna("").astype(str)
    ec = out.get("ec_number", pd.Series([""] * len(out))).fillna("").astype(str)

    key = np.where(uni != "", uni, seq.str.slice(0, 40))
    out["group_id"] = [f"{a}|{b}" for a, b in zip(key, lig)]
    # family: EC 优先; 无 EC 时用序列(前 40 aa)占位, 但要能看出是占位
    out["family_id"] = np.where(ec != "", "EC:" + ec,
                                "SEQ:" + seq.str.slice(0, 40))
    out["family_is_ec"] = (ec != "")
    return out


def three_layer_isolation(df: pd.DataFrame, seq_identity: float = 0.40
                          ) -> dict[str, Any]:
    """三层隔离: **EC 号 > CATH > 40% 序列同一性**。

    返回每层的统计与最终 `isolation_key`。
    实现说明:
    * 第 1 层 (EC): 有 EC 号的样本按 EC 分组 -> 保证**不跨 EC**;
    * 第 2 层 (CATH): 本机**没有** CATH 注释（V2.0 实测覆盖 0.0%）,
      因此该层**跳过并记录**, 不假装做了;
    * 第 3 层 (序列): 对**没有 EC**的样本用 40% 同一性贪心聚类兜底。
    """
    out = df.copy()
    ec = out.get("ec_number", pd.Series([""] * len(out))).fillna("").astype(str)
    has_ec = (ec != "")
    iso = np.array(["UNASSIGNED"] * len(out), dtype=object)
    iso[has_ec.to_numpy()] = ("EC:" + ec[has_ec]).to_numpy()

    n_cath = 0
    cath_note = ("本机无 CATH 注释 (V2.0 实测覆盖率 0.0%) -> **第二层跳过**, "
                 "不伪装已执行")

    # 第 3 层: 序列聚类兜底
    no_ec = ~has_ec.to_numpy()
    if no_ec.any() and "sequence" in out.columns:
        seqs = out.loc[no_ec, "sequence"].fillna("").astype(str).to_numpy()
        cluster = greedy_cluster_kmers(seqs, identity=seq_identity)
        for i, c in zip(np.where(no_ec)[0], cluster):
            iso[i] = f"SEQCLUST:{c}"
    out["isolation_key"] = iso
    return {
        "n_rows": int(len(out)),
        "n_ec_assigned": int(has_ec.sum()),
        "n_seqclust_assigned": int((~has_ec).sum() and no_ec.sum()),
        "n_cath_used": n_cath,
        "cath_note": cath_note,
        "n_unique_isolation_key": int(len(set(iso))),
    }


def greedy_cluster_kmers(seqs: np.ndarray, identity: float = 0.40,
                         k: int = 3) -> list[int]:
    """3-mer Jaccard 贪心聚类 (与 V2.0 的家族隔离同口径)。

    输入顺序即优先级; 返回每个序列的簇号。O(n·c)。
    """
    def kmers(s: str) -> set:
        return {s[i:i + k] for i in range(max(0, len(s) - k + 1))}

    reps: list[set] = []
    labels: list[int] = []
    for s in seqs:
        ks = kmers(s)
        if not ks:
            labels.append(-1)
            continue
        best, best_j = -1, 0.0
        for ci, r in enumerate(reps):
            inter = len(ks & r)
            if inter == 0:
                continue
            j = inter / len(ks | r)
            if j > best_j:
                best_j, best = j, ci
        if best >= 0 and best_j >= identity:
            labels.append(best)
        else:
            reps.append(ks)
            labels.append(len(reps) - 1)
    return labels


# ----------------------------------------------------------- 下载骨架 (不执行)
def fetch_sabio_rk(dry_run: bool = True, **kw) -> pd.DataFrame:
    """SABIO-RK 接入骨架。

    **重要**: 只取 `Ki` / `Kd` 记录; `Km` / `kcat` 必须显式走另一条通路。

    真实实现 (未执行, 需网络):
        GET https://sabiork.h-its.org/sabioRestWebServices/searchKineticLaws/
        参数: q=<query>, format=json
    许可: 学术免费; **原始数据不得再分发** -> 只落地 `features/v2.5_*`
    """
    print("[fetch_sabio_rk] DRY-RUN —— 不联网、不下载。")
    print("  计划: 查询 Ki/Kd 子集 -> 映射 Unified schema -> 丢弃 Km/kcat")
    print("  许可: 学术免费, 原始数据禁止再分发")
    print("  注意: 必须带 uniprot_id 与 ec_number, 否则无法做三层隔离")
    return pd.DataFrame(columns=UNIFIED_SCHEMA)


def fetch_catpred_db(dry_run: bool = True, **kw) -> pd.DataFrame:
    """CatPred-DB 接入骨架（~1.2 万条 Ki）。

    ⚠️ **先核实许可与数据血缘**: 若其 Ki 派生自 BRENDA, 则受 BRENDA 许可限制。
    本项目**不预先下载**。
    """
    print("[fetch_catpred_db] DRY-RUN —— 不联网、不下载。")
    print("  前置条件: 确认发行渠道与许可 (若派生自 BRENDA -> 不可用)")
    return pd.DataFrame(columns=UNIFIED_SCHEMA)


def fetch_v22_enzyme_subset(dry_run: bool = True, **kw) -> pd.DataFrame:
    """**优先方案**: 不下载任何新数据, 从现有 V2.2 表里筛酶条目。

    这一步是**唯一零许可风险**的通路, 也是推荐的起点:
    用 UniProt/EC 注释 (可用 UniProt REST 批量取) 判断哪些 PDB 是酶。
    """
    print("[fetch_v22_enzyme_subset] DRY-RUN —— 只读本地 V2.2 表, 不联网。")
    v22 = ROOT / "features" / "v2.2_proteinea_features.csv"
    print(f"  计划: 读 `{v22.name}` -> 按 UniProt/EC 筛酶条目 -> 映射 Unified schema")
    print("  优点: 零许可风险; 复用已有特征与划分")
    if dry_run or not v22.exists():
        return pd.DataFrame(columns=UNIFIED_SCHEMA)
    return pd.read_csv(v22, low_memory=False)


# ----------------------------------------------------------- self-test
def self_test() -> int:
    print("=== V2.5 任务 6.2 骨架自检 ===")
    print()
    ok = True

    # 1. schema 完整性
    need = {"source", "label_type", "value", "unit", "real_dG", "group_id",
            "family_id", "ec_number", "uniprot_id", "temperature_K"}
    miss = need - set(UNIFIED_SCHEMA)
    print(f"1) schema 字段数 {len(UNIFIED_SCHEMA)}; 缺失关键字段: {sorted(miss) or '无'}")
    ok &= not miss

    # 2. Kd -> dG 往返
    for k in (1e-3, 1e-6, 1e-9, 1e-12):
        dg = kd_to_dG(k)
        back = dG_to_kd(dg)
        rel = abs(back - k) / k
        print(f"2) Kd={k:.0e} -> ΔG={dg:+.4f} -> Kd={back:.3e} (相对误差 {rel:.1e})")
        ok &= rel < 1e-9

    # 3. Km 必须被拒绝
    bad = pd.DataFrame({"label_type": ["Kd", "Km"], "value": [1e-9, 1e-6]})
    try:
        assert_affinity_only(bad)
        print("3) !! 断言**没有**拦住 Km —— 这是 bug")
        ok = False
    except ValueError as e:
        print(f"3) ✅ 断言正确拦住了 Km: {str(e)[:60]}...")

    good = pd.DataFrame({"label_type": ["Kd", "Ki"], "value": [1e-9, 1e-6]})
    try:
        assert_affinity_only(good)
        print("4) ✅ Kd/Ki 通过断言")
    except ValueError as e:
        print(f"4) !! Kd/Ki 被误拦: {e}")
        ok = False

    # 5. group/family 构造
    demo = pd.DataFrame({
        "uniprot_id": ["P001", "P001", "P002"],
        "ligand_id": ["L1", "L1", "L2"],
        "sequence": ["ACDEFGHIKL", "ACDEFGHIKL", "MNPQRSTVWY"],
        "ec_number": ["3.4.21.4", "3.4.21.4", ""],
    })
    d2 = build_group_family_ids(demo)
    print(f"5) group_id: {d2['group_id'].tolist()}")
    print(f"   family_id: {d2['family_id'].tolist()}")
    ok &= (d2.loc[0, "group_id"] == d2.loc[1, "group_id"])
    ok &= (d2.loc[0, "family_id"] == d2.loc[1, "family_id"])
    ok &= (d2.loc[2, "family_id"] != d2.loc[0, "family_id"])

    # 6. 三层隔离 (含 CATH 缺失的诚实记录)
    iso = three_layer_isolation(demo)
    print(f"6) 隔离统计: {iso}")
    ok &= iso["n_ec_assigned"] == 2
    ok &= iso["n_cath_used"] == 0
    print(f"   CATH 说明: {iso['cath_note']}")

    # 7. 序列聚类
    lab = greedy_cluster_kmers(np.array(["ACDEFGHIKL", "ACDEFGHIKM", "MNPQRSTVWY"]),
                               identity=0.5)
    print(f"7) 3-mer 聚类标签: {lab} (前两条应同簇)")
    ok &= lab[0] == lab[1]
    ok &= lab[2] != lab[0]

    # 8. 下载函数在 dry_run 下不联网
    for fn in (fetch_sabio_rk, fetch_catpred_db, fetch_v22_enzyme_subset):
        r = fn(dry_run=True)
        ok &= isinstance(r, pd.DataFrame) and len(r) == 0

    print()
    print(f"=== 自检结果: {'全部通过 ✅' if ok else '存在失败 ❌'} ===")
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser(description="V2.5 酶数据接入骨架 (默认不下载)")
    ap.add_argument("--self-test", action="store_true", help="只跑骨架自检")
    ap.add_argument("--dry-run", action="store_true", default=True,
                    help="默认 True —— 不做任何网络请求")
    ap.add_argument("--execute", action="store_true",
                    help="**显式**关闭 dry_run (仍不会自动下载; 需自行改实现)")
    args = ap.parse_args()

    if args.self_test or not args.execute:
        return self_test()
    print("真实下载需要先实现各 fetch_* 的网络部分; 当前**未实现、未执行**。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
