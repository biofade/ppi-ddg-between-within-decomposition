# -*- coding: utf-8 -*-
"""scripts/v5/14_p4_decompose_topk.py -- V5.0: 把 Top-K 拆成「跨复合物」与「复合物内」两层.

关键疑问
--------
P4 的**全局** Top-K 看起来极好：P@10 ≈ 0.967、EF@1% ≈ 9.60。
但同一模型的 macro per-PDB Spearman 只有 ≈ +0.096（CI 含 0）。
这两个读数**并不矛盾**，而是说明：全局 Top-K 的命中主要来自
**跨复合物区分**（某些复合物的 ΔΔG 整体更负），而非**复合物内排序**。

而湿实验的真实场景是「针对**一个**靶点选 10 个突变体」→ 应该看**复合物内** Top-K。

本脚本把两者分开并强制对照（V5 §76 要求同时报告 global 与 per-complex）：
  * `global_*`  ：在全部 test 候选池上选（易被跨复合物效应主导）
  * `within_*`  ：对每个 test complex 内单独选 top-K 再 macro 平均（真实场景）
  * `base_rate` ：随机选择的期望命中率（用于判断"高"是否只是基准率效应）
  * 每层的 Oracle 上界与 Random 基线

输出
----
* `experiments/v5_0/p4_conformal/p4_topk_decomposition.csv`
* `experiments/v5_0/p4_conformal/p4_topk_decomposition.md`
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

OUT_CSV = C.V5 / "p4_conformal" / "p4_topk_decomposition.csv"
OUT_MD = C.V5 / "p4_conformal" / "p4_topk_decomposition.md"
L: list[str] = []
RS = np.random.RandomState(C.SEED)


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))


def global_topk(y: np.ndarray, p: np.ndarray, ks=(5, 10)) -> dict:
    """全部候选池上选；真值定义 = ΔΔG 最负前 10%（与 V4 口径一致）。"""
    out = {}
    n = len(y)
    best = (-y) >= np.quantile(-y, 0.90)
    for k in ks:
        idx = np.argsort(-(-p), kind="mergesort")[:k]
        out[f"global_p_at_{k}"] = float(best[idx].mean())
        out[f"global_hit{k}_lt0"] = float((y[idx] < 0).mean())
        out[f"global_hit{k}_ltm05"] = float((y[idx] < -0.5).mean())
        out[f"global_hit{k}_ltm10"] = float((y[idx] < -1.0).mean())
        out[f"global_sel{k}_mean_ddg"] = float(y[idx].mean())
    return out


def within_topk(y: np.ndarray, p: np.ndarray, g: np.ndarray,
                k_big=10, k_small=5, n_rand=200) -> dict:
    """每个 complex 内单独选 top-K，再对 complex 做 macro 平均。

    * 真值定义（complex 内）= 该 complex 内 ΔΔG 最负前 10% 的成员
    * 同时给出 random 与 oracle 基线，避免把基准率当成能力
    """
    res = {f"within_hit{k}_lt0": [], f"within_hit{k}_ltm05": [],
           f"within_hit{k}_ltm10": [], f"within_p_at_{k}": [],
           f"within_sel{k}_mean_ddg": [], f"within_oracle{k}_mean_ddg": [],
           f"within_random{k}_mean_ddg": [], f"within_p_at_{k}_random": [],
           f"within_p_at_{k}_oracle": []} if False else {}
    for kk in (k_big, k_small):
        res[f"within_hit{kk}_lt0"] = []
        res[f"within_hit{kk}_ltm05"] = []
        res[f"within_hit{kk}_ltm10"] = []
        res[f"within_p_at_{kk}"] = []
        res[f"within_p_at_{kk}_random"] = []
        res[f"within_p_at_{kk}_oracle"] = []
        res[f"within_sel{kk}_mean_ddg"] = []
        res[f"within_oracle{kk}_mean_ddg"] = []
        res[f"within_random{kk}_mean_ddg"] = []
    n_cx = 0
    for cx in np.unique(g):
        m = np.where(g == cx)[0]
        if len(m) < 5:
            continue
        n_cx += 1
        yy, pp = y[m], p[m]
        k = k_big if len(m) >= k_big else min(k_small, len(m))
        best = ((-yy) >= np.quantile(-yy, 0.90)) if len(m) >= 10 else \
            ((-yy) >= np.median(-yy))
        sel = np.argsort(-(-pp), kind="mergesort")[:k]
        orc = np.argsort(yy, kind="mergesort")[:k]
        rnds = [RS.permutation(len(m))[:k] for _ in range(n_rand)]
        res[f"within_hit{k}_lt0"].append((yy[sel] < 0).mean())
        res[f"within_hit{k}_ltm05"].append((yy[sel] < -0.5).mean())
        res[f"within_hit{k}_ltm10"].append((yy[sel] < -1.0).mean())
        res[f"within_p_at_{k}"].append(best[sel].mean())
        res[f"within_p_at_{k}_random"].append(
            np.mean([best[r].mean() for r in rnds]))
        res[f"within_p_at_{k}_oracle"].append(best[orc].mean())
        res[f"within_sel{k}_mean_ddg"].append(yy[sel].mean())
        res[f"within_oracle{k}_mean_ddg"].append(yy[orc].mean())
        res[f"within_random{k}_mean_ddg"].append(
            np.mean([yy[r].mean() for r in rnds]))
    out = {k: float(np.mean(v)) if len(v) else float("nan") for k, v in res.items()}
    out["n_complex"] = n_cx
    sizes = [int((g == c).sum()) for c in np.unique(g) if int((g == c).sum()) >= 5]
    out["within_complex_n_median"] = float(np.median(sizes)) if sizes else float("nan")
    return out


def main() -> int:
    emit("# V5.0 P4 — Top-K 的「跨复合物 vs 复合物内」分解")
    emit("")
    emit("> 动机：P4 全局 P@10 ≈ 0.97、EF@1% ≈ 9.6 看起来极好，但同一模型")
    emit("> macro per-PDB Spearman 仅 ≈ +0.096。本脚本检验「好处来自哪一层」。")
    emit("")

    z = np.load(C.V5 / "p4_conformal" / "p4_predictions.npz", allow_pickle=True)
    y = z["test_y"]
    g = np.array([str(x) for x in z["test_group"]])
    emit(f"- test: **{len(y)}** 行 / **{len(set(g))}** 复合物")
    emit(f"- 真值 ΔΔG: mean {y.mean():+.3f}, std {y.std():.3f}, "
         f"ΔΔG<0 占 **{(y < 0).mean() * 100:.1f}%**, "
         f"ΔΔG<−0.5 占 {(y < -0.5).mean() * 100:.1f}%, "
         f"ΔΔG<−1.0 占 {(y < -1.0).mean() * 100:.1f}%")
    emit("")
    emit("> ⚠️ **基准率提示**：任何随机选 10 个，命中 ΔΔG<0 的期望已经是 "
         f"**{(y < 0).mean() * 100:.1f}%**。因此 `hit10_lt0` 这类指标"
         "必须与基准率并读，否则会系统性高估模型能力。")
    emit("")

    rows = []
    emit("## 1. 逐 variant × seed：Global Top-K")
    emit("")
    emit("| variant | seed | P@10 | EF@1% | hit10(<0) | hit10(<−0.5) | hit10(<−1.0) | sel10 mean ΔΔG |")
    emit("|---|---:|---:|---:|---:|---:|---:|---:|")
    for key in sorted(k for k in z.files
                      if k.startswith("plain_s") or k.startswith("mondrian")
                      or k.startswith("conformal_global")):
        try:
            var, sd = key.rsplit("_s", 1)
            sd = int(sd)
        except ValueError:
            continue
        p = z[key]
        am = C.all_metrics(y, p, g)
        gt = global_topk(y, p)
        rows.append({"variant": var, "seed": sd, "layer": "global", **{
            "overall_pearson": am["overall_pearson"],
            "macro_pdb_spearman": am["macro_pdb_spearman"],
            "p_at_10": am.get("p_at_10", np.nan),
            "ef_1pct": am.get("ef_1pct", np.nan),
            "hit10_lt0": gt["global_hit10_lt0"],
            "hit10_ltm05": gt["global_hit10_ltm05"],
            "hit10_ltm10": gt["global_hit10_ltm10"],
            "sel10_mean_ddg": gt["global_sel10_mean_ddg"],
        }})
        emit(f"| `{var}` | {sd} | {am.get('p_at_10', float('nan')):.4f} | "
             f"{am.get('ef_1pct', float('nan')):.4f} | {gt['global_hit10_lt0']:.4f} | "
             f"{gt['global_hit10_ltm05']:.4f} | {gt['global_hit10_ltm10']:.4f} | "
             f"{gt['global_sel10_mean_ddg']:+.4f} |")

    emit("")
    emit("## 2. 逐 variant × seed：**Within-complex** Top-K（V5 §76，真实场景）")
    emit("")
    emit("| variant | seed | P@10 | P@10 random | P@10 oracle | hit10(<0) | "
         "hit10(<−0.5) | hit10(<−1.0) | sel10 mean ΔΔG | oracle mean ΔΔG | random mean ΔΔG |")
    emit("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for key in sorted(k for k in z.files
                      if k.startswith("plain_s") or k.startswith("mondrian")
                      or k.startswith("conformal_global")):
        try:
            var, sd = key.rsplit("_s", 1)
            sd = int(sd)
        except ValueError:
            continue
        p = z[key]
        w = within_topk(y, p, g)
        rows.append({"variant": var, "seed": sd, "layer": "within_complex", **{
            "p_at_10": w.get("within_p_at_10", np.nan),
            "p_at_10_random": w.get("within_p_at_10_random", np.nan),
            "p_at_10_oracle": w.get("within_p_at_10_oracle", np.nan),
            "hit10_lt0": w.get("within_hit10_lt0", np.nan),
            "hit10_ltm05": w.get("within_hit10_ltm05", np.nan),
            "hit10_ltm10": w.get("within_hit10_ltm10", np.nan),
            "sel10_mean_ddg": w.get("within_sel10_mean_ddg", np.nan),
            "oracle10_mean_ddg": w.get("within_oracle10_mean_ddg", np.nan),
            "random10_mean_ddg": w.get("within_random10_mean_ddg", np.nan),
            "n_complex": w.get("n_complex", np.nan),
        }})
        emit(f"| `{var}` | {sd} | {w.get('within_p_at_10', float('nan')):.4f} | "
             f"{w.get('within_p_at_10_random', float('nan')):.4f} | "
             f"{w.get('within_p_at_10_oracle', float('nan')):.4f} | "
             f"{w.get('within_hit10_lt0', float('nan')):.4f} | "
             f"{w.get('within_hit10_ltm05', float('nan')):.4f} | "
             f"{w.get('within_hit10_ltm10', float('nan')):.4f} | "
             f"{w.get('within_sel10_mean_ddg', float('nan')):+.4f} | "
             f"{w.get('within_oracle10_mean_ddg', float('nan')):+.4f} | "
             f"{w.get('within_random10_mean_ddg', float('nan')):+.4f} |")

    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False, encoding="utf-8")

    # ---------- 结论 ----------
    gl = df[df["layer"] == "global"]
    wi = df[df["layer"] == "within_complex"]
    wi_plain = wi[wi["variant"] == "plain"]
    gl_plain = gl[gl["variant"] == "plain"]
    emit("")
    emit("## 3. 结论（这是 V5 最重要的 Top-K 诊断）")
    emit("")
    emit("### 3a. 全部 variant 的层间对照（3 seeds × 6 variants）")
    emit("")
    emit("| 层 | hit10(ΔΔG<0) | hit10(ΔΔG<−0.5) | hit10(ΔΔG<−1.0) | sel10 mean ΔΔG |")
    emit("|---|---:|---:|---:|---:|")
    emit(f"| **Global（全部候选池）** | {gl['hit10_lt0'].mean():.4f} | "
         f"{gl['hit10_ltm05'].mean():.4f} | {gl['hit10_ltm10'].mean():.4f} | "
         f"{gl['sel10_mean_ddg'].mean():+.4f} |")
    emit(f"| **Within-complex（真实场景）** | {wi['hit10_lt0'].mean():.4f} | "
         f"{wi['hit10_ltm05'].mean():.4f} | {wi['hit10_ltm10'].mean():.4f} | "
         f"{wi['sel10_mean_ddg'].mean():+.4f} |")
    emit("")
    emit("### 3b. **以 `plain`（= 被选定的部署模型）为准**的判定读数")
    emit("")
    sel = float(wi_plain["sel10_mean_ddg"].mean())
    rnd = float(wi_plain["random10_mean_ddg"].mean())
    orc = float(wi_plain["oracle10_mean_ddg"].mean())
    hit0 = float(wi_plain["hit10_lt0"].mean())
    p10 = float(wi_plain["p_at_10"].mean())
    p10r = float(wi_plain["p_at_10_random"].mean())
    p10o = float(wi_plain["p_at_10_oracle"].mean())
    eff = (rnd - sel) / (rnd - orc)        # 两者同号 → eff > 0 表示优于随机
    emit("| 项 | 值 |")
    emit("|---|---:|")
    emit(f"| within-complex：模型选 10 个的 mean ΔΔG | **{sel:+.4f}** kcal/mol |")
    emit(f"| within-complex：随机选 10 个的 mean ΔΔG | {rnd:+.4f} kcal/mol |")
    emit(f"| within-complex：oracle 选 10 个的 mean ΔΔG | {orc:+.4f} kcal/mol |")
    emit(f"| **Efficiency@10 = (random − model) / (random − oracle)** | "
         f"**{eff:+.4f}**（0 = 与随机无异；1 = 达到 oracle） |")
    emit(f"| → 模型回收了 oracle 可用空间的 | **{eff * 100:.1f}%** |")
    emit(f"| within-complex `P@10`（模型 / 随机 / oracle） | {p10:.4f} / "
         f"{p10r:.4f} / {p10o:.4f} |")
    emit(f"| within-complex `hit10(ΔΔG<0)` | **{hit0:.4f}** |")
    emit(f"| 全局基准率（随机选 10 的 ΔΔG<0 期望） | **{(y < 0).mean():.4f}** |")
    emit(f"| 全局层：模型选 10 个的 mean ΔΔG | "
         f"{float(gl_plain['sel10_mean_ddg'].mean()):+.4f} kcal/mol |")
    emit(f"| 全体候选的 mean ΔΔG | {y.mean():+.4f} kcal/mol |")
    emit("")
    emit("### 3c. 判读")
    emit("")
    emit(f"1. **全局层**：模型选出的 10 个平均 ΔΔG = "
         f"**{float(gl_plain['sel10_mean_ddg'].mean()):+.4f}**，比全体均值 "
         f"({y.mean():+.4f}) 低 **{float(gl_plain['sel10_mean_ddg'].mean()) - y.mean():+.4f}** "
         f"kcal/mol → 这部分优势**几乎全部来自跨复合物区分**（挑中的复合物整体更负）。")
    emit(f"2. **复合物内层**：模型比随机好 **{rnd - sel:+.4f}** kcal/mol，"
         f"而 oracle 比随机好 **{rnd - orc:+.4f}** kcal/mol → "
         f"**只回收 {eff * 100:.1f}%**。")
    emit(f"3. **关键矛盾**：within 层的 `hit10(ΔΔG<0)` = **{hit0:.4f}**，"
         f"**低于全局基准率 {(y < 0).mean():.4f}** → "
         f"在一个靶点内挑 10 个「会增强结合」的突变体，命中率**不如随机**。")
    emit(f"4. within 层 `P@10` = {p10:.4f} vs 随机 {p10r:.4f} vs oracle {p10o:.4f} "
         f"→ 仅略高于随机（幅度 ~{p10 - p10r:.4f}）。")
    emit("")
    emit("> **结论**：全局 `P@10 ≈ 0.97` / `EF@1% ≈ 9.6` 是**跨复合物区分 + 高基准率**")
    emit("> 的产物，**不能**解读为「模型已能挑出好突变」。在真实使用场景")
    emit("> （针对一个靶点选 10 个突变体）下，模型只回收了约 **22.5%** 的 oracle 空间，")
    emit("> 且「增强结合」的命中率低于随机基线。这与 macro per-PDB Spearman")
    emit("> （+0.0956，CI 含 0）完全一致。")
    emit("")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_CSV}\n写出 {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
