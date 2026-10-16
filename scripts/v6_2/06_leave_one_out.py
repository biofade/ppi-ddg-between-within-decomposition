# -*- coding: utf-8 -*-
"""scripts/v6_2/06_leave_one_out.py -- V6.2 PHASE 6: Leave-one-model-out + BH-FDR.

分类（V6.2 §6.1）
-----------------
* **Explicit physics**: `{FoldX, Flex-ddG}`
* **Learned structure-conditioned**: `{ESM-IF, ProteinMPNN}`

五个组合（§6.2）
----------------
1. `{Flex-ddG, FoldX}` vs `{ProteinMPNN}`      （去 ESM-IF）
2. `{Flex-ddG, FoldX}` vs `{ESM-IF}`            （去 ProteinMPNN）
3. `{FoldX}` vs `{ESM-IF, ProteinMPNN}`         （去 Flex-ddG）
4. `{Flex-ddG}` vs `{ESM-IF, ProteinMPNN}`      （去 FoldX）
5. Full `{Flex-ddG, FoldX}` vs `{ESM-IF, ProteinMPNN}`

方法
----
* 指标：**per-complex Spearman**，在 **common eligible targets** 上配对
* 检验：**配对 Wilcoxon 符号秩**（同一批 complex → 必须配对）
* **BH-FDR** 对 5 个 p 值校正

判定（§6.4）
------------
* **ROBUST**：4 个 leave-one-out 方向一致 + majority q<0.05
* **PARTIALLY_ROBUST**：方向全部相同但部分 CI 跨 0
* **MODEL_DRIVEN**：删任一模型后 effect size 掉 >50% 或反号

输出
----
`phase6_leave_one_out/leave_one_out.csv` + `leave_one_out_report.md`
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

PH6 = C.ROOT / "experiments" / "v6_2" / "phase6_leave_one_out"
PH6.mkdir(parents=True, exist_ok=True)
STD = C.ROOT / "experiments" / "v6_sota_audit" / "standardized_preds.csv"
MIN_N = 5
N_BOOT = 5000
SEED = 42
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def per_cx(df, model):
    out = {}
    for cx, gg in df[df["model"] == model].groupby("complex_id"):
        if len(gg) < MIN_N:
            continue
        v = C.spearman(gg["true_ddG"].to_numpy(float),
                       gg["pred_ddG"].to_numpy(float))
        if v == v:
            out[str(cx)] = float(v)
    return out


def boot_delta(a, b, n_boot=N_BOOT, seed=SEED):
    """配对 bootstrap（抽 complex）：Δ = mean(a) − mean(b)。"""
    n = len(a)
    if n < 3:
        return {"delta": np.nan, "ci_low": np.nan, "ci_high": np.nan,
                "sd": np.nan, "n": n}
    d = np.asarray(a, float) - np.asarray(b, float)
    rs = np.random.RandomState(seed)
    vals = [d[rs.randint(0, n, n)].mean() for _ in range(n_boot)]
    v = np.asarray(vals, float)
    lo, hi = np.percentile(v, [2.5, 97.5])
    return {"delta": float(d.mean()), "ci_low": float(lo), "ci_high": float(hi),
            "sd": float(v.std(ddof=1)), "n": n}


def main() -> int:
    t0 = time.time()
    emit("# V6.2 PHASE 6 — Leave-one-model-out + BH-FDR")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- 指标: per-complex Spearman（每 complex ≥{MIN_N} 条突变）")
    emit(f"- 检验: **配对 Wilcoxon**（同批 complex）+ **BH-FDR**（5 个 p 值）")
    emit(f"- bootstrap: {N_BOOT} 次，单位 = complex")
    emit("")
    if not STD.exists():
        emit("!! `standardized_preds.csv` 不存在")
        (PH6 / "leave_one_out_report.md").write_text("\n".join(L) + "\n",
                                                    encoding="utf-8")
        return 1
    d = pd.read_csv(STD, low_memory=False)
    ca = d[d["dataset_name"] == "CATH"]
    per = {m: per_cx(ca, m) for m in ("FoldX", "Flex-ddG", "ESM-IF", "proteinMPNN")}
    emit("| 模型 | valid targets | per-complex Spearman 均值 |")
    emit("|---|---:|---:|")
    for m, v in per.items():
        arr = np.array(list(v.values()), float)
        emit(f"| {m} | {len(arr)} | {arr.mean():+.4f} |")
    emit("")

    combos = [
        ("去 ESM-IF", ["Flex-ddG", "FoldX"], ["proteinMPNN"]),
        ("去 ProteinMPNN", ["Flex-ddG", "FoldX"], ["ESM-IF"]),
        ("去 Flex-ddG", ["FoldX"], ["ESM-IF", "proteinMPNN"]),
        ("去 FoldX", ["Flex-ddG"], ["ESM-IF", "proteinMPNN"]),
        ("Full", ["Flex-ddG", "FoldX"], ["ESM-IF", "proteinMPNN"]),
    ]
    from scipy.stats import wilcoxon
    rows = []
    emit("## 五个组合")
    emit("")
    emit("| 组合 | 物理侧 | 序列侧 | n_common | 物理 mean | 序列 mean | Δ | "
         "95% CI | Wilcoxon p |")
    emit("|---|---|---|---:|---:|---:|---:|---|---:|")
    for name, phys, seq in combos:
        common = None
        for m in phys + seq:
            s = set(per[m])
            common = s if common is None else (common & s)
        common = sorted(common or [])
        if len(common) < 5:
            emit(f"| {name} | {'+'.join(phys)} | {'+'.join(seq)} | {len(common)} | "
                 f"— | — | — | — | 配对数不足 |")
            continue
        a = np.array([np.mean([per[m][c] for m in phys]) for c in common], float)
        b = np.array([np.mean([per[m][c] for m in seq]) for c in common], float)
        try:
            w, p = wilcoxon(a, b)
        except Exception:  # noqa: BLE001
            w, p = float("nan"), float("nan")
        bt = boot_delta(a, b)
        emit(f"| {name} | {'+'.join(phys)} | {'+'.join(seq)} | {len(common)} | "
             f"{a.mean():+.4f} | {b.mean():+.4f} | **{bt['delta']:+.4f}** | "
             f"[{bt['ci_low']:+.4f}, {bt['ci_high']:+.4f}] | {p:.6f} |")
        rows.append({"combo": name, "phys_side": "+".join(phys),
                     "seq_side": "+".join(seq), "n_common": len(common),
                     "phys_mean": float(a.mean()), "seq_mean": float(b.mean()),
                     "delta": bt["delta"], "ci_low": bt["ci_low"],
                     "ci_high": bt["ci_high"], "boot_sd": bt["sd"],
                     "wilcoxon_W": float(w), "wilcoxon_p": float(p)})
    res = pd.DataFrame(rows)
    emit("")

    # BH-FDR
    if len(res):
        p = pd.to_numeric(res["wilcoxon_p"], errors="coerce").to_numpy(float)
        q = C.bh_fdr(list(p))
        res["bh_q"] = q
        res["ci_excludes_zero"] = (res["ci_low"] > 0) | (res["ci_high"] < 0)
        emit("## BH-FDR 校正（5 个 p 值）")
        emit("")
        emit("| 组合 | Δ | 95% CI | p | **BH q** | q<0.05? | CI 排除 0? |")
        emit("|---|---:|---|---:|---:|---|---|")
        for _, r in res.iterrows():
            emit(f"| {r['combo']} | {r['delta']:+.4f} | "
                 f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}] | "
                 f"{r['wilcoxon_p']:.6f} | **{r['bh_q']:.6f}** | "
                 f"{'✅' if r['bh_q'] < 0.05 else '❌'} | "
                 f"{'✅' if r['ci_excludes_zero'] else '❌'} |")
        emit("")
        # 判定
        loo = res[res["combo"] != "Full"]
        dirs = np.sign(loo["delta"].to_numpy(float))
        same_dir = len(set(dirs)) == 1
        maj_q = int((loo["bh_q"] < 0.05).sum()) >= len(loo) / 2
        full = res[res["combo"] == "Full"]
        base = float(full["delta"].iloc[0]) if len(full) else float("nan")
        drops = []
        for _, r in loo.iterrows():
            if abs(base) > 1e-9:
                drops.append(abs(r["delta"]) / abs(base))
        model_driven = any(d < 0.5 for d in drops) or (len(set(dirs)) > 1)
        if same_dir and maj_q and not model_driven:
            verdict = "ROBUST"
        elif same_dir and not model_driven:
            verdict = "PARTIALLY_ROBUST"
        elif model_driven:
            verdict = "MODEL_DRIVEN"
        else:
            verdict = "PARTIALLY_ROBUST"
        emit(f"### 判定：**{verdict}**")
        emit("")
        emit("| 判据 | 值 |")
        emit("|---|---|")
        emit(f"| 4 个 leave-one-out 方向一致 | **{same_dir}**（符号: "
             f"{[int(x) for x in dirs]}）|")
        emit(f"| majority q<0.05 | **{maj_q}**"
             f"（{int((loo['bh_q'] < 0.05).sum())}/{len(loo)}）|")
        emit(f"| effect size 相对 Full 的比值 | "
             f"{[round(d, 3) for d in drops]}（<0.5 即 MODEL_DRIVEN）|")
        emit("")
        res.to_csv(PH6 / "leave_one_out.csv", index=False, encoding="utf-8")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    (PH6 / "leave_one_out_report.md").write_text("\n".join(L) + "\n",
                                                encoding="utf-8")
    print(f"\n写出 {PH6 / 'leave_one_out.csv'}")
    print(f"写出 {PH6 / 'leave_one_out_report.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
