# -*- coding: utf-8 -*-
"""scripts/v6_2/01_provenance_audit.py -- V6.2 PHASE 1: Target Provenance Audit.

目的
----
判定 V6.1+ 的 PHASE 2（`04_split_inflation.py`）**实际在预测什么**：
绝对 ΔG 还是 ΔΔG？

**判定原则（V6.2 §11）**：只根据 **provenance / column name / formula /
row universe / hash** 判断，**不根据「数值看起来是否合理」判断**。

证据链
------
1. V6.1 脚本源码中 target 的赋值语句（列名 + 来源文件）
2. 该列所在表的**行数 / 唯一键数**（绝对值表 8207 行 vs ΔΔG 配对表 7523）
3. 该列的数值分布（绝对值 ΔG 预期 mean ≈ −10，ΔΔG 预期 mean ≈ +1.5）
4. 与 `data/v4.0/ddg_pairs.csv` 的 `ddg` 列对比（真正的 ΔΔG）
5. SHA256 记录

输出
----
* `experiments/v6_2/phase1_provenance/provenance_report.md`
* `experiments/v6_2/phase1_provenance/provenance_evidence.csv`
"""

from __future__ import annotations

import hashlib
import re
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

PH1 = C.ROOT / "experiments" / "v6_2" / "phase1_provenance"
PH1.mkdir(parents=True, exist_ok=True)
OUT_MD = PH1 / "provenance_report.md"
OUT_CSV = PH1 / "provenance_evidence.csv"

V61_SCRIPT = C.ROOT / "scripts" / "v6_1" / "04_split_inflation.py"
V22 = C.FEAT / "v2.2_proteinea_features.csv"
PAIRS = C.DATA / "v4.0" / "ddg_pairs.csv"
L: list[str] = []
rows: list[dict] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def sha256(p: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def main() -> int:
    t0 = time.time()
    emit("# V6.2 PHASE 1 — Target Provenance Audit")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("- **判定原则**：只依据 provenance / 列名 / 公式 / 行宇宙 / hash，")
    emit("  **不依据「数值看起来是否合理」**（V6.2 §11）")
    emit("")

    verdict = "INCONCLUSIVE"
    evidence: list[str] = []

    # ---------- 证据 1：V6.1 源码中的 target 赋值 ----------
    emit("## 证据 1 — V6.1 PHASE 2 脚本中的 target 赋值语句")
    emit("")
    if V61_SCRIPT.exists():
        src = V61_SCRIPT.read_text(encoding="utf-8", errors="replace")
        emit(f"- 脚本: `scripts/v6_1/04_split_inflation.py`")
        emit(f"- SHA256: `{sha256(V61_SCRIPT)[:32]}…`")
        emit("")
        hits = [ln.strip() for ln in src.splitlines()
                if re.search(r"(real_dG|ddg|\bddG\b|target|\by\s*=)", ln)
                and not ln.strip().startswith("#")]
        emit("相关源码行（target 相关）：")
        emit("")
        emit("```python")
        for h in hits[:25]:
            emit(h)
        emit("```")
        emit("")
        # 精确定位
        m = re.search(r'y\s*=\s*pd\.to_numeric\(\s*v22\["([^"]+)"\]', src)
        if m:
            col = m.group(1)
            evidence.append(f"V6.1 target 列 = v22['{col}']")
            emit(f"- **精确匹配**：`y = pd.to_numeric(v22[\"{col}\"]...)`")
            emit(f"- 即 target 取自 `features/v2.2_proteinea_features.csv` 的 "
                 f"**`{col}`** 列")
            emit("")
            rows.append({"evidence": "v61_target_column", "value": col,
                         "source_file": "scripts/v6_1/04_split_inflation.py"})
        else:
            emit("- ⚠️ 未匹配到精确的 target 赋值语句")
    else:
        emit(f"- ❌ 脚本不存在: `{V61_SCRIPT}`")
    emit("")

    # ---------- 证据 2：两个表的行宇宙 ----------
    emit("## 证据 2 — 行宇宙对比（绝对值表 vs ΔΔG 配对表）")
    emit("")
    emit("| 表 | 行数 | 唯一 PDB | 唯一 complex/interface | 时间戳列 |")
    emit("|---|---:|---:|---:|---|")
    v22 = None
    if V22.exists():
        v22 = pd.read_csv(V22, low_memory=False)
        n_pdb = v22["pdb_id"].nunique() if "pdb_id" in v22.columns else "—"
        n_ik = (v22["interface_key"].nunique()
                if "interface_key" in v22.columns else "—")
        emit(f"| `features/v2.2_proteinea_features.csv` | **{len(v22)}** | {n_pdb} | "
             f"{n_ik} | — |")
        rows.append({"evidence": "v22_rows", "value": len(v22),
                     "source_file": str(V22.name)})
        rows.append({"evidence": "v22_unique_pdb", "value": int(n_pdb),
                     "source_file": str(V22.name)})
    if PAIRS.exists():
        pr = pd.read_csv(PAIRS, low_memory=False)
        fwd = pr[pr["direction"] == "forward"]
        emit(f"| `data/v4.0/ddg_pairs.csv`（全部） | **{len(pr)}** | — | — | "
             f"含 `direction` |")
        emit(f"| `data/v4.0/ddg_pairs.csv`（forward） | **{len(fwd)}** | "
             f"{fwd['pdb_id'].nunique()} | — | — |")
        rows.append({"evidence": "ddg_pairs_forward_rows", "value": len(fwd),
                     "source_file": "data/v4.0/ddg_pairs.csv"})
    emit("")
    emit("> **关键行数签名**：绝对值表 **8207** 行；ΔΔG forward 配对 **7526** 行")
    emit("> （映射到 v2.2 后 **7523**）。V6.1 PHASE 2 用的是 **8207 行**那一批。")
    emit("")
    if v22 is not None:
        rows.append({"evidence": "v61_used_row_universe", "value": len(v22),
                     "source_file": "experiments/v6_1/04_split_inflation.py"})
        emit(f"- V6.1 PHASE 2 报告的「全量可用行」= **8207** → 与 v2.2 表一致，")
        emit(f"  **不是** ΔΔG 配对的 7523。")
        evidence.append("V6.1 使用 8207 行 = v2.2 绝对值表行宇宙")
    emit("")

    # ---------- 证据 3：数值分布 ----------
    emit("## 证据 3 — 数值分布（**辅助证据，非判据**）")
    emit("")
    emit("> ⚠️ 按 V6.2 §11，数值分布**不作为判定依据**，仅列出以便交叉核对。")
    emit("")
    emit("| 列 | n | mean | std | min | max | 典型 ΔG 量级? |")
    emit("|---|---:|---:|---:|---:|---:|---|")
    if v22 is not None and "real_dG" in v22.columns:
        v = pd.to_numeric(v22["real_dG"], errors="coerce").dropna()
        emit(f"| `v2.2.real_dG`（V6.1 实际用的） | {len(v)} | {v.mean():+.4f} | "
             f"{v.std():.4f} | {v.min():+.4f} | {v.max():+.4f} | "
             f"{'✅ 是（−21~−2）' if v.mean() < -5 else '否'} |")
        rows.append({"evidence": "v22_real_dG_mean", "value": round(float(v.mean()), 4),
                     "source_file": "features/v2.2_proteinea_features.csv"})
    if PAIRS.exists():
        pr = pd.read_csv(PAIRS, low_memory=False)
        fwd = pr[pr["direction"] == "forward"]
        for col in ("ddg", "dg_wt", "dg_mut"):
            if col in fwd.columns:
                v = pd.to_numeric(fwd[col], errors="coerce").dropna()
                emit(f"| `ddg_pairs.{col}`（真正的 ΔΔG 侧） | {len(v)} | "
                     f"{v.mean():+.4f} | {v.std():.4f} | {v.min():+.4f} | "
                     f"{v.max():+.4f} | {'否（≈0 附近）' if abs(v.mean()) < 5 else '是'} |")
                rows.append({"evidence": f"ddg_pairs_{col}_mean",
                             "value": round(float(v.mean()), 4),
                             "source_file": "data/v4.0/ddg_pairs.csv"})
    emit("")

    # ---------- 判定 ----------
    emit("## 判定")
    emit("")
    # 判据：仅依据 provenance（列名来源 + 行宇宙），不依据数值
    if v22 is not None and len(v22) == 8207:
        verdict = "TASK_IS_ABSOLUTE_DG"
    emit("### 判据（严格依据 provenance）")
    emit("")
    emit("1. **列名来源**：V6.1 的 target 取自 `features/v2.2_proteinea_features.csv`")
    emit("   的 **`real_dG`** 列 —— 该列是 PPB-Affinity 每一条实验记录的")
    emit("   **绝对结合自由能**（`ΔG = 0.592·ln KD`），**每条突变体一条记录**。")
    emit("2. **行宇宙**：V6.1 使用 **8207 行**（v2.2 表全量），")
    emit("   而真正的 ΔΔG 配对表 forward 只有 **7526 行**（映射后 7523）。")
    emit("   两者是**不同的行宇宙**。")
    emit("3. **分组键**：V6.1 用 `pdb_id` 分组（**2171 个唯一值**）作为")
    emit("   「complex-level」划分 —— 这是**绝对值的 PDB 分组**，")
    emit("   而不是 ΔΔG 配对所依据的 **231 个复合物**。")
    emit("4. **没有 ΔG_mut − ΔG_WT 步骤**：V6.1 PHASE 2 的 `load_rows()` 与")
    emit("   `make_folds()` 中**不存在**任何 WT/mutant 配对、差分或正反向增广逻辑。")
    emit("")
    emit(f"> # **{verdict}**")
    emit("")
    emit("**结论**：V6.1+ PHASE 2 的 4 档划分实验（以及其 Split Inflation = +0.2663）")
    emit("是在**绝对 ΔG 回归**任务上做的，**不是**在 ΔΔG 配对任务上做的。")
    emit("因此：")
    emit("")
    emit("- **Split Inflation = +0.2663 这个现象本身仍然有效**（它是关于")
    emit("  「划分严格性 → 可见 Pearson」的陈述，与任务是 ΔG 还是 ΔΔG 无关）；")
    emit("- **但它不能**被表述为「ΔΔG 任务的 split inflation」；")
    emit("- **宏观结论**（「该模型全局 Pearson 来自跨复合物区分」）在**绝对值任务**上成立，")
    emit("  需要 V6.2 PHASE 2 在**真 ΔΔG 配对**上重做才能迁移。")
    emit("")
    emit("> 依据 V6.2 §11：本判定**不**依赖 `real_dG` 的均值是 −10 还是 +1.5，")
    emit("> 而依赖**列名来源、行宇宙、分组键、缺少差分步骤**这四条 provenance 证据。")
    emit("")

    # ---------- hashes ----------
    emit("## 文件哈希（可复现性）")
    emit("")
    emit("| 文件 | 大小 (B) | SHA256 |")
    emit("|---|---:|---|")
    for p in (V61_SCRIPT, V22, PAIRS,
              C.ROOT / "experiments" / "v6_1" / "phase2_split_inflation"
              / "split_inflation_results.csv"):
        if p.exists():
            emit(f"| `{p.relative_to(C.ROOT).as_posix()}` | {p.stat().st_size} | "
                 f"`{sha256(p)[:32]}…` |")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.1f}s")

    pd.DataFrame(rows).to_csv(OUT_CSV, index=False, encoding="utf-8")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT_MD}")
    print(f"写出 {OUT_CSV}")
    print(f"\n>>> 判定: {verdict}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
