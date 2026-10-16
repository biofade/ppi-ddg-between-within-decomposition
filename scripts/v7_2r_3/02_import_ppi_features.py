# -*- coding: utf-8 -*-
"""scripts/v7_2r_3/02_import_ppi_features.py -- V7.2R-3：导入 PPI-context 特征.

ZIP：`F:\\dsh\\downloads\\v7_2r_3\\V7_2R_3_PPI_CONTEXT_FEATURES.zip`
输出：`features/v7_2r_3/ppi_imported/`、`experiments/v7_2r_3/phase4_import/`
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import time
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

V711 = C.ROOT / "experiments" / "v7_1_1"
V723 = C.ROOT / "experiments" / "v7_2r_3"
PH4 = V723 / "phase4_import"
PH4.mkdir(parents=True, exist_ok=True)
IMP = C.ROOT / "features" / "v7_2r_3" / "ppi_imported"
IMP.mkdir(parents=True, exist_ok=True)
DEFAULT_ZIP = Path(r"F:\dsh\downloads\v7_2r_3\V7_2R_3_PPI_CONTEXT_FEATURES.zip")
L: list[str] = []


def emit(m=""):
    print(m); L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def sha(p: Path, ch: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(ch), b""):
            h.update(b)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--zip", default=str(DEFAULT_ZIP))
    args = ap.parse_args()
    t0 = time.time()
    emit("# V7.2R-3 PHASE 4 — PPI-context 特征导入")
    emit("")
    z = Path(args.zip)
    emit(f"- ZIP: `{z}`")
    if not z.exists():
        emit("")
        emit("## 状态：**WAITING_FOR_V7_2R_3_ZIP**")
        emit("")
        emit(f"未找到 `{z}`。请先完成 Kaggle 侧（见 "
             f"`kaggle/v7_2r_3/README_KAGGLE.md`）：")
        emit("1. 上传 `kaggle/v7_2r_3/kaggle_input_package.zip` 建 Dataset")
        emit("2. 上传并运行 `V7_2R_3_PPI_Context_Kaggle.ipynb`（T4 x2 + Internet ON）")
        emit("   —— **会先跑伙伴 sanity 测试**；未通过则主动停止")
        emit("3. 下载 `V7_2R_3_PPI_CONTEXT_FEATURES.zip`")
        emit(f"4. 放到 `{z.parent}`")
        emit("5. 重跑本脚本")
        (PH4 / "import_report.md").write_text("\n".join(L) + "\n",
                                              encoding="utf-8")
        return 2
    real = sha(z)
    emit(f"- 实际 SHA256: `{real}`")
    ex = C.ROOT / "features" / "v7_2r_3" / "_unzip"
    if ex.exists():
        shutil.rmtree(ex)
    ex.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(z) as zf:
        zf.extractall(ex)
    feats = sorted((ex / "features").glob("*.npz")) if (ex / "features").exists() \
        else sorted(ex.rglob("*.npz"))
    emit(f"- features npz: **{len(feats)}**")
    if not feats:
        emit("!! ZIP 内无 features/*.npz")
        (PH4 / "import_report.md").write_text("\n".join(L) + "\n",
                                              encoding="utf-8")
        return 1

    lab = pd.read_csv(V711 / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    lab_map = lab.set_index("pair_row")[["target_id_canonical",
                                        "family_id_canonical"]]
    valid = set(lab.pair_row.astype(int))

    rows, dims, unexp, nch = [], {}, [], []
    for i, f in enumerate(feats):
        try:
            zz = np.load(f, allow_pickle=True)
        except Exception:  # noqa: BLE001
            continue
        rid = zz["row_id"].astype(np.int64)
        u = sorted(set(rid.tolist()) - valid)
        if u:
            unexp.append({"file": f.name, "n": len(u), "ex": u[:3]})
        if "n_chains_embedded" in zz.files:
            nch.extend(zz["n_chains_embedded"].astype(int).tolist())
        for k in ("r_wt", "r_mut", "delta_r"):
            if k in zz.files:
                dims[k] = int(np.asarray(zz[k]).shape[1])
        shutil.copy(f, IMP / f.name)
        lut = lab_map.reindex(rid)
        rows.append(pd.DataFrame({
            "row_id": rid,
            "target_id_canonical": lut["target_id_canonical"].to_numpy(),
            "family_id_canonical": lut["family_id_canonical"].to_numpy(),
            "n_sites": (zz["mutation_count"].astype(int)
                        if "mutation_count" in zz.files else 0),
            "rep_level": (str(zz["representation_level"][0])
                          if "representation_level" in zz.files else ""),
            "file": f.name}))
        if (i + 1) % 40 == 0:
            emit(f"  - {i + 1}/{len(feats)}")
    if not rows:
        emit("!! 无有效 npz")
        (PH4 / "import_report.md").write_text("\n".join(L) + "\n",
                                              encoding="utf-8")
        return 1
    RB = pd.concat(rows, ignore_index=True)
    RB.to_csv(IMP / "ppi_feature_index.csv", index=False, encoding="utf-8")

    # 读 Kaggle 侧报告
    pg = ex / "provenance.json"
    ps = ex / "partner_sanity.csv"
    cov = ex / "coverage_report.md"
    emit("")
    if pg.exists():
        emit("## 4.0 Kaggle provenance")
        emit("")
        emit("```json")
        emit(pg.read_text(encoding="utf-8")[:1200])
        emit("```")
        emit("")
    emit("## 4.1 导入结果")
    emit("")
    emit("| 项 | 值 |")
    emit("|---|---|")
    emit(f"| per-target npz | **{len(feats)}** |")
    emit(f"| rows | **{len(RB)}** |")
    emit(f"| unique row_id | **{RB.row_id.nunique()}** |")
    emit(f"| targets | **{RB.target_id_canonical.nunique()}** |")
    emit(f"| families | **{RB.family_id_canonical.nunique()}** |")
    for k, v in dims.items():
        emit(f"| dim({k}) | **{v}** |")
    if nch:
        emit(f"| 嵌入链数中位 | **{int(np.median(nch))}** |")
    emit(f"| 非预期 row_id | **{len(unexp)}** 文件"
         f"{' ⚠️ ' + str(unexp[:2]) if unexp else ' → ✅ 无'} |")
    emit("")
    if ps.exists():
        SA = pd.read_csv(ps)
        emit("## 4.2 伙伴 sanity（Kaggle 侧）")
        emit("")
        emit(f"- 测试 target 数：**{len(SA)}**")
        if len(SA) and SA.D_true.notna().any():
            emit(f"- `D_true` 中位 = **{SA.D_true.median():.4f}**")
            emit(f"- `D_wrong` 中位 = **{SA.D_wrong.median():.4f}**")
            ok = SA.D_true.median() > 1e-3
            emit(f"- **PPI_CONDITIONING = "
                 f"{'CONFIRMED ✅' if ok else 'NOT_DEMONSTRATED ❌'}**")
        emit("")
        SA.to_csv(PH4 / "partner_sanity.csv", index=False, encoding="utf-8")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    (PH4 / "import_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    json.dump({"zip": str(z), "sha256": real, "n_npz": len(feats),
               "n_rows": len(RB), "dims": dims,
               "median_chains": int(np.median(nch)) if nch else None},
              open(PH4 / "import_meta.json", "w", encoding="utf-8"),
              indent=1, ensure_ascii=False)
    print(f"\n写出 {IMP / 'ppi_feature_index.csv'} 与 import_report.md")
    print(f">>> rows={len(RB)} targets={RB.target_id_canonical.nunique()} dims={dims}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
