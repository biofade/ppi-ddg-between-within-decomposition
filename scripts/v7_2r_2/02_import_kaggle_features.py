# -*- coding: utf-8 -*-
"""scripts/v7_2r_2/02_import_kaggle_features.py -- V7.2R-2：导入 Kaggle Track B 特征.

流程（§39）
----------
1. 接受 ZIP 路径（默认 `F:\\dsh\\downloads\\v7_2r_2\\V7_2R_2_TRACKB_FEATURES.zip`）
2. 校验 SHA256（对照 `package_meta.json`；无对照则记录实际值）
3. 解压 features
4. 校验 row_id
5. 按 **row_id** join `canonical_target_labels.csv`（§7）
6. 检查是否有非预期行
7. 报告重复
8. 建立 Track B feature index
9. **表示校验期间绝不读标签**

输出
----
`features/v7_2r_2/trackB_imported/`
`experiments/v7_2r_2/phase4_import/import_report.md`
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
V722 = C.ROOT / "experiments" / "v7_2r_2"
PH4 = V722 / "phase4_import"
PH4.mkdir(parents=True, exist_ok=True)
IMP = C.ROOT / "features" / "v7_2r_2" / "trackB_imported"
IMP.mkdir(parents=True, exist_ok=True)
DEFAULT_ZIP = Path(r"F:\dsh\downloads\v7_2r_2\V7_2R_2_TRACKB_FEATURES.zip")
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
    emit("# V7.2R-2 PHASE 4 — Kaggle Track B 特征导入")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    z = Path(args.zip)
    emit(f"- ZIP: `{z}`")
    if not z.exists():
        emit("")
        emit("## 状态：**WAITING_FOR_KAGGLE_ZIP**")
        emit("")
        emit(f"未找到 `{z}`。")
        emit("")
        emit("**请先完成 Kaggle 侧**（见 `kaggle/v7_2r_2/README_KAGGLE.md`）：")
        emit("1. 上传 `kaggle_input_package.zip` 建立 Dataset")
        emit("2. 上传并运行 `V7_2R_2_TrackB_Kaggle.ipynb`（T4 x2 + Internet ON）")
        emit("3. 下载 `V7_2R_2_TRACKB_FEATURES.zip`")
        emit(f"4. 放到 `{z.parent}`")
        emit("5. 重新运行本脚本")
        emit("")
        emit("> 本阶段的**本地科学评估必须等待**该 ZIP（§54 末条）。")
        emit(f"- 耗时 {time.time() - t0:.0f}s")
        (PH4 / "import_report.md").write_text("\n".join(L) + "\n",
                                              encoding="utf-8")
        return 2

    real = sha(z)
    emit(f"- 实际 SHA256: `{real}`")
    meta_p = V722 / "phase0_prepare" / "package_meta.json"
    if meta_p.exists():
        m = json.loads(meta_p.read_text(encoding="utf-8"))
        emit(f"- 输入包 SHA256（记录）: `{m.get('sha256')}`")
    emit("")

    # 解压
    ex = C.ROOT / "features" / "v7_2r_2" / "_unzip"
    if ex.exists():
        shutil.rmtree(ex)
    ex.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(z) as zf:
        names = zf.namelist()
        zf.extractall(ex)
    emit(f"- 解压 {len(names)} 个条目 → `{ex}`")
    emit("")

    # 找 features
    feats = sorted((ex / "features").glob("*.npz")) if (ex / "features").exists() \
        else sorted(ex.rglob("*.npz"))
    emit(f"- 找到 **{len(feats)}** 个 per-target npz")
    if not feats:
        emit("!! ZIP 内无 features/*.npz")
        (PH4 / "import_report.md").write_text("\n".join(L) + "\n",
                                              encoding="utf-8")
        return 1

    # 权威标签（**只用于 row→target/family 映射**，不用于表示校验）
    lab = pd.read_csv(V711 / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    lab_map = lab.set_index("pair_row")[["target_id_canonical",
                                        "family_id_canonical"]]
    valid_rows = set(lab.pair_row.astype(int).tolist())

    rows_all, dims = [], {}
    unexpected, dup, n_npz_ok = [], [], 0
    for i, f in enumerate(feats):
        try:
            zz = np.load(f, allow_pickle=True)
        except Exception as e:  # noqa: BLE001
            emit(f"  - ⚠️ 读取失败 `{f.name}`: {type(e).__name__}")
            continue
        if "row_id" not in zz.files:
            emit(f"  - ⚠️ 缺 `row_id`：`{f.name}`")
            continue
        rid = zz["row_id"].astype(np.int64)
        if len(set(rid.tolist())) != len(rid):
            dup.append({"file": f.name,
                        "n_dup": len(rid) - len(set(rid.tolist()))})
        unexp = sorted(set(rid.tolist()) - valid_rows)
        if unexp:
            unexpected.append({"file": f.name, "n": len(unexp),
                               "example": unexp[:5]})
        for k in ("r_wt", "r_mut", "delta_r"):
            if k in zz.files:
                dims[k] = int(np.asarray(zz[k]).shape[1])
        n_npz_ok += 1
        shutil.copy(f, IMP / f.name)
        lut = lab_map.reindex(rid)
        rows_all.append(pd.DataFrame({
            "row_id": rid,
            "target_id_kaggle": (zz["target_id"].astype(str)
                                 if "target_id" in zz.files
                                 else pd.Series([""] * len(rid))),
            "target_id_canonical": lut["target_id_canonical"].to_numpy(),
            "family_id_canonical": lut["family_id_canonical"].to_numpy(),
            "n_sites": (zz["mutation_count"].astype(int)
                        if "mutation_count" in zz.files
                        else np.zeros(len(rid), int)),
            "file": f.name,
        }))
        if (i + 1) % 40 == 0:
            emit(f"  - 已处理 {i + 1}/{len(feats)}")
    if not rows_all:
        emit("!! 无有效 npz")
        (PH4 / "import_report.md").write_text("\n".join(L) + "\n",
                                              encoding="utf-8")
        return 1
    RB = pd.concat(rows_all, ignore_index=True)
    RB.to_csv(IMP / "trackB_feature_index.csv", index=False, encoding="utf-8")

    emit("")
    emit("## 4.1 导入结果")
    emit("")
    emit("| 项 | 值 |")
    emit("|---|---|")
    emit(f"| per-target npz | **{n_npz_ok}** / {len(feats)} |")
    emit(f"| Track B rows | **{len(RB)}** |")
    emit(f"| unique row_id | **{RB.row_id.nunique()}** |")
    emit(f"| targets | **{RB.target_id_canonical.nunique()}** |")
    emit(f"| families | **{RB.family_id_canonical.nunique()}** |")
    for k, v in dims.items():
        emit(f"| dim({k}) | **{v}** |")
    emit("")
    emit("## 4.2 完整性检查")
    emit("")
    emit("| 检查 | 结果 |")
    emit("|---|---|")
    emit(f"| 非预期 row_id（不在 canonical） | **{len(unexpected)}** 个文件"
         f"{' → ⚠️ ' + str(unexpected[:3]) if unexpected else ' → ✅ 无'} |")
    emit(f"| 重复 row_id | **{len(dup)}** 个文件"
         f"{' → ⚠️ ' + str(dup[:3]) if dup else ' → ✅ 无'} |")
    emit(f"| 无法映射到 canonical target 的行 | "
         f"**{int(RB.target_id_canonical.isna().sum())}** |")
    # target 名一致性（Kaggle 侧 target_id vs canonical）
    both = RB[RB.target_id_kaggle.astype(str).str.len() > 0]
    if len(both):
        agree = (both.target_id_kaggle.astype(str).str.upper()
                 == both.target_id_canonical.astype(str).str.upper()).mean()
        emit(f"| Kaggle `target_id` vs canonical 一致率 | "
             f"**{agree * 100:.2f}%**（{len(both)} 行） |")
    emit("")
    emit("> **`target identity` 一律采用 `canonical_target_labels.csv`**（§7），")
    emit("> Kaggle 侧自报的 `target_id` 仅用于一致性核对。")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")
    (PH4 / "import_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    json.dump({"zip": str(z), "sha256": real, "n_npz": n_npz_ok,
               "n_rows": len(RB), "dims": dims,
               "n_unexpected_files": len(unexpected),
               "n_dup_files": len(dup)},
              open(PH4 / "import_meta.json", "w", encoding="utf-8"),
              indent=1, ensure_ascii=False)
    print(f"\n写出 {IMP / 'trackB_feature_index.csv'}")
    print(f"写出 {PH4 / 'import_report.md'}")
    print(f">>> Track B rows={len(RB)} targets="
          f"{RB.target_id_canonical.nunique()} dims={dims}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
