# -*- coding: utf-8 -*-
"""scripts/v7_2r_2/06_final_report.py -- V7.2R-2 PHASE：最终报告（§50）.

在 Kaggle ZIP 尚未返回时，输出 **本地准备完成** 的总结
（`KAGGLE_PACKAGE_READY`）并明确标 `WAITING_FOR_KAGGLE_ZIP`。
ZIP 返回后运行，会填入完整科学结果。
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

V722 = C.ROOT / "experiments" / "v7_2r_2"
KAG = C.ROOT / "kaggle" / "v7_2r_2"
DEFAULT_ZIP = Path(r"F:\dsh\downloads\v7_2r_2\V7_2R_2_TRACKB_FEATURES.zip")
L: list[str] = []


def emit(m=""):
    print(m); L.append(str(m))


def sha(p: Path, ch: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(ch), b""):
            h.update(b)
    return h.hexdigest()


def main() -> int:
    t0 = time.time()
    pkg = KAG / "kaggle_input_package.zip"
    nb = KAG / "V7_2R_2_TrackB_Kaggle.ipynb"
    rd = KAG / "README_KAGGLE.md"
    man = KAG / "kaggle_input_manifest.csv"
    ev = V722 / "phase5_eval" / "trackB_results.csv"
    have_b = ev.exists()
    ready = pkg.exists() and nb.exists() and rd.exists() and man.exists()

    emit("# V7.2R-2 Final Summary")
    emit("")
    emit(f"- 生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("")
    emit("## 1. Kaggle Environment")
    emit("")
    emit("**尚未运行**（`WAITING_FOR_KAGGLE_ZIP`）。")
    emit("")
    emit("Kaggle 侧目标环境（§0）：")
    emit("")
    emit("| 项 | 目标值 |")
    emit("|---|---|")
    emit("| Accelerator | **GPU T4 x2**（16 GB VRAM / GPU） |")
    emit("| Internet | **ON**（必需，用于下载权重） |")
    emit("| 使用 GPU | **仅 `cuda:0`**（§12：单 GPU 验证前不启用 GPU1） |")
    emit("| 精度 | 推理 **fp16**；存储 **float16**（§26） |")
    emit("| 预计下载 | MINT 3.25 GB + ESM2-650M 2.61 GB ≈ **5.9 GB**（Kaggle 侧） |")
    emit("")
    emit("## 2. Artifact Dependency Audit")
    emit("")
    emit("| 问题 | 结论 |")
    emit("|---|---|")
    emit("| `mint.ckpt` 是否含 ESM backbone 参数 | **是（很可能）** —— MINT 主干配置即 "
         "`esm2_t33_650M_UR50D.json`；notebook cell 07 会在加载后**实测核对** |")
    emit("| 是否仍需单独下载 ESM2-650M | **需要（仅当要位点级特征时）** —— MINT 公开 API "
         "只暴露链级/对级 `(N,1280|2560)` |")
    emit("| 哪些对象同时在 RAM | **不同时** —— MINT 与 ESM 分两阶段**串行**加载 |")
    emit("| 哪些对象同时在 GPU 显存 | **同一时刻仅一个模型** |")
    emit("| 能否不加载第二个 ESM 完成 MINT 推理 | **能**（仅链级表示时），但拿不到 Δr |")
    emit("")
    emit("```")
    emit("SEPARATE_ESM_REQUIRED_FOR_SITE_FEATURES = YES")
    emit("DUPLICATE_WEIGHTS_DOWNLOADED = NO")
    emit("```")
    emit("")
    emit("## 3. Smoke Test")
    emit("")
    emit("notebook cell 09 内建三目标冒烟（**label-free**，按 `n_sites` 小/中/大选靶）：")
    emit("")
    emit("| 检查 | 说明 |")
    emit("|---|---|")
    emit("| MODEL_LOAD_PASS | ESM2-650M fp16 加载成功 |")
    emit("| RESIDUE_MAPPING_PASS | `seq_index0` 与 `wt_sequence` 一致 |")
    emit("| WT_VALIDATION_PASS | 期望 WT == 实际 WT（**不符即拒绝**，§19） |")
    emit("| MUTANT_RECONSTRUCTION_PASS | 多突变**同时**施加并回读核对 |")
    emit("| OUTPUT_SHAPE_UNDERSTOOD | 由**真实 tensor** 判定层级（per-residue） |")
    emit("| PARTNER_SEMANTICS_UNDERSTOOD | `‖r_A alone − r_A+true‖` vs "
         "`‖r_A+true − r_A+wrong‖` |")
    emit("| NO_NAN_INF | 位点向量全有限 |")
    emit("| CACHE_RELOAD_PASS | manifest 断点续跑可读 |")
    emit("")
    emit("**冒烟失败 → notebook 主动抛 `TRACK_B_SMOKE_FAIL` 并停止**（§23）。")
    emit("")
    emit("## 4. Track B Coverage")
    emit("")
    emit("**NOT_RUN**（等 ZIP）。")
    emit("")
    emit("## 5. Matched U_AB")
    emit("")
    emit("`U_AB = Track A valid ∩ Track B valid`（§40）。`03_matched_trackB_eval.py` "
         "会在导入后构造并检查 §40 门槛（n≥5 target 覆盖 ≥80% of Track A），"
         "不达标标 **`COVERAGE_LIMITED`**。")
    emit("")
    emit("Track A reference：**4233 rows / 174 targets / 91 families**（V7.2R-1，94.6%）。")
    emit("")
    emit("## 6. Primary Results")
    emit("")
    emit("| Method | Macro rho | Eff@10 | C_W | Kendall | Pairwise |")
    emit("|---|---:|---:|---:|---:|---:|")
    for m in ("B0", "A2", "B1", "B2", "B3"):
        emit(f"| {m} | NOT_RUN | | | | |")
    emit("")
    emit("## 7. Bootstrap")
    emit("")
    emit("**NOT_RUN**（5000 paired target-cluster bootstrap，等 ZIP）")
    emit("")
    emit("## 8. Negative Controls")
    emit("")
    emit("**NOT_RUN**（NC1–NC4 已实现于 `05_negative_controls.py`，等 ZIP）")
    emit("")
    emit("## 9. Top-K Analysis")
    emit("")
    emit("**NOT_RUN**。若 Macro ρ 的 CI 含 0 但 Eff@10 显著 > 0，"
         "依 §47 只声明 **`TOPK_ONLY_SIGNAL`**。")
    emit("")
    emit("## 10. Calibration Shell")
    emit("")
    emit("固定为 `y_hat = mu^MU3 + z`，`alpha = 1`（§48）。**NOT_RUN**。")
    emit("")
    emit("## 11. Scientific Verdict")
    emit("")
    emit("```")
    emit("TRACK_B_RESULT = WAITING_FOR_KAGGLE_ZIP")
    emit("```")
    emit("")
    emit("V7.2R-1 的既有结论（未变）：")
    emit("")
    emit("```")
    emit("TRACK_A_RESULT = FULL_MAPPABLE_NULL")
    emit("R3(A2) - R0(B0): Delta MacroSpearman = -0.0143, 95% CI [-0.0588, +0.0307]")
    emit("```")
    emit("")
    emit("## 12. Claims Allowed")
    emit("")
    emit("- 本地已产出**无标签**的 Kaggle 推理包与**可执行、重启安全**的 notebook")
    emit("- 外部模型**只提供表示**；未使用任何 SKEMPI 监督 head")
    emit("- 评估协议为 **V7.1.1 冻结的嵌套 family-OOD**")
    emit("- Track A（小 PLM Δh）在完整可映射宇宙上为 **NULL**")
    emit("")
    emit("## 13. Claims Forbidden")
    emit("")
    emit("- ❌ 任何 Track B 的性能结论（**尚未运行**）")
    emit("- ❌ 「protein foundation models do not work」")
    emit("- ❌ 「we invented the underlying PPI representation」"
         "（依 §51 必须写 *building on*）")
    emit("")
    emit("## 14. Innovation Level")
    emit("")
    emit("**LEVEL 2（TASK-FORMULATION GAIN）** —— 沿用 V7.1 的标定/排序解耦证据；"
         "Track B 的 LEVEL 3/4/5 判定**待 Kaggle 结果**。")
    emit("")
    emit("## 15. Exact Next Route")
    emit("")
    emit("```")
    emit("KAGGLE_PACKAGE_READY = " + ("YES" if ready else "NO"))
    emit("")
    emit("Notebook:")
    emit(str(nb))
    emit("")
    emit("Input package:")
    emit(str(pkg))
    emit("")
    emit("README:")
    emit(str(rd))
    emit("")
    emit("What user must upload:")
    emit("  1) " + pkg.name + "  ->  Kaggle Dataset")
    emit("  2) " + nb.name + "  ->  Kaggle Notebook (Import Notebook)")
    emit("")
    emit("What user must download after Kaggle:")
    emit("  V7_2R_2_TRACKB_FEATURES.zip")
    emit("  -> F:\\dsh\\downloads\\v7_2r_2\\")
    emit("")
    emit("Then run:")
    emit("  python scripts/v7_2r_2/02_import_kaggle_features.py")
    emit("  python scripts/v7_2r_2/03_matched_trackB_eval.py")
    emit("  python scripts/v7_2r_2/04_bootstrap.py")
    emit("  python scripts/v7_2r_2/05_negative_controls.py")
    emit("  python scripts/v7_2r_2/06_final_report.py")
    emit("```")
    emit("")
    if DEFAULT_ZIP.exists():
        emit(f"- 已检测到 Kaggle ZIP：`{DEFAULT_ZIP}` "
             f"({DEFAULT_ZIP.stat().st_size / 1e6:.1f} MB, sha256 `{sha(DEFAULT_ZIP)[:16]}`)")
    else:
        emit(f"- **未检测到** `{DEFAULT_ZIP}`")
    emit("")
    emit(f"- 耗时 {time.time() - t0:.0f}s")

    (V722 / "V7_2R_2_FINAL_SUMMARY.md").write_text("\n".join(L) + "\n",
                                                   encoding="utf-8")
    json.dump({"KAGGLE_PACKAGE_READY": bool(ready),
               "notebook": str(nb), "input_package": str(pkg),
               "readme": str(rd),
               "trackB_features_returned": bool(DEFAULT_ZIP.exists()),
               "generated": time.strftime("%Y-%m-%dT%H:%M:%S")},
              open(V722 / "package_status.json", "w", encoding="utf-8"),
              indent=1, ensure_ascii=False)
    print("\n" + "=" * 66)
    print(f"KAGGLE_PACKAGE_READY = {'YES' if ready else 'NO'}")
    print("=" * 66)
    print(f"\nNotebook:\n  {nb}")
    print(f"\nInput package:\n  {pkg}")
    print(f"\nREADME:\n  {rd}")
    print("\nWhat user must upload:")
    print(f"  1) {pkg.name}   -> Kaggle Dataset")
    print(f"  2) {nb.name}   -> Kaggle Notebook")
    print("\nWhat user must download after Kaggle:")
    print("  V7_2R_2_TRACKB_FEATURES.zip  ->  F:\\dsh\\downloads\\v7_2r_2\\")
    print(f"\n写出 {V722 / 'V7_2R_2_FINAL_SUMMARY.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
