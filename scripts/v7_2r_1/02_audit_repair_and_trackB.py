# -*- coding: utf-8 -*-
"""scripts/v7_2r_1/02_audit_repair_and_trackB.py -- V7.2R-1 PHASE 4-5.

PHASE 4：**修复外部审计错误**（Prompt-DDG 许可；ReCLIP 三对象分离；MINT 新审计）
PHASE 5：**Track B 可行性 + 计算可行性**评估（**本轮不下载 >2 GB 产物**）
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R1 = C.ROOT / "experiments" / "v7_2r_1"
PH4, PH5 = R1 / "phase4_audit_repair", R1 / "phase5_trackB"
for d in (PH4, PH5):
    d.mkdir(parents=True, exist_ok=True)
L4: list[str] = []
L5: list[str] = []


def e4(m=""):
    print(m); L4.append(str(m))


def e5(m=""):
    print(m); L5.append(str(m))


# 实测的 artifact 大小（HTTP HEAD，未经下载）
ARTIFACTS = [
    ("MINT checkpoint（mint.ckpt）",
     "https://huggingface.co/varunullanat2012/mint/resolve/main/mint.ckpt",
     "MIT", 3253.8, "MINT 主干权重（基于 ESM2-650M）"),
    ("ESM2-650M（model.safetensors）",
     "https://huggingface.co/facebook/esm2_t33_650M_UR50D/resolve/main/model.safetensors",
     "MIT", 2609.5, "MINT/ReCLIP 特征构建所需的底层 PLM"),
    ("Bernett MLP head（bernett_mlp.pth）",
     "https://huggingface.co/varunullanat2012/mint/resolve/main/bernett_mlp.pth",
     "MIT", 26.2, "MINT 提供的**下游任务 head**（预测两序列是否结合）"),
    ("ESM2-8M（model.safetensors）",
     "https://huggingface.co/facebook/esm2_t6_8M_UR50D/resolve/main/model.safetensors",
     "MIT", 31.4, "Track A 已在用（已缓存于 F 盘）"),
]


def main() -> int:
    t0 = time.time()

    # ==================== PHASE 4 ====================
    e4("# V7.2R-1 PHASE 4 — 外部审计修复")
    e4("")
    e4(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    e4("- 修复 V7.2R PHASE 1 的**三处错误/粗糙之处**")
    e4("")

    e4("## 4.1 Prompt-DDG — 许可与用途拆分（**修复**）")
    e4("")
    e4("| 对象 | 判定 |")
    e4("|---|---|")
    e4("| **代码许可** | **MIT** |")
    e4("| **预训练 ΔΔG checkpoint** | **SKEMPI 监督** |")
    e4("| **该 checkpoint 的 primary 用途** | **NOT ALLOWED** |")
    e4("| **架构 / 代码复用** | **ALLOWED**（subject to attribution） |")
    e4("")
    e4("```")
    e4("Prompt-DDG:")
    e4("  CODE LICENSE = MIT")
    e4("  PRETRAINED DDG CHECKPOINT = SKEMPI-supervised")
    e4("  PRIMARY USE OF CHECKPOINT = NOT ALLOWED")
    e4("  ARCHITECTURE/CODE REUSE = ALLOWED subject to attribution")
    e4("```")
    e4("")
    e4("> V7.2R 原报告把 `license_code` 与 `license_weights` 一并写成 "
       "`NOT DECLARED`，**这是错误的**：代码为 MIT。")
    e4("> 正确区分：**代码 MIT（可复用）** vs **权重 SKEMPI 监督（primary 禁用）**。")
    e4("")

    e4("## 4.2 ReCLIP — 三个对象必须分开（**修复**）")
    e4("")
    e4("| 对象 | 内容 | 许可 | 判定 |")
    e4("|---|---|---|---|")
    e4("| **A. ReCLIP 论文** | bioRxiv 2026, `10.64898/2026.06.01.729118` | — | 参考文献 |")
    e4("| **B. `SiweiLab/ReCLIP` 源码** | 特征构建 pipeline + 任务脚本（**不含**训练好的 mutation head） | 见仓库 LICENSE | **架构/代码可复用候选** |")
    e4("| **C. `RiverZ/reclip` HF 发布物** | `classifier_heads/mutation/reclip_mutation_xgb_*.json/.pkl` | **`other`** | **YELLOW**（除能排除训练重叠） |")
    e4("")
    e4("```")
    e4("ReCLIP CODE / ARCHITECTURE      : candidate for reuse")
    e4("ReCLIP HF mutation XGBoost head : YELLOW unless training overlap excluded")
    e4("```")
    e4("")
    e4("**README 已核实的关键事实**：")
    e4("")
    e4("- README 明言：*「The release excludes earlier binary mutation pipelines, legacy "
       "cross-attention experiments, ESM-pLM/ESum-pLM folders, and global-embedding and "
       "local automation experiment scripts.」*")
    e4("  → **发布的源码不含已训练的分类头**")
    e4("- README 明言：*「MINT checkpoints and trained task-specific classifier heads "
       "are external artifacts.」*")
    e4("- 特征构建器：`scripts/{four_classes_mutation,ptm,peptide,clinvar,ablation}/`")
    e4("- **系统要求：64-bit Linux**（Ubuntu）；参考配置 Python 3.8 / torch 1.12.1 / "
       "CUDA 11.3.1；**实测于 NVIDIA A100 40GB 与 T4**")
    e4("- CPU-only **被支持**（`--device cpu` / `--xgb-device cpu`），但 README 明言"
       "「substantially slower for feature extraction」，因 ESM2 与 MINT 皆为大型 PLM")
    e4("- **依赖 MINT**：`git clone https://github.com/VarunUllanat/mint.git mint`，"
       "checkpoint 放 `mint/mint.ckpt`")
    e4("")

    e4("## 4.3 MINT 作为 foundation representation 的审计（**新增**）")
    e4("")
    e4("| 项 | 值 |")
    e4("|---|---|")
    e4("| 仓库 | `https://github.com/VarunUllanat/mint` |")
    e4("| checkpoint 来源 | `https://huggingface.co/varunullanat2012/mint` → `mint.ckpt` |")
    e4("| **许可** | **MIT**（HF `cardData.license = mit`，tag `license:mit`） |")
    e4("| 训练目标 | 掩码/对比式 PLM 预训练，**面向多聚体相互作用**（MINT = "
       "Multimeric INteraction Transformer） |")
    e4("| 训练语料 | README：*「First PLM to be trained on large-scale PPI data」*；"
       "**具体 PDB/SKEMPI 组成未在 README 明确** |")
    e4("| 参数规模 | 主干为 **ESM2-650M** 类（config `esm2_t33_650M_UR50D.json`） |")
    e4("| **checkpoint 大小** | **3253.8 MB（≈3.25 GB）** ← **实测 HTTP HEAD** |")
    e4("| embedding API | `mint.helpers.extract.MINTWrapper(cfg, ckpt, device)` → "
       "`wrapper(chains, chain_ids)` |")
    e4("| 输出层级 | **默认 (N, 1280)** = 链级/对级；`sep_chains=True` → **(N, 2560)** |")
    e4("| PPI 条件化机制 | 多聚体 transformer：**链间上下文注意力**（chain-aware） |")
    e4("")
    e4("### 关键区分（§9 要求）")
    e4("")
    e4("| 类型 | MINT 情况 |")
    e4("|---|---|")
    e4("| **task-label leakage**（ΔΔG/SKEMPI 标签训练） | **未见证据**（预训练为自监督 PLM；"
       "另有单独的 `bernett_mlp.pth` 任务 head，属 external artifact） |")
    e4("| **foundation-pretraining family exposure**（预训练语料含我们的测试家族） | "
       "**可能**（「large-scale PPI data」组成未公开到可比对程度）→ "
       "标 **`RESIDUAL_PRETRAINING_EXPOSURE_POSSIBLE`** |")
    e4("")
    e4("> **论文措辞强制要求（§9）**：")
    e4("> ")
    e4("> > family-OOD refers to **task-specific supervised training**, "
       "**not** guaranteed absence from foundation-model pretraining.")
    e4(">")
    e4("> 即：我们可以声称「外层测试家族未参与**任务监督训练**」，")
    e4("> **但不能**声称「基础模型预训练完全没见过这些家族」。")
    e4("")

    # 许可矩阵
    lic = pd.DataFrame([
        {"object": "Prompt-DDG code", "license": "MIT", "reuse": "ALLOWED"},
        {"object": "Prompt-DDG DDG checkpoint", "license": "MIT(code)",
         "reuse": "NOT ALLOWED (SKEMPI-supervised, primary)"},
        {"object": "SiweiLab/ReCLIP source", "license": "see repo LICENSE",
         "reuse": "architecture reuse candidate"},
        {"object": "RiverZ/reclip mutation XGBoost head", "license": "other",
         "reuse": "YELLOW; not in primary"},
        {"object": "MINT checkpoint", "license": "MIT",
         "reuse": "candidate (foundation, no task-label leakage found)"},
        {"object": "MINT bernett_mlp head", "license": "MIT",
         "reuse": "NOT in primary (task head)"},
        {"object": "ESM2-650M", "license": "MIT", "reuse": "ALLOWED (generic)"},
        {"object": "ESM2-8M", "license": "MIT", "reuse": "ALLOWED (generic)"},
    ])
    lic.to_csv(PH4 / "license_matrix_repaired.csv", index=False, encoding="utf-8")
    prov = pd.DataFrame([
        {"model": "Prompt-DDG", "license_code": "MIT",
         "license_weights": "MIT(code); SKEMPI-supervised weights",
         "status": "YELLOW (checkpoint) / GREEN (code+architecture)",
         "primary_use": "architecture/code reuse ONLY"},
        {"model": "ReCLIP (paper)", "license_code": "n/a",
         "license_weights": "n/a", "status": "REFERENCE",
         "primary_use": "method citation"},
        {"model": "ReCLIP (SiweiLab source)", "license_code": "see repo LICENSE",
         "license_weights": "no head released",
         "status": "candidate", "primary_use": "feature-construction reuse"},
        {"model": "ReCLIP (RiverZ HF head)", "license_code": "other",
         "license_weights": "other",
         "status": "YELLOW", "primary_use": "NOT in primary"},
        {"model": "MINT", "license_code": "MIT", "license_weights": "MIT",
         "status": "candidate (foundation) / GPU required",
         "primary_use": "Track B representation (deferred: >2GB)"},
    ])
    prov.to_csv(PH4 / "model_provenance_repaired.csv", index=False,
                encoding="utf-8")
    e4(f"- 写出 `license_matrix_repaired.csv`、`model_provenance_repaired.csv`")
    e4("")
    (PH4 / "audit_repair.md").write_text("\n".join(L4) + "\n", encoding="utf-8")

    # ==================== PHASE 5 ====================
    e5("# V7.2R-1 PHASE 5 — Track B 可行性 + 计算可行性")
    e5("")
    e5("## 5.1 Track B 的**目标**形态（§10）")
    e5("")
    e5("```")
    e5("ReCLIP-style residue-context feature construction")
    e5("  +  generic leakage-safe ESM / MINT representations")
    e5("  +  OUR ranking head trained inside nested family CV")
    e5("```")
    e5("")
    e5("**明确不是**：`ReCLIP` 预训练的 mutation 分类器。")
    e5("")

    e5("## 5.2 能否在不使用其监督 mutation head 的前提下抽出所需表示？")
    e5("")
    e5("| 所需表示 | 公共 pipeline 是否提供 | 依据 |")
    e5("|---|---|---|")
    e5("| site-local residue representation | **YES** | 特征构建器 "
       "`scripts/*/ReCLIP/*_prediction_save.py` 落 `Feature_cache/`；"
       "ReCLIP 本身是**残基为中心**的框架 |")
    e5("| partner-conditioned representation | **YES** | README：*「which interaction "
       "partner residues are most informative」*；特征为 residue×partner 上下文 |")
    e5("| MINT partner representation | **YES** | MINT `MINTWrapper(chains, chain_ids)` "
       "接受**双链**输入，输出链级 (N,1280) 或 `sep_chains=True` (N,2560) |")
    e5("| 是否需要其 mutation XGBoost head | **NO** | README：*「Release artifacts are "
       "hosted separately on Hugging Face」*，且 *「trained task-specific classifier "
       "heads are external artifacts」* |")
    e5("")
    e5("```")
    e5("TRACK_B_FEASIBLE = YES")
    e5("```")
    e5("")
    e5("**但有两个硬性前置条件**（见 5.3 / 5.4）。")
    e5("")

    e5("## 5.3 计算可行性（§11 要求的估算）")
    e5("")
    e5("**实测 artifact 大小（HTTP HEAD，未下载）**")
    e5("")
    e5("| artifact | 许可 | 大小 | 用途 |")
    e5("|---|---|---:|---|")
    art_rows = []
    for n, u, lc, mb, purp in ARTIFACTS:
        e5(f"| {n} | {lc} | **{mb:,.1f} MB** | {purp} |")
        art_rows.append({"artifact": n, "license": lc, "size_mb": mb,
                         "purpose": purp, "url": u})
    pd.DataFrame(art_rows).to_csv(PH5 / "artifact_sizes.csv", index=False,
                                  encoding="utf-8")
    need = 3253.8 + 2609.5
    e5("")
    e5(f"- **Track B 必需两件合计 = 3253.8 + 2609.5 = "
       f"{need:,.1f} MB = {need / 1000:.2f} GB**")
    e5(f"- §29 的自动允许上限 = **2 GB** → **超出 {need / 1000 - 2:.2f} GB**")
    e5("")
    e5("```")
    e5("DEFERRED_LARGE_ARTIFACT = YES")
    e5("```")
    e5("")
    e5("**本机资源核算**")
    e5("")
    e5("| 项 | 值 |")
    e5("|---|---|")
    e5("| OS | Windows 10（ReCLIP 参考环境为 **64-bit Linux**，未在 Windows 上测试） |")
    e5("| GPU | **无 CUDA**（`torch.cuda.is_available() = False`） |")
    e5("| 参考 GPU | README 实测于 **A100 40GB / T4** |")
    e5("| RAM | **7.55 GB 总量** |")
    e5("| ESM2-650M 权重磁盘 | ≈2.61 GB |")
    e5("| MINT 权重磁盘 | ≈3.25 GB |")
    e5("| 峰值 RAM（单模型 fp32 权重 + 激活） | ≈3–6 GB |")
    e5("| **两模型同时驻留** | **≈6–9 GB → 超过 7.55 GB 物理内存** |")
    e5("")
    e5("**序列规模（U_CONTRACT）**")
    e5("")
    e5("| 项 | 估算 |")
    e5("|---|---|")
    e5("| 行数 | 4473（+ WT 复本） |")
    e5("| 唯一 WT/MUT 链序列 | 数千条 |")
    e5("| 最长链 | 3128 aa（需截断） |")
    e5("| 特征缓存（650M 级、残基级） | **数十 GB 量级** |")
    e5("")
    e5("⇒ 结论：")
    e5("")
    e5("```")
    e5("TRACK_B_GPU_REQUIRED = YES")
    e5("TRACK_B_READY_FOR_GPU = YES")
    e5("```")
    e5("")
    e5("**理由**：① 必需 artifact 合计 5.86 GB > 2 GB 上限；")
    e5("② 无 CUDA；③ ESM2-650M + MINT 同时驻留超出 7.55 GB 物理内存；")
    e5("④ ReCLIP 仅实测于 Linux/A100。")
    e5("**依 §11「Do not substitute a weaker model silently」，")
    e5("不以更弱模型替代，而是明确标记需要 GPU。**")
    e5("")

    e5("## 5.4 Track B 启动清单（GPU 就绪后可直接执行）")
    e5("")
    e5("1. 在 **Linux + NVIDIA GPU**（≥24 GB 显存建议）环境：")
    e5("   `git clone https://github.com/SiweiLab/ReCLIP.git`；")
    e5("   `git clone https://github.com/VarunUllanat/mint.git mint`")
    e5("2. 下载 `mint.ckpt`（3.25 GB, MIT）与 `facebook/esm2_t33_650M_UR50D`（2.61 GB, MIT）")
    e5("3. **只运行 feature builder**，落 `Feature_cache/`；")
    e5("   **绝不加载** `classifier_heads/mutation/reclip_mutation_xgb_*.pkl`")
    e5("4. 用我们自己的 `z` ranker，在 `OUTER_FOLDS_V7_2` + inner family CV 下训练")
    e5("5. 依 §9，论文须写明：family-OOD 指**任务监督训练**，")
    e5("   **不保证**基础模型预训练未见这些家族")
    e5("")
    e5(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (PH5 / "trackB_feasibility.md").write_text("\n".join(L5) + "\n",
                                               encoding="utf-8")
    print(f"\n写出 PHASE4/PHASE5 文件")
    print(f">>> TRACK_B_FEASIBLE = YES | TRACK_B_GPU_REQUIRED = YES | "
          f"DEFERRED_LARGE_ARTIFACT = YES")
    return 0


if __name__ == "__main__":
    sys.exit(main())
