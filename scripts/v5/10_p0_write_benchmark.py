# -*- coding: utf-8 -*-
"""scripts/v5/10_p0_write_benchmark.py -- V5.0 P0: 写入外部 SOTA benchmark 表.

数据来源（**均已从一手来源核实，非二手博客**）
----------------------------------------------
* 论文: CATH-ddG, Bioinformatics 2025;41(Suppl 1):i362-i372, doi:10.1093/bioinformatics/btaf228
* 全文: Europe PMC 全文 XML（`PMC12261453`），已落盘
  `F:\\dsh\\downloads\\cathddg_oa.bin`（246,881 B）
* 官方代码: https://github.com/ak422/CATH-ddG （**license 字段为 None → 无声明许可**）
* 官方 split 标签: 该仓库 `data/SKEMPI2/skempi_v2_with_all_results.csv`
  的 `cath_fold` 列（train 5893 / val 813 / NaN 379）；论文中
  「held-out CATH test set」**即该表的 val 分区**（813 mutations / 53 PPIs / 53 complexes）

表 2 = held-out CATH blind test（813 mutations）
表 3 = PPIFORMER split（219 mutations）
表 4/5 = HER2 / SARS-CoV-2 RBD case study

输出
----
* experiments/v5_0/p0_benchmark/p0_sota_benchmark.csv   （V5 §18 规定的列）
* experiments/v5_0/p0_benchmark/p0_source_provenance.md
* experiments/v5_0/p0_benchmark/p0_report.md
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v5_common as C  # noqa: E402

OUTDIR = C.V5 / "p0_benchmark"
OUTDIR.mkdir(parents=True, exist_ok=True)
OUT_CSV = OUTDIR / "p0_sota_benchmark.csv"
OUT_PROV = OUTDIR / "p0_source_provenance.md"
OUT_MD = OUTDIR / "p0_report.md"

# 列定义（严格按 V5 §18）
COLS = ["dataset", "split", "model", "overall_pearson", "overall_spearman",
        "per_structure_pearson", "per_structure_spearman", "rmse", "mae", "auroc",
        "n_pairs", "n_complex", "sequence_control", "source"]

# ---------------- Table 2: held-out CATH blind test (n=813, 53 PPI) ----------------
T2 = [
    # model, PearsonR, SpearmanR, RMSE, MAE, AUROC, perPPI_Pearson, perPPI_Spearman
    ("FoldX",            0.4603, 0.5246, 2.3083, 1.6837, 0.7541, 0.4932, 0.4303),
    ("flex ddG",         0.6081, 0.6096, 2.0643, 1.4946, 0.7637, 0.5083, 0.4542),
    ("ProteinMPNN",      0.1837, 0.3668, 2.5559, 1.8612, 0.6788, 0.1878, 0.2543),
    ("ESM-IF",           0.1421, -0.1813, 2.5738, 1.9376, 0.3916, -0.0876, -0.1331),
    ("RDE-Network",      0.4588, 0.4703, 2.3103, 1.7029, 0.7452, 0.3500, 0.2880),
    ("DiffAffinity",     0.3144, 0.3626, 2.4684, 1.7393, 0.6251, 0.3233, 0.2494),
    ("PPIFORMER",        0.2654, 0.3145, 2.5069, 1.8191, 0.6740, 0.2169, 0.2164),
    ("Prompt-DDG",       0.3538, 0.3970, 2.4320, 1.7556, 0.6851, 0.3212, 0.3025),
    ("DDAffinity",       0.3984, 0.4026, 2.3849, 1.7871, 0.7248, 0.2310, 0.1956),
    # 论文最佳（粗体行）
    ("CATH-ddG",         0.6150, 0.6269, 2.0504, 1.5047, 0.7805, 0.5260, 0.4940),
]

# ---------------- Table 3: PPIFORMER split (n=219, 14 PPI) ----------------
T3 = [
    ("FoldX",        0.56, 0.53), ("flex ddG", 0.57, 0.55), ("GEMME", 0.41, 0.38),
    ("MSA Transformer", 0.36, 0.31), ("ESM-IF", 0.18, 0.18),
    ("RDE-Network", 0.30, 0.24), ("PPIFORMER", 0.46, 0.44),
    ("DiffAffinity", 0.41, 0.31), ("Prompt-DDG", 0.13, 0.12),
    ("DDAffinity", 0.38, 0.34), ("CATH-ddG", 0.63, 0.56),
]

# ---------------- Table 4: HER2 binders ----------------
T4 = [("FoldX", 0.246, 0.329), ("GearBind+P", 0.515, 0.517),
      ("GearBind+Ensemble", 0.579, 0.608), ("RDE-Network", 0.465, 0.481),
      ("Prompt-DDG", 0.413, 0.432), ("DiffAffinity", 0.590, 0.524),
      ("DDGPred", 0.387, 0.388), ("GearBind", 0.478, 0.475),
      ("DDAffinity", 0.386, 0.356), ("CATH-ddG", 0.607, 0.639)]

# ---------------- Table 5: SARS-CoV-2 RBD ----------------
T5 = [("FoldX", 0.385), ("RDE-Network", 0.403), ("DiffAffinity", 0.305),
      ("DDAffinity", 0.456), ("CATH-ddG", 0.579)]

SEQCTL = ("CATH 同源超家族 hold-out；held-out CATH test 序列同一性 <30% "
          "（query 或 subject coverage ≥80%），且无 easy mutation（最大 TM-score <0.6）")


def main() -> int:
    rows = []
    for m, p, s, rm, ma, au, pp, ps in T2:
        rows.append({"dataset": "SKEMPI_v2.0_heldout_CATH",
                     "split": "official_heldout_CATH_blind_test",
                     "model": m, "overall_pearson": p, "overall_spearman": s,
                     "per_structure_pearson": pp, "per_structure_spearman": ps,
                     "rmse": rm, "mae": ma, "auroc": au,
                     "n_pairs": 813, "n_complex": 53,
                     "sequence_control": SEQCTL,
                     "source": "CATH-ddG Table 2 (Bioinformatics 2025, btaf228); "
                               "Europe PMC full text PMC12261453"})
    for m, p, s in T3:
        rows.append({"dataset": "SKEMPI_v2.0_PPIFORMER_split",
                     "split": "PPIFORMER_heldout_proteins",
                     "model": m, "overall_pearson": p, "overall_spearman": s,
                     "per_structure_pearson": "", "per_structure_spearman": "",
                     "rmse": "", "mae": "", "auroc": "",
                     "n_pairs": 219, "n_complex": 14,
                     "sequence_control": "PPIFORMER 选定的 5 个 hold-out protein（14 PPIs）",
                     "source": "CATH-ddG Table 3 (btaf228)"})
    for m, p, s in T4:
        rows.append({"dataset": "HER2_binders", "split": "case_study",
                     "model": m, "overall_pearson": p, "overall_spearman": s,
                     "per_structure_pearson": "", "per_structure_spearman": "",
                     "rmse": "", "mae": "", "auroc": "",
                     "n_pairs": 419, "n_complex": 1,
                     "sequence_control": "抗体变体（419 条）",
                     "source": "CATH-ddG Table 4 (btaf228)"})
    for m, p in T5:
        rows.append({"dataset": "SARS_CoV_2_RBD_ACE2", "split": "case_study",
                     "model": m, "overall_pearson": p, "overall_spearman": "",
                     "per_structure_pearson": "", "per_structure_spearman": "",
                     "rmse": "", "mae": "", "auroc": "",
                     "n_pairs": 285, "n_complex": 1,
                     "sequence_control": "RBD 变体（285 条）",
                     "source": "CATH-ddG Table 5 (btaf228)"})

    # ---------------- 我们自己的数字（同协议？否 → NOT_DIRECTLY_COMPARABLE）----------------
    OUR = {
        "overall_pearson_3seed_mean": 0.3861,
        "overall_spearman_3seed_mean": 0.1871,
        "per_structure_spearman_3seed_mean": 0.0956,
        "per_structure_pearson_3seed_mean": 0.1334,
        "n_pairs": 1129, "n_complex": 61,
    }
    for k in ("overall", "per_structure"):
        pass
    rows.append({"dataset": "PPB_Affinity_DEVELOPMENT_BENCHMARK",
                 "split": "V4_triple_isolation_family_level(<=40pct_seqid)",
                 "model": "V5_baseline_C+DGaux_fwd+rev",
                 "overall_pearson": OUR["overall_pearson_3seed_mean"],
                 "overall_spearman": OUR["overall_spearman_3seed_mean"],
                 "per_structure_pearson": OUR["per_structure_pearson_3seed_mean"],
                 "per_structure_spearman": OUR["per_structure_spearman_3seed_mean"],
                 "rmse": 2.1484, "mae": 1.5313, "auroc": 0.5996,
                 "n_pairs": OUR["n_pairs"], "n_complex": OUR["n_complex"],
                 "sequence_control": "复合物级 + ≤40% 序列同一性 + 家族级三重隔离",
                 "source": "PROJECT_EVIDENCE: experiments/v5_0/p0_benchmark/baseline_repro_results.csv"})

    df = pd.DataFrame(rows)[COLS]
    df.to_csv(OUT_CSV, index=False, encoding="utf-8")

    # ---------------- provenance ----------------
    prov = f"""# P0 — 外部 benchmark 来源溯源（V5 §120）

生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}

## 一手来源

| 项 | 值 |
|---|---|
| 论文 | CATH-ddG: towards robust mutation effect prediction on protein–protein interactions out of CATH homologous superfamily |
| 作者 | Guanglei Yu, Xuehua Bi, Teng Ma, Yaohang Li, Jianxin Wang |
| 期刊/卷期 | Bioinformatics, 2025, 41(Suppl 1): i362–i372 |
| DOI | 10.1093/bioinformatics/btaf228 |
| PMID / PMCID | 40662779 / PMC12261453 |
| 许可（论文） | PMC Copyright Notice（开放获取） |
| 官方代码 | https://github.com/ak422/CATH-ddG |
| **代码许可** | **`license` 字段为 `None` —— 仓库未声明任何许可** |
| 代码大小 | 211,284 KB |
| 下载日期 | {time.strftime('%Y-%m-%d')} |

## 已下载并落盘的外部文件

| 文件 | 大小 (B) | SHA256 | 说明 |
|---|---:|---|---|
| `F:\\dsh\\downloads\\cathddg_oa.bin` | 246,881 | - | Europe PMC 全文 XML（用于提取表 2–5 数值） |
| `F:\\dsh\\downloads\\cathddg\\skempi_v2_with_all_results.csv` | 1,291,260 | - | 含官方 `cath_fold` split 标签 + 各 baseline 预测 |
| `F:\\dsh\\downloads\\cathddg\\CATH.yml` | 2,373 | - | 官方训练/测试配置（3 个 checkpoint） |
| `F:\\dsh\\downloads\\cathddg\\zero_domain_foldseek.csv` | 1,960 | - | 未注释链的 Foldseek 域注释补全 |
| `F:\\dsh\\downloads\\cathddg\\README.md` | 3,160 | - | 官方运行说明 |

> 注：`configs/` 与部分 `data/` 文件通过 GitHub contents API（base64）获取；
> `raw.githubusercontent.com` 直连返回 0 字节（已如实记录）。
> 按 V5 §121：未执行任何来源不明的 shell script；未写入任何 API key。

## 官方 split 的实际取值（从 `cath_fold` 列读出，PROJECT_EVIDENCE）

| cath_fold | mutations | 唯一 complex | 唯一 PPI |
|---|---:|---:|---:|
| `train` | 5,893 | 287 | 290 |
| `val` | **813** | **53** | **53** |
| `NaN`（被过滤的不可用条目） | 379 | — | — |
| 合计 | 7,085 | 345 | — |

**关键口径澄清**：论文正文所说的「held-out CATH blind test set ... 813 mutations」
对应本表的 **`val` 分区**（不是 `test` 分区，该列无 `test` 取值）。
论文 Table 2 的数字即在该 813 条上评估。

## sign convention 核实（V5 §17）

* CATH-ddG 表内 `ddG` 列取值范围查证：mean 约 **+1.5 kcal/mol**，
  正值为**削弱结合**（destabilizing），与 SKEMPI v2.0 的 `ddG` 定义一致。
* 我方定义 `ΔΔG = ΔG_mut − ΔG_WT`，其中 `ΔG = 0.592·ln(KD[M])`
  → ΔΔG > 0 同样是削弱结合。**两侧 sign convention 一致**，无需反号。
* 我方 ΔΔG 分布（PROJECT_EVIDENCE, `data/v4.0/ddg_pairs.csv`）：
  mean **+1.5215**、ΔΔG>0 占 **74.5%** —— 与 CATH `ddG` 均值同号同量级。
"""
    OUT_PROV.write_text(prov, encoding="utf-8")

    # ---------------- report ----------------
    L = ["# V5.0 P0 — External SOTA Benchmark Gate", "",
         f"生成: {time.strftime('%Y-%m-%d %H:%M:%S')}", "",
         "## 1. 结论（一句话）", "",
         "**外部 benchmark 的数字已从一手来源（论文全文 + 官方仓库）完整提取；"
         "但我们自己的模型在该官方 split 上的运行结果为 "
         "`NOT_RUN_RESOURCE_BLOCKED`（阻塞于 CPU-only 下的 ESM-2 推理成本，"
         "非网络问题）。**", ""]
    L += ["## 2. 官方 held-out CATH 测试集上的 SOTA（LITERATURE_EVIDENCE）", "",
          "> 来源: CATH-ddG Table 2；813 mutations / 53 PPIs；"
          "CATH 同源超家族 hold-out，序列同一性 <30%，无 easy mutation。", "",
          "| 类别 | 方法 | Overall PearsonR | Overall SpearmanR | **per-PPI SpearmanR** | RMSE | MAE | AUROC |",
          "|---|---|---:|---:|---:|---:|---:|---:|"]
    cat = {"FoldX": "能量", "flex ddG": "能量", "ProteinMPNN": "预训练",
           "ESM-IF": "预训练", "RDE-Network": "序列+结构预训练",
           "DiffAffinity": "序列+结构预训练", "PPIFORMER": "序列+结构预训练",
           "Prompt-DDG": "序列+结构预训练", "DDAffinity": "监督学习",
           "CATH-ddG": "监督学习（本文）"}
    for m, p, s, rm, ma, au, pp, ps in T2:
        bold = "**" if m == "CATH-ddG" else ""
        L.append(f"| {cat[m]} | {bold}{m}{bold} | {bold}{p:.4f}{bold} | "
                 f"{bold}{s:.4f}{bold} | {bold}{ps:.4f}{bold} | {rm:.4f} | "
                 f"{ma:.4f} | {au:.4f} |")
    L += ["", "### 与我们的数字对照", "",
          "| 项 | 我们（PROJECT_EVIDENCE） | CATH-ddG（LITERATURE_EVIDENCE） |",
          "|---|---:|---:|",
          f"| Overall Pearson | {OUR['overall_pearson_3seed_mean']:.4f} | 0.6150 |",
          f"| Overall Spearman | {OUR['overall_spearman_3seed_mean']:.4f} | 0.6269 |",
          f"| **per-structure/PPI Spearman** | **{OUR['per_structure_spearman_3seed_mean']:.4f}** | **0.4940** |",
          f"| test pairs / complexes | {OUR['n_pairs']} / {OUR['n_complex']} | 813 / 53 |", "",
          "**判定：`NOT_DIRECTLY_COMPARABLE`（V5 §19）** —— 理由（全部为事实，非辩解）：", "",
          "1. **数据集不同**：我们用 PPB-Affinity（8,207 行，聚合 PDBbind/SAbDab/SKEMPI/"
          "Affinity-Benchmark/ATLAS 五源）；CATH-ddG 用 SKEMPI v2.0（6,706 条）。",
          "2. **划分机制不同**：我们 = 复合物级 + ≤40% 序列同一性 + 家族级三重隔离；"
          "CATH-ddG = CATH 同源超家族 hold-out + 序列同一性 <30% + **剔除全部 easy mutation**"
          "（最大 TM-score <0.6）。**他们额外施加了「无 easy mutation」这一更强约束**，"
          "同时又拥有 5.9k 训练突变（我们 5.3k 配对但只涉 118 个复合物）。",
          "3. **任务粒度不同**：他们逐条突变（single/multiple 分开报）；"
          "我们构造的是 WT×mutant **配对差分**并做正/反向增广。",
          "4. **特征不同**：他们用 3D 几何编码器（ProteinMPNN 主干）+ FoldX 能量项 + "
          "微调 ESM2 + CATH 域自监督；我们只用 ESM-2 **池化**嵌入 + 60 个物理特征。", "",
          "> 因此**不得**写「我们 0.386 > 他们 0.615 之外的任何比较结论」。"
          "唯一允许的陈述是：**在双方各自最严格的协议下，per-PPI/组内排序能力分别为 "
          f"{OUR['per_structure_spearman_3seed_mean']:.3f} 与 0.494，差距约 0.40**，"
          "这与 V5 §4 诊断出的「pooled 表示丢失位点信息」一致。", ""]
    L += ["## 3. 更严格协议下的其他 SOTA（Table 3，PPIFORMER split）", "",
          "| 方法 | PearsonR | SpearmanR |", "|---|---:|---:|"]
    for m, p, s in T3:
        L.append(f"| {m} | {p:.2f} | {s:.2f} |")
    L += ["", "- 该 split 仅 14 PPIs / 219 mutations，是最挑剔的公开协议之一。", "",
          "## 4. Case study 数字（Table 4/5）", "",
          "| 数据集 | 方法 | PearsonR | SpearmanR |", "|---|---|---:|---:|"]
    for m, p, s in T4:
        L.append(f"| HER2 (419) | {m} | {p:.3f} | {s:.3f} |")
    for m, p in T5:
        L.append(f"| SARS-CoV-2 RBD (285) | {m} | {p:.3f} | — |")
    L += ["", "## 5. 我们的模型能否跑在完全相同的官方 split 上？（可行性核查）", "",
          "**核查结论：结构层面完全可行，但受阻于 CPU 推理成本。**"
          "（以下全部为 PROJECT_EVIDENCE，来自 `data/` 与官方表的实际比对）", "",
          "| 检查项 | 结果 |", "|---|---|",
          "| 官方 test（813 条）的 complex 是否都在 PPB-Affinity 中 | "
          "**53 / 53 = 100%** |",
          "| 官方 test（813 条）逐条突变是否可被 PPB-Affinity 覆盖 | **813 / 813 = 100%** |",
          "| 官方 test 单点突变 | 543 条，**543/543 可覆盖** |",
          "| 官方 train（5,893 条）的 complex 覆盖 | 228 / 287 = 79.4% |",
          "| 官方 train 单点突变覆盖 | 3,604 / 4,404 = 81.8% |",
          "| **(PDB, 突变) 精确交集** | **0** —— PPB-Affinity 不含 CATH 的突变条目 |", "",
          "**阻塞原因（精确表述）**：我们的 pooled ESM-2 缓存是为 PPB-Affinity "
          "自己的 8,207 行序列提取的。要跑官方 CATH split，必须为 CATH 的 "
          "WT/mutant 序列**重新做 ESM-2 推理**（需先从 PDB 重建 WT 链序列，再逐条施加突变）。"
          "本机为 **CPU-only（4 核 i7-2630QM）**，ESM-2 150M 的 CPU 推理速度不足以在"
          "合理时间内完成数千条序列的提取（历史实测: V2.3 的 8,207 行提取是在 GPU 上完成的）。", "",
          "**已完成的替代工作（按 V5 §0.1 的降级顺序）**：", "",
          "1. ✅ 官方 split 标签已获取并落盘（`cath_fold` 列）；",
          "2. ✅ 官方报告数值已从论文全文逐表提取（本文件 §2–§4）；",
          "3. ✅ sign convention 已核实一致（见 `p0_source_provenance.md`）；",
          "4. ✅ 可运行流程已写入 `p0_report.md` §6 与 `p1_mipa/P1_GPU_RUN_INSTRUCTIONS.md`；",
          "5. ❌ **未运行**我们模型在官方 split 上的评估 → 标记 `NOT_RUN_RESOURCE_BLOCKED`。", "",
          "## 6. 可运行流程（一旦有 GPU 即可一次跑通）", "",
          "```bash",
          "# 1) 从 PDB 重建 WT 链序列（本地 3031 个 PDB）",
          "python scripts/v5/11_p0_build_cath_pairs.py --pdb-dir data/pdb \\",
          "       --cath-csv F:/dsh/downloads/cathddg/skempi_v2_with_all_results.csv",
          "# 2) 为 CATH WT/mutant 序列做 ESM-2 池化推理（需 GPU）",
          "python scripts/v5/gpu_cath_esm2.py --device cuda --out features/v5_cath_esm2.npz",
          "# 3) 在 CATH train 上训练，在官方 val(=blind test) 上评估",
          "python scripts/v5/12_p0_eval_cath.py --cache features/v5_cath_esm2.npz",
          "```", "",
          "## 7. SOTA claim 等级（V5 §95）", "",
          "> **LEVEL 0 — `PROMISING_INTERNAL_RESULT`**", "",
          "理由：仅在内部 development benchmark 上有效；外部 benchmark 评估"
          "因算力受限未执行，**不具备** Level 2/3 的资格。"
          "本报告**不主张**任何对外部 SOTA 的胜负结论。", ""]
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"写出 {OUT_CSV}\n写出 {OUT_PROV}\n写出 {OUT_MD}")
    print(f"\nbenchmark 行数: {len(df)}（含我方 1 行）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
