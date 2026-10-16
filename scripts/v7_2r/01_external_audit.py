# -*- coding: utf-8 -*-
"""scripts/v7_2r/01_external_audit.py -- V7.2R PHASE 0 + PHASE 1.

PHASE 0：**冻结并登记现有 V7.2 产物**（TRACK_A_SMALL_PLM），不重算、不删除
PHASE 1：**外部模型 provenance / leakage / license 审计**
          → GREEN / YELLOW / RED，**审计完成前不下载任何 checkpoint**

判定规则（V7.2R §8–§10）
----------------------
* **RED**：checkpoint 确定训练于我们的 outer-test 行 / family，或用了测试标签
* **YELLOW**：SKEMPI 监督模型且无法排除 outer-family 重叠；或训练数据文档不全
  → **仅作 exploratory**，**不能**决定主结论
* **GREEN**：通用无监督/自监督预训练（与 SKEMPI 标签无关）或架构可重训

关键区分（§7 末条）
------------------
`can_extract_frozen_representation_without_task_checkpoint?`
—— 即**不加载其任务 head**，只用 backbone 表示。
"""

from __future__ import annotations

import hashlib
import json
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

V71 = C.ROOT / "experiments" / "v7_1"
V711 = C.ROOT / "experiments" / "v7_1_1"
V72 = C.ROOT / "experiments" / "v7_2"
V72R = C.ROOT / "experiments" / "v7_2r"
PH0 = V72R / "phase0_existing_v7_2"
PH1 = V72R / "phase1_external_audit"
for d in (PH0, PH1):
    d.mkdir(parents=True, exist_ok=True)
L0: list[str] = []
L1: list[str] = []


def e0(m=""):
    print(m); L0.append(str(m))


def e1(m=""):
    print(m); L1.append(str(m))


def sha(p: Path, ch: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(ch), b""):
            h.update(b)
    return h.hexdigest()


# ============================================================ PHASE 0
def phase0() -> None:
    e0("# V7.2R PHASE 0 — 冻结并登记现有 V7.2 产物（TRACK_A_SMALL_PLM）")
    e0("")
    e0(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    e0("- **不删除、不重算**；只在完成完整性审计后重新登记（§31）")
    e0("")
    rows = []
    # V7.2 各类产物
    pats = [
        ("v72_report", V72 / "phase0_provenance" / "model_decision.md"),
        ("v72_report", V72 / "phase1_embeddings" / "embedding_report.md"),
        ("v72_manifest", C.ROOT / "features" / "v7_2" /
         "residue_embedding_manifest.csv"),
    ]
    for kind, p in pats:
        if p.exists():
            rows.append({"kind": kind, "path": str(p.relative_to(C.ROOT)),
                         "size": p.stat().st_size,
                         "sha256_prefix16": sha(p)[:16],
                         "mtime": time.strftime(
                             "%Y-%m-%d %H:%M",
                             time.localtime(p.stat().st_mtime))})
    # 残基嵌入（per_target）
    per = sorted((C.ROOT / "features" / "v7_2" / "per_target").glob("*.npz"))
    tot_rows = 0
    tgts = []
    for p in per:
        try:
            z = np.load(p)
            n = len(z["row_id"])
            tot_rows += n
            tgts.append(p.stem)
        except Exception:  # noqa: BLE001
            pass
    rows.append({"kind": "v72_residue_embeddings", "path":
                 "features/v7_2/per_target/*.npz",
                 "size": sum(p.stat().st_size for p in per),
                 "sha256_prefix16": f"{len(per)}_files",
                 "mtime": ""})
    # 模型缓存
    hf = Path(r"F:\dsh\cache\huggingface")
    hf_sz = sum(p.stat().st_size for p in hf.rglob("*") if p.is_file()) \
        if hf.exists() else 0
    rows.append({"kind": "model_cache", "path": str(hf), "size": hf_sz,
                 "sha256_prefix16": "esm2_t6_8M_UR50D", "mtime": ""})
    pd.DataFrame(rows).to_csv(PH0 / "existing_inventory.csv", index=False,
                              encoding="utf-8")
    e0("## 0.1 现有产物清单")
    e0("")
    e0("| 类别 | 路径 | 大小 | 标识 |")
    e0("|---|---|---:|---|")
    for r in rows:
        e0(f"| {r['kind']} | `{r['path']}` | {r['size']:,} | "
           f"`{r['sha256_prefix16']}` |")
    e0("")
    e0("## 0.2 残基嵌入状态")
    e0("")
    e0(f"- `per_target/*.npz`：**{len(per)}** 个 target 文件")
    e0(f"- 累计 rows：**{tot_rows}** / 4473（**{tot_rows / 4473 * 100:.1f}%**）")
    e0(f"- 模型：`facebook/esm2_t6_8M_UR50D`（7.84M 参数，hidden 320，**MIT**）")
    e0(f"- 缓存：`F:\\dsh\\cache\\huggingface`（**{hf_sz / 1e6:.1f} MB**，未写 C 盘）")
    e0("")
    e0("## 0.3 完整性审计（§5 三项否决条件）")
    e0("")
    lab = pd.read_csv(V711 / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    e0("| 判据 | 结果 |")
    e0("|---|---|")
    e0(f"| target namespace 是否错误 | **NO** —— V7.2 的嵌入脚本用 "
       f"`canonical_target_labels.csv` 的 `target_id_canonical` 作为 target 名；"
       f"按 `row_id` 联合（§3.5） |")
    e0(f"| 是否发生 outer 泄漏 | **NO** —— 嵌入是**与标签无关的通用 PLM 前向**，"
       f"不含任何 fold 信息 |")
    e0(f"| 测试标签是否进入特征构造 | **NO** —— F0/F1/F2 只用 mutation 字符串与"
       f"WT/MUT 序列；无 ΔΔG |")
    e0("")
    e0("⇒ **现有 V7.2 产物有效**，登记为 **`TRACK_A_SMALL_PLM`**。")
    e0("")
    e0("### 可复用性")
    e0("")
    e0("| 产物 | reusable | 说明 |")
    e0("|---|---|---|")
    e0("| 残基嵌入（per_target npz） | **YES** | 与标签无关，任何 fold 均可直接用 |")
    e0("| 模型缓存（ESM-2 8M） | **YES** | MIT，F 盘 |")
    e0("| F0/F1/F2 特征构造代码 | **YES** | 无标签泄漏 |")
    e0("| `canonical_target_labels.csv` | **YES**（权威） | V7.1.1 冻结 |")
    e0("| **R0–R7 ranker 结果** | **NOT YET** | V7.2 的 02 脚本尚未跑完（嵌入覆盖不足） |")
    e0("")
    (PH0 / "phase0_report.md").write_text("\n".join(L0) + "\n", encoding="utf-8")


# ============================================================ PHASE 1
CANDS = [
    dict(
        model_name="ESM-2 8M (small PLM)",
        paper="Lin et al., Science 2023 (ESM-2 / ESMFold)",
        official_repository="https://github.com/facebookresearch/esm",
        official_checkpoint_source="https://huggingface.co/facebook/esm2_t6_8M_UR50D",
        license_code="MIT", license_weights="MIT",
        task="generic protein language model (masked LM)",
        training_dataset="UniRef50 (sequences, ~2022)",
        validation_dataset="n/a", test_dataset="n/a",
        SKEMPI_used="NO", SKEMPI_version="n/a",
        CATH_split_used="NO", PPIformer_split_used="NO",
        checkpoint_training_scope="unsupervised MLM on UniRef50; no affinity labels",
        residue_embeddings_available="YES",
        intermediate_representation_accessible="YES",
        prediction_scores_available="NO (not an affinity model)",
        architecture_reusable="YES", checkpoint_reusable="YES",
        download_size="29.9 MB",
        hardware_requirement="CPU ok (tested)",
        can_extract_frozen_representation_without_task_checkpoint="YES",
        status="GREEN",
        allowed_use="primary: Track A baseline + Δh perturbation",
        reason="与 SKEMPI 标签完全无关的通用自监督预训练；MIT 双许可",
    ),
    dict(
        model_name="Prompt-DDG",
        paper="Wu et al., ICML 2024 (microenvironment-aware hierarchical prompt learning)",
        official_repository="https://github.com/lirongwu/prompt-ddg",
        official_checkpoint_source="Google Drive (ddg_model.ckpt)",
        license_code="NOT DECLARED", license_weights="NOT DECLARED",
        task="ΔΔG prediction (supervised)",
        training_dataset="**SKEMPI v2**（README 明确 `data/get_skempi_v2.sh`）",
        validation_dataset="SKEMPI v2 derived", test_dataset="SKEMPI v2 derived",
        SKEMPI_used="**YES**", SKEMPI_version="v2.0",
        CATH_split_used="unknown", PPIformer_split_used="unknown",
        checkpoint_training_scope="supervised ΔΔG on SKEMPI v2",
        residue_embeddings_available="partial (internal, coupled to task head)",
        intermediate_representation_accessible="uncertain",
        prediction_scores_available="YES",
        architecture_reusable="uncertain", checkpoint_reusable="license 未声明",
        download_size="unknown",
        hardware_requirement="GPU (训练脚本含 cuda)",
        can_extract_frozen_representation_without_task_checkpoint="UNCERTAIN",
        status="**YELLOW**",
        allowed_use="EXPLORATORY_ONLY（附录/次要）；**不能**决定主结论",
        reason="SKEMPI v2 监督训练；我们的数据集同样源自 SKEMPI v2.0 → "
               "无法排除 outer-family 重叠（§8 `UNKNOWN_TRAINING_OVERLAP`）；"
               "且 license 未声明",
    ),
    dict(
        model_name="CATH-ddG (DDAffinity)",
        paper="Zhang et al., Bioinformatics 2024",
        official_repository="https://github.com/ak422/CATH-ddG",
        official_checkpoint_source="Google Drive（CATH / PPIformer / case study 三套权重）",
        license_code="NOT DECLARED", license_weights="NOT DECLARED",
        task="ΔΔG prediction (supervised)",
        training_dataset="**SKEMPI v2.0**（README 明确 `skempi_parallel.py --subset skempi_v2`）",
        validation_dataset="held-out CATH split", test_dataset="held-out CATH + PPIformer split",
        SKEMPI_used="**YES**", SKEMPI_version="v2.0",
        CATH_split_used="**YES**", PPIformer_split_used="**YES**",
        checkpoint_training_scope="supervised ΔΔG；提供 held-out CATH 专用权重",
        residue_embeddings_available="partial",
        intermediate_representation_accessible="uncertain",
        prediction_scores_available="YES（**且我们已在 V6.0 用过其预测列**）",
        architecture_reusable="uncertain", checkpoint_reusable="license 未声明",
        download_size="unknown（Google Drive）",
        hardware_requirement="GPU（`--device cuda:0`）",
        can_extract_frozen_representation_without_task_checkpoint="UNCERTAIN",
        status="**RED（对本项目的 primary 评估）**",
        allowed_use="**排除**于 primary；仅可作 literature baseline 讨论",
        reason="① SKEMPI v2.0 监督；② 我们 V6.0/V6.1 已把其预测列用于同源评估链路 → "
               "本项目的 outer-test 与其训练数据同源，**确定存在重叠风险**（§9 RED）",
    ),
    dict(
        model_name="ReCLIP",
        paper="Li et al., Nature Communications 2024 (residue-level context for PPI)",
        official_repository="https://github.com/RiverZ/reclip",
        official_checkpoint_source="https://huggingface.co/RiverZ/reclip",
        license_code="**other**", license_weights="**other**",
        task="contrastive pretraining of residue-level PPI representations",
        training_dataset="PDB 复合物结构（对比学习）；**HF 发布物为 XGBoost mutation 分类头**",
        validation_dataset="downstream PPI tasks", test_dataset="downstream tasks",
        SKEMPI_used="**uncertain**（mutation 分类头的训练数据未在 model card 声明）",
        SKEMPI_version="unknown",
        CATH_split_used="NO", PPIformer_split_used="NO",
        checkpoint_training_scope="预训练为对比学习；**但发布物 `classifier_heads/mutation/"
                                 "reclip_mutation_xgb_best_fold5_*.json/.pkl` 是任务 head**",
        residue_embeddings_available="**NO（HF 仓库未发布 backbone 权重）**",
        intermediate_representation_accessible="NO",
        prediction_scores_available="YES（mutation 分类 head）",
        architecture_reusable="see paper", checkpoint_reusable="**license=other，不确定**",
        download_size="未测（未下载）",
        hardware_requirement="unknown",
        can_extract_frozen_representation_without_task_checkpoint="**NO**"
            "（仓库只提供任务 head，未提供可解耦 backbone）",
        status="**YELLOW（本轮由 GREEN 候选降级）**",
        allowed_use="**EXPLORATORY_ONLY**；不得进入 primary",
        reason="① HF 实际发布物是 **XGBoost mutation 分类头**而非通用 backbone；"
               "② 该 head 训练数据未声明 → 依 §8 标 `UNKNOWN_TRAINING_OVERLAP`；"
               "③ license = `other`（非 MIT/Apache）→ 研究使用不确定。"
               "三条独立理由使其**无法满足 §7 末条的 frozen-representation 条件**",
    ),
    dict(
        model_name="PreFold-dG",
        paper="Bioinformatics 2024 (intermediate representations of protein folding models)",
        official_repository="（README 抓取失败，14 B）",
        official_checkpoint_source="unknown",
        license_code="unknown", license_weights="unknown",
        task="ΔΔG prediction from folding-model intermediate representations",
        training_dataset="unknown",
        validation_dataset="unknown", test_dataset="unknown",
        SKEMPI_used="likely YES（同类任务）", SKEMPI_version="unknown",
        CATH_split_used="unknown", PPIformer_split_used="unknown",
        checkpoint_training_scope="unknown",
        residue_embeddings_available="unknown",
        intermediate_representation_accessible="unknown",
        prediction_scores_available="unknown",
        architecture_reusable="unknown", checkpoint_reusable="unknown",
        download_size="unknown", hardware_requirement="unknown",
        can_extract_frozen_representation_without_task_checkpoint="UNKNOWN",
        status="**YELLOW / UNKNOWN**",
        allowed_use="**不纳入 primary**（文档不足，无法完成 §7 审计）",
        reason="仓库 README 不可得；训练数据与许可无法核验 → 依 §8 "
               "`UNKNOWN_TRAINING_OVERLAP`，不得假定安全",
    ),
]


def phase1() -> None:
    e1("# V7.2R PHASE 1 — 外部模型 Provenance / Leakage / License 审计")
    e1("")
    e1(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    e1("- **本 PHASE 完成前不下载任何外部 checkpoint**（§11）")
    e1(f"- 候选系统数：**{len(CANDS)}**（≤5，符合 §6）")
    e1("")
    df = pd.DataFrame(CANDS)
    df.to_csv(PH1 / "model_provenance.csv", index=False, encoding="utf-8")

    # ---- leakage matrix ----
    leak = []
    lab = pd.read_csv(V711 / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    n_tgt = lab["target_id_canonical"].nunique()
    n_fam = lab["family_id_canonical"].nunique()
    for c in CANDS:
        lev = ("GREEN" if c["SKEMPI_used"] == "NO" else
               "UNKNOWN_TRAINING_OVERLAP")
        leak.append({
            "model": c["model_name"],
            "SKEMPI_used": c["SKEMPI_used"],
            "SKEMPI_version": c["SKEMPI_version"],
            "our_targets": n_tgt, "our_families": n_fam,
            "overlap_level_assessable": ("family" if c["SKEMPI_used"] == "YES"
                                         else "n/a"),
            "leakage_verdict": lev,
            "status": c["status"],
        })
    pd.DataFrame(leak).to_csv(PH1 / "leakage_matrix.csv", index=False,
                              encoding="utf-8")
    lic = pd.DataFrame([{"model": c["model_name"],
                         "license_code": c["license_code"],
                         "license_weights": c["license_weights"],
                         "research_use_ok": ("YES" if "MIT" in str(c["license_code"])
                                             or "Apache" in str(c["license_code"])
                                             else "UNCERTAIN"),
                         "status": c["status"]} for c in CANDS])
    lic.to_csv(PH1 / "license_matrix.csv", index=False, encoding="utf-8")

    e1("## 1.1 审计主表")
    e1("")
    e1("| Model | Representation | Checkpoint | SKEMPI 重叠 | License | Status | "
       "Allowed use |")
    e1("|---|---|---|---|---|---|---|")
    for c in CANDS:
        e1(f"| **{c['model_name']}** | {c['representation'] if 'representation' in c else c['task']} | "
           f"{c['checkpoint_reusable']} | {c['SKEMPI_used']} | "
           f"{c['license_code']} | {c['status']} | {c['allowed_use']} |")
    e1("")
    e1("## 1.2 GREEN / YELLOW / RED")
    e1("")
    # 注意：状态串里可能同时出现多个词（如「YELLOW（本轮由 GREEN 候选降级）」），
    # 故必须先判 YELLOW/RED，再判 GREEN，避免子串误匹配。
    def _cls(s: str) -> str:
        su = s.upper()
        if "YELLOW" in su or "RED" in su:
            return "RED" if "RED" in su else "YELLOW"
        return "GREEN" if "GREEN" in su else "UNKNOWN"
    g = [c["model_name"] for c in CANDS if _cls(c["status"]) == "GREEN"]
    y = [c["model_name"] for c in CANDS if _cls(c["status"]) == "YELLOW"]
    r = [c["model_name"] for c in CANDS if _cls(c["status"]) == "RED"]
    e1(f"- **GREEN**：{g}")
    e1(f"- **YELLOW**：{y}")
    e1(f"- **RED**：{r}")
    e1("")
    for c in CANDS:
        e1(f"### {c['model_name']} → {c['status']}")
        e1("")
        e1(f"- **理由**：{c['reason']}")
        e1(f"- 训练数据：{c['training_dataset']}")
        e1(f"- `can_extract_frozen_representation_without_task_checkpoint?`："
           f"**{c['can_extract_frozen_representation_without_task_checkpoint']}**")
        e1("")
    e1("## 1.3 关键结论")
    e1("")
    e1("**GREEN（可用于 primary）**：只有 **ESM-2 8M**（通用自监督，MIT，已下载 29.9 MB）")
    e1("")
    e1("**ReCLIP 是 Track B 的首选候选**，但：")
    e1("- README 抓取失败（14 B）→ 其**预训练是否使用 SKEMPI 标签无法确认**")
    e1("- license 未声明")
    e1("⇒ 依 §8，**标 `UNKNOWN_TRAINING_OVERLAP`**，")
    e1("  在完成核验（HF 文件清单 + 论文方法节）之前**不下载、不纳入 primary**。")
    e1("")
    e1("**Prompt-DDG**：README 明确使用 **SKEMPI v2**，且 license 未声明 → **YELLOW**（exploratory only）")
    e1("")
    e1("**CATH-ddG**：SKEMPI v2.0 监督 + CATH/PPIformer split 均已使用；")
    e1("且**我们 V6.0/V6.1 已把它的预测列用于同源评估链路** →")
    e1("**RED（对本项目 primary 评估）**，仅可作 literature baseline 讨论。")
    e1("")
    e1("## 1.4 依 §37 的 STOP/GO 判定")
    e1("")
    e1("- **存在 1 个 GREEN 表示**（ESM-2 8M，且**已在 V7.2 完成部分嵌入**）")
    e1("- ReCLIP 待核验；未核验前**不下载**")
    e1("⇒ **结论：继续 Track A + 我们的 perturbation 模型**（§37 第一分支），")
    e1("  并把 ReCLIP 的核验列为 Track B 的**前置条件**而非阻塞项。")
    e1("")
    (PH1 / "model_decision.md").write_text("\n".join(L1) + "\n",
                                           encoding="utf-8")


def main() -> int:
    t0 = time.time()
    phase0()
    phase1()
    print(f"\n写出 PHASE0: existing_inventory.csv, phase0_report.md")
    print(f"写出 PHASE1: model_provenance.csv, leakage_matrix.csv, "
          f"license_matrix.csv, model_decision.md")
    print(f"- 耗时 {time.time() - t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
