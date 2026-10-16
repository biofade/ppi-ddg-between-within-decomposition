# -*- coding: utf-8 -*-
"""scripts/v8/02_literature_matrix.py -- V8 PHASE 2 / §24/§25：文献基准矩阵 + 新颖性矩阵。

**资源来源是项目内已有的 `papers/*.txt`（13 篇已抽取文本）**，不依赖网络。

诚实规则（§4/§25/§69）：
  * 机器只做**正则候选抽取**，逐条标注 `extraction="MACHINE_CANDIDATE"`，
    并给出命中原文片段；**不把候选当作已核实事实**。
  * 未能抽取的一律 `NOT_FOUND_IN_LOCAL_TEXT`，**不臆造**。
  * 新颖性矩阵按 §25 如实标注 YES/NO/PARTIAL，**不为制造新颖性而软化既有工作**。
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
OUT = R / "experiments" / "v8"
LIT = OUT / "literature"
LIT.mkdir(parents=True, exist_ok=True)
PAPERS = R / "papers"

METHODS = {
    "Pythia-PPI": "Pythia-PPI",
    "SAMAffinity": "SAMAffinity",
    "DSSA-PPI": "DSSA",
    "BindPred": "BindPred",
    "PPAP": "PPAP",
    "AINN-P1": "AINN",
    "CatPred": "CatPred",
    "SELFprot": "SELFprot",
    "PLM": "PLM",
    "TwinPeaks": "Twin Peaks",
    "MultimodalPPI": "Multimodal",
    "EnzymeKinetics": "enzyme kinetic",
    "EvoEF2": None,
}

FEATURES = {
    "random_split": r"random(?:ly)?\s+split|random\s+division",
    "protein_split": r"protein[- ](?:level\s+)?split|protein[-\s]?wise",
    "family_split": r"family[-\s]?(?:level\s+)?split|family[-\s]?wise|family[-\s]?out",
    "CATH_split": r"CATH",
    "homology_control": r"homolog|sequence\s+identity|identity\s+cluster|"
                        r"30%\s*identity|40%\s*identity|sequence\s+similarity",
    "overall_Pearson": r"(?:overall|global|test)\s+Pearson|Pearson\s*\(?r\)?\s*=|"
                       r"Pearson\s+correlation\s+coefficient",
    "overall_Spearman": r"(?:overall|global|test)\s+Spearman|Spearman\s*\(?\s*[ρr]",
    "RMSE": r"\bRMSE\b",
    "MAE": r"\bMAE\b",
    "per_structure_Pearson": r"per[-\s]?structure|per[-\s]?complex|per[-\s]?pdb",
    "per_structure_Spearman": r"per[-\s]?structure\s+Spearman|per[-\s]?complex\s+Spearman",
    "macro_ranking": r"macro[-\s]?(?:average|ranking)|macro\s+Spearman",
    "Kendall": r"Kendall",
    "pairwise_accuracy": r"pairwise\s+accuracy|pair[-\s]?wise\s+classification",
    "HitK": r"Hit@|Hit\s*@|Hits@|top-?K\s+hit",
    "PrecisionK": r"Precision@|Precision\s*@",
    "NDCG": r"NDCG",
    "enrichment": r"enrichment\s+factor|\bEF@|\bEF\b",
    "explicit_calibration_ranking_decomposition":
        r"calibration\s+(?:and|vs\.?)\s+ranking|inter[-\s]?target\s+.{0,25}intra[-\s]?target|"
        r"between[-\s]?target\s+.{0,25}within[-\s]?target|C_B\s+.{0,10}C_W",
    "target_weighting_declared":
        r"target[-\s]?weight|weight(?:ed|ing)\s+by\s+(?:the\s+)?number\s+of\s+mutations|"
        r"equal[-\s]?target|per[-\s]?target\s+weight",
    "target_weighting_sensitivity": r"weight(?:ing)?\s+sensitiv|robust(?:ness)?\s+to\s+weight",
    "matched_universe_comparison":
        r"matched\s+(?:row|universe|set)|common\s+(?:row|subset)|intersection\s+of\s+rows",
    "candidate_set_calibration":
        r"candidate[-\s]?set|set[-\s]?level\s+calibration|library\s+calibration",
    "rank_reversal_analysis":
        r"rank\s+reversal|ranking\s+reversal|model\s+ranking\s+(?:change|disagree)",
    "public_row_predictions":
        r"predictions?\s+(?:are\s+)?(?:available|provided)|supplementary\s+(?:data|table)|"
        r"github\.com/[\w\-]+/[\w\-]+",
    "official_code": r"github\.com/[\w\-]+/[\w\-]+",
    "checkpoint": r"checkpoint|pretrained\s+(?:model|weights)|\.pth\b|\.ckpt\b",
    "processed_data": r"processed\s+data|preprocessed\s+data|data\s+available",
    "training_provenance": r"training\s+(?:set|data|split)|SKEMPI|trained\s+on",
}


def load_texts():
    t = {}
    if not PAPERS.exists():
        return t
    for p in sorted(PAPERS.glob("*.txt")):
        try:
            t[p.name] = p.read_text(encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            t[p.name] = ""
    return t


def snippet(text, rx, width=90):
    m = re.search(rx, text, re.I)
    if not m:
        return ""
    s = max(0, m.start() - width // 3)
    return re.sub(r"\s+", " ", text[s:m.end() + width])[:220]


def main() -> int:
    texts = load_texts()
    print(f"本地论文文本 = {len(texts)} 篇: {sorted(texts)}")
    rows = []
    for mname, token in METHODS.items():
        # 找匹配该方法的文本
        keys = []
        if token:
            for fn, tx in texts.items():
                if re.search(token, fn, re.I) or re.search(token, tx[:4000], re.I):
                    keys.append(fn)
        if not keys:
            keys = list(texts.keys())
        blob = "\n".join(texts[k] for k in keys)
        row = {"method": mname, "local_paper_texts": len(keys),
               "source_files": "; ".join(sorted(keys)[:3]),
               "extraction": "MACHINE_CANDIDATE"}
        for feat, rx in FEATURES.items():
            sn = snippet(blob, rx)
            row[feat] = "YES_CANDIDATE" if sn else "NOT_FOUND_IN_LOCAL_TEXT"
            row[feat + "__evidence"] = sn
        rows.append(row)
    D = pd.DataFrame(rows)
    D.to_csv(LIT / "benchmark_landscape.csv", index=False, encoding="utf-8")

    print(f"\n完成 {len(D)} 个方法的候选抽取")
    core = ["random_split", "family_split", "homology_control", "overall_Pearson",
            "per_structure_Spearman", "macro_ranking", "NDCG",
            "explicit_calibration_ranking_decomposition",
            "target_weighting_sensitivity", "matched_universe_comparison",
            "rank_reversal_analysis", "candidate_set_calibration"]
    print(f"\n{'method':<16}" + "".join(f"{c[:11]:>12}" for c in core))
    for r in D.itertuples():
        line = f"{r.method:<16}"
        for c in core:
            v = getattr(r, c)
            line += f"{('YES' if v.startswith('YES') else '-'):>12}"
        print(line)

    # ---- §25 新颖性矩阵 ----
    NOV_ROWS = ["PreFold-dG", "CATH-ddG", "PPIformer", "HERMES", "Prompt-DDG",
                "RDE-Network", "BA-DDG", "EBM-DDG", "Light-DDG", "Pythia-PPI",
                "H3-DDG", "MARS-DDG", "USP-ddG", "ReCLIP", "CR-PPI-ddG (ours)"]
    NOV_COLS = ["strict task-specific family-OOD", "overall metrics",
                "per-target ranking", "explicit inter/within decomposition",
                "calibration axis", "selection utility axis",
                "target weighting audit", "matched-universe comparison",
                "rank reversal analysis", "candidate-set calibration",
                "benchmark-guided remedy case study"]
    nd = pd.DataFrame(index=NOV_ROWS, columns=NOV_COLS, dtype=object)
    for m in NOV_ROWS:
        for c in NOV_COLS:
            nd.loc[m, c] = "UNKNOWN"
    # 依据我们自己的证据如实填充 ours 行
    ours = {
        "strict task-specific family-OOD": "YES",
        "overall metrics": "YES",
        "per-target ranking": "YES",
        "explicit inter/within decomposition": "YES",
        "calibration axis": "YES",
        "selection utility axis": "YES (但区分度不足, 见 DECISION_ALIGNMENT_REPORT)",
        "target weighting audit": "YES",
        "matched-universe comparison": "PARTIAL (row_id 严格匹配; 外部未审计)",
        "rank reversal analysis": "YES",
        "candidate-set calibration": "NOT_RUN",
        "benchmark-guided remedy case study": "NOT_RUN",
    }
    for k, v in ours.items():
        nd.loc["CR-PPI-ddG (ours)", k] = v
    nd.to_csv(LIT / "novelty_matrix.csv", encoding="utf-8")

    lines = ["# BENCHMARK NOVELTY MATRIX（§25）\n",
             "> **诚实规则**：所有非 ours 行的单元格为 `UNKNOWN` —— ",
             "> 因为**外部仓库与论文的完整审计在本轮未执行**（见 `EXTERNAL_AUDIT.md`）。",
             "> 依 §25，**不为制造新颖性而软化既有工作**；在缺乏核实的情况下填 YES 或 NO 都是错的。\n",
             "> `papers/*.txt`（13 篇）的机器候选抽取见 `benchmark_landscape.csv`，",
             "> 仅为**候选**，逐条附命中原文，需人工核实。\n",
             "| method | " + " | ".join(NOV_COLS) + " |",
             "|---|" + "---|" * len(NOV_COLS)]
    for m in NOV_ROWS:
        lines.append(f"| `{m}` | " + " | ".join(str(nd.loc[m, c]) for c in NOV_COLS) + " |")
    lines += [
        "\n## 已可确认的「新颖性已被占据」项（避免 §4 的禁用声明）\n",
        "以下概念在本项目**自己的历史工作里就已出现**，因此 V8 **不得**声称 first：\n",
        "- per-target / per-structure Spearman（V6.1/V6.2 已用）",
        "- family-OOD 划分（OUTER_FOLDS_V7_2 已用）",
        "- Top-K prioritization（V6.2 / V7.3 Eff@10 已用）",
        "- baseline/mutation 分解（V7.1 `C_B`/`C_W` 已用）",
        "\n## V8 可能仍然独有的组合（**待外部核实后才能定稿**）\n",
        "- 严格 task-specific family-OOD **+** 精确匹配宇宙 **+** 三轴能力分解",
        "  **+** 目标加权敏感性审计 **+** 模型秩反转 **+** decision alignment",
        "  **+** benchmark-guided remedy 案例",
        "\n> ⚠️ 这只是一个**候选组合**。在外部审计完成之前，",
        "> **不得**在论文中据此声明任何优先级。",
    ]
    (OUT / "BENCHMARK_NOVELTY_MATRIX.md").write_text("\n".join(lines) + "\n",
                                                     encoding="utf-8")
    meta = {"n_papers_local": len(texts),
            "extraction": "MACHINE_CANDIDATE (regex + evidence snippet)",
            "novelty_matrix_status": "UNKNOWN for all external methods (not audited)",
            "forbidden_claims": ["first family-OOD ddG benchmark",
                                 "first per-target evaluation",
                                 "first Top-K benchmark",
                                 "all existing benchmarks are invalid",
                                 "Overall Pearson is useless"]}
    (LIT / "literature_meta.json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n写出 benchmark_landscape.csv / novelty_matrix.csv / "
          f"BENCHMARK_NOVELTY_MATRIX.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
