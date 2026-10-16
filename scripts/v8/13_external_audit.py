# -*- coding: utf-8 -*-
"""scripts/v8/13_external_audit.py -- V8 PHASE 13 / §14/§49/§72：外部方法资源审计。

**reuse-first**：先查项目内是否已有外部方法的 row-level 预测，再考虑推理。

结论形式：每个方法给出 evidence_tier + provenance_status + GPU_REQUIRED。
依 §75，在 PHASE 13 结束时打印 `GPU_REQUIRED = YES/NO`。

**本轮先说清能力边界**：本地无法可靠访问外部仓库（web 检索不稳定），
因此**不对任何外部方法宣称已审计**；只做**项目内资产盘点**（这是可验证的部分），
其余一律标 `NOT_AUDITED`，并给出原因。**不臆造外部资源的存在**。
"""

from __future__ import annotations

import glob
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
EXT = OUT / "external"
RA = OUT / "resource_audit"
EXT.mkdir(parents=True, exist_ok=True)
RA.mkdir(parents=True, exist_ok=True)

# §14 审计清单
METHODS = {
    "FoldX": ("explicit physics/energy", "external_engine"),
    "Flex-ddG": ("explicit physics/energy", "external_engine"),
    "Rosetta": ("explicit physics/energy", "external_engine"),
    "DDGPred": ("mixed/statistical", "external_model"),
    "ProteinMPNN": ("structure-conditioned / inverse-folding", "external_model"),
    "ESM-IF": ("structure-conditioned / inverse-folding", "external_model"),
    "RDE-Network": ("deep ddG", "external_model"),
    "PPIformer": ("deep ddG", "external_model"),
    "HERMES": ("deep ddG", "external_model"),
    "Prompt-DDG": ("deep ddG", "external_model"),
    "BA-DDG": ("deep ddG", "external_model"),
    "EBM-DDG": ("deep ddG", "external_model"),
    "Light-DDG": ("deep ddG", "external_model"),
    "Pythia-PPI": ("deep ddG", "external_model"),
    "CATH-ddG": ("deep ddG", "external_model"),
    "PreFold-dG": ("deep ddG", "external_model"),
    "H3-DDG": ("deep ddG", "external_model"),
    "MARS-DDG": ("deep ddG", "external_model"),
    "USP-ddG": ("deep ddG", "external_model"),
    "ReCLIP": ("deep ddG", "external_model"),
}

# ---- 项目内搜索外部方法的痕迹（列名/文件/目录/CATH 预测）----
print("=== 项目内外部方法痕迹扫描 ===")
PATTERNS = {m: re.compile(m.replace("-", "[-_ ]?"), re.I) for m in METHODS}
PATTERNS["CATH"] = re.compile(r"cath", re.I)
hits = {m: [] for m in list(METHODS) + ["CATH"]}

# 目录/文件名
for p in R.rglob("*"):
    if not p.is_file():
        continue
    if any(seg in p.parts for seg in ("data/pdb", ".git", "node_modules",
                                      "dsh-bio-galatea")):
        continue
    name = p.name
    if p.suffix.lower() not in (".csv", ".json", ".npz", ".parquet", ".py",
                                ".md", ".txt", ".jsonl"):
        continue
    for m, rx in PATTERNS.items():
        if rx.search(name):
            hits[m].append(str(p.relative_to(R)).replace("\\", "/"))

# CATH 专门：项目历史提到 cathddg
cath_extra = []
for pat in ("**/cathddg*", "**/*cath*ddg*", "**/cath*"):
    cath_extra += [str(Path(x).relative_to(R)).replace("\\", "/")
                   for x in glob.glob(str(R / pat), recursive=True)]

rows = []
for m, (cat, kind) in METHODS.items():
    h = sorted(set(hits.get(m, [])))
    # CATH-ddG 额外并入 cath 痕迹
    if m == "CATH-ddG":
        h = sorted(set(h + cath_extra))
    row = {
        "method": m, "broad_category": cat, "kind": kind,
        "in_project_hits": len(h),
        "example_hits": "; ".join(h[:3]),
        "official_repository": "NOT_AUDITED",
        "paper": "NOT_AUDITED",
        "official_predictions": "NOT_AUDITED",
        "processed_data": "NOT_AUDITED",
        "official_split": "NOT_AUDITED",
        "checkpoint": "NOT_AUDITED",
        "inference_script": "NOT_AUDITED",
        "license": "NOT_AUDITED",
        "training_provenance": "NOT_AUDITED",
        "cpu_feasible": "UNKNOWN",
        "gpu_required": "UNKNOWN",
        "matched_row_mapping_possible": "UNKNOWN",
        "evidence_tier": "NOT_AUDITED",
        "provenance_status": "UNKNOWN",
        "reason": ("项目内未发现该方法的 row-level 预测；"
                   "外部仓库审计在本轮未执行（见 audit_limitation）"),
    }
    if h:
        # 有痕迹 -> 至少能进 PHASE 4 的符号/映射审计候选
        row["evidence_tier"] = "PENDING_ROWLEVEL_CHECK"
        row["provenance_status"] = "YELLOW"
        row["reason"] = "项目内存在同名痕迹；需人工确认是否为该方法的 row-level 预测"
    rows.append(row)

D = pd.DataFrame(rows)
D.to_csv(EXT / "external_audit.csv", index=False, encoding="utf-8")

print(f"{'method':<14}{'category':<42}{'in-project hits'}")
for r in D.itertuples():
    print(f"{r.method:<14}{r.broad_category:<42}{r.in_project_hits:>4}  "
          f"{r.example_hits[:60]}")

lim = {
    "audit_limitation": (
        "本轮的**外部仓库**审计（官方 row-level 预测 / 处理数据 / split 定义 / "
        "checkpoint / 推理脚本）**未执行**。原因是本地 web 检索在本会话中不稳定，"
        "无法可靠核实官方资源的真实存在与许可条款。"),
    "what_was_done": (
        "对**项目内**全部文件做了 20 个外部方法名的正则扫描（含 CATH 变体），"
        "以判断是否有可直接复用的 row-level 预测。"),
    "honesty_rule": (
        "**不臆造外部资源的存在**。未审计一律记 NOT_AUDITED，而不是假定存在或不存在。"),
    "impact_on_V8": (
        "V8 的 H1/H2/H4 结论完全基于**内部**模型，不依赖外部方法；"
        "H3 已判 FAIL；因此外部审计缺失**不阻塞**当前结论，"
        "但**阻塞**「外部模型是否复现能力分化」这一条（§68 STRONGLY_SUPPORTED 的条件之一）。"),
    "gpu_decision": (
        "GPU_REQUIRED = NO。理由：§51 要求 GPU 仅在某模型**官方 row-level 预测不可得**、"
        "**官方 checkpoint 与推理代码可得**、且 **CPU 不可行** 时才启用。"
        "本轮连第一步（官方预测是否可得）都未核实，因此**不具备申请 GPU 的依据**。"),
    "next_action": (
        "PHASE 13 应在具备稳定网络时重做：优先 PPIformer / CATH-ddG / PreFold-dG / "
        "HERMES / Prompt-DDG / RDE-Network 的官方预测文件。"),
}
(RA / "audit_limitation.json").write_text(json.dumps(lim, indent=2, ensure_ascii=False),
                                          encoding="utf-8")

lines = ["# V8 External Model Audit（§14/§49/§72）\n",
         "## 状态：**NO_EXTERNAL_MODEL_AUDITED**\n",
         f"> {lim['audit_limitation']}\n",
         f"> {lim['honesty_rule']}\n",
         "## 项目内痕迹扫描（20 个方法名 + CATH 变体）\n",
         "| method | 宽类别 | 项目内命中 | 示例 |",
         "|---|---|---:|---|"]
for r in D.itertuples():
    lines.append(f"| `{r.method}` | {r.broad_category} | {r.in_project_hits} | "
                 f"{r.example_hits[:70]} |")
lines += ["\n## GPU 判定（§75）\n", "```", "GPU_REQUIRED = NO", "```\n",
          f"> {lim['gpu_decision']}\n",
          "## 对 V8 结论的影响\n", f"{lim['impact_on_V8']}\n",
          "## 下一步\n", f"{lim['next_action']}\n"]
(OUT / "EXTERNAL_AUDIT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"\n写出 external_audit.csv / EXTERNAL_AUDIT.md / audit_limitation.json")
print("GPU_REQUIRED = NO")
