# -*- coding: utf-8 -*-
"""scripts/v6_1/07_search_external_models.py -- V6.1 PHASE 7.2: 寻找更多外部模型.

对每个目标模型，逐一检查 **GitHub 仓库** 与 **Zenodo / Figshare**，
寻找**已跑好的预测结果文件**（本任务禁止跑模型推理，V6.1 全局约束 3）。

判定标签
--------
* `PREDICTION_FOUND` —— 找到逐条预测文件（并记录路径）
* `CODE_ONLY_NO_PRED` —— 仓库存在但只有代码/权重
* `DATA_ONLY_NO_PRED` —— 有数据集/真值但无预测列
* `REPO_NOT_FOUND` —— 仓库不存在
* `SEARCH_BLOCKED` —— 检索被网络/限流阻断（如实记录）

输出
----
* `experiments/v6_1/phase7_toolkit/external_model_search.md`
* `experiments/v6_1/phase7_toolkit/external_model_search.json`
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import time
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
OUT = ROOT / "experiments" / "v6_1" / "phase7_toolkit"
OUT.mkdir(parents=True, exist_ok=True)
UA = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/124.0 Safari/537.36")

TARGETS = {
    "SAMAffinity": {"repo": "iobio-zjut/SAMAffinity",
                    "note": "V6.0 已下载 58 个 test_fold CSV（仅含 ddG 真值）"},
    "Pythia-PPI": {"repo": "Wublab/pythia_ppi", "note": "V6.0 未找到预测文件"},
    "MINT": {"repo": "VarunUllanat/mint", "note": "V6.0 未找到预测文件"},
    "DeepPPAPred": {"repo": None, "note": "需搜索"},
    "MuToN": {"repo": None, "note": "需搜索"},
    "ProBASS": {"repo": None, "note": "需搜索"},
    "CATH-ddG": {"repo": "ak422/CATH-ddG", "note": "V6.0 已获取基线预测列"},
    "DSSA-PPI": {"repo": "Agitw/DSSA-PPI", "note": "V6.0 已获取 8 模型预测列"},
}

PRED_PAT = re.compile(r"(predict|pred[s_]?|result|metric|score|eval)", re.I)
EXCL = re.compile(r"\.(pt|pth|ckpt|bin|safetensors|h5|pkl|tar|gz|zip|npz|npy)$", re.I)
DATA_EXT = (".csv", ".tsv", ".txt", ".json", ".xlsx")

L: list[str] = []
results: dict[str, dict] = {}


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def get(url: str) -> tuple[object | None, str]:
    try:
        r = subprocess.run(["curl.exe", "-sL", "-A", UA, "--max-time", "50", url],
                           capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=80)
        t = (r.stdout or "").strip()
        return (json.loads(t) if t.startswith(("{", "[")) else None), t
    except Exception as exc:  # noqa: BLE001
        return None, f"<error {exc}>"


def gh_search(q: str, n: int = 5) -> list[dict]:
    import urllib.parse
    u = (f"https://api.github.com/search/repositories?q="
         f"{urllib.parse.quote(q)}&per_page={n}")
    j, _ = get(u)
    if not isinstance(j, dict):
        return []
    return [{"full_name": it["full_name"],
             "stars": it.get("stargazers_count"),
             "pushed_at": (it.get("pushed_at") or "")[:10],
             "license": (it.get("license") or {}).get("spdx_id"),
             "size_kb": it.get("size"),
             "desc": (it.get("description") or "")[:100]}
            for it in j.get("items", [])]


def inspect_repo(slug: str) -> dict:
    meta, raw = get(f"https://api.github.com/repos/{slug}")
    if not isinstance(meta, dict) or "full_name" not in meta:
        return {"exists": False, "reason": (str(meta)[:120] if meta else raw[:120])}
    br = meta.get("default_branch", "main")
    tree, _ = get(f"https://api.github.com/repos/{slug}/git/trees/{br}?recursive=1")
    blobs = []
    if isinstance(tree, dict) and "tree" in tree:
        blobs = [x for x in tree["tree"] if x.get("type") == "blob"]
    data_files = [b for b in blobs
                  if str(b.get("path", "")).lower().endswith(DATA_EXT)
                  and not EXCL.search(str(b.get("path", "")))]
    predish = [b for b in data_files if PRED_PAT.search(str(b.get("path", "")))]
    large_data = [b for b in data_files if (b.get("size") or 0) > 200_000]
    return {"exists": True, "license": (meta.get("license") or {}).get("spdx_id"),
            "size_kb": meta.get("size"), "default_branch": br,
            "pushed_at": meta.get("pushed_at"), "n_files": len(blobs),
            "n_data_files": len(data_files),
            "pred_candidates": [{"path": b["path"], "size": b.get("size")}
                                for b in sorted(predish, key=lambda z: z.get("size", 0))][:12],
            "large_data_files": [{"path": b["path"], "size": b.get("size")}
                                 for b in sorted(large_data,
                                                 key=lambda z: -(z.get("size") or 0))][:6]}


def main() -> int:
    emit("# V6.1 PHASE 7.2 — 外部模型预测文件搜索")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit("- 原则: **不跑模型推理**（V6.1 约束 3）；只找已跑好的预测文件；"
         "只下 <2 MB 的小文件")
    emit("- 检索: GitHub Repository Search + 各仓库 git tree；"
         "Zenodo / Figshare 见文末说明")
    emit("")

    for name, info in TARGETS.items():
        emit(f"## {name}")
        emit("")
        rec: dict = {"note": info.get("note", "")}
        # 若没有已知 repo，先搜索
        slug = info.get("repo")
        if not slug:
            cands = gh_search(name)
            rec["search_candidates"] = cands
            emit(f"- GitHub 搜索 `{name}` → {len(cands)} 个候选：")
            for c in cands:
                emit(f"  - `{c['full_name']}` (★{c['stars']}, {c['size_kb']} KB, "
                     f"license={c['license']}, pushed {c['pushed_at']}) — {c['desc']}")
            if cands:
                slug = cands[0]["full_name"]
                emit(f"- 取最相关的一个做深入检查：`{slug}`")
            else:
                rec["status"] = "REPO_NOT_FOUND"
                emit("- **未找到 GitHub 仓库**")
        if slug:
            emit("")
            info2 = inspect_repo(slug)
            rec["repo"] = slug
            rec.update(info2)
            if not info2.get("exists"):
                rec["status"] = "REPO_NOT_FOUND"
                emit(f"- ❌ 仓库不可访问：`{info2.get('reason', '')}`")
            else:
                emit(f"- license: `{info2.get('license')}` | size {info2.get('size_kb')} KB "
                     f"| 文件 {info2.get('n_files')} | 数据文件 {info2.get('n_data_files')} "
                     f"| pushed {info2.get('pushed_at')}")
                if info2["pred_candidates"]:
                    emit("")
                    emit("| 候选预测文件 | size (B) |")
                    emit("|---|---:|")
                    for c in info2["pred_candidates"]:
                        emit(f"| `{c['path']}` | {c['size']} |")
                    rec["status"] = "PREDICTION_FOUND"
                else:
                    emit("- **无候选预测文件**")
                    if info2.get("large_data_files"):
                        emit("- 但有较大的数据文件（可能是真值/数据集，无预测列）：")
                        for c in info2["large_data_files"]:
                            emit(f"  - `{c['path']}` ({c['size']} B)")
                        rec["status"] = "DATA_ONLY_NO_PRED"
                    else:
                        rec["status"] = "CODE_ONLY_NO_PRED"
        emit("")
        results[name] = rec

    # ---------- Zenodo / Figshare ----------
    emit("## Zenodo / Figshare 检索")
    emit("")
    emit("> V6.1 要求同时查 Zenodo / Figshare。以下为检索尝试结果（如实记录）。")
    emit("")
    for api, label in (("https://zenodo.org/api/records?q=", "Zenodo"),
                       ("https://api.figshare.com/v2/articles?search_for=", "Figshare")):
        for term in ("PPI ddG prediction", "SKEMPI prediction results"):
            import urllib.parse
            u = api + urllib.parse.quote(term) + ("&size=3" if "zenodo" in api else "")
            j, raw = get(u)
            if j is None:
                emit(f"- `{label}` 查询 `{term}` → **检索失败/被阻断** "
                     f"（{str(raw)[:80]}）")
                continue
            try:
                items = j.get("hits", {}).get("hits", []) if isinstance(j, dict) else j
                emit(f"- `{label}` 查询 `{term}` → {len(items)} 条结果")
                for it in items[:3]:
                    t = it.get("metadata", {}).get("title") or it.get("title") or "?"
                    emit(f"  - {str(t)[:90]}")
            except Exception as exc:  # noqa: BLE001
                emit(f"- `{label}` 查询 `{term}` → 解析失败 ({type(exc).__name__})")
    emit("")

    # ---------- 汇总 ----------
    emit("## 汇总")
    emit("")
    emit("| 模型 | 状态 | repo | license | 说明 |")
    emit("|---|---|---|---|---|")
    for name, rec in results.items():
        emit(f"| {name} | `{rec.get('status', 'UNKNOWN')}` | "
             f"`{rec.get('repo', '—')}` | {rec.get('license', '—')} | "
             f"{rec.get('note', '')[:60]} |")
    emit("")

    (OUT / "external_model_search.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=1), encoding="utf-8")
    (OUT / "external_model_search.md").write_text("\n".join(L) + "\n",
                                                  encoding="utf-8")
    print(f"\n写出 {OUT / 'external_model_search.md'}")
    print(f"写出 {OUT / 'external_model_search.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
