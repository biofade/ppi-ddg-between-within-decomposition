# -*- coding: utf-8 -*-
"""scripts/v6/01_scavenge_sota_preds.py -- V6.0 PHASE 1: 搜刮外部 SOTA 预测文件.

原则（V6 §0.1）
---------------
**不跑任何大模型推理**。本脚本只做「数据考古」：
用 GitHub API 列出目标仓库的文件树，寻找**已跑好的预测结果文件**
（`*pred*.csv` / `*test*result*` / `prediction` 目录等），
只下载**小文件**（< 2 MB，远低于 §0.1 的 100 MB 上限）。

输出
----
* `experiments/v6_sota_audit/data_scavenge_raw.json`  （每仓库的文件清单与判定）
* `experiments/v6_sota_audit/data_scavenge_summary.md`（人类可读）
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
OUT = ROOT / "experiments" / "v6_sota_audit"
OUT.mkdir(parents=True, exist_ok=True)
DL = Path(r"F:\dsh\downloads")
UA = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/124.0 Safari/537.36")
MAX_BYTES = 2 * 1024 * 1024          # 只下载小文件

REPOS = {
    "SAMAffinity": "iobio-zjut/SAMAffinity",
    "MINT": "VarunUllanat/mint",
    "DSSA-PPI": "Agitw/DSSA-PPI",
    "Pythia-PPI": "Wublab/pythia_ppi",
}

# 预测文件特征（大小写不敏感）
PRED_PAT = re.compile(
    r"(predict|pred[s_]?|result|output|metric|score|eval|test.*(csv|tsv|txt)|"
    r"skempi.*(csv|tsv)|ddg.*(csv|tsv))", re.I)
# 明确排除的大文件 / 权重
EXCL_PAT = re.compile(r"\.(pt|pth|ckpt|bin|safetensors|h5|pkl|tar|gz|zip|npz|npy|pdb|cif)$",
                      re.I)

L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def api(url: str) -> tuple[object | None, str]:
    """调 GitHub API；返回 (parsed_json, raw_text)。"""
    try:
        r = subprocess.run(["curl.exe", "-sL", "-A", UA, "--max-time", "60", url],
                           capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=90)
        txt = (r.stdout or "").strip()
        return (json.loads(txt) if txt.startswith(("{", "[")) else None), txt
    except Exception as exc:  # noqa: BLE001
        return None, f"<error: {exc}>"


def main() -> int:
    emit("# V6.0 PHASE 1 — SOTA 预测文件搜刮")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- 原则: **不跑任何模型推理**，只找已跑好的预测文件；"
         f"只下载 < {MAX_BYTES // 1024 // 1024} MB 的文件")
    emit("")
    emit("| 仓库 | 存在 | license | size KB | 文件数 | 候选预测文件 |")
    emit("|---|---|---|---:|---:|---:|")

    results: dict[str, dict] = {}
    for name, slug in REPOS.items():
        meta, raw = api(f"https://api.github.com/repos/{slug}")
        rec: dict = {"slug": slug, "url": f"https://github.com/{slug}"}
        if not isinstance(meta, dict) or "full_name" not in meta:
            rec["exists"] = False
            rec["reason"] = f"repo_not_found: {str(meta)[:120] if meta else raw[:120]}"
            emit(f"| {name} | ❌ | — | — | — | — |")
            results[name] = rec
            continue
        rec["exists"] = True
        rec["license"] = (meta.get("license") or {}).get("spdx_id") or "None"
        rec["size_kb"] = meta.get("size")
        rec["default_branch"] = meta.get("default_branch", "main")
        rec["pushed_at"] = meta.get("pushed_at")
        tree, _ = api(f"https://api.github.com/repos/{slug}/git/trees/"
                      f"{rec['default_branch']}?recursive=1")
        blobs: list[dict] = []
        if isinstance(tree, dict) and "tree" in tree:
            blobs = [x for x in tree["tree"] if x.get("type") == "blob"]
        rec["n_files"] = len(blobs)
        cands = []
        for b in blobs:
            p = b.get("path", "")
            if EXCL_PAT.search(p):
                continue
            if PRED_PAT.search(p):
                cands.append({"path": p, "size": b.get("size", 0)})
        cands.sort(key=lambda x: x["size"])
        rec["candidates"] = cands
        emit(f"| {name} | ✅ | {rec['license']} | {rec['size_kb']} | "
             f"{rec['n_files']} | **{len(cands)}** |")
        results[name] = rec

    emit("")
    emit("## 各仓库的候选预测文件（按大小升序，只列前 15）")
    emit("")
    for name, rec in results.items():
        emit(f"### {name} — `{rec['slug']}`")
        if not rec.get("exists"):
            emit("")
            emit(f"- ❌ 仓库不可访问：`{rec.get('reason', '')}`")
            emit("")
            continue
        emit("")
        emit(f"- license: `{rec.get('license')}` | size {rec.get('size_kb')} KB | "
             f"文件 {rec.get('n_files')} | pushed {rec.get('pushed_at')}")
        if not rec["candidates"]:
            emit("- **无候选预测文件** → `PREDICTION_FILE_NOT_FOUND`")
        else:
            emit("")
            emit("| path | size (B) |")
            emit("|---|---:|")
            for c in rec["candidates"][:15]:
                emit(f"| `{c['path']}` | {c['size']} |")
        emit("")

        # 下载小文件（< MAX_BYTES）
        d = DL / re.sub(r"[^a-z0-9]+", "_", name.lower())
        d.mkdir(parents=True, exist_ok=True)
        got = []
        for c in rec["candidates"]:
            if c["size"] > MAX_BYTES:
                continue
            fn = d / Path(c["path"]).name
            if fn.exists() and fn.stat().st_size > 0:
                got.append((c["path"], fn.stat().st_size, "cached"))
                continue
            url = (f"https://raw.githubusercontent.com/{rec['slug']}/"
                   f"{rec.get('default_branch', 'main')}/{c['path']}")
            try:
                r = subprocess.run(["curl.exe", "-sL", "-A", UA, "--max-time", "60",
                                    "-o", str(fn), url],
                                   capture_output=True, text=True, timeout=90)
                sz = fn.stat().st_size if fn.exists() else 0
                if sz > 0:
                    got.append((c["path"], sz, "downloaded"))
                else:
                    got.append((c["path"], 0, "empty_or_blocked"))
            except Exception as exc:  # noqa: BLE001
                got.append((c["path"], 0, f"error:{exc}"))
        rec["downloaded"] = got
        if got:
            emit(f"- 尝试下载 {len(got)} 个小文件到 `{d}`:")
            for p, sz, st in got:
                emit(f"  - `{p}` → {sz} B ({st})")
        emit("")

    (OUT / "data_scavenge_raw.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=1), encoding="utf-8")
    (OUT / "data_scavenge_summary.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\n写出 {OUT / 'data_scavenge_raw.json'}")
    print(f"写出 {OUT / 'data_scavenge_summary.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
