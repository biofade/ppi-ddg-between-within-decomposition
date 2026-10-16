# -*- coding: utf-8 -*-
"""scripts/v3.1/00_extract_pdfs.py -- V3.1 任务 1.1: PDF 文本提取 (带页码标记).

为什么必须带页码标记
--------------------
任务硬约束 5: **每条从 PDF 提取的信息必须标注页码**。
因此提取时**逐页**写出, 并在每页前插入 `[[page=N]]` 标记。
后续任何引用都必须能追溯到具体页。

提取失败的处理
--------------
* 若某页提取为空 -> 记录;
* 若整篇提取字符数 < 阈值 (默认 500) -> 判定为**扫描版**, 写
  `logs/v3.1/pdf_failures.log` 并**跳过**(不产出 .txt, 也不假装读过);
* 绝不对失败的 PDF 编造内容。

输出
----
* `papers/<name>.txt`   (逐页文本, 含 `[[page=N]]` 标记)
* `logs/v3.1/pdf_failures.log`
* `features/v3.1_pdf_extract_stats.csv`  (每篇的页数/字符数/状态)
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(r"F:\DSH_Galatea\Affinity_Intelligence")
PAPERS = ROOT / "papers"
LOG_DIR = ROOT / "logs" / "v3.1"
LOG_DIR.mkdir(parents=True, exist_ok=True)
FAIL_LOG = LOG_DIR / "pdf_failures.log"
STATS_CSV = ROOT / "features" / "v3.1_pdf_extract_stats.csv"

MIN_CHARS = 500          # 低于此判定为扫描版


def clean(t: str) -> str:
    """轻度清洗: 去掉过多空行, 但**不改动文字内容**。"""
    t = t.replace("\r\n", "\n").replace("\r", "\n")
    t = re.sub(r"[ \t]+\n", "\n", t)
    t = re.sub(r"\n{4,}", "\n\n\n", t)
    return t.strip()


def extract_one(pdf: Path, out_txt: Path) -> dict:
    from pypdf import PdfReader
    rec = {"file": pdf.name, "n_pages": 0, "n_chars": 0,
           "empty_pages": 0, "status": "", "note": ""}
    try:
        reader = PdfReader(str(pdf))
        rec["n_pages"] = len(reader.pages)
        parts: list[str] = []
        for i, page in enumerate(reader.pages, start=1):
            try:
                txt = page.extract_text() or ""
            except Exception as exc:  # noqa: BLE001
                txt = ""
                rec["note"] += f"p{i}:{type(exc).__name__}; "
            txt = clean(txt)
            if not txt:
                rec["empty_pages"] += 1
            parts.append(f"[[page={i}]]\n{txt}")
        body = "\n\n".join(parts)
        rec["n_chars"] = len(body)
        if rec["n_chars"] < MIN_CHARS:
            rec["status"] = "SCANNED_OR_EMPTY"
            rec["note"] += f"chars={rec['n_chars']} < {MIN_CHARS}"
        else:
            rec["status"] = "OK"
            out_txt.write_text(body, encoding="utf-8")
    except Exception as exc:  # noqa: BLE001
        rec["status"] = "FAILED"
        rec["note"] += f"{type(exc).__name__}: {exc}"
    return rec


def main() -> int:
    ap = argparse.ArgumentParser(description="V3.1 PDF 文本提取")
    ap.add_argument("--papers", type=Path, default=PAPERS)
    args = ap.parse_args()

    t0 = time.time()
    if not args.papers.exists():
        print(f"!! 找不到 {args.papers}")
        return 2
    pdfs = sorted(p for p in args.papers.glob("*.pdf"))
    print(f"=== V3.1 任务 1.1: PDF 文本提取 ===")
    print(f"目录: {args.papers}")
    print(f"发现 {len(pdfs)} 个 PDF")
    print()

    rows = []
    fails = []
    for i, pdf in enumerate(pdfs, 1):
        out_txt = args.papers / (pdf.stem + ".txt")
        rec = extract_one(pdf, out_txt)
        rows.append(rec)
        flag = "OK " if rec["status"] == "OK" else rec["status"]
        print(f"  [{i:2d}/{len(pdfs)}] {flag:18s} p={rec['n_pages']:3d} "
              f"chars={rec['n_chars']:8d} empty_p={rec['empty_pages']:3d}  {pdf.name}")
        if rec["status"] != "OK":
            fails.append(rec)

    import csv
    STATS_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(STATS_CSV, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    if fails:
        with open(FAIL_LOG, "w", encoding="utf-8") as fh:
            fh.write("# V3.1 PDF 提取失败/跳过记录\n\n")
            for r in fails:
                fh.write(f"- `{r['file']}`: status={r['status']}, "
                         f"pages={r['n_pages']}, chars={r['n_chars']}, "
                         f"note={r['note']}\n")

    n_ok = sum(1 for r in rows if r["status"] == "OK")
    print()
    print(f"=== 完成: {n_ok}/{len(rows)} 成功, {len(fails)} 失败/跳过 ===")
    if fails:
        print("失败清单 (已写入 logs/v3.1/pdf_failures.log):")
        for r in fails:
            print(f"  - {r['file']}: {r['status']} ({r['note'][:70]})")
    print(f"写出 {STATS_CSV.name}")
    print(f"耗时 {time.time() - t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
