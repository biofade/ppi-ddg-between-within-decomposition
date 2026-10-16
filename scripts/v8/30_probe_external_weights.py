# -*- coding: utf-8 -*-
"""scripts/v8/30_probe_external_weights.py -- 探测外部权重的**大小与可下载性**。

目的：在花时间建环境之前，先确认
  1) Google Drive 的 CATH-ddG 权重到底多大（本机可用 RAM 仅 ~1.1 GB free / 7.5 GB）
  2) RDE-Network 的权重是否可直接下载
  3) 仓库结构（configs / 数据脚本）是否支撑 CPU 推理

**只探测，不下大文件**（超过 THRESHOLD_MB 就只记录不下载）。
"""

from __future__ import annotations

import json
import re
import sys
import urllib.request
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
OUT = R / "experiments" / "v8" / "external"
OUT.mkdir(parents=True, exist_ok=True)
THRESHOLD_MB = 400
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}

CATH_IDS = {
    "CATH_model_0": "1Hs2h_9VBNzFSZ37vDObZ8pu47gI-cwgq",
    "CATH_model_1": "10HEIJM_NU5Cz-vEcqxYNCbYFdvbzZdv7",
    "CATH_model_2": "13tAPMEnSLB9uvIWJpWv4JZoSoZvs3xZ_",
    "PPIFORMER_split": "12XqQ3ucelemtaUB9MSARaKd7g49NZmi1",
    "case_study": "1oOJjA7Lp7xXQqN1ChNScXQ3YS7wfnmvd",
}
RAW = {
    "cath_config_CATH.yml": "https://raw.githubusercontent.com/ak422/CATH-ddG/main/configs/train/CATH.yml",
    "cath_config_PPIformer.yml": "https://raw.githubusercontent.com/ak422/CATH-ddG/main/configs/train/PPIformer.yml",
    "cath_test_script": "https://raw.githubusercontent.com/ak422/CATH-ddG/main/test_DDAffinity.py",
    "rde_readme": "https://raw.githubusercontent.com/luost26/RDE-PPI/main/README.md",
    "rde_test_script": "https://raw.githubusercontent.com/luost26/RDE-PPI/main/test_rde_network_skempi.py",
    "rde_get_skempi": "https://raw.githubusercontent.com/luost26/RDE-PPI/main/data/get_skempi_v2.sh",
    "cath_get_skempi": "https://raw.githubusercontent.com/ak422/CATH-ddG/main/skempi_parallel.py",
}


def head_size(url: str):
    """返回 (status, content_length, content_type, content_disposition, snippet)。"""
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=60) as r:
            cl = r.headers.get("Content-Length")
            ct = r.headers.get("Content-Type")
            cd = r.headers.get("Content-Disposition")
            body = r.read(512)
            return (r.status, cl, ct, cd,
                    body[:120].decode("utf-8", "replace"))
    except Exception as e:  # noqa: BLE001
        return (None, None, None, None, f"ERR {type(e).__name__}: {e}")


print("=" * 72)
print("1) Google Drive: CATH-ddG 权重大小探测")
print("=" * 72)
g = []
for name, fid in CATH_IDS.items():
    url = f"https://drive.google.com/uc?id={fid}&export=download"
    st, cl, ct, cd, sn = head_size(url)
    mb = (int(cl) / 1048576) if (cl and str(cl).isdigit()) else None
    print(f"\n  {name}")
    print(f"    status={st} content-length={cl} "
          f"({mb:.1f} MB)" if mb else f"    status={st} content-length={cl}")
    print(f"    content-type={ct}")
    print(f"    disposition={cd}")
    print(f"    snippet={sn[:100]!r}")
    g.append({"name": name, "file_id": fid, "status": st,
              "content_length": cl, "mb": mb, "content_type": ct,
              "disposition": cd})
print()
print("  判读：")
big = [x for x in g if x["mb"] and x["mb"] > THRESHOLD_MB]
unknown = [x for x in g if not x["mb"]]
print(f"    >{THRESHOLD_MB} MB 的文件数 = {len(big)}")
print(f"    大小未知（Drive 需确认页）的文件数 = {len(unknown)}")
if unknown:
    print("    ⇒ Google Drive 对大文件返回 HTML 确认页，需 gdown 或确认 token。")

print("\n" + "=" * 72)
print("2) 仓库结构（raw 文件）")
print("=" * 72)
raw_out = {}
for name, url in RAW.items():
    st, cl, ct, cd, sn = head_size(url)
    print(f"  {name:<26} status={st} len={cl} ct={ct}")
    raw_out[name] = {"url": url, "status": st, "len": cl}

print("\n" + "=" * 72)
print("3) 关键：CATH-ddG 的 SKEMPI 处理脚本是否可离线跑")
print("=" * 72)
p = Path("F:/dsh/downloads/v8_probe/skempi_parallel.py")
if raw_out.get("cath_get_skempi", {}).get("status") == 200:
    try:
        req = urllib.request.Request(RAW["cath_get_skempi"], headers=UA)
        with urllib.request.urlopen(req, timeout=60) as r:
            txt = r.read().decode("utf-8", "replace")
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(txt, encoding="utf-8")
        print(f"  已保存 {p} ({len(txt)} 字符)")
        for kw in ("skempi", "HER2", "S285", "cath", "CATH", "split", "http"):
            hits = [l.strip() for l in txt.splitlines() if kw in l][:3]
            if hits:
                print(f"    [{kw}] " + " | ".join(h[:90] for h in hits))
    except Exception as e:  # noqa: BLE001
        print("  读取失败:", e)
else:
    print("  未取到")

res = {"google_drive": g, "raw_files": raw_out,
       "threshold_mb": THRESHOLD_MB}
(OUT / "external_weight_probe.json").write_text(
    json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"\n写出 {OUT / 'external_weight_probe.json'}")
