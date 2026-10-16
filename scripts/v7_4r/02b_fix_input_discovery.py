# -*- coding: utf-8 -*-
"""scripts/v7_4r/02b_fix_input_discovery.py -- cell 04 自动解压 ZIP（Colab 上传的是 zip，不是目录）.

## 问题（用户实测）

```
FileNotFoundError: 未找到 V7.4R 输入包。请把 colab_input_package.zip 解压到 /content 后再跑。
```

用户在 Colab 里上传的是 `colab_input_package.zip`（**文件**），
而 cell 04 只找 **已解压的目录**。这一步不该让用户手动做。

## 修复

cell 04 现在：
1. 先在常见位置找**已解压**的输入目录（原有逻辑）
2. 找不到就**自动搜索并解压** `colab_input_package.zip`
   （搜索 `/content`、`/content/drive/MyDrive`、`/content/sample_data` 等）
3. 解压到 `WORK / "input"`（`/content/v7_4r/input`）
4. 若 ZIP 里含单层目录则自动下钻
5. 仍找不到才报错，并**列出已看到的候选文件**便于诊断
"""

from __future__ import annotations

import ast
import json
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
NB = R / "colab" / "v7_4r" / "V7_4R_OpenMM_Relaxation_COLAB.ipynb"
nb = json.loads(NB.read_text(encoding="utf-8"))
cells = nb["cells"]


def find(marker):
    for i, c in enumerate(cells):
        if c["cell_type"] == "code" and marker in "".join(c["source"]):
            return i
    raise KeyError(marker)


NEW04 = r'''
# ============================================================ 04 Input Discovery (auto-unzip)
import zipfile as _zip

REQ = ["variants.csv", "mutation_sites.csv", "target_chains.csv", "canonical_rows.csv"]

def _ok_dir(d: Path) -> bool:
    try:
        return all((d / r).exists() for r in REQ)
    except Exception:
        return False

def _search_dirs():
    roots = [Path("/content"), Path("/content/drive/MyDrive"),
             Path("/content/sample_data"), Path("/root"), Path.cwd()]
    seen, hits = set(), []
    for root in roots:
        try:
            if not root.exists():
                continue
        except Exception:
            continue
        for f in root.rglob("variants.csv"):
            d = f.parent
            if str(d) in seen:
                continue
            seen.add(str(d))
            if _ok_dir(d):
                hits.append(d)
    return hits

def _autounzip():
    """在常见位置找 colab_input_package.zip（或任何含 variants.csv 的 zip）并解压。"""
    dest = WORK / "input"
    dest.mkdir(parents=True, exist_ok=True)
    cands = []
    for root in (Path("/content"), Path("/content/drive/MyDrive"),
                 Path("/content/sample_data"), Path.cwd()):
        try:
            if not root.exists():
                continue
        except Exception:
            continue
        cands += list(root.glob("*.zip"))
        cands += list(root.glob("**/colab_input_package.zip"))
    cands = sorted(set(cands), key=str)
    print("找到 ZIP 候选:", [str(c) for c in cands] or "无")
    for z in cands:
        try:
            with _zip.ZipFile(z) as zf:
                names = zf.namelist()
                if not any(n.endswith("variants.csv") for n in names):
                    continue
                zf.extractall(dest)
            print("已解压:", z, "->", dest)
            # 处理 ZIP 内的单层目录
            sub = [p for p in dest.iterdir() if p.is_dir()]
            for s in sub:
                if _ok_dir(s):
                    print("下钻单层目录:", s)
                    return s
            if _ok_dir(dest):
                return dest
        except Exception as e:
            print("解压失败", z, type(e).__name__, str(e)[:120])
    return None

PKG = None
hits = _search_dirs()
if hits:
    PKG = sorted(hits, key=str)[0]
    print("发现已解压输入目录:", PKG)
else:
    print("未发现已解压输入目录 -> 尝试自动解压 colab_input_package.zip ...")
    PKG = _autounzip()

if PKG is None:
    dbg = []
    for root in (Path("/content"), Path("/content/drive/MyDrive")):
        try:
            if root.exists():
                dbg += [str(p) for p in list(root.glob("*.zip"))[:10]]
                dbg += [str(p) for p in list(root.glob("*.csv"))[:10]]
        except Exception:
            pass
    print("当前可见的 zip/csv:", dbg or "无")
    raise FileNotFoundError(
        "仍未找到输入包。请确认 colab_input_package.zip 已上传到 /content "
        "（左侧 Files 面板应能看到它）。本 cell 会自动解压，无需手动 unzip。")

import pandas as pd
var = pd.read_csv(PKG / "variants.csv")
sit = pd.read_csv(PKG / "mutation_sites.csv")
tch = pd.read_csv(PKG / "target_chains.csv")
canon = pd.read_csv(PKG / "canonical_rows.csv")
print("选定输入目录:", PKG)
print("variants:", var.shape, "| sites:", sit.shape, "| chains:", tch.shape)
print("targets:", var.target_id.nunique(), "| rows:", len(var))

labtok = ("ddg", "delta_g", "affinity", "kd(", "label", "y_true", "pred")
bad = [c for df in (var, sit, tch) for c in df.columns
       if any(t in c.lower() for t in labtok)]
print("疑似标签列:", bad if bad else "无 ✅")
assert not bad, "输入包含疑似标签列 -> 停止"

sites_by_row = {int(k): g for k, g in sit.groupby("row_id")}
TARGETS = sorted(var.target_id.astype(str).unique())
if MAX_VARIANTS:
    TARGETS = TARGETS[:int(MAX_VARIANTS)]
print("处理 target 数:", len(TARGETS))

THREE = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
         "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
         "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
         "Y": "TYR", "V": "VAL"}

import requests
def fetch_pdb(pdb_id):
    dst = PDB_DIR / f"{pdb_id}.pdb"
    if dst.exists() and dst.stat().st_size > 1000:
        return dst
    r = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=60)
    r.raise_for_status()
    dst.write_bytes(r.content)
    return dst
'''

i04 = find("# ============================================================ 04 Input Discovery")
cells[i04]["source"] = NEW04.strip("\n").splitlines(True)
NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding="utf-8")

nb2 = json.loads(NB.read_text(encoding="utf-8"))
code2 = [c for c in nb2["cells"] if c["cell_type"] == "code"]
bad = []
for i, c in enumerate(code2):
    try:
        ast.parse("".join(c["source"]))
    except SyntaxError as e:
        bad.append((i, e.lineno, e.msg))
t = "\n".join("".join(c["source"]) for c in code2)
print(f"cells={len(nb2['cells'])} code={len(code2)} 语法错误={len(bad)}")
for b in bad:
    print("   ", b)
for k, v in (("自动解压函数", "_autounzip" in t),
             ("搜索多个 root", "/content/drive/MyDrive" in t),
             ("单层目录下钻", "下钻单层目录" in t),
             ("失败时列出可见 zip/csv", "当前可见的 zip/csv" in t),
             ("无需手动 unzip 的提示", "无需手动 unzip" in t)):
    print(f"  {'✅' if v else '❌'} {k}")
