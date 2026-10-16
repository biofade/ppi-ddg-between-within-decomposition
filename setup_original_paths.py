#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""setup_original_paths.py -- 让本仓库在任意 clone 位置可运行，且**不改一行代码**。

背景
----
仓库中的部分脚本（主要在 `scripts/v8_1/`、`scripts/v6_*/`、`scripts/v7_*`）
硬编码了原始工作路径与第三方数据路径。直接 clone 到别处运行会找不到文件。

本脚本在**原路径不存在**时，建立符号链接 / junction 指回你 clone 的位置，
从而无需修改任何脚本即可复现。

需要建立的路径
--------------
    <原项目根>            -> 本仓库根
    <第三方数据目录>      -> 你放置 skempi_v2.csv 等外部数据的目录

默认映射（可用环境变量覆盖）
--------------------------
    PPIDDG_ORIG_ROOT   默认 F:\\DSH_Galatea\\Affinity_Intelligence
    PPIDDG_EXT_DIR     默认 F:\\dsh\\downloads
    本仓库根           = 本脚本所在目录

用法
----
    python setup_original_paths.py            # 建立链接（需要权限）
    python setup_original_paths.py --check    # 只检查当前是否已可运行
    python setup_original_paths.py --remove   # 删除本脚本建立的链接

Windows 说明
-----------
junction 不需要管理员权限（但需要 NTFS）；
符号链接需要管理员权限或开启「开发者模式」。
本脚本在 Windows 上优先用 junction。
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ORIG_ROOT = Path(os.environ.get("PPIDDG_ORIG_ROOT",
                                r"F:\DSH_Galatea\Affinity_Intelligence"))
EXT_DIR = Path(os.environ.get("PPIDDG_EXT_DIR", r"F:\dsh\downloads"))

# 运行完整管线所必需的第三方文件（存在其一即视为数据目录可用）
REQUIRED_EXT = [
    "cathddg/skempi_v2.csv",
]


def is_link(p: Path) -> bool:
    return p.is_symlink() or (hasattr(p, "is_junction") and p.is_junction())


def make_link(link: Path, target: Path) -> tuple[bool, str]:
    if os.name == "nt":
        # junction 不需要管理员权限
        r = subprocess.run(["cmd", "/c", "mklink", "/J", str(link), str(target)],
                           capture_output=True, text=True)
        if r.returncode == 0:
            return True, "junction"
        r2 = subprocess.run(["cmd", "/c", "mklink", "/D", str(link), str(target)],
                            capture_output=True, text=True)
        if r2.returncode == 0:
            return True, "symlink(dir)"
        return False, (r.stderr or r.stdout or "mklink failed").strip()[:160]
    try:
        link.symlink_to(target, target_is_directory=True)
        return True, "symlink"
    except OSError as exc:
        return False, str(exc)[:160]


def remove_link(link: Path) -> tuple[bool, str]:
    if not is_link(link):
        return False, "not a link created by this script"
    try:
        if os.name == "nt":
            r = subprocess.run(["cmd", "/c", "rmdir", str(link)],
                               capture_output=True, text=True)
            if r.returncode != 0:
                link.unlink()
        else:
            link.unlink()
        return True, "removed"
    except Exception as exc:  # noqa: BLE001
        return False, str(exc)[:160]


def check() -> int:
    ok = True
    print("  check:")
    for label, p in (("project root", ORIG_ROOT), ("third-party", EXT_DIR)):
        if p.exists():
            print(f"    [ok]   {label:<14} {p}")
        else:
            print(f"    [MISS] {label:<14} {p}")
            ok = False
    for rel in REQUIRED_EXT:
        p = EXT_DIR / rel
        if p.exists():
            print(f"    [ok]   required file  {p}")
        else:
            print(f"    [MISS] required file  {p}  (run fetch_third_party.py)")
            ok = False
    print()
    print("  -> ready to run" if ok else
          "  -> not ready; run `python setup_original_paths.py` first")
    return 0 if ok else 1


def main() -> int:
    args = set(sys.argv[1:])
    if "--check" in args:
        return check()

    print(f"  clone root : {HERE}")
    print(f"  orig root  : {ORIG_ROOT}")
    print(f"  ext dir    : {EXT_DIR}")
    print()

    pairs = [(ORIG_ROOT, HERE), (EXT_DIR, EXT_DIR)]
    if "--remove" in args:
        for link, _ in pairs[:1]:
            ok, msg = remove_link(link)
            print(f"    {'removed' if ok else 'skip'}: {link}  ({msg})")
        return 0

    rc = 0
    for link, target in pairs:
        if link == target:
            print(f"    [skip] {link} is already the clone root")
            continue
        if link.exists():
            if is_link(link):
                print(f"    [ok]   {link} is already a link")
            else:
                print(f"    [warn] {link} exists and is a real directory; "
                      f"refusing to touch it")
                rc = 1
            continue
        link.parent.mkdir(parents=True, exist_ok=True)
        ok, kind = make_link(link, target if link != EXT_DIR else EXT_DIR)
        print(f"    {'[made]' if ok else '[FAIL]'} {link} -> {target}  ({kind})")
        if not ok:
            rc = 1
    print()
    if rc == 0:
        print("  done. Now run:  python validate_final_manuscript.py")
    else:
        print("  some links could not be created. On Windows, enable Developer Mode")
        print("  or run this script from an elevated shell. Alternatively copy the")
        print("  repository to the original path instead.")
    return rc


if __name__ == "__main__":
    sys.exit(main())
