#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Download third-party inputs and verify their SHA-256.

The SKEMPI-derived annotation table distributed with CATH-ddG is NOT redistributed
in this repository. Fill in the exact URL below (or edit THIRD_PARTY_INPUTS.csv),
then run this script; it verifies the SHA-256 before accepting the file.
"""
import hashlib
import sys
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
TARGET_DIR = HERE / "third_party"


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    TARGET_DIR.mkdir(exist_ok=True)
    man = pd.read_csv(HERE / "THIRD_PARTY_INPUTS.csv")
    ok = True
    for r in man.itertuples():
        name = Path(str(r.path)).name
        dest = TARGET_DIR / name
        url = str(r.source_url)
        if "<" in url:
            print(f"  [SKIP] {name}: source URL is a placeholder; edit "
                  f"THIRD_PARTY_INPUTS.csv with the real URL")
            ok = False
            continue
        print(f"  downloading {name} ...")
        try:
            urllib.request.urlretrieve(url, dest)
        except Exception as exc:
            print(f"  [FAIL] {name}: {exc}")
            ok = False
            continue
        got = sha256(dest)
        exp = str(r.sha256)
        if exp and exp != "nan" and got != exp:
            print(f"  [FAIL] {name}: sha256 mismatch\n    expected {exp}\n    got      {got}")
            ok = False
        else:
            print(f"  [OK] {name}  sha256={got[:16]}...")
    print()
    print("  all third-party inputs verified" if ok else
          "  some inputs were not verified - see messages above")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
