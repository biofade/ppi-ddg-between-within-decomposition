# -*- coding: utf-8 -*-
"""scripts/v8_1/155_repo_files.py -- 生成仓库的 README / LICENSE / .gitignore / fetch 脚本。

这些文件写入 FINAL_SUBMISSION/repo_staging/，供作者直接推送到 GitHub。
"""
import hashlib
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(".")
STAGE = R / "FINAL_SUBMISSION/repo_staging"
assert STAGE.exists(), "请先运行 153_build_repo_package.py"

README = r"""# Disentangling Between-Target and Within-Target Signal in PPI dG Benchmarking

Reproducibility package for the manuscript *Disentangling Between-Target and Within-Target
Signal in PPI $\Delta\Delta G$ Benchmarking*.

**Author:** Le Yu (Independent Researcher) - ORCID 0009-0000-8807-3973 - biofader@sjtu.edu.cn

---

## Contents

| Path | What it is |
|---|---|
| `latex_v6/` | LaTeX sources + **compiled PDFs** (main text, supplementary, bibliography, build script) |
| `figures_v5/` | All figures as vector PDF + PNG, plus the plotting script that generates them |
| `v8_1/` | **Frozen artifacts: the machine-readable source of every headline number** |
| `experiments/` | Frozen inputs (predictions, split assignments, scoreboard) needed to reproduce |
| `data/v5_0/mutation_sites.csv` | Mutation-site mapping used to build the annotation namespace |
| `scripts/v8_1/` | Analysis, audit, document-generation and validation scripts |
| `scripts/v6_2/` | Four-level controlled split and the exact covariance decomposition |
| `scripts/v7_1/` | Rankers and nested cross-validation |
| `validate_final_manuscript.py` | Automated consistency checker (runs against the compiled PDFs) |
| `VALIDATION_REPORT.md` | Output of the checker |
| `REPO_MANIFEST.csv` | Every file with size and SHA-256 |
| `THIRD_PARTY_INPUTS.csv` | Third-party data: source URL + SHA-256 (not redistributed) |

> The directory layout is deliberate: `validate_final_manuscript.py` resolves every input
> relative to its own location, so a fresh clone can run the full consistency check with no
> configuration.

---

## Reproduction

```bash
pip install numpy pandas scipy scikit-learn lightgbm matplotlib pymupdf

# fetch third-party inputs (not redistributed here)
python fetch_third_party.py

# rebuild the PDFs
cd manuscript/latex_v6
xelatex -interaction=nonstopmode manuscript_v6.tex
bibtex manuscript_v6
python ../../scripts/v8_1/56_reorder_bbl.py
xelatex -interaction=nonstopmode manuscript_v6.tex
xelatex -interaction=nonstopmode manuscript_v6.tex
xelatex -interaction=nonstopmode supplementary_v6.tex
xelatex -interaction=nonstopmode supplementary_v6.tex

# automated consistency checker (exits non-zero on any FAIL)
cd ../..
python validate_final_manuscript.py
```

The checker recomputes every headline number from the frozen artifacts and compares it with the
value printed in the compiled PDF, including a text scan of the rendered figures.

---

## The two results that carry the paper

**1. An exact decomposition.** Pooled Pearson correlation splits exactly into a between-target
and a within-target covariance contribution, `r_pooled = C_B + C_W`. The identity holds for every
individual configuration (maximum absolute residual 4.4e-16). In a controlled comparison that
changes *only* the split regime - identical 7,523 mutation rows, identical model configurations -
pooled Pearson falls from 0.6844 (random split) to 0.1508 (family-out-of-distribution), and
**81.1%** of the summed absolute component change is carried by the between-target component.

**2. The geometric constraint on fixed-K evaluation.** Under the canonical rank-linear reference
the discriminable interval `D(n,K)` peaks at `n* = K + sqrt(K^2 + K)`. For `n <= K` every model
selects all candidates, so fixed-K shortlist metrics carry **zero** information; here 53 of 130
targets with `n >= 5` (40.8%) fall in that regime. See Supplementary Note 4 for the derivation
and Methods for its stated scope.

---

## Third-party data

The five external prediction columns are annotation columns of the SKEMPI-derived table
distributed with CATH-ddG. We do **not** redistribute that table; `fetch_third_party.py`
downloads it and verifies the SHA-256 recorded in `THIRD_PARTY_INPUTS.csv`. All upstream
predictions and datasets remain subject to their original licences.

---

## Known limitations

* **Key-level grouping.** The controlled split unions rows sharing an observation key, so no
  identical mutation key straddles folds; verified empirically for all four strictness levels and
  all three seeds (`reports/RECOMPUTE_FEASIBILITY.md`).
* **External score orientation.** The published external columns document no sign convention. The
  raw analysis is primary; a sign-orientation sensitivity (explicitly label-informed) is reported
  alongside it. The two orientations lead to opposite conclusions about cross-axis agreement, so
  no directional conclusion is drawn from that five-method panel.
* **RMSE** is reported only for outputs on a kcal/mol-like scale; ESM-IF and ProteinMPNN return
  dimensionless mutation-effect scores.
* **External panel power.** Five methods is too few to support a directional claim either way.

---

## Licence

Code: MIT (`LICENSE`). Text and figures: CC-BY-4.0 unless stated otherwise. Third-party data
retains its original licence.
"""

LICENSE = """MIT License

Copyright (c) 2026 Le Yu

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Note: the manuscript text and figures are released under CC-BY-4.0. Third-party
datasets and predictions retain their original licences and are not redistributed
here (see THIRD_PARTY_INPUTS.csv and fetch_third_party.py).
"""

GITIGNORE = """# LaTeX build products
*.aux
*.log
*.out
*.bbl
*.blg
*.synctex.gz
*.fls
*.fdb_latexmk

# Python
__pycache__/
*.py[cod]
.venv/
venv/

# OS
.DS_Store
Thumbs.db

# large regenerable outputs
*.zip
"""

# fetch 脚本：用 THIRD_PARTY_INPUTS.csv 里的 URL 与 sha256 下载并校验
FETCH = '''#!/usr/bin/env python
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
            print(f"  [FAIL] {name}: sha256 mismatch\\n    expected {exp}\\n    got      {got}")
            ok = False
        else:
            print(f"  [OK] {name}  sha256={got[:16]}...")
    print()
    print("  all third-party inputs verified" if ok else
          "  some inputs were not verified - see messages above")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
'''

(STAGE / "README.md").write_text(README, encoding="utf-8")
(STAGE / "LICENSE").write_text(LICENSE, encoding="utf-8")
(STAGE / ".gitignore").write_text(GITIGNORE, encoding="utf-8")
(STAGE / "fetch_third_party.py").write_text(FETCH, encoding="utf-8")
(STAGE / "metadata" / "AUTHOR_METADATA.json").write_text(
    (R / "FINAL_SUBMISSION/AUTHOR_METADATA.json").read_text(encoding="utf-8"),
    encoding="utf-8")

print("  已写入 README.md / LICENSE / .gitignore / fetch_third_party.py / metadata")
print()
print("  仓库顶层:")
for p in sorted(STAGE.iterdir()):
    kind = "dir " if p.is_dir() else "file"
    print(f"    [{kind}] {p.name}")
