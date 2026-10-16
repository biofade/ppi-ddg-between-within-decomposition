# Disentangling Between-Target and Within-Target Signal in PPI dG Benchmarking

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
| `PATH_SETUP.md` / `setup_original_paths.py` | One-time path setup, needed by **some** analysis scripts only |

> The directory layout is deliberate: `validate_final_manuscript.py` resolves every input
> relative to its own location, so a fresh clone can run the full consistency check with no
> configuration.
>
> **First-time setup.** The validator and the LaTeX build need no configuration. A subset of the
> historical analysis scripts (about 48 files) hard-codes the original working paths; to re-run
> those, see `PATH_SETUP.md` — one command, no code changes.

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
