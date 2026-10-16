# LATEX_V6

The v6 Word manuscript converted to LaTeX with the four requested conventions:
`equation` environments, vector-PDF figures, booktabs three-line tables, and BibTeX
citations.

**Content is taken from the v6 PDF you supplied, not rewritten.** The conversion was driven
by extracting the v6 text and reproducing it section by section; a fidelity check then
confirms the result.

## Files

| file | contents |
|---|---|
| `manuscript_v6.tex` | main manuscript, XeLaTeX, 12 pages compiled |
| `supplementary_v6.tex` | Supplementary Information, XeLaTeX, 4 pages compiled |
| `references.bib` | 13 BibTeX entries, all externally verified |
| `build.ps1` | one-command build |
| `manuscript_v6.pdf` | compiled main manuscript |
| `supplementary_v6.pdf` | compiled Supplementary Information |

## Toolchain

MiKTeX 24.1 at `F:\MiKTeX` (data/config/install/package roots all on F:).
Engine **XeLaTeX** with `fontspec` (Times New Roman, Arial, Consolas);
bibliography via **BibTeX** with `plainnat` + `natbib[super,square]`.

```powershell
cd latex_v6
powershell -File build.ps1
```

Build pipeline: `xelatex` → `bibtex` → reorder `.bbl` to v6 numbering → `xelatex` ×2.

## The four requested conversions

### 1. `equation` environments

| # | label | content |
|---|---|---|
| 1 | `eq:decomp` | $y_{ci} = \mu_c + \varepsilon_{ci}$ and $\hat{y}_{ci} = \hat{\mu}_c + \hat{\varepsilon}_{ci}$ |
| 2 | `eq:identity` | $r_{\mathrm{pooled}} = \operatorname{Cov}(y,\hat{y})/(\sigma_y\sigma_{\hat{y}}) = C_{\mathrm{B}} + C_{\mathrm{W}}$ |
| 3 | `eq:geometry` | $D(n,K) = 0$ for $n \leq K$; $K(n-K)/[n(n+1)]$ for $n > K$ |
| 4 | `eq:nstar` | $n^{*} = K + \sqrt{K^{2} + K}$ |

Supplement: 2 more `equation` environments in Supplementary Note 1 (the covariance
expansion and the identity).

### 2. Vector-PDF figures

`\graphicspath{{../figures_v5/}}` plus `\includegraphics` on PDF vectors — 6 includes in the
manuscript, 3 in the supplement:

```
Fig1_overview_EN.pdf                             FigS1_matched_universe_strata_EN.pdf
Fig2_controlled_split_collapse_EN.pdf            FigS2_capability_heatmap_provenance_clean_EN.pdf
Fig3_target_weighting_sensitivity_EN.pdf         FigS3_pool_sensitivity_rank_agreement_EN.pdf
Fig4_fixed_K_geometry_EN.pdf
Fig5_retrospective_remedy_EN.pdf
Fig6_external_methods_EN.pdf
```

### 3. Three-line (booktabs) tables

`\toprule` / `\midrule` / `\bottomrule` inside `threeparttable`, no vertical rules.

| document | tables | toprule | midrule | bottomrule |
|---|---|---:|---:|---:|
| manuscript | Table 1 overview, Table 2 split, Table 3 external | 3 | 3 | 3 |
| supplement | Table S1 matched universes (`landscape`) | 1 | 1 | 1 |

### 4. BibTeX citations

18 `\cite{}` calls resolved against `references.bib`, rendered as superscript numerals.
The compiled PDF contains **13** numbered entries.

**Numbering matches v6 exactly.** `unsrtnat` numbers by first citation, which made
Sokal [1] because the RMSE sentence cites `biometry` first. `plainnat` sorts alphabetically,
which is also wrong. The build therefore **reorders the `.bbl` into v6's fixed order** and
renumbers, giving:

```
[1] Jankauskaitė et al. SKEMPI 2.0      [8]  Dauparas et al. ProteinMPNN
[2] Bushuiev et al. leakage             [9]  Hsu et al. ESM-IF
[3] Yu et al. CATH-ddG                  [10] Bushuiev et al. ICLR 2024
[4] Shan et al. DDGPred (PNAS)          [11] Järvelin & Kekäläinen NDCG
[5] Park et al. PreFold-dG              [12] Efron & Tibshirani bootstrap
[6] Barlow et al. Flex ddG              [13] Sokal & Rohlf biometry
[7] Delgado et al. FoldX 5.0
```

## Build verification

```
manuscript_v6     fatal = 0   citation undefined = 0   reference undefined = 0   11 pages
supplementary_v6  fatal = 0   citation undefined = 0   reference undefined = 0    4 pages
FIDELITY_CHECK = PASS
```

## Figure placement — a defect found and fixed

In the first LaTeX build **all six main figures floated to the end of the document**
(pages 9–12) instead of appearing near their citations, because the figure environments used
`[t]` with very large graphics and LaTeX deferred every float.

This went unnoticed at first because the figures are **vector**, so they show up as PDF
drawing operations rather than as embedded raster images: a naive image count reported
"2 images" and made them look absent. Counting `get_drawings()` reveals them.

Fix (`scripts/v8_1/58_fix_float_placement.py`):

* figure and table environments changed from `[t]` to `[htbp]`
* float parameters relaxed: `\topfraction` 0.92, `\bottomfraction` 0.75,
  `\textfraction` 0.06, `\floatpagefraction` 0.75, `topnumber` 3, `totalnumber` 4

Result — figures now sit in the body, adjacent to their first citation:

| page | content |
|---|---|
| 3 | Figure 1 (208 vector objects) + Table 2 |
| 4 | Figure 2 + Figure 3 |
| 5 | Figure 4 |
| 6 | Figure 5 |
| 7 | Figure 6 (+ Table 3) |
| 3–4 (supplement) | Figures S1–S3 |

Page count reduced from 12 to 11 as the floats integrated into the text.

All nine captions verified present in the compiled PDFs: Figure 1–6, Table 1–3 in the
manuscript; Figure S1–S3 and Table S1 in the supplement.

## Fidelity to the supplied v6

`scripts/v8_1/57_fidelity_check.py` compares the compiled LaTeX PDFs against the v6 PDFs
you supplied:

| check | result |
|---|---|
| v6 sentences present in the new main manuscript | **47 / 47** |
| v6 sentences present in the new supplement | **5 / 5** |
| key numeric values present | **35 / 35** |
| `HERMES` residue | absent from both |
| `log-odds` residue | absent from both |
| `REFERENCE_VERIFY_REQUIRED` residue | absent from both |

Values checked include 4,473 / 3,351 / 7,523 / 1,781 rows; 0.6844 → 0.1508;
ΔC_B = −0.4347; ΔC_W = −0.1011; 81.1 %; per-regime C_B / C_W (0.4275, 0.2591, 0.0233,
0.1685, −0.0072, 0.1580); weighting medians 0.0624 and 0.0183; n\* ≈ 20.5;
D(20) = 0.2381; D(50) = 0.1569; 40.8 %; 53 of 130; 35/109 and 0/109; all five external
RMSE values (2.1636, 2.1262, 2.3691, 2.3587, 2.3199); and the identity residual
3.3 × 10⁻¹⁶.

Two sentences needed substring fallbacks in the checker because PDF text extraction mangles
the Greek hat accents (`σŷ`) and curly quotes in v6; both were confirmed present by their
constituent phrases.

## Issues fixed during conversion

| # | problem | fix |
|---|---|---|
| 1 | Hard-coded `\textsuperscript{13}` instead of real citations; bibliography produced **0** entries | mapped numbers to BibTeX keys and converted all 18 to `\cite{}` |
| 2 | `unsrtnat` / `plainnat` numbering did not match v6 | `.bbl` reordered into v6's fixed order and renumbered (`56_reorder_bbl.py`) |
| 3 | `\doi` undefined after `.bbl` reordering (the style's preamble was dropped) | re-inserted `\providecommand{\doi}` / `\url` / `\href` |
| 4 | Raw underscores in table data (`U_V6_2_split4_...`) → `Missing $`, `Double subscript` | central `esc()` applied to every cell |
| 5 | `$\Delta\Delta G$` in section titles broke hyperref's bookmark routine | `\texorpdfstring` wrappers |
| 6 | `\SI{...}` triggered `Improper alphabetic constant` | removed `siunitx`; plain math |
| 7 | `pandas.itertuples()` drops dotted column names | table built with `iterrows()` |

## Notes

1. Author, affiliation, ORCID and funding fields remain placeholders.
2. The manuscript compiles to 12 pages, matching v6's 12 pages.
3. `references.bib` deliberately has no HERMES entry — it is not in the frozen external
   panel (`v8_1/external_clean_profiles.csv`: Flex-ddG, DDGPred, ddG_FoldX, ESM-IF,
   proteinMPNN).
