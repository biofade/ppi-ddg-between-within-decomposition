# VALIDATION_REPORT

Automated consistency check of the final manuscript and supplementary PDF.

```
VALIDATION:  PASS = 67    WARNING = 1    FAIL = 0
```

This report is regenerated on every run of `validate_final_manuscript.py`. Every check
either recomputes a headline number from the frozen artifacts in `v8_1/`, or reads the
**compiled** PDFs in `latex_v6/` (so figure-internal text and table values are checked as
printed, not as intended).

## Open items

- **WARN** `R. Data/code availability truthful` - Manuscript states frozen hashes + machine-readable source + consistency checker are provided with the Supplementary Information, and that a public repository will be linked on acceptance. Repository URL not yet supplied (author action pending), so the availability statement currently relies on the SI package.

## All checks

| # | Check | Status | Detail |
|---:|---|---|---|
| 1 | A. per-configuration r=C_B+C_W identity | PASS | V6_2_split4 n=36, max\|residual\|=4.441e-16 (tol 1e-12) |
| 2 | B. Table 2 vs frozen source | PASS | all 4 regimes x 5 metrics reproduced from all_models_by_rowset.csv |
| 3 | C. Table 2 in .tex == machine source | PASS | parsed 4 rows from Table 2 |
| 4 | D. Figure S2 / 1D vs ablation text | PASS | F1/F2/F4/F8 four-metric values present in FigS2 vector text and ablation paragraph |
| 5 | E. Figure 4 geometry values | PASS | D(20)=0.2381 D(50)=0.1569 n*=20.5 n<=K=40.8% (53/130) |
| 6 | F. Table 3 external values (sign-oriented) | PASS | parsed 5 rows; all equal mu_c x raw from external_clean_profiles.csv |
| 7 | F2. external coverage claim | PASS | coverage equal=False (counts {'Flex-ddG': 1755, 'DDGPred': 1756, 'ddG_FoldX': 1776, 'ESM-IF': 1776, 'proteinMPNN': 1776}); 'same mutations' claim present=False |
| 8 | G. rank agreement recomputed | PASS | rho=0.9000 exact perm p=0.0833 (n=5, 120 permutations); text contains 0.90=True |
| 9 | G2. Flex-ddG/FoldX rank-swap text | PASS | text states a single adjacent swap; old 'no rank reversal' claim absent=True |
| 10 | H. Table S1 vs Figure S1 | PASS | 12 strata; all row counts and hashes found in FigS1 vector text |
| 11 | I. fingerprints / row-id sets | PASS | mas=7523 rows, canonical=4473 rows, subset relation holds |
| 12 | J. LaTeX build cleanly | PASS | [source=log] main fatal=0 cite_undef=0 ref_undef=0 pages=15 \| supp fatal=0 pages=12 |
| 13 | K. 4,473 vs 7,523 universe explained | PASS | canonical subset of mas = True; mas\canon=3050 rows; Methods section + Supplementary Note 2 present |
| 14 | L. median-vs-exact wording | PASS | Table 2 note + Figure 2 caption state marginal medians and per-configuration identity |
| 15 | M. 81.1% interpretation | PASS | marginal-median attribution=0.8113; paired=0.8104; 95% CI=[0.7097198883539126, 0.8448723244877555] |
| 16 | N. Supplementary figure/table numbering | PASS | S-prefix counters set = True; pool table label = tab:s4pools |
| 17 | O. D(n,K) derived with E_max/E_rand | PASS | Supplementary Note 4 defines E_max, E_rand, derives D(n,K) and n* |
| 18 | P. LA/LC/* defined + pool table | PASS | Figure S3 caption defines LA/LC; Table S4 lists pool composition |
| 19 | Q. weighted decomposition derived | PASS | Supplementary Note 5 states a_c schemes and proves r_w = C_B,w + C_W,w |
| 20 | R. Data/code availability truthful | WARN | Manuscript states frozen hashes + machine-readable source + consistency checker are provided with the Supplementary Information, and that a public repository will be linked on acceptance. Repository URL not yet supplied  |
| 21 | S. author metadata filled | PASS | real author name, affiliation, corresponding email, ORCID and funding statement are present in the manuscript; no placeholder remains |
| 22 | T. no fabricated values | PASS | Every headline number traced to v8_1 artifacts; Figure 3 histogram regenerated from the real 109 configuration values (previous hardcoded bins did not match the data and were replaced) |
| 23 | U. weighted identity r_w = C_B,w + C_W,w | PASS | 109 configs x 5 schemes; max\|residual\| = 5.551e-17; by scheme {'ROW': 5.551115123125783e-17, 'EQUAL': 5.551115123125783e-17, 'SQRT': 5.551115123125783e-17, 'LOG': 5.551115123125783e-17, 'CAP30': 5.551115123125783e-17} |
| 24 | V. Note 5 self-consistency (C_W not scheme-invariant) | PASS | forbidden 'identical across all five schemes' present = False; C_W cross-scheme range median = 0.0624 |
| 25 | W. weighting audit values use final (T) decomposition | PASS | primary=T_target_total_mass w=a_c / n_c; range_C_B median=0.3052 range_C_W median=0.1240; text has 40/109 and 34/109; stale 35/109=False 0.0183=False |
| 26 | X. D(n,K) states rank-linear utility reference | PASS | rank-linear reference declared = True; removed universal-utility claim = True |
| 27 | Y. duplicate-key sensitivity analysis | PASS | Note 6 present = True; decomposition verified to 3.3e-16; attribution 81.0-85.1% (direction preserved = True) |
| 28 | Z. supplementary caption numbering | PASS | doubled prefixes found = 0; S-counter set = True; Table S4 caption clean = True |
| 29 | AA. Figure 6 caption wording | PASS | stale 'matched unambiguous mapping' wording present = False |
| 30 | AB. Abstract novelty and attribution wording | PASS | formulate-wording present = True; 'summed absolute component change' in Abstract/S = True |
| 31 | AC. final PDF has no stale Figure 3 strings | PASS | stale strings found in compiled PDF = none |
| 32 | AD. final PDF Figure 3 shows updated values | PASS | expected in PDF: ['0.3052', '0.1240', '34/109']; missing = none |
| 33 | AE. final PDF Figure 3 C_W flip count matches audit | PASS | C_W flip count 5/109 present in PDF = True |
| 34 | AF. text and Figure 3 consistent on C_B | PASS | range_C_B median 0.3052 appears in text = True |
| 35 | AG. D(n,K) qualified in Abstract and Results | PASS | abstract qualifier = True; results qualifier = True |
| 36 | AH. supplementary table cross-references | PASS | Figure S3 -> Table S2 = True; stale 'Table~S4 below' gone = True |
| 37 | AI. final PDF embeds current figure sources | PASS | 9 figures compared against figures_v5 vector sources using the compiled PDF text; failures = none |
| 38 | AJ. final PDF Figure 6 labels updated | PASS | stale labels = none; required labels present = 3/3 |
| 39 | AK. in-figure attribution wording consistent | PASS | stale 'absolute component decline' in PDF = False; 'summed absolute' present = True |
| 40 | AL. external common-intersection sensitivity performed | PASS | Note 7 present = True; intersection rows = 1752; pooled ranking identical = True; macro ranking identical = True |
| 41 | AM. no unverified coverage claim, no 'de-replicated' | PASS | stale claim present = False; 'de-replicated' present = False |
| 42 | AN. Discussion does not generalise the surrogate | PASS | over-general sentence present = False; qualified wording present = True |
| 43 | AO. no RMSE for dimensionless external scores | PASS | dimensionless scores = ['ESM-IF', 'proteinMPNN']; stale RMSE still in PDF = none |
| 44 | AP. external provenance table and Note 8 present | PASS | Table S3 (provenance) = True; Note 8 = True |
| 45 | AQ. units caveat stated in main text | PASS | dimensionless caveat present in manuscript = True |
| 46 | AR. Table 3 marks undefined RMSE | PASS | dash marker present in manuscript PDF = True |
| 47 | AS. no Markdown emphasis residue in final PDFs | PASS | bold residue = none; italic residue = none |
| 48 | AT. DeltaDeltaG sign convention defined in Methods | PASS | sign convention stated = True; 'no sign inversion' present = True |
| 49 | AU. 'positive affine rescaling' wording | PASS | positive-affine wording = True; unqualified wording gone = True |
| 50 | AV. no absolute-value DDG notation for signed predictors | PASS | absolute-value notation present = False |
| 51 | AW. external source stated precisely | PASS | CATH-ddG named = True; 'official per-mutation' removed = True; skempi_v2.csv named in supplement = True |
| 52 | AX. Table S3 caption not self-duplicated | PASS | self-duplicating caption present = False |
| 53 | AY. ESM-IF sign-oriented in main table | PASS | oriented ESM-IF values (+0.2250, -0.0213) in manuscript PDF = True |
| 54 | AZ. rank agreement reported with orientation stated | PASS | raw rho 0.90 present=True; oriented rho 0.30 present=True; p=0.683 present=True |
| 55 | BB. no RMSE for dimensionless scores in supplement | PASS | dimensionless rows with an RMSE value = none |
| 56 | BC. provenance table records raw orientation and multiplier | PASS | orientation column present = True |
| 57 | BD. no incorrect 'same targets' claim | PASS | stale phrase present = False |
| 58 | BE. orientation is declared label-informed, not label-free | PASS | stale 'uses no evaluation labels' = False; label-informed declaration = True |
| 59 | BF. both raw and oriented results reported | PASS | raw rho present = True; oriented rho present = True |
| 60 | BG. no contradictory 'no sign flip' claim | PASS | main=False supplement=False |
| 61 | BH. supplementary table reference is S3 not S4 | PASS | Table S3 present = True; 'Table S4' in supplement = False |
| 62 | BI. Figure 6 Panel B title says sign-oriented | PASS | sign-oriented title = True; stale raw title = False |
| 63 | BJ. energy-type RMSE retained in common-intersection table | PASS | Flex/FoldX/DDGPred RMSE present = True |
| 64 | BK. orientation multiplier uses o_m, not mu_c | PASS | o_m present in both = True; mu_c reused for orientation = False |
| 65 | BL. Note 7 updated to oriented analysis | PASS | oriented wording in Note 7 = True |
| 66 | BM. seq40 is no longer identical to complex | PASS | complex r=0.1931 vs seq40 r=0.1961 (identical=False) |
| 67 | BN. no 'identical' claim for seq40 in final PDF | PASS | stale wording found = none |
| 68 | BO. key-level leakage eliminated across all levels | PASS | levels with keys crossing folds = none |

---

`FAIL = 0` with the remaining `WARNING` items listed above. The manuscript is not
declared submission-ready while any open item requires author action.
