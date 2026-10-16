# BENCHMARK NOVELTY MATRIX — V8.1 FINAL FREEZE

> Synced to the V8.1 frozen state (`FROZEN_RESULTS.json`). Supersedes the V8 draft matrix.

## 0. Honest coverage of prior work (audited in V8 PHASE 1)

| Claim | Prior work | Our position |
|---|---|---|
| Leakage-aware / family-OOD PPI ΔΔG evaluation | Bushuiev et al. 2024 ([arXiv:2404.10457](https://ar5iv.labs.arxiv.org/html/2404.10457)) quantifies structural leakage (PDB-code 65%, sequence-30% 30%, SKEMPI v2.0 PPI-code 56%, expert grouping 0%) | **NOT ours. Do not claim.** |
| CATH-superfamily (family-OOD) splits + hard-mutation definition | CATH-ddG, Bioinformatics 2025 ([btaf228](https://doi.org/10.1093/bioinformatics/btaf228)) | **NOT ours. Do not claim.** |
| Per-target / per-structure Spearman reporting | already present in prior work | **NOT ours. Do not claim.** |
| Top-K / NDCG-style shortlist reporting | already present in prior work | **NOT ours. Do not claim.** |
| Overall leakage-aware benchmarking as a concept | prior work | **NOT ours. Do not claim.** |

## 1. External methods actually audited (V8.1 state)

| Method | Row-level predictions | Mapping used | Provenance | Evidence status |
|---|---|---|---|---|
| `DDGPred` | yes (CATH-ddG table) | `EXACT_UNIQUE` only | literature-only provenance | TIER 2 |
| `ddG_FoldX` | yes (CATH-ddG table) | `EXACT_UNIQUE` only | proprietary executable; predictions analysed | TIER 2 |
| `Flex-ddG` | yes (CATH-ddG table) | `EXACT_UNIQUE` only | literature-only provenance | TIER 2 |
| `proteinMPNN` | yes (HERMES repo + CATH-ddG table) | exact key | inverse-folding, not ΔΔG-supervised | TIER 2 |
| `ESM-IF` | yes (CATH-ddG table) | `EXACT_UNIQUE` only | inverse-folding, not ΔΔG-supervised | TIER 2 |
| `HERMES` | yes (official CSV) | exact 1:1 | cDNA117k stability fine-tune; **not** SKEMPI-binding-supervised; RMSE/MAE = NOT_APPLICABLE | TIER 2 |
| `PPIformer` | no official row-level predictions | — | upstream CPU blocker | NOT_RUN |
| `RDE-Network` | no official row-level predictions | — | needs self-run inference | NOT_RUN |
| `Prompt-DDG` | no official row-level predictions | — | no scorer CLI | NOT_RUN |
| `CATH-ddG` (itself) | **not released** | — | 88% of rows are its own training set | NOT_RUN |
| `PreFold-dG` | **not released** | No weights, no scorer | — | NOT_RUN |

## 2. Capability-resolved components — coverage

| Component | Covered by prior work? | Our status |
|---|---|---|
| Exact **between-target / within-target covariance-signal decomposition** `r = C_B + C_W` with identity QC | not identified | **contribution** (residual ≤ 3.3e-16) |
| **Controlled split-strictness** comparison of ΔC_B vs ΔC_W | not identified | **contribution** (between-fraction = 0.811) |
| **Target-weighting sensitivity** of the decomposition | not identified | **contribution** (C_B near-zero instability; C_W sign-stable) |
| **Matched-universe rank-reversal** analysis | not identified | **contribution** (range −0.27 … +0.82 across 10 strata) |
| **Fixed-K candidate-set geometry** audit, `D(n,K)=K(n−K)/[n(n+1)]` | not identified | **contribution** (40.8% of n≥5 targets have n ≤ K) |
| **Retrospective benchmark-guided remedy** attribution | not identified | **contribution** (retrospective only) |
| Integration into one capability-resolved evaluation framework | not identified | **contribution** |

## 3. Forbidden claims (V8.1, §19)

* first leakage-aware / first family-OOD / first per-target / first Top-K benchmark
* all existing benchmarks are invalid
* Overall Pearson is useless
* universal calibration–ranking orthogonality
* all PLMs fail; physics generally beats learning
* target size causes prediction difficulty
* HERMES log-odds is calibrated kcal/mol ΔΔG (**RMSE/MAE = NOT_APPLICABLE**)
* candidate-set floor proves SKEMPI ΔΔG values are intrinsically too similar
* benchmark prospectively predicted the successful remedy
* prospective wet-lab validity; universal model recommender
* `STRONGLY_SUPPORTED` merely because external models were added
* "recommend n ≥ 50 because fixed-K discrimination increases with n"
* any external capability claim computed on chain-ambiguous mappings

## 4. Positioning statement

> To our knowledge we did not identify prior work that **jointly** reports an exact
> between-/within-target covariance decomposition under matched, provenance-clean,
> family-OOD evaluation together with split-strictness, target-weighting, rank-reversal,
> fixed-K geometry and retrospective remedy analyses for PPI ΔΔG prediction.
> We make **no "first" claim** for leakage-aware benchmarking, family-OOD splitting,
> per-target reporting or Top-K evaluation.
