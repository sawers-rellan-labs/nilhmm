# nilHMM R package — S9.4 validation (Task 5)

Result of validating the R+Rcpp reimplementation against the frozen Python
baselines (`tests/fixtures/baseline_pre_refactor/`). Per REFACTOR_R_PACKAGE.md
S9.4 the bar is **distributional equivalence** (call concordance, segment-size
KS), not bit-identical — EM init / Viterbi tie-breaking and scipy-vs-lgamma
BetaBinomial rounding differ across languages.

Harness: `agent/nilhmm_validate_full_cohort.R` (a zealtiger consumer script; owns
the data paths and per-taxon parameters — the package stays data-agnostic).
Run: `Rscript agent/nilhmm_validate_full_cohort.R [skim|brb|rtiger|all]`.

## Count caller (`nnil`, fixed means) — the strict gate: PASS

Fixed emission means reproduce the Python count caller, so these are effectively
bit-identical (the residual ~1e-5 disagreement is a handful of near-tie markers
flipping on floating-point differences — within S9.4 tolerance).

| Cohort | Design | params | segments (R vs Py) | donor frac (R vs Py) | block-size KS D | marker concordance |
|---|---|---|---|---|---|---|
| **BRB** (341 NILs) | BC2S3 | r=1e-8, err=.01, conc=20 | 8394 vs 8394 | 0.0415 vs 0.0415 | 0.0011 | **0.99999** (7.39M markers) |
| **skim** (1403 NILs) | BC2S2 | per-taxon r (Zd/Zx 7e-6, Zv/Zl 3e-6, Zh 3e-5) | 27630 vs 27630 | 0.0674 vs 0.0674 | 0.0008 | **0.99999** (68.7M markers) |

Conclusion: the R count caller reproduces the frozen nilHMM baselines on both
regimes. **The refactor's core acceptance gate is met.**

## Rigidity (`rtiger`) vs RTIGER's `calls_taxa_r5.csv` — distributional, UNCALIBRATED

RTIGER EM-fits its emissions in Julia and uses its own transition; our rigidity
mode here runs with **fixed means**, `rigidity=5`, `p_switch=1e-3` — not tuned to
RTIGER. So this is a regime check, not a match.

| metric (40-sample subset) | R rigidity | RTIGER r5 |
|---|---|---|
| donor genome fraction | 0.058 | 0.102 |
| donor-block median (Mb) | 3.9 | 8.2 |
| block-size KS D | 0.223 | — |

The gap is the expected fixed-means under-calling (cf. `BRB_run_findings.md`):
RTIGER recovers ~2x more donor via EM-fit emissions and produces longer blocks.
**Matching RTIGER is a calibration task** — `fit_means=TRUE` plus a KS-calibrated
`(rigidity, p_switch)` — not part of this validation. The rigidity *mechanism* is
correct (min-run enforced; r=1 == geometric; see `tests/testthat/test-rigidity.R`).

## Performance (arm64-native, head-to-head)

Benchmark: count caller on Zd (255 samples x 49002 markers = 12.5M cells), HMM
core only (VCF I/O excluded), best of 3. **Both native arm64** — the `nilhmm`
conda env is x86_64 (whole anaconda3 is Intel), so it runs under Rosetta and was
NOT used; the Python figure is from a native-arm64 venv (numpy 2.5 / scipy 1.18).

| implementation | time | throughput |
|---|---|---|
| Python (original, numpy-batched Viterbi + emission memoization) | 1.35s | 9.3M cells/s |
| R + Rcpp (after emission memoization + split-once + Rcpp RLE) | 4.3s | 2.9M cells/s |

R is **~3x slower head-to-head** — and it is *not* the Rcpp hot loops (Viterbi +
emission are ~0.7s / 18% of the time). The gap is R-level orchestration: the
caller loops over 2550 (sample, chromosome) sequences, each paying per-call
`unique`/`match`, matrix allocation, `split`/`order`, and GC, where the Python
runs ONE numpy-vectorized Viterbi across all samples per chromosome. Emission
memoization (BetaBinomial per distinct (n,a) pair, not per cell) landed its
expected win (5.85s -> 4.22s); the residual is the per-sequence loop.

Practical impact: combined with an O(n) split-by-sample fix (was O(samples x
rows)), the full 1403-sample skim validation dropped from ~260s to ~25s. Closing
the remaining 3x to per-op parity needs **batching the Viterbi across samples in
Rcpp** (a `viterbi_batch_cpp` over a samples x markers x states array, mirroring
the Python) — a larger restructure, not yet done.

## Status

Count caller (geometric, both designs) validated at cohort scale. Rigidity mode
implemented and mechanism-tested; distributional calibration to RTIGER is future
work. `gt`/`dosage` emissions remain stubs (skimbin, MolBreeding gt path).
