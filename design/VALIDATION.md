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

## Rigidity (`rtiger`) vs RTIGER's `calls_taxa_r5.csv` — CALIBRATED

RTIGER EM-fits its emissions in Julia and uses its own transition, so this is a
distributional match (not bit-identical). Calibrated on Zx via a
`(fit_means, p_switch, rigidity)` sweep against the actual RTIGER calls
(`agent/calibrate_rigidity.R`), scored on donor fraction, donor-block-size KS,
and per-marker state concordance with RTIGER.

**Calibrated setting: `fit_means = TRUE, p_switch = 2e-3, r = 5`.**

| Zx (target: donor 0.103, median 8.2 Mb) | donor frac | median block | block-KS | concordance vs RTIGER |
|---|---|---|---|---|
| fixed means, `p_switch=1e-3` (uncalibrated) | 0.059 | 4.9 Mb | 0.19 | 0.906 |
| **`fit_means`, `p_switch=2e-3`** | **0.108** | 10.4 Mb | **0.133** | **0.949** |

Findings:
- **`fit_means=TRUE` is what matches RTIGER** — it lifts donor fraction from ~0.06
  to ~0.108 (RTIGER 0.103) and raises per-marker concordance 0.906 -> **0.949**
  (95% of calls agree with RTIGER's independent Julia implementation).
- **Rigidity `r` (3/5/8) is ~irrelevant on skim density** — a few-marker minimum
  run is tiny vs marker spacing; `p_switch` is the effective length knob.
  `p_switch >= ~5e-3` over-calls (donor fraction inflates, concordance drops).
- **Residual:** median block runs ~25% long (10.4 vs 8.2 Mb) and can't shrink
  without inflating donor fraction — because we EM-fit only the emission *means*,
  while RTIGER jointly EM-fits emissions *and* transitions. Matching both block
  size and donor fraction simultaneously would need transition fitting (not done).
- **Per-taxon:** RTIGER's median block varies 4–15 Mb across taxa (real
  donor-specific structure), so block-size matching wants a per-taxon `p_switch`;
  donor fraction and ~95% concordance hold with the single calibrated setting.

The rigidity *mechanism* is also unit-tested (min-run enforced; r=1 == geometric;
`tests/testthat/test-rigidity.R`).

## Performance (arm64-native, head-to-head)

Benchmark: count caller on Zd (255 samples x 49002 markers = 12.5M cells), HMM
core only (VCF I/O excluded), best of 3. **Both native arm64** — the `nilhmm`
conda env is x86_64 (whole anaconda3 is Intel), so it runs under Rosetta and was
NOT used; the Python figure is from a native-arm64 venv (numpy 2.5 / scipy 1.18).

| implementation | time | throughput |
|---|---|---|
| Python (original, numpy-batched Viterbi + emission memoization) | 1.35s | 9.3M cells/s |
| R + Rcpp (naive: per-(sample,chr) loop) | 5.85s | 2.1M cells/s |
| R + Rcpp (optimized: see below) | **2.5s** | 5.0M cells/s |

Optimizations applied (each correctness-preserving — golden slice stays
bit-identical):
1. **Emission memoization** — BetaBinomial per distinct (n,a) pair, not per cell
   (5.85 -> 4.22s).
2. **O(n) split-by-sample** — was O(samples x rows) re-scan (full 1403-sample
   skim validation ~260s -> ~15s).
3. **`viterbi_batch_cpp`** — one C++ Viterbi across all samples of a chromosome,
   reading memoized emission through an index matrix (mirrors numpy batching);
   plus `rle_segments_batch_cpp` (segment RLE for the whole sample matrix in C++).

Serial result: R is **~1.9x** of Python (2.5s vs 1.35s). Profiling shows the
**compiled kernel alone (~1.2s) matches Python's *total* runtime** — the residual
gap is R-level pivot overhead (`unique`/`match`/`split` for the long->wide
reshape, ~0.6s) plus numpy's vectorized inner Viterbi.

4. **Threaded decode** (`call_ancestry(parallel = TRUE)`, opt-in) — RcppParallel
   /TBB splits the independent sample axis across cores (`viterbi_batch_par_cpp`).
   Bit-identical to serial.

| threads (Zd) | time | vs Python |
|---|---|---|
| serial | 2.50s | 1.9x |
| 4 | 1.30s | 0.96x |
| **8** | **1.24s** | **0.92x (faster)** |

At 8 threads R **beats the original Python** (1.24s vs 1.35s) while staying
bit-identical. Scaling is ~2x (not 8x) because the serial long->wide pivot
(~0.6s, Amdahl) and the slower efficiency cores bound it. The batched path
activates for the fixed-means count caller on a rectangular cohort; fit_means /
rigidity / ragged inputs fall back to the (correct, slower) per-sample loop. NEON
SIMD (xsimd) was considered but is only 2-wide for doubles on arm64 and the
argmax branch blocks auto-vectorization (`-O3 -mcpu=native` gave 0 gain) —
threading is the higher-ROI lever here. `parallel = FALSE` remains the default so
the validated path is deterministic without thread setup.

## Status

Count caller (geometric, both designs) validated at cohort scale. Rigidity mode
implemented and mechanism-tested; distributional calibration to RTIGER is future
work. `gt`/`dosage` emissions remain stubs (skimbin, MolBreeding gt path).
