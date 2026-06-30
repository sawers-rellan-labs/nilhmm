# Golden slice — fast dev-loop fixture (chr1 only)

A tiny, deterministic (input -> expected-output) slice for the **inner development
loop**: it decodes in seconds so engine changes can be regression-checked on
every `devtools::test()`. This is the cheap counterpart to the full-cohort
**distributional** gate in `../../../tests/fixtures/baseline_pre_refactor/`
(§9.4) — use this while coding, that for acceptance.

Lives under `tests/testthat/fixtures/` (NOT the Rbuildignored top-level
`tests/fixtures/`) so it ships in the built package and runs under `R CMD check`.
The HMM decodes each `(sample, chromosome)` independently, so a chr1 subset of
the frozen calls is a **valid** expected output for chr1 inputs.

## Contents

| Source | Samples | Input | Expected |
|---|---|---|---|
| `brb/`  | PN1_SID25, PN1_SID32, PN1_SID41 (Zd, share a ~297.1–297.46 Mb HET), PN1_SID11 (heavy, multi-segment) | `counts/<s>.chr1.tsv` (~3.4k markers each) | `expected_calls_chr1.csv` (12 segs, 4 non-REF) |
| `skim/` | PN14_SID1259 (Zd), PN3_SID235 (Zx), PN17_SID1630 (Zx) | `counts/<s>.chr1.tsv` (7389 chr1 biallelic markers) | `expected_calls_chr1.csv` (25 segs, 11 non-REF) |

Count TSV columns (no header): `chr  pos  ref_base  n_ref  alt_base  n_alt`.
Expected-calls schema: `source,donor,name,chr,start_bp,end_bp,state` (state 0/1/2).

## Provenance (inputs not otherwise vendored)
- BRB counts: `results/sim_calibration/brbseq_ks_wideseq/counts/` (BC2S3, ~2.8×).
- Skim counts: extracted from `data/nilhmm/{Zd,Zx}_counts.vcf.gz` (BC2S2, ~0.4×) —
  the **same** source the baseline consumed (biallelic FORMAT/AD, chr1 only, missing→0),
  so input↔expected is exactly paired. nilHMM count-caller recipe (from
  `agent/stage_nilhmm_calls.py`): `err=0.01, conc=20, f_1=0.0625, f_2=0.0938`,
  per-taxon `r` from `calibration/calibrated_params_all_taxa.csv`
  (Zd/Zx = 7e-6, Zv/Zl = 3e-6, Zh = 3e-5).
- Expected calls: chr1 subset of the frozen `baseline_pre_refactor/{brb,skim}_nilhmm`
  Python outputs. `shasum -c SHA256SUMS` to verify.

## Caveat (carried from the full baseline)
On **BRB**, the Python collapses ALT->HET (fixed `theta_ALT=0.99` vs ref-biased
~0.5 alt-fraction); the R engine fixes this with fittable emission means (§10).
So treat the BRB slice's ALT calls as *expected to differ* — use it for
HET-tail / run-length sanity, and the **skim** slice (fixed means adequate) as
the strict concordance target in the dev loop.
