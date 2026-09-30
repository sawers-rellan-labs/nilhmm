# TODO (user task): binhmm emission Gaussian -> BetaBinomial

Owner: Fausto (user task, noted 2026-09-30). Not started.

## Now
`binhmm` (R/binhmm.R) bins the genome (default 1 Mb), reduces each bin to one number,
`alt_freq = alt reads / total reads`, and decodes REF/HET/ALT per bin with a 3-state
**Gaussian**-emission HMM (`.binhmm_gauss_states`, the default `cluster_method = "gauss"`):
REF anchored at the empirical baseline floor, REF width floored, HET/ALT means seeded from
the data and refined by Viterbi hard-EM, collapse guards and a `min_run` de-speckle. The
legacy backends (`gmm` / `kmeans` / `rebmix` + the fixed-confusion `.BINHMM_EMISS` Viterbi)
cluster the same `alt_freq`.

Reducing a bin to a frequency throws away its depth: a bin with 3 of 10 reads ALT and one
with 300 of 1,000 get the same observation and the same emission, and the Gaussian width is
a pooled guess rather than a consequence of the counts.

## Task
Give `binhmm` a **count emission**: per bin, the summed `ref` / `alt` read depths (`n = ref + alt`,
`k = alt`), scored with a **BetaBinomial(n, a_i, b_i)** per state, as the count callers
already do (`emission_count()` in R/emissions.R: BetaBinomial over ref/alt read depths with
`err` and `conc`; used by `bbnil` and `rtiger`, CLAUDE.md caller table). Depth then sets how
sharp each bin's evidence is, and overdispersion (`conc`) replaces the Gaussian width.

Points to settle:
- **State means on the diluted scale**: binhmm's REF floor is a small positive alt fraction,
  not 0 (the reason for the anchoring). Map REF/HET/ALT to BetaBinomial means that keep that
  floor (fit or anchor the REF mean like `.binhmm_gauss_states` does; `fit_means` in
  `emission_count()` may already cover it).
- **Reuse, not a new kernel**: build the log-emission matrix with the existing BetaBinomial
  code (`RcppExports`: "Log BetaBinomial emission matrix over REF/HET/ALT states") and decode
  with the package's Viterbi (see TODO_viterbi_kernel_consolidation.md / `viterbi_sweep_cpp`
  on branch feat/viterbi-sweep-kernel).
- **What stays**: binning, sticky transitions from the design priors, collapse guards, `min_run`.
  Keep `cluster_method = "gauss"` as an option until the new emission is validated.
- **Validation**: against the Gaussian backend and the other count callers on the benchmark
  sets used in VALIDATION.md; watch the HET over-call and high-coverage fragmentation that
  the Gaussian anchoring was built to fix.
- Update CLAUDE.md's caller table (`binhmm` row: emission) and NEWS.md.
