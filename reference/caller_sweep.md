# Sweep a caller's segmentation parameter with an amortized fit

Fit the emission once (rtiger: the joint EM; bbnil: per-sample means, or
nothing when `fit_means = FALSE`) and sweep `values`, fanning the
decodes over `threads`. `refit` controls the accuracy/speed trade (see
details).

## Usage

``` r
caller_sweep(
  data,
  caller = c("rtiger", "bbnil", "lbimpute"),
  values,
  refit = c("none", "cold"),
  design = NULL,
  f_1 = NULL,
  f_2 = NULL,
  threads = 1L,
  ref = NULL,
  min_reads = 1L,
  err = 0.01,
  conc = 20,
  fit_means = FALSE,
  seed = 1L,
  postprocess = TRUE,
  unit = c("bp", "cm"),
  genotypeerr = 0.05,
  drp = FALSE,
  source = "nilHMM",
  donor = NA_character_
)
```

## Arguments

- data:

  Common input: `name, chr, pos, n_ref, n_alt` (+ optional `donor`;
  `lbimpute` with `unit = "cm"` also needs a `cm` map-position column).

- caller:

  `"rtiger"` (sweeps `rigidity`), `"bbnil"` (sweeps `rrate`), or
  `"lbimpute"` (sweeps `recombdist`).

- values:

  Parameter grid to sweep.

- refit:

  `"none"` (fit once at `ref` and reuse – exact at `ref`, a close
  approximation elsewhere; recommended for calibration, as it isolates
  the segmentation prior) or `"cold"` (fit each value from scratch –
  exact per value, the baseline). For bbnil with `fit_means = FALSE` the
  emission is `rrate`-independent, so both are identical (and exact).
  This sweep *finds* the best value; for the exact final calls, refit
  once with
  [`call_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md)
  at the chosen value. **Ignored for `lbimpute`**: `recombdist` touches
  only the transition, never the emission, so every swept value is
  already EXACT – identical to a cold
  `call_ancestry(caller = "lbimpute", recombdist = v)` (the emission is
  computed once per run and only the Viterbi transition is re-run over
  the grid, batched in C++).

- design, f_1, f_2:

  Population priors (a design name, or explicit `f_1`,`f_2`). For
  `lbimpute` these only seed the start distribution (flat if absent).

- threads:

  Fan-out width
  ([`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) on
  unix; serial otherwise).

- ref:

  Reference value for the shared fit (default `median(values)`, rounded
  for rtiger). Unused for `lbimpute` (no shared fit).

- min_reads:

  Minimum read depth to keep a marker before decoding (default `1L`);
  `0L` keeps all. **No-op for `lbimpute`**, which keeps zero-read
  markers (flat emission) so the distance transition sees true marker
  spacing.

- err, conc, fit_means:

  bbnil count-emission parameters (fixed across the grid). `err` is also
  LB-Impute's per-read error (`readerr`).

- seed, postprocess:

  rtiger fit seed and border post-processing.

- unit, genotypeerr, drp:

  `lbimpute` only: `unit` is the transition coordinate (`"bp"` physical
  / `"cm"` genetic map; output stays bp), `genotypeerr` the emission
  floor/ceiling, `drp` the single-vs-double homozygous-switch cost. Same
  validation/mismatch warnings as
  [`call_states()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_states.md).

- source, donor:

  Output labels.

## Value

A common-schema segment table with an added column named for the swept
parameter (`rigidity`, `rrate`, or `recombdist`) tagging each value's
calls.
