# Batched LB-Impute Viterbi over a recombdist grid (shared emission)

For calibration: `recombdist` affects only the transition, never the
emission, so the (memoized) emission matrix is computed once per run and
only the Viterbi transition is re-run for each grid value – inside C++,
avoiding per-value R\<-\>C++ marshaling. Column `k` is **bit-identical**
to
[`lb_viterbi_cpp()`](https://sawers-rellan-labs.github.io/nilhmm/reference/lb_viterbi_cpp.md)
called with `recombdist = recombdists[k]`, so a sweep value equals a
cold `call_ancestry(caller = "lbimpute", recombdist = v)`.

## Usage

``` r
lb_viterbi_sweep_cpp(log_init, log_emit, tpos, recombdists, drp)
```

## Arguments

- log_init:

  Length-3 log initial-state probabilities (REF/HET/ALT).

- log_emit:

  T x 3 log emissions (from
  [`lb_emission_loglik_cpp()`](https://sawers-rellan-labs.github.io/nilhmm/reference/lb_emission_loglik_cpp.md)).

- tpos:

  Numeric length-T transition coordinate, non-decreasing (bp or cM).

- recombdists:

  Numeric grid of recombdist values (same units as `tpos`).

- drp:

  If `TRUE`, a homozygous-\>homozygous switch costs a single
  recombination rather than a double event.

## Value

A T x length(recombdists) integer matrix of state paths (0/1/2); column
k is the decode at `recombdists[k]`.
