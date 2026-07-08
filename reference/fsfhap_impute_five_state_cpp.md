# FSFHap stage 3: 5-state Viterbi-training EM imputation (imputeUsingViterbiFiveState)

The real FSFHap imputation step (via `ViterbiAlgorithmPlugin`). Faithful
port of `imputeUsingViterbiFiveState`: per-taxon 5-state Viterbi on
non-missing parent calls with a distance-scaled transition, EM
re-estimating emission + transition from state-count matrices until the
emission-count matrix stabilizes.

## Usage

``` r
fsfhap_impute_five_state_cpp(G, pos, phet, max_iter = 50L)
```

## Arguments

- G:

  Integer matrix, taxa x sites, parent-called `g` in 0 A-hom, 1 het, 2
  C-hom, 3 missing (stage-1b output); one family, one chromosome,
  sorted.

- pos:

  Integer marker positions (bp), length = ncol(G).

- phet:

  Design-derived expected heterozygosity (`(1-F)/2`); sets the initial
  state distribution `{phom, .25 phet, .5 phet, .25 phet, phom}`.

- max_iter:

  EM iteration cap (TASSEL 50).

## Value

List: `imputed` (taxa x sites, `0`=A / `1`=het / `2`=C / `3`=missing on
undecoded sites), `iters` (EM iterations run), `emission` (final 5x3).
