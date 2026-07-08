# FSFHap stage 3: 5-state EM imputation for one family x chromosome

Thin wrapper over
[`fsfhap_impute_five_state_cpp()`](https://sawers-rellan-labs.github.io/nilhmm/reference/fsfhap_impute_five_state_cpp.md)
— the real FSFHap imputation (via TASSEL's `ViterbiAlgorithmPlugin` →
`imputeUsingViterbiFiveState`), then
[`fsfhap_fill_gaps_cpp()`](https://sawers-rellan-labs.github.io/nilhmm/reference/fsfhap_fill_gaps_cpp.md)
(`fillGapsInAlignment`) if `fill_gaps`. Runs on the parent-called
A/het/C frame from
[`.fsfhap_call_parents_bc()`](https://sawers-rellan-labs.github.io/nilhmm/reference/dot-fsfhap_call_parents_bc.md).

## Usage

``` r
.fsfhap_impute(G, pos, phet, max_iter = 50L, fill_gaps = TRUE)
```

## Arguments

- G:

  Integer matrix, taxa x sites, parent-called `g` in 0,1,2,3.

- pos:

  Integer marker positions (bp), length = ncol(G).

- phet:

  Expected heterozygosity. **Design-derived, no default** — the caller
  supplies `(1 - F)/2` from the pedigree inbreeding coefficient (BC1S4
  F=0.9375 → 0.03125), mirroring TASSEL's `ViterbiAlgorithmPlugin`.

- max_iter:

  EM iteration cap (TASSEL 50).

- fill_gaps:

  Forward-fill missing runs between equal flanking calls.

## Value

The imputed matrix (taxa x sites, `0`/`1`/`2`/`3`), with attributes
`iters` and `emission` (final 5x3) from the EM.
