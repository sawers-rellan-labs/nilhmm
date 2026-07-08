# Design-derived expected heterozygosity from the pedigree inbreeding coefficient

`phet = (1 - F)/2` when `F` in `[0, 1]` (TASSEL's
`ViterbiAlgorithmPlugin`), else the `default` fallback (TASSEL
`probHeterozygous` = 0.07). Keeps `phet` design-derived, not a magic
constant (no hardcoded design parameters).

## Usage

``` r
.fsfhap_phet(F, default = 0.07)
```

## Arguments

- F:

  Inbreeding coefficient (pedigree column 7 / `1 - 0.5^nself`).

- default:

  Fallback when `F` is outside `[0, 1]`.

## Value

The expected-heterozygosity scalar for
[`.fsfhap_impute()`](https://sawers-rellan-labs.github.io/nilhmm/reference/dot-fsfhap_impute.md).
