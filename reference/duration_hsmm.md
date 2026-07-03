# Explicit-duration (HSMM) sojourn – reserved

Explicit sojourn distribution. Reserved for a genuine duration model;
**must not** be used to bake in the design's fitted Gamma (S7
circularity trap).

## Usage

``` r
duration_hsmm(sojourn = NULL)
```

## Arguments

- sojourn:

  Sojourn distribution spec (placeholder).

## Value

A duration spec for
[`fit()`](https://sawers-rellan-labs.github.io/nilhmm/reference/fit.md).

## Examples

``` r
duration_hsmm()   # reserved placeholder spec
#> $type
#> [1] "hsmm"
#> 
#> $sojourn
#> NULL
#> 
#> attr(,"class")
#> [1] "nilHMM_duration_hsmm" "nilHMM_duration"     
```
