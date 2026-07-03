# Count (BetaBinomial) emission

Emission on `(n_ref, n_alt)` read counts. State means `theta`; depth-0
markers emit flat. **Means are FITTABLE** (EM-fit /
reference-bias-corrected), not the fixed `c(err, 0.5, 1 - err)` – the
BRB run showed fixed `theta_ALT = 1 - err` collapses ALT-\>HET under
reference bias (S10, BRB_run_findings.md).

## Usage

``` r
emission_count(err = 0.01, conc = 20, fit_means = FALSE)
```

## Arguments

- err:

  Baseline genotyping/sequencing error (initialises `theta`).

- conc:

  BetaBinomial concentration (overdispersion); near-no-op on BRB (S10)
  but retained for the regime axis.

- fit_means:

  If `TRUE`, EM-fit the state means; if `FALSE` (default) use the fixed
  `c(err, 0.5, 1 - err)` that reproduces the Python baseline.
  Reference-biased data (RNA / BRB) needs `TRUE` (S10).

## Value

An emission spec for
[`fit()`](https://sawers-rellan-labs.github.io/nilhmm/reference/fit.md).

## Examples

``` r
emission_count(err = 0.01, conc = 20)               # fixed means (Python baseline)
#> $type
#> [1] "count"
#> 
#> $err
#> [1] 0.01
#> 
#> $conc
#> [1] 20
#> 
#> $fit_means
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "nilHMM_emission_count" "nilHMM_emission"      
emission_count(err = 0.01, fit_means = TRUE)         # EM-fit means (RNA / BRB)
#> $type
#> [1] "count"
#> 
#> $err
#> [1] 0.01
#> 
#> $conc
#> [1] 20
#> 
#> $fit_means
#> [1] TRUE
#> 
#> attr(,"class")
#> [1] "nilHMM_emission_count" "nilHMM_emission"      
```
