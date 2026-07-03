# Geometric (memoryless) duration

Self-transition rate `r`; the nilHMM transition. `r` is a resolution
hyperparameter, **not** an MLE – calibrate by KS-vs-sim (memory
`rigidity-not-mle`,
[`calibrate_r()`](https://sawers-rellan-labs.github.io/nilhmm/reference/calibrate_r.md)).

## Usage

``` r
duration_geometric(r = 0.01)
```

## Arguments

- r:

  Self-transition / recombination rate between adjacent markers.

## Value

A duration spec for
[`fit()`](https://sawers-rellan-labs.github.io/nilhmm/reference/fit.md).

## Examples

``` r
duration_geometric(r = 1e-4)
#> $type
#> [1] "geometric"
#> 
#> $r
#> [1] 1e-04
#> 
#> attr(,"class")
#> [1] "nilHMM_duration_geometric" "nilHMM_duration"          
```
