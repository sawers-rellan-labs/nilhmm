# Rigidity duration (RTIGER reimplementation)

Enforces a hard minimum run length of `rigidity` markers by expanding
each state into a chain of `rigidity` sub-states (phase-type / Erlang).
Reimplements the subset of RTIGER we use (BetaBinomial emission +
rigidity + Viterbi + KS-calibrated `r`); `optimize_R` is intentionally
NOT ported (it over-rigidifies, memory
`optimize-R-overrigidifies-nil-sim`). Validate against
`tests/fixtures/baseline_pre_refactor/rtiger_rigidity_ref/`.

## Usage

``` r
duration_rigidity(rigidity = 5L, xrate = 0.01)
```

## Arguments

- rigidity:

  Minimum run length in markers (integer \>= 1).

- xrate:

  Per-marker switch / exit probability at the free (post-minimum) state
  – the geometric tail beyond the enforced minimum run. `rigidity` of 1
  reduces the expansion to a plain geometric transition with this exit
  rate.

## Value

A duration spec for
[`fit()`](https://sawers-rellan-labs.github.io/nilhmm/reference/fit.md).

## Examples

``` r
duration_rigidity(rigidity = 5, xrate = 2e-3)
#> $type
#> [1] "rigidity"
#> 
#> $r
#> [1] 5
#> 
#> $p_switch
#> [1] 0.002
#> 
#> attr(,"class")
#> [1] "nilHMM_duration_rigidity" "nilHMM_duration"         
```
