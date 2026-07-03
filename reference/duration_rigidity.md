# Rigidity duration (RTIGER reimplementation)

Enforces a hard minimum run length of `r` markers by expanding each
state into a chain of `r` sub-states (phase-type / Erlang machinery).
Reimplements the subset of RTIGER we use (BetaBinomial emission +
rigidity + Viterbi + KS-calibrated `r`); `optimize_R` is intentionally
NOT ported (it over-rigidifies, memory
`optimize-R-overrigidifies-nil-sim`). Validate against
`tests/fixtures/baseline_pre_refactor/rtiger_rigidity_ref/`.

## Usage

``` r
duration_rigidity(r = 5L, p_switch = 0.01)
```

## Arguments

- r:

  Minimum run length in markers (the rigidity, integer \>= 1).

- p_switch:

  Per-marker switch probability at the free (post-minimum) state – the
  geometric tail beyond the enforced minimum run. Rigidity of 1 reduces
  the expansion to a plain geometric transition with this switch rate.

## Value

A duration spec for
[`fit()`](https://sawers-rellan-labs.github.io/nilhmm/reference/fit.md).

## Examples

``` r
duration_rigidity(r = 5, p_switch = 2e-3)
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
