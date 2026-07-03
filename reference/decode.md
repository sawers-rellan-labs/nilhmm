# Decode the most-likely state path (Viterbi)

Decode the most-likely state path (Viterbi)

## Usage

``` r
decode(model, obs)
```

## Arguments

- model:

  A fitted model from
  [`fit()`](https://sawers-rellan-labs.github.io/nilhmm/reference/fit.md).

- obs:

  A per-(sample, chromosome) observation table.

## Value

Integer macro-state path over markers (0 = REF, 1 = HET, 2 = ALT); for
rigidity the expanded sub-states are mapped back to macro-states.

## Examples

``` r
obs <- data.frame(n = c(10, 9, 11, 8, 12), a = c(0, 0, 1, 4, 6))
model <- fit(obs, emission_count(), duration_geometric(1e-4),
             priors = design_priors("BC2S2"))
decode(model, obs)
#> [1] 0 0 0 1 1
```
