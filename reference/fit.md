# Fit HMM emission/transition parameters

For a fixed-mean emission (the default, reproducing the Python count
caller) this resolves `theta` and the transition and is otherwise a
no-op. For a fittable-mean emission (`fit_means = TRUE`; required for
RNA / reference-biased data per S10) it Baum-Welch EM-fits the emission
means. Viterbi/EM hot loops are in Rcpp.

## Usage

``` r
fit(obs, emission, duration, priors, control = list())
```

## Arguments

- obs:

  A per-(sample, chromosome) observation table from an emission's reader
  (for `count`: columns `n`, `a`).

- emission:

  Emission spec
  ([`emission_count()`](https://sawers-rellan-labs.github.io/nilhmm/reference/emission_count.md)
  /
  [`emission_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/emission_gt.md)).

- duration:

  Duration spec
  ([`duration_geometric()`](https://sawers-rellan-labs.github.io/nilhmm/reference/duration_geometric.md)
  /
  [`duration_rigidity()`](https://sawers-rellan-labs.github.io/nilhmm/reference/duration_rigidity.md)
  /
  [`duration_hsmm()`](https://sawers-rellan-labs.github.io/nilhmm/reference/duration_hsmm.md)).

- priors:

  Single-locus genotype-frequency priors `list(f_1, f_2)`.

- control:

  EM control (`max_iter`, `tol`).

## Value

A fitted model: `list(theta, log_start, log_trans, emission, duration)`.

## Examples

``` r
# Per-(sample, chromosome) count obs: n = depth, a = alt count.
obs <- data.frame(n = c(10, 9, 11, 8, 12), a = c(0, 0, 1, 4, 6))
model <- fit(obs, emission_count(), duration_geometric(1e-4),
             priors = design_priors("BC2S2"))
```
