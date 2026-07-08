# Degrade simulated truth to observed allelic read counts

Turn a
[`simulate_nil()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_nil.md)
truth table into observed data under a generic sequencing model: per
marker draw depth `~ Poisson(depth)`, mask a fraction `p_missing` (and
zero-depth markers) to missing, and draw ALT reads
`~ Binomial(depth, p)` where `p` is the state's donor fraction (`0` /
`0.5` / `1`) flipped by the per-read `error`. Returns the
[`read_counts()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_counts.md)-style
columns plus a hard genotype `g` (masked to `3` where missing), so the
result drives both the count callers (`n_ref`/`n_alt`) and the genotype
callers (`g`). The named experiment-specific coverage regimes (skim /
BRB / target) live in the pipeline, not here.

## Usage

``` r
simulate_counts(truth, depth = 1, error = 0.01, p_missing = 0, seed = NULL)
```

## Arguments

- truth:

  A
  [`simulate_nil()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_nil.md)
  table (needs `name, donor, chr, pos, state`).

- depth:

  Mean per-marker read depth (Poisson mean).

- error:

  Per-read error rate (homozygote miscall probability).

- p_missing:

  Extra fraction of markers set missing (on top of zero-depth).

- seed:

  Optional RNG seed.

## Value

data.frame `name, donor, chr, pos, n_ref, n_alt, g` (`g` in
`{0,1,2,3}`).

## See also

[`simulate_nil()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_nil.md),
[`read_counts()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_counts.md).

## Examples

``` r
if (requireNamespace("simcross", quietly = TRUE)) {
  truth <- simulate_nil("BC2S2", n = 2, n_chr = 2, n_markers = 100, seed = 1)
  obs <- simulate_counts(truth, depth = 6, seed = 1)
  head(obs)
}
#> Error in simulate_nil("BC2S2", n = 2, n_chr = 2, n_markers = 100, seed = 1): unused argument (n_chr = 2)
```
