# Fit the RTIGER emission once (to reuse across per-chromosome decodes)

Fits the faithful RTIGER EM (emission alpha/beta + transition) on ALL
chains in `data` (samples x chromosomes) and returns the fit, to be
passed as `call_states(..., caller = "rtiger", rtiger_fit = <fit>)`.
This keeps the joint per-family emission while letting each chromosome
be decoded in a separate, low-memory call – e.g. a genome-scale sweep
that would otherwise hold the whole panel per worker. Decoding is
per-chromosome regardless, so results are identical to a single
whole-family call.

## Usage

``` r
fit_rtiger(data, rigidity, threads = 1L, seed = 1L)
```

## Arguments

- data:

  Long observation table (`name, chr, pos, n_ref, n_alt`), as for
  [`call_states()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_states.md).
  Filter to covered markers upstream if using `min_reads`.

- rigidity:

  Integer RTIGER minimum run length.

- threads, seed:

  Passed to the fit (RcppParallel E-step; seeded init).

## Value

An RTIGER fit `list(A, pi, alpha, beta, iterations)`.

## See also

[`call_states()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_states.md)
