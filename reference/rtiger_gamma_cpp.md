# RTIGER gamma (state posteriors from zeta)

Literal port of the fork's `gamma`. Works in linear space (zeta is
already exponentiated): the first r columns are the normalized
alpha\*beta at the window edge; interior columns add the windowed
cumulative sum of the off-diagonal zeta to the diagonal term; the tail
repeats; columns are then normalized to sum to 1.

## Usage

``` r
rtiger_gamma_cpp(zeta, logalpha, logbeta, r)
```

## Arguments

- zeta:

  T x s x s array from rtiger_zeta_cpp.

## Value

s x T matrix of state posteriors.
