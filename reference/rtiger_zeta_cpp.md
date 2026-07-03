# RTIGER zeta (pairwise posteriors over the rigidity window)

Literal port of the fork's `zeta`: build the log values, then normalize
by the global max and exponentiate (so exp(-Inf - PO) = 0). Returned as
a T x s x s array (column-major), non-log, as in the Julia.

## Usage

``` r
rtiger_zeta_cpp(logalpha, logbeta, logA, logPSI, logpsi, r)
```

## Value

T x s x s numeric array.
