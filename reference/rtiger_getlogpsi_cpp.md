# RTIGER emission log-probabilities (getlogpsi)

logpsi(i,t) = logpdf(BetaBinomial(n_t, a_i, b_i), k_t). Memoized over
distinct (k,n) pairs (as in the fork's getlogpsi cache; bit-identical
values).

## Usage

``` r
rtiger_getlogpsi_cpp(k, n, a, b)
```

## Arguments

- k:

  Integer vector of ref-allele counts (length T).

- n:

  Integer vector of totals (length T).

- a, b:

  Per-state BetaBinomial shape vectors (length s).

## Value

s x T matrix of log emission probabilities.
