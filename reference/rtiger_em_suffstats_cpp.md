# RTIGER one-EM-iteration sufficient statistics (E-step + fold, all in C++)

Runs the per-chain E-step (getlogpsi..gamma) and accumulates the pooled
sufficient statistics, replacing the slow R-level fold
(apply()/tapply()). Mirrors the fork's EM accumulation: transition
band-sum of zeta, start = Sum of gamma column 1, and per-state emission
weights grouped by distinct (k,n).

## Usage

``` r
rtiger_em_suffstats_cpp(ks_list, ns_list, logPI, logA, alpha, beta, r, nstates)
```

## Arguments

- ks_list, ns_list:

  Lists of per-chain integer (k, n) vectors.

- logPI, logA:

  log start / log transition.

- alpha, beta:

  current BetaBinomial shape vectors.

- r, nstates:

  Rigidity and number of states.

## Value

list(sumZeta, startAcc, nOb, kvals, nvals, wmat, sumk, sumn).
