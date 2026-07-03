# RTIGER full EM fit (E-step + M-step + convergence loop, all in C++)

The entire fit (port of the fork's `fit`/`EM`): per-chain rigidity
E-step, pooled M-steps (transition, start, emission via C++ Brent),
iterated until max(\|Δα\|,\|Δβ\|) \<= eps or max.iter. Deterministic
init (generate_params forms, randomize off). Returns the fitted
parameters and iteration count.

## Usage

``` r
rtiger_fit_cpp(
  ks_list,
  ns_list,
  r,
  nstates,
  eps,
  max_iter,
  threads,
  init_alpha,
  init_beta
)
```

## Arguments

- ks_list, ns_list:

  Lists of per-chain integer (k, n) vectors.

- r, nstates:

  Rigidity and number of states.

- eps, max_iter:

  Convergence tolerance and iteration cap.

## Value

list(A, pi, alpha, beta, iterations).
