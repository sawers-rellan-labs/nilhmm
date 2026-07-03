# RTIGER backward pass (rigidity-aware Baum-Welch backward)

Literal port of the fork's `backward` (the time-reversed mirror of
forward).

## Usage

``` r
rtiger_backward_cpp(logPSI, logA, logpsi, r)
```

## Arguments

- logPSI:

  s x (T+r) windowed-emission matrix (rtiger_productpsi_cpp).

- logA:

  s x s log transition matrix.

- logpsi:

  s x T log emission matrix (rtiger_getlogpsi_cpp).

- r:

  Rigidity.

## Value

s x T log backward matrix.
