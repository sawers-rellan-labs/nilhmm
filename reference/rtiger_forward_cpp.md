# RTIGER forward pass (rigidity-aware Baum-Welch forward)

Literal port of the fork's `forward`: at each step a state either
"stays" (diagonal transition) or was "entered" from another state r
positions back.

## Usage

``` r
rtiger_forward_cpp(logPI, logPSI, logA, logpsi, r)
```

## Arguments

- logPI:

  length-s log start probabilities.

- logPSI:

  s x (T+r) windowed-emission matrix (rtiger_productpsi_cpp).

- logA:

  s x s log transition matrix.

- logpsi:

  s x T log emission matrix (rtiger_getlogpsi_cpp).

- r:

  Rigidity.

## Value

s x T log forward matrix (cols 1 and \>T-r+1 stay -Inf).
