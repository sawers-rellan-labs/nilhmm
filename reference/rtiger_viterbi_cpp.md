# RTIGER Viterbi (rigidity max-product + rigid backtrace)

Literal port of the fork's `viterbi`. Each step is a max over: "stay" in
the target state (diagonal), or "enter" it from another state r
positions back (using the windowed PSI). The backtrace fills a whole
r-block with the current state on a switch, giving the r-rigid path.
Tie-break is first state (seed at state 1, replace only on strict \>),
matching the fork's findmax.

## Usage

``` r
rtiger_viterbi_cpp(PI, PSI, psi, A, r)
```

## Arguments

- PI:

  length-s log start probabilities.

- PSI:

  s x (T+r) windowed-emission matrix.

- psi:

  s x T log emission matrix.

- A:

  s x s log transition matrix.

- r:

  Rigidity.

## Value

length-T integer state path (1-based states).
