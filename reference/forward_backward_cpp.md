# Forward-backward state posteriors (log-space; time-homogeneous transitions)

Forward-backward state posteriors (log-space; time-homogeneous
transitions)

## Usage

``` r
forward_backward_cpp(log_init, log_trans, log_emit)
```

## Arguments

- log_init:

  Length-K log initial-state probabilities.

- log_trans:

  K x K log transition matrix (row = from, col = to).

- log_emit:

  T x K log emission matrix.

## Value

T x K matrix of posterior state probabilities (rows sum to 1).
