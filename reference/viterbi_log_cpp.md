# Generic log-space Viterbi decode (time-homogeneous transitions)

Generic log-space Viterbi decode (time-homogeneous transitions)

## Usage

``` r
viterbi_log_cpp(log_init, log_trans, log_emit)
```

## Arguments

- log_init:

  Length-K vector of log initial-state probabilities.

- log_trans:

  K x K matrix of log transition probabilities (row = from-state, column
  = to-state).

- log_emit:

  T x K matrix of log emission probabilities per marker/state.

## Value

Integer length-T most-likely state path, 0-indexed.
