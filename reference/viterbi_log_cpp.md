# Generic log-space Viterbi decode (time-homogeneous transitions)

Generic log-space Viterbi decode (time-homogeneous transitions)

## Usage

``` r
viterbi_log_cpp(log_init, log_trans, log_emit, tie_break = 0L)
```

## Arguments

- log_init:

  Length-K vector of log initial-state probabilities.

- log_trans:

  K x K matrix of log transition probabilities (row = from-state, column
  = to-state).

- log_emit:

  T x K matrix of log emission probabilities per marker/state.

- tie_break:

  Backpointer tie policy for exactly equal predecessor scores: `0`
  (default) keeps the FIRST (lowest-index) predecessor; `1`
  ("incumbent") keeps the LAST (highest-index) predecessor. With
  categorical/GT emissions, emission-degenerate positions (missing
  genotype -\> all states equal; het -\> REF and donor equal) make the
  segment boundary a genuine tie; because the introgression states
  outrank REF, `tie_break = 1` keeps the introgression and delays the
  switch back to REF, matching the intent of Holland's `hmmlearn` nNIL
  caller. Count (BetaBinomial) emissions essentially never tie, so this
  is a no-op for them. Only the transition backpointer honours this; the
  terminal argmax keeps the first max (as `numpy.argmax`).

## Value

Integer length-T most-likely state path, 0-indexed.
