# Threaded batched log-space Viterbi (RcppParallel)

Identical result to
[viterbi_batch_cpp](https://sawers-rellan-labs.github.io/nilhmm/reference/viterbi_batch_cpp.md);
splits the sample axis across cores.

## Usage

``` r
viterbi_batch_par_cpp(log_init, log_trans, em_uniq, inv)
```

## Arguments

- log_init:

  Length-K log initial-state probabilities.

- log_trans:

  K x K log transition matrix (row = from, col = to).

- em_uniq:

  U x K log emission table, one row per distinct observation.

- inv:

  T x S integer matrix; inv(t,s) is the 0-based row of `em_uniq` giving
  the emission for sample s at marker t.

## Value

T x S integer matrix of most-likely state paths (0-based).
