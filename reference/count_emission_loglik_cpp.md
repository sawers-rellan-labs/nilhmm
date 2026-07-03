# Log BetaBinomial emission matrix over REF/HET/ALT states

Log BetaBinomial emission matrix over REF/HET/ALT states

## Usage

``` r
count_emission_loglik_cpp(n, a, theta, conc)
```

## Arguments

- n:

  Integer vector of total depths (n_ref + n_alt) per marker.

- a:

  Integer vector of alt counts (n_alt) per marker.

- theta:

  Length-K vector of expected alt-fractions per state (e.g. c(err, 0.5,
  1 - err)).

- conc:

  BetaBinomial concentration (alpha + beta); larger -\> closer to
  Binomial.

## Value

T x K matrix of log emission probabilities.
