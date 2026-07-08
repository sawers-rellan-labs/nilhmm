# Pairwise marker relatedness matrix (r2 / MI / VI)

Build a symmetric markers x markers relatedness matrix from a markers x
samples dosage matrix, under a pluggable measure, for feeding
`fast_indep_cpp`. Computed over the samples where both markers are
non-missing (any non-finite entry is missing).

## Usage

``` r
pairwise_distance_cpp(geno, method, base)
```

## Arguments

- geno:

  Numeric matrix, rows = markers, cols = samples; dosages (typically
  0/1/2). Missing entries are `NA`/`NaN`.

- method:

  0 = r2 (squared Pearson correlation, a similarity), 1 = mi
  (Miller-Madow mutual information, a similarity), 2 = vi (variation of
  information, a distance / metric).

- base:

  For mi/vi only: 0 = nats, 1 = bits. Ignored for r2.

## Value

Symmetric numeric matrix (markers x markers). Diagonal: r2 = 1, mi =
H(X), vi = 0.
