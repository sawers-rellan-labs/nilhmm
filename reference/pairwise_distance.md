# Pairwise marker relatedness (r2, mutual information, or variation of information)

Compute a symmetric markers x markers relatedness matrix from a markers
x samples dosage matrix, under a pluggable measure. The measure is
computed over the samples where both markers are non-missing (any
`NA`/`NaN`). The result is handed to
[`select_independent()`](https://sawers-rellan-labs.github.io/nilhmm/reference/select_independent.md)
for maximal-independent-set marker thinning; the measure is thresholded
in its own units (no r2\<-\>MI conversion).

## Usage

``` r
pairwise_distance(
  geno,
  method = c("r2", "mi", "vi"),
  base = c("nats", "bits"),
  max_markers = 7000L
)
```

## Arguments

- geno:

  Numeric matrix, rows = markers, cols = samples; dosages (typically
  `0`/`1`/`2`). Missing entries are `NA`/`NaN`. Row names (marker ids)
  are carried onto both dimensions of the result.

- method:

  One of `"r2"` (default), `"mi"`, `"vi"`.

- base:

  For `mi`/`vi`: information unit, `"nats"` (default) or `"bits"`.
  Ignored for `r2`.

- max_markers:

  Safety guard on the marker count (default `7000L`, the validated
  scale). The output is an O(n^2) markers x markers matrix, so a request
  above `max_markers` errors cleanly (before any large allocation)
  rather than risking a `malloc` failure; large-but-allowed sizes warn.
  Raise it (up to the ceiling in `getOption("nilHMM.marker_hard_cap")`,
  default `30000`) only if the memory is available; the intended
  workflow prunes per chromosome, well under the default.

## Value

Symmetric numeric matrix (markers x markers) with `attr(, "kind")`
`"similarity"` (r2, mi) or `"distance"` (vi) and `attr(, "method")` the
measure. Dimnames from `rownames(geno)`.

## Details

- `r2` (default, a similarity):

  `cor(g_i, g_j)^2`, the squared Pearson correlation of the dosage rows
  – the PLINK `--indep-pairwise` convention. High = related. A constant
  (zero-variance) marker contributes `0`.

- `mi` (a similarity):

  plug-in mutual information from the joint genotype contingency table,
  with the Miller-Madow bias correction
  `I_MM = I_plugin + (m_X + m_Y - m_XY - 1) / (2N)`. High = related.

- `vi` (a distance, a true metric):

  variation of information `VI(X,Y) = H(X) + H(Y) - 2*I(X,Y)`, on
  MM-corrected entropies. LOW = related; `VI` is a metric (obeys the
  triangle inequality), unlike r2/MI.

For near-linearly-dependent genotypes (e.g. RIL ancestry dosages) `r2`
and `mi` rank pairs almost identically; `mi`/`vi` earn their keep on
non-monotone or multiallelic data and for the metric/embedding
direction.

## See also

[`select_independent()`](https://sawers-rellan-labs.github.io/nilhmm/reference/select_independent.md)

## Examples

``` r
set.seed(1)
geno <- matrix(sample(0:2, 5 * 20, replace = TRUE), nrow = 5,
               dimnames = list(paste0("m", 1:5), NULL))
d <- pairwise_distance(geno, "r2")
attr(d, "kind")   # "similarity"
#> [1] "similarity"
```
