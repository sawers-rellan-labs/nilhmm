# Pairwise marker distance from map/physical coordinates

Build a symmetric markers x markers **distance** matrix of the absolute
coordinate difference `|pos_i - pos_j|` from a vector of marker
positions (genetic cM or physical bp). The result is handed to
[`select_independent()`](https://sawers-rellan-labs.github.io/nilhmm/reference/select_independent.md)
with `sense = "distance"`: thresholding at `t` links markers `<= t`
apart, so the selected independent set is a maximal set of markers
pairwise `> t` apart – distance-based thinning on the same machinery as
r2/MI/VI relatedness.

## Usage

``` r
position_distance(pos, chr = NULL, method = "cm", max_markers = 7000L)
```

## Arguments

- pos:

  Numeric vector of marker coordinates (cM or bp). `names(pos)`, if
  present, are carried onto both dimensions of the result.

- chr:

  Optional chromosome label per marker (same length as `pos`; any type
  comparable with `!=`). When given, cross-chromosome pairs are `Inf`.

- method:

  Label recorded in `attr(, "method")` (default `"cm"`); use `"pos"` or
  `"bp"` for physical coordinates. Does not change the arithmetic.

- max_markers:

  Safety guard on the marker count (default `7000L`, the validated
  scale). The output is an O(n^2) markers x markers matrix, so a request
  above `max_markers` errors cleanly (before any large allocation)
  rather than risking a `malloc` failure; large-but-allowed sizes warn.
  Raise it (up to the ceiling in `getOption("nilHMM.marker_hard_cap")`,
  default `30000`) only if the memory is available; the intended
  workflow prunes per chromosome, well under the default.

## Value

Symmetric numeric matrix (markers x markers) of `|pos_i - pos_j|`, zero
diagonal, with `attr(, "kind") = "distance"` and `attr(, "method")` the
label. Dimnames from `names(pos)` when present. Cross-chromosome pairs
are `Inf` when `chr` is supplied.

## Details

Distance sense: an edge exists when `|pos_i - pos_j| <= threshold`, so
the independent set has all pairwise gaps `> threshold`. For a cM map,
`threshold` is in cM (e.g. `0.1` ~ Chen's 0.1-cM thinning). On a
one-dimensional coordinate the maximal independent set coincides with a
greedy min-spacing sweep (greedy is optimal on a line); the matrix path
exists to give a uniform interface across r2/MI/VI/cM rather than to
beat the sweep.

If `chr` is supplied, cross-chromosome pairs are set to `Inf` (never
within any finite `threshold`, hence never linked), so a genome-wide
matrix still prunes correctly. Per-chromosome use is preferred and
cheaper (the matrix is O(n^2)); the `chr` argument is a convenience for
a single pooled call.

## See also

[`select_independent()`](https://sawers-rellan-labs.github.io/nilhmm/reference/select_independent.md),
[`pairwise_distance()`](https://sawers-rellan-labs.github.io/nilhmm/reference/pairwise_distance.md)

## Examples

``` r
cm <- c(m1 = 0, m2 = 0.05, m3 = 0.2, m4 = 0.9)
d <- position_distance(cm)
attr(d, "kind")   # "distance"
#> [1] "distance"
# Prune to markers pairwise > 0.1 cM apart:
select_independent(d, threshold = 0.1, sense = "distance")
#> [1] "m1" "m3" "m4"
#> attr(,"sets")
#> attr(,"sets")[[1]]
#> [1] "m1" "m3" "m4"
#> 
#> attr(,"size_dist")
#> 3 
#> 1 
#> attr(,"kind")
#> [1] "distance"
```
