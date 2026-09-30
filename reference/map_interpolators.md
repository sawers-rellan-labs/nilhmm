# Map interpolators: physical \<-\> genetic position

Build a per-chromosome monotone Hyman spline (a Marey map) through a
consensus map's `(bp, cm)` pairs, clamped to each chromosome's observed
range, and return a vectorized `function(chr, x)`. Duplicate coordinates
are collapsed by mean. The default map is the bundled B73 v5 consensus
map; pass any data frame with `chr`, `bp`, `cm` columns to use your own
(e.g. a population's native map).

## Usage

``` r
bp_to_cm(map = load_map())

cm_to_bp(map = load_map())
```

## Arguments

- map:

  A consensus map with columns `chr`, `bp`, `cm`. Defaults to
  [`load_map()`](https://sawers-rellan-labs.github.io/nilhmm/reference/load_map.md)
  (bundled B73 v5).

## Value

`bp_to_cm`: `function(chr, bp) -> cm`; `cm_to_bp`:
`function(chr, cm) -> bp`.

## See also

[`cm_to_mb()`](https://sawers-rellan-labs.github.io/nilhmm/reference/cm_to_mb.md),
[`load_map()`](https://sawers-rellan-labs.github.io/nilhmm/reference/load_map.md)

## Examples

``` r
to_cm <- bp_to_cm()          # bundled map
to_cm(1L, 1e6)
#> [1] 0.3921432
```
