# Build a physical marker grid from a genetic map

Spread `n_markers` markers across the chromosomes of `map` in proportion
to physical span, interpolating each marker's cM from the map (monotone
Hyman spline, bp -\> cM). This is how
[`simulate_nil()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_nil.md)
samples the bundled consensus map
([`load_map()`](https://sawers-rellan-labs.github.io/nilhmm/reference/load_map.md))
at a chosen density.

## Usage

``` r
build_marker_grid(map, n_markers = 2000L)
```

## Arguments

- map:

  A map with columns `chr`, `cm`, `bp` (e.g. from
  [`load_map()`](https://sawers-rellan-labs.github.io/nilhmm/reference/load_map.md)).

- n_markers:

  Approximate total marker count across all chromosomes.

## Value

data.frame `chr`, `pos` (bp), `cm`, sorted by `chr, pos`.

## See also

[`load_map()`](https://sawers-rellan-labs.github.io/nilhmm/reference/load_map.md),
[`simulate_nil()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_nil.md)

## Examples

``` r
g <- build_marker_grid(load_map(), n_markers = 500)
table(g$chr)
#> 
#>  1  2  3  4  5  6  7  8  9 10 
#> 73 57 56 59 53 41 44 43 38 36 
```
