# Load the bundled consensus map

Returns the version-namespaced bundled consensus map (marker -\> chr,
cM, bp). Only `"v5"` (B73 v5,
[maize_map_v5](https://sawers-rellan-labs.github.io/nilhmm/reference/maize_map_v5.md))
is bundled; the map is overridable everywhere it is used (e.g. pass your
own `map` to
[`simulate_nil()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_nil.md)).

## Usage

``` r
load_map(version = "v5")
```

## Arguments

- version:

  Map version key (default `"v5"`).

## Value

The map data.frame `locus, chr, cm, bp` (with an `"assembly"` attr).

## See also

[maize_map_v5](https://sawers-rellan-labs.github.io/nilhmm/reference/maize_map_v5.md),
[`simulate_nil()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_nil.md),
[`build_marker_grid()`](https://sawers-rellan-labs.github.io/nilhmm/reference/build_marker_grid.md)

## Examples

``` r
map <- load_map("v5")
head(map)
#>             locus chr   cm       bp
#> 1    LOC103644366   1 0.00  37410.5
#> 2       pco082477   1 0.01  43988.0
#> 3 Zm00001eb000050   1 0.04 111468.0
#> 4         bnlg149   1 0.09 189070.0
#> 5 Zm00001eb000070   1 0.09 194512.0
#> 6           gbss2   1 0.10 201827.5
```
