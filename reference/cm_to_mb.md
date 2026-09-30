# Project segment cM coordinates to physical Mb

Converts segment coordinates in genetic (cM) space to physical Mb via
the inverse Marey spline
[`cm_to_bp()`](https://sawers-rellan-labs.github.io/nilhmm/reference/map_interpolators.md).
cM coordinates are assembly-robust; the bp/Mb they map to are tied to
the map's assembly (bundled = B73 v5).

## Usage

``` r
cm_to_mb(seg, map = load_map())
```

## Arguments

- seg:

  Segments with columns `chr`, `start_cm`, `end_cm`.

- map:

  A consensus map (`chr`, `bp`, `cm`); defaults to
  [`load_map()`](https://sawers-rellan-labs.github.io/nilhmm/reference/load_map.md).

## Value

`seg` with `start_mb` and `end_mb` added.

## See also

[`cm_to_bp()`](https://sawers-rellan-labs.github.io/nilhmm/reference/map_interpolators.md),
[`bp_to_cm()`](https://sawers-rellan-labs.github.io/nilhmm/reference/map_interpolators.md)

## Examples

``` r
cm_to_mb(data.frame(chr = 1L, start_cm = 0, end_cm = 10))
#>   chr start_cm end_cm  start_mb   end_mb
#> 1   1        0     10 0.0374105 4.346492
```
