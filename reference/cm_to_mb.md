# Convert segment sizes from cM to Mb using a map

cM-space Gammas are assembly-robust; anything bp/Mb is tied to B73 v5.

## Usage

``` r
cm_to_mb(seg, map)
```

## Arguments

- seg:

  Segments with cM coordinates.

- map:

  A consensus map (see
  [`load_map()`](https://sawers-rellan-labs.github.io/nilhmm/reference/load_map.md)).

## Value

`seg` with Mb coordinates added.

## Examples

``` r
if (FALSE) { # \dontrun{
# Planned (Task 4): project cM segment coordinates to Mb via a consensus map.
cm_to_mb(seg, load_map("v5"))
} # }
```
