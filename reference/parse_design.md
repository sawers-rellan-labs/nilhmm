# Parse a `"BCnSm"` breeding-design string into its generation counts

Parse a `"BCnSm"` breeding-design string into its generation counts

## Usage

``` r
parse_design(design)
```

## Arguments

- design:

  Design key of the form `"BC<n>S<m>"` (e.g. `"BC2S2"`, `"BC1S4"`).

## Value

`list(n_bc, n_self)` – backcross and selfing generation counts.

## Examples

``` r
parse_design("BC2S2")   # list(n_bc = 2, n_self = 2)
#> $n_bc
#> [1] 2
#> 
#> $n_self
#> [1] 2
#> 
```
