# Expected fragment-size distribution (the Null / KS target)

Returns the fitted Gamma `(k, lambda)` for a design plus density/CDF
closures, for plotting the grey Null and as the KS-vs-sim calibration
target. Never feeds the engine transition.

## Usage

``` r
expected_fragment_dist(design)
```

## Arguments

- design:

  Design key.

## Value

`list(k, lambda, density, cdf)`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Planned (Task 4): the grey Null / KS target for a design.
nd <- expected_fragment_dist("BC2S3")
} # }
```
