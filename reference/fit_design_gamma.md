# Fit a design Gamma from simulated segments

Extensible entry point used by `data-raw/make_breeding_designs.R` to
derive the bundled `(k, lambda)` from simcross output.

## Usage

``` r
fit_design_gamma(sim_segments)
```

## Arguments

- sim_segments:

  data.frame of simulated segment sizes (cM).

## Value

`list(k, lambda)`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Planned (Task 4): fit the design Gamma from simcross segment sizes.
fit_design_gamma(sim_segments)
} # }
```
