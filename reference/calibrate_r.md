# Calibrate the duration hyperparameter `r` by KS-vs-sim

Sweeps `r`, calls on data, and returns the `r` minimising the KS
distance between called donor-block sizes and the simulated-truth
distribution.

## Usage

``` r
calibrate_r(data, sim_segments, r_grid, ...)
```

## Arguments

- data:

  Marker-level input for a taxon/group.

- sim_segments:

  Simulated-truth segments for the matching design (the KS target; see
  [`expected_fragment_dist()`](https://sawers-rellan-labs.github.io/nilhmm/reference/expected_fragment_dist.md)
  and the frozen `sim_truth/` fixtures).

- r_grid:

  Candidate `r` values to sweep.

- ...:

  Forwarded to
  [`call_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md).

## Value

data.frame of `(r, D, ...)` plus the argmin; warns if `at_grid_edge`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Planned (Task 4): pick r by KS-vs-sim minimum distance (rigidity-not-mle).
sweep <- calibrate_r(data, sim_segments,
                     r_grid = c(1e-6, 1e-5, 1e-4),
                     caller = "bbnil", design = "BC2S2")
} # }
```
