# Plot called fragment sizes against the expected Null

Plot called fragment sizes against the expected Null

## Usage

``` r
plot_fragment_sizes(calls, design, space = c("Mb", "cM"))
```

## Arguments

- calls:

  Segment calls (common schema).

- design:

  Design key for the grey Null (see
  [`expected_fragment_dist()`](https://sawers-rellan-labs.github.io/nilhmm/reference/expected_fragment_dist.md)).

- space:

  One of `"cM"`, `"Mb"`.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Planned (Task 4): called fragment sizes vs the expected Null.
calls <- call_ancestry(read_counts("counts/"), caller = "nnil", design = "BC2S2")
plot_fragment_sizes(calls, design = "BC2S2", space = "Mb")
} # }
```
