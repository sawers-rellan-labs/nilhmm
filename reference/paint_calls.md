# Chromosome-paint called ancestry segments

Paint the common-schema segments as coloured rectangles along each
chromosome — REF / HET / ALT (`state` 0/1/2). Each **sample** (`name`)
is a facet row and each chromosome a facet column. When `track` is
supplied, its levels are stacked as horizontal bands *within* every
cell, so several callers (or data sources) can be overlaid on identical
axes to check whether they agree on where the donor blocks land — the
"sanity paint" comparison. Build the multi-track input by
[`rbind()`](https://rdrr.io/r/base/cbind.html)-ing each caller's
[`call_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md)
output with a column naming the track, then pass that column name as
`track`.

## Usage

``` r
paint_calls(
  calls,
  track = NULL,
  samples = NULL,
  palette = c(REF = "gold", HET = "springgreen4", ALT = "purple4")
)
```

## Arguments

- calls:

  Segment calls in the common schema: columns `name`, `chr`, `start_bp`,
  `end_bp`, `state` (`0`/`1`/`2`). Extra columns are ignored.

- track:

  Optional column name whose levels are stacked as bands within each
  `name` x `chr` cell (e.g. `"method"` for a caller comparison, or
  `"donor"`). `NULL` (default) paints one band per sample. A factor
  `track` keeps its level order (top band = first level).

- samples:

  Optional character vector restricting to a subset of `name` values
  (keeps busy figures legible).

- palette:

  Named length-3 fill for `c(REF, HET, ALT)`.

## Value

A ggplot2 object (requires the suggested ggplot2).

## Examples

``` r
calls <- data.frame(
  name = "NIL01", donor = "B", chr = 1L,
  start_bp = c(1e6, 3e6, 5e6), end_bp = c(3e6, 5e6, 9e6),
  state = c(0L, 2L, 0L))
if (requireNamespace("ggplot2", quietly = TRUE)) paint_calls(calls)
```
