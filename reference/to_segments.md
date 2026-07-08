# Collapse per-unit ancestry states into genomic segments

The coordinate-reinstatement step: run-length-collapses the equal-state
runs within each `(name, chr)` from
[`call_states()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_states.md)
into segments. **Point** units (markers) use `pos` for both boundaries;
**interval** units (bins) carry their own per-unit `start_bp`/`end_bp`
(as the `binhmm` caller emits). This is the only place genomic
coordinates re-enter — the decode itself is coordinate-free.

## Usage

``` r
to_segments(states)
```

## Arguments

- states:

  A per-unit state table from
  [`call_states()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_states.md):
  `name, chr, pos, state` (+ optional `source, donor`; + optional
  `start_bp, end_bp` for interval/bin units).

## Value

data.frame in the common schema
(`source, donor, name, chr, start_bp, end_bp, state`).

## Examples

``` r
st <- data.frame(name = "NIL1", chr = 1L, pos = (1:6) * 1e5L,
                 state = c(0L, 0L, 2L, 2L, 0L, 0L))
to_segments(st)
#>   source donor name chr start_bp end_bp state
#> 1 nilHMM  <NA> NIL1   1   100000 200000     0
#> 2 nilHMM  <NA> NIL1   1   300000 400000     2
#> 3 nilHMM  <NA> NIL1   1   500000 600000     0
```
