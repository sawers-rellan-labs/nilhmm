# Top-level ancestry-calling API (decode + segment)

Convenience wrapper that chains the two stages: coordinate-free decode
([`call_states()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_states.md))
followed by coordinate reinstatement
([`to_segments()`](https://sawers-rellan-labs.github.io/nilhmm/reference/to_segments.md)).
`call_ancestry(...)` is exactly `to_segments(call_states(...))`. See
[`call_states()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_states.md)
for the full argument list.

## Usage

``` r
call_ancestry(data, ...)
```

## Arguments

- data, ...:

  Passed to
  [`call_states()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_states.md).

## Value

data.frame in the common schema
(`source, donor, name, chr, start_bp, end_bp, state`).

## Examples

``` r
toy <- data.frame(name = "NIL1", chr = 1L, pos = (1:6) * 1e5L,
                  n_ref = c(9, 8, 4, 5, 9, 8), n_alt = c(0, 0, 4, 5, 0, 0))
call_ancestry(toy, caller = "bbnil", design = "BC2S2", rrate = 1e-4, err = 0.01)
#>   source donor name chr start_bp end_bp state
#> 1 nilHMM  <NA> NIL1   1   100000 600000     0
```
