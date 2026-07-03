# Write segment calls in the common schema

Columns: `source, donor, name, chr, start_bp, end_bp, state` (0/1/2).
`name` is the NIL id (pedigree string) per the doc terminology
convention.

## Usage

``` r
write_common_schema(calls, path)
```

## Arguments

- calls:

  A data.frame of segment calls.

- path:

  Output CSV path.

## Value

`path`, invisibly.

## Examples

``` r
calls <- data.frame(source = "nilHMM", donor = NA, name = "NIL1", chr = 1L,
                    start_bp = 1L, end_bp = 5e6L, state = 2L)
write_common_schema(calls, tempfile(fileext = ".csv"))
```
