# Run-length-encode a state path into (start_bp, end_bp, state) segments

Run-length-encode a state path into (start_bp, end_bp, state) segments

## Usage

``` r
rle_segments_cpp(path, pos)
```

## Arguments

- path:

  Integer state path (0/1/2), length T.

- pos:

  Integer marker positions, length T (same order as path).

## Value

List of equal-length integer vectors: start_bp, end_bp, state.
