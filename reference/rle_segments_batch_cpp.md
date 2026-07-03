# Run-length-encode a batch of state paths (one per column) into segments

Run-length-encode a batch of state paths (one per column) into segments

## Usage

``` r
rle_segments_batch_cpp(paths, pos)
```

## Arguments

- paths:

  T x S integer matrix of state paths (column = sample).

- pos:

  Integer marker positions, length T (shared across samples).

## Value

List of equal-length vectors: sample (1-based column), start_bp, end_bp,
state.
