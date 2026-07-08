# FSFHap stage 3: forward-fill gaps (fillGapsInAlignment)

Per taxon (row), across sites in order: when two non-missing calls
flanking a run of missing (`3`) are EQUAL, fill the run with that value;
a differing non-missing call resets the anchor (no fill). Faithful to
`fillGapsInAlignment`.

## Usage

``` r
fsfhap_fill_gaps_cpp(G)
```

## Arguments

- G:

  Integer matrix, taxa x sites, `0`/`1`/`2`/`3` (`3` = missing).

## Value

`G` with eligible missing runs forward-filled.
