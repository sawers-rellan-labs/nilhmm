# Read allelic counts into the engine's observation table

Reads CollectAllelicCounts-style read counts (the GATK-table readcount
standard, memory `gatk-table-readcount-standard`) into a per-sample,
per-marker `(n_ref, n_alt, chr, pos)` table for the `count` emission. A
TSV-\>counts adapter covers the wideseq-thinned BRB inputs (cf. the
Python `agent/brb_nilhmm_counts.py`).

## Usage

``` r
read_counts(path, format = c("tsv", "gatk_table", "vcf_ad"), name = NULL)
```

## Arguments

- path:

  A count TSV file, or a directory of them.

- format:

  One of `"tsv"` (the `chr pos ref n_ref alt n_alt` headerless layout
  used by the skim/BRB counts), `"gatk_table"`, `"vcf_ad"`.

- name:

  Optional sample name for a single file; defaults to the file's
  basename with extensions stripped.

## Value

A long observation table: `name, chr, pos, n_ref, n_alt`.

## Examples

``` r
# Headerless "chr pos ref n_ref alt n_alt" TSV, as produced upstream.
f <- tempfile(fileext = ".tsv")
write.table(
  data.frame(chr = "chr1", pos = c(1e5, 2e5, 3e5),
             ref = "A", n_ref = c(8, 5, 0), alt = "T", n_alt = c(0, 3, 7)),
  f, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
read_counts(f, name = "NIL1")
#>   name chr    pos n_ref n_alt
#> 1 NIL1   1 100000     8     0
#> 2 NIL1   1 200000     5     3
#> 3 NIL1   1 300000     0     7
```
