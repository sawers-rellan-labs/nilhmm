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

  A count TSV file (or directory of them) for `format = "tsv"`; a single
  `.vcf`/`.vcf.gz` for `format = "vcf_ad"`.

- format:

  One of `"tsv"` (the `chr pos ref n_ref alt n_alt` headerless layout
  used by the skim/BRB counts), `"vcf_ad"` (per-sample allelic depths
  from a biallelic VCF's `AD` FORMAT field – e.g. the LB-Impute example
  data and GATK/bcftools output), or `"gatk_table"` (not yet
  implemented).

- name:

  Optional sample name for a single `tsv` file; defaults to the file's
  basename with extensions stripped. Ignored for `vcf_ad` (sample names
  come from the VCF `#CHROM` header).

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
#> Warning: Attempt to override column 2 of inherent type 'float64' down to 'int32' ignored. Only overrides to a higher type are currently supported. If this was intended, please coerce to the lower type afterwards.
#>   name chr    pos n_ref n_alt
#> 1 NIL1   1 100000     8     0
#> 2 NIL1   1 200000     5     3
#> 3 NIL1   1 300000     0     7

# Biallelic VCF with a GT:AD FORMAT: AD = "n_ref,n_alt" per sample.
v <- tempfile(fileext = ".vcf")
writeLines(c(
  "##fileformat=VCFv4.2",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNIL1\tNIL2",
  "1\t100000\t.\tA\tT\t.\t.\t.\tGT:AD\t0/0:8,0\t0/1:4,5",
  "1\t200000\t.\tC\tG\t.\t.\t.\tGT:AD\t1/1:0,7\t./.:."), v)
read_counts(v, format = "vcf_ad")
#>   name chr    pos n_ref n_alt
#> 1 NIL1   1 100000     8     0
#> 2 NIL1   1 200000     0     7
#> 3 NIL2   1 100000     4     5
#> 4 NIL2   1 200000     0     0
```
