# Read a TASSEL HapMap into the engine's observation table

Adapter for FSFHap's native input: a TASSEL **HapMap** (11 fixed columns
`rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode`
then one column per taxon) into the canonical `(name, chr, pos, g)`
table, `g` in `{0,1,2,3}` = allele0-hom / het / allele1-hom / missing,
where allele0/allele1 are the two alleles from the per-site `alleles`
field (e.g. `A/C`). Genotype cells may be single-character IUPAC
(`A`/`C`/`M`/`N`) or two-character diploid (`AA`/`AC`/`CC`/`NN`). Feed
to `call_ancestry(..., caller = "fsfhap"|"nnil")` (a `g`-only input).

## Usage

``` r
read_hapmap(path, samples = NULL)
```

## Arguments

- path:

  A HapMap file (`.hmp.txt`, optionally `.gz`).

- samples:

  Optional character vector restricting to a subset of taxa.

## Value

A long observation table: `name, chr, pos, g`.

## See also

[`read_pedigree()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_pedigree.md)
for the matching family/priors,
[`read_vcf_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_vcf_gt.md).

## Examples

``` r
f <- tempfile(fileext = ".hmp.txt")
writeLines(c(
  "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\tL1\tL2",
  "m1\tA/C\t1\t100\t+\tNA\tNA\tNA\tNA\tNA\tNA\tA\tC",
  "m2\tA/C\t1\t200\t+\tNA\tNA\tNA\tNA\tNA\tNA\tM\tN"), f)
read_hapmap(f)          # g in {0 A-hom, 1 het, 2 C-hom, 3 missing}
#>   name chr pos g
#> 1   L1   1 100 0
#> 2   L1   1 200 1
#> 3   L2   1 100 2
#> 4   L2   1 200 3
```
