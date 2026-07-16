# Read hard genotype calls (VCF GT) into the engine's observation table

Extracts the per-sample `GT` field of a **biallelic diploid** VCF into a
`(name, chr, pos, g)` table for the categorical `gt` emission –
Holland's nNIL genotype path, for the saturated-depth regime (e.g.
MolBreeding target sequencing at ~20x or more) where the caller's
*called* genotype is trustworthy. `g` is the alt-allele dosage: `0`
REF-hom, `1` het, `2` ALT-hom, `3` missing. Unlike
[`read_counts()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_counts.md)
(which feeds the count/BetaBinomial emission from `AD` read depths),
this reads the called genotype directly, so it never touches `AD`.

## Usage

``` r
read_vcf_gt(path, samples = NULL)
```

## Arguments

- path:

  A `.vcf` or `.vcf.gz` file. `GT` must be the first FORMAT field (the
  VCF spec requirement when GT is present).

- samples:

  Optional character vector to restrict to a subset of samples.

## Value

A long observation table: `name, chr, pos, g`.

## Details

Run with
`call_ancestry(read_vcf_gt(path), caller = "nnil", design = ...)`
(`nnil`/`catiger` are the categorical gt callers; a `g`-only input feeds
their gt emission directly).

## Examples

``` r
# A minimal biallelic-diploid VCF with a GT field.
f <- tempfile(fileext = ".vcf")
writeLines(c(
  "##fileformat=VCFv4.2",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNIL1",
  "1\t100000\t.\tA\tT\t.\t.\t.\tGT\t0/0",
  "1\t200000\t.\tA\tT\t.\t.\t.\tGT\t0/1",
  "1\t300000\t.\tA\tT\t.\t.\t.\tGT\t1/1"), f)
read_vcf_gt(f)   # g in {0 REF-hom, 1 het, 2 ALT-hom, 3 missing}
#>   name chr    pos g
#> 1 NIL1   1 100000 0
#> 2 NIL1   1 200000 1
#> 3 NIL1   1 300000 2
```
