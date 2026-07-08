# Read a PLINK binary genotype (.bed/.bim/.fam) into the engine's table

Adapter for PLINK 1 binary genotypes into the canonical
`(name, chr, pos, g)` table. `g` is the **A1-allele dosage** (A1 =
`.bim` column 5, PLINK's minor allele by default): `2` = A1-hom, `1` =
het, `0` = A2-hom, `3` = missing. SNP-major `.bed` only (the PLINK
default).

## Usage

``` r
read_plink(prefix)
```

## Arguments

- prefix:

  Path prefix (the `.bed`, `.bim`, `.fam` trio); a trailing `.bed` is
  stripped if given.

## Value

A long observation table: `name, chr, pos, g`.

## See also

[`read_pedigree()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_pedigree.md)
(`format = "fam"`) for the family grouping.
