# Read a pedigree (FSFHap TSV or PLINK .fam) into a family/priors table

The IO adapter for FSFHap family grouping + design priors (keeps the
caller data-agnostic: the algorithm never reads files). Returns one row
per taxon.

## Usage

``` r
read_pedigree(path, format = c("fsfhap", "fam"))
```

## Arguments

- path:

  A pedigree file.

- format:

  `"fsfhap"` (default) or `"fam"`.

## Value

data.frame
`taxon, family, parent1, parent2, contribution1, contribution2, F`.

## Details

- `format = "fsfhap"`: TASSEL's 7-column, header-first FSFHap pedigree
  (`family taxon parent1 parent2 contribution1 contribution2 F`).

- `format = "fam"`: PLINK `.fam` (`FID IID PID MID sex pheno`) -\>
  `family` = FID, `taxon` = IID; contribution/`F` are `NA` (not in
  `.fam`).

Use it to attach the `family` key and derive the design for
`caller = "fsfhap"`:
`ped <- read_pedigree(p); data$family <- ped$family[match(data$name, ped$taxon)]`;
then pass `design`/`phet` (phet = `(1 - F)/2`, see
[`.fsfhap_phet()`](https://sawers-rellan-labs.github.io/nilhmm/reference/dot-fsfhap_phet.md)).

## Examples

``` r
p <- tempfile()
writeLines(c("family\ttaxon\tparent1\tparent2\tp1\tp2\tF",
             "SIM\tL1\tA\tC\t0.75\t0.25\t0.9375"), p)
read_pedigree(p)
#>   taxon family parent1 parent2 contribution1 contribution2      F
#> 1    L1    SIM       A       C          0.75          0.25 0.9375
```
