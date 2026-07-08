# B73 v5 maize consensus genetic map

A consensus genetic map for maize: each marker with its chromosome,
genetic (cM) and physical (bp, B73 v5 assembly) position. cM-space is
assembly-robust; the `bp` column is tied to B73 v5 (see
[`load_map()`](https://sawers-rellan-labs.github.io/nilhmm/reference/load_map.md),
which warns on an assembly override). Used as the default map for
[`simulate_nil()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_nil.md)
and for cM\<-\>Mb position conversion.

## Usage

``` r
maize_map_v5
```

## Format

A data.frame with 19,486 rows (loci) and 4 columns, plus an `"assembly"`
attribute (`"B73v5"`):

- locus:

  marker id (character)

- chr:

  chromosome, 1-10 (integer)

- cm:

  genetic position within the chromosome (cM)

- bp:

  physical position, B73 v5 (bp)

## Source

Cleaned B73 v5 consensus map assembled in the companion zealhmm
pipeline; regenerate with `data-raw/make_maize_map.R`.

## See also

[`load_map()`](https://sawers-rellan-labs.github.io/nilhmm/reference/load_map.md),
[`simulate_nil()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_nil.md)
