# Simulate a tracked forest of BCnSm sibling families (shared ancestry)

Unlike
[`simulate_nil()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_nil.md)
(independent lines with no relatives), this draws ONE simcross meiosis
per family and reads off `sibs` terminal siblings from that single draw,
so the sibs genuinely share the IBD blocks inherited through the
(latent, ungenotyped) founder -\> selfing chain. This is the
ground-truth generator for validating a pedigree-aware refinement
([`refine_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/refine_ancestry.md)):
the coupling refine exploits only exists when relatives share ancestry.

## Usage

``` r
simulate_family(
  design = "BC2S3",
  families = 10L,
  sibs = 10L,
  map = NULL,
  n_markers = 1000L,
  chr = NULL,
  m = 10L,
  p = 0,
  seed = NULL,
  donor = "B",
  prefix = "fam"
)
```

## Arguments

- design:

  "BCnSm" (default "BC2S3"); needs \>= 1 selfing generation.

- families:

  Number of independent families (independent founders).

- sibs:

  Genotyped siblings per family (they share the founder + chain).

- map, n_markers, chr, m, p, seed, donor:

  As in
  [`simulate_nil()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_nil.md).

- prefix:

  Family-id prefix for taxon/family names (default "fam").

## Value

`list(truth = <per-marker table: source, donor, name, family, chr, pos, cm, state>, pedigree = <data.frame taxon, family, parent1, parent2>)`.

## Details

Returns both the per-marker `truth` for the genotyped sibs (a
[`simulate_nil()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_nil.md)
table, so it feeds
[`simulate_counts()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_counts.md)
/
[`to_segments()`](https://sawers-rellan-labs.github.io/nilhmm/reference/to_segments.md))
and a
[`read_pedigree()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_pedigree.md)-shaped
`pedigree`: one row per taxon, with the founder and the intermediate
selfing nodes present as **parent-only latent rows** (no genotype), the
sibs as leaves. Family `f` is a star: founder `g0` -\> `g1` -\> ... -\>
`g{m-1}` -\> `sibs` leaves; families are independent (independent
founders).

## See also

[`simulate_nil()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_nil.md),
[`simulate_counts()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_counts.md),
[`read_pedigree()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_pedigree.md),
[`refine_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/refine_ancestry.md).

## Examples

``` r
if (requireNamespace("simcross", quietly = TRUE)) {
  fam <- simulate_family("BC2S3", families = 2, sibs = 3, chr = 1,
                         n_markers = 100, seed = 1)
  head(fam$pedigree)
  table(fam$truth$name)
}
#> 
#> fam01_L01 fam01_L02 fam01_L03 fam02_L01 fam02_L02 fam02_L03 
#>       100       100       100       100       100       100 
```
