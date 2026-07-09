# FSFHap stage 2b: BiparentalHaplotypeFinder parent-calling (non-BC1 route)

The non-backcross route:
[`fsfhap_prefilter_sites_cpp()`](https://sawers-rellan-labs.github.io/nilhmm/reference/fsfhap_prefilter_sites_cpp.md)
→ reconstruct two parental haplotypes
([`fsfhap_biparental_alleles_cpp()`](https://sawers-rellan-labs.github.io/nilhmm/reference/fsfhap_biparental_alleles_cpp.md))
→ recode each genotype into the parent-origin A/het/C frame using the
per-site `alleleA`/`alleleC` → coverage filter. Returns the same shape
as
[`.fsfhap_call_parents_bc()`](https://sawers-rellan-labs.github.io/nilhmm/reference/dot-fsfhap_call_parents_bc.md)
so the shared stage-3 imputation can consume it.

## Usage

``` r
.fsfhap_biparental_call(
  G,
  pos,
  min_maf = 0.05,
  min_coverage = 0.2,
  max_het_deviation = 5,
  min_r2 = 0.2,
  min_gametes = 200L
)
```

## Arguments

- G:

  Integer matrix, taxa x sites, canonical `g` in `{0,1,2,3}`; one
  family, one chromosome, sites sorted.

- pos:

  Integer marker positions (bp), length = ncol(G).

- min_maf, min_coverage, max_het_deviation, min_r2:

  preFilterSites knobs (TASSEL 0.05 / 0.2 / 5 / 0.2).

- min_gametes:

  Coverage floor (TASSEL 200).

## Value

List `G` (recoded, kept taxa x kept sites), `keep_sites`, `keep_taxa`,
`pos`; empty if preFilterSites leaves too few sites or no seed is found.
