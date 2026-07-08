# FSFHap stage 2b: BiparentalHaplotypeFinder.assignHaplotyes (faithful port)

Reconstruct two parental haplotypes across a bidirectional window scan
on the preFiltered genotype matrix, producing per-site parent alleles.
Seed: first window with exactly two clusters (3rd `< minClusterSize`)
whose majority haplotypes differ by `>= 2*window-4`; then extend forward
and backward, matching each window's candidate haplotypes to the running
parents and writing the non-overlap alleles.

## Usage

``` r
fsfhap_biparental_alleles_cpp(Gf)
```

## Arguments

- Gf:

  Integer matrix, taxa x sites, canonical `g` in 0,1,2,3; the ALREADY
  preFiltered chromosome (see
  [`fsfhap_prefilter_sites_cpp()`](https://sawers-rellan-labs.github.io/nilhmm/reference/fsfhap_prefilter_sites_cpp.md)).

## Value

List: `alleleA`, `alleleC` (length ncol(Gf); g-frame `0`=REF-hom /
`2`=ALT-hom / `3`=NN) parent-of-origin alleles per site; `seeded`
(logical).
