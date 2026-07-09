# FSFHap preFilterSites (full): filterSnpsByTag -\> het-deviation -\> biallelic -\> LD

Faithful port of `BiparentalHaplotypeFinder.preFilterSites`. Step 1 =
[`fsfhap_filter_snps_by_tag_cpp()`](https://sawers-rellan-labs.github.io/nilhmm/reference/fsfhap_filter_snps_by_tag_cpp.md)
with `(minMaf, 1 - minCoverage, 1.0)`. On the result: (2) drop sites
whose het fraction exceeds `mean + maxHetDeviation * sd` (sample sd over
the filtered sites); (3) drop monomorphic (non-biallelic) sites; (4) if
`minR2 > 0`, drop sites whose average `r^2` over a ±50-filtered-site
window (hets→missing, `calc_rsqr`) is `< minR2` (NaN pairs skipped;
all-NaN average is not a rejection, matching TASSEL's `NaN < minR2` =
false).

## Usage

``` r
fsfhap_prefilter_sites_cpp(
  G,
  pos,
  min_maf,
  min_coverage,
  max_het_deviation,
  min_r2
)
```

## Arguments

- G:

  Integer matrix, taxa x sites, canonical `g` in `{0,1,2,3}`; one
  chromosome.

- pos:

  Integer marker positions (bp), length = ncol(G).

- min_maf, min_coverage, max_het_deviation, min_r2:

  TASSEL fields (0.05 / 0.2 / 5 / 0.2).

## Value

Logical vector, length = ncol(G): TRUE = site kept by preFilterSites.
