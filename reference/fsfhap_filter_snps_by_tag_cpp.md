# FSFHap preFilterSites step 1: filterSnpsByTag (faithful port)

Thins SNPs from the same GBS tag and applies per-site quality gates.
Finds the first site passing MAF/missing/het thresholds (the head), then
keeps each later site that (is \>= 64 bp from the head OR has
presence-correlation `< 0.7` to it) AND passes the gates; the head
advances to each kept site. Faithful to TASSEL's
`filterSnpsByTag(a, minMaf, maxMissing, maxHet)`.

## Usage

``` r
fsfhap_filter_snps_by_tag_cpp(G, pos, min_maf, max_missing, max_het)
```

## Arguments

- G:

  Integer matrix, taxa x sites, canonical `g` in 0,1,2,3; one
  chromosome, sites sorted by position.

- pos:

  Integer marker positions (bp), length = ncol(G).

- min_maf, max_missing, max_het:

  Per-site gates (minor-allele freq, missing proportion, heterozygous
  proportion). preFilterSites calls this with
  `(minMaf, 1 - minCoverage, 1.0)`.

## Value

Logical vector, length = ncol(G): TRUE = kept site.
