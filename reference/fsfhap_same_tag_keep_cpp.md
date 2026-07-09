# FSFHap same-tag SNP filter (whichSnpsAreFromSameTag), faithful port

Within a 64 bp window, consecutive SNPs whose presence-based R^2 to the
window's anchor is \>= `min_rsq` are treated as coming from the same GBS
tag and dropped, keeping only the anchor; a SNP that is \>= 64 bp away,
or has R^2 \< `min_rsq` (or NaN), becomes the next kept anchor. **One
chromosome per call** (positions only; the chromosome-equality guard is
implicit), markers sorted by position.

## Usage

``` r
fsfhap_same_tag_keep_cpp(G, pos, major_is_ref, min_rsq)
```

## Arguments

- G:

  Integer matrix, taxa x sites, canonical `g` in `{0,1,2,3}`.

- pos:

  Integer marker positions (bp), length = ncol(G), sorted ascending.

- major_is_ref:

  Logical, length = ncol(G): TRUE if REF is the major allele at that
  site (so allele presence is defined per the site's own major/minor, as
  TASSEL's `majorAllele`/`minorAllele` do).

- min_rsq:

  R^2 threshold for "same tag" (TASSEL default 0.8).

## Value

Logical vector, length = ncol(G): TRUE = keep (a distinct tag).
