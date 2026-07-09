# FSFHap stage 1b: backcross parent-allele calling for one family x chromosome

Faithful orchestration of TASSEL's
`callParentAllelesByWindowForBackcrosses`: keep sites that both
segregate as backcross (stage 1a) and survive the same-tag filter
([`fsfhap_same_tag_keep_cpp()`](https://sawers-rellan-labs.github.io/nilhmm/reference/fsfhap_same_tag_keep_cpp.md));
assign A = major / C = minor per kept site; drop low-coverage taxa
(`<= min_gametes` non-missing gametes across kept sites); recode
genotypes into the parent-origin A/C frame (`0` = A-hom / `1` = het /
`2` = C-hom / `3` = missing), where A is the site's major allele. The
result is the observation matrix the 5-state smoother consumes next.

## Usage

``` r
.fsfhap_call_parents_bc(
  G,
  pos,
  max_missing = 1,
  min_rsq = 0.8,
  min_gametes = 200L,
  min_r = 0
)
```

## Arguments

- G:

  Integer matrix, taxa x sites, canonical `g` in `{0,1,2,3}`; one
  family, one chromosome, sites sorted by position.

- pos:

  Integer marker positions (bp), length = ncol(G).

- max_missing:

  Max missing proportion for the segregating-site test (TASSEL default
  0.9; TeoNAM run 1.0).

- min_rsq:

  Same-tag R^2 threshold (TASSEL 0.8).

- min_gametes:

  Coverage floor: taxa with `> min_gametes` non-missing gametes (2 x
  non-missing sites) are kept (TASSEL 200).

- min_r:

  LD-filter threshold; `> 0` enables TASSEL's `ldfilter`, **not yet
  ported** (warns and ignores). The TeoNAM run uses `0`.

## Value

List: `G` (recoded, kept taxa x kept sites), `keep_sites`/`keep_taxa`
(indices into the input), `pos` (kept-site positions), `major_is_ref`
(per kept site), and counts `n_seg`/`n_sametag`/`n_kept_sites`.
