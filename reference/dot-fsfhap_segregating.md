# FSFHap stage 1a: segregating-site test for one family x chromosome

Thin wrapper over
[`fsfhap_segregating_sites_cpp()`](https://sawers-rellan-labs.github.io/nilhmm/reference/fsfhap_segregating_sites_cpp.md)
— a faithful port of TASSEL's `whichSitesSegregateCorrectly`. Given one
family's genotype matrix on one chromosome, flags the sites whose
minor-allele count fits the family's segregation model (backcross vs F2)
better than a monomorphic/error model.

## Usage

``` r
.fsfhap_segregating(G, max_missing = 0.9, ratio)
```

## Arguments

- G:

  Integer matrix, taxa x sites, values 0/1/2/3 = REF-hom / het / ALT-hom
  / missing (engine canonical `g`); sites sorted by position.

- max_missing:

  Max missing-genotype proportion for a site to be tested (TASSEL
  default 0.9; the TeoNAM run passes `1.0`).

- ratio:

  Expected minor-allele frequency selecting the segregation model:
  `0.25`/`0.75` (backcross) or `0.5` (F2). **Design-derived, no
  default** — it must be supplied by the caller from the breeding design
  (the routing dispatcher, mirroring TASSEL's `contribution1` → route →
  ratio), never guessed or defaulted here.

## Value

A list: `seg` (logical, kept sites), `Mj`/`Mn` (major/minor allele
counts), `p_missing`, and model probabilities
`pmono`/`pquarter`/`phalf`.

## Details

`ratio` selects the branch, matching TASSEL: `0.25`/`0.75` → backcross
(keep iff `pquarter > phalf && pquarter > pmono`), anything else → F2
(`phalf / (pmono + pquarter) > 2`). Note TASSEL's backcross branch
always tests the **0.25** model regardless of backcross depth, so deeper
backcrosses (e.g. a BC2 donor fraction ~0.125) are tested against 0.25
too; the port reproduces this deliberately for parity.
