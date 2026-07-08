# FSFHap segregating-site test (whichSitesSegregateCorrectly), faithful port

Per site, keep it if its minor-allele count fits the design's
segregation model better than the alternatives. Backcross (`ratio` 0.25
or 0.75): keep iff `pquarter > phalf && pquarter > pmono`. F2 (`ratio`
0.5): keep iff `phalf / (pmono + pquarter) > 2`. Sites that are
monomorphic (single allele) or exceed `max_missing` are dropped.

## Usage

``` r
fsfhap_segregating_sites_cpp(G, max_missing, ratio)
```

## Arguments

- G:

  Integer matrix, taxa x sites, values 0/1/2/3 = REF-hom / het / ALT-hom
  / missing (engine canonical `g`). One family, one chromosome, sites
  sorted by position.

- max_missing:

  Max missing-genotype proportion for a site to be tested.

- ratio:

  Expected minor-allele frequency: 0.25/0.75 backcross, 0.5 F2.

## Value

List: `seg` (logical, kept sites), `Mj`/`Mn` (major/minor allele
counts), `p_missing`, and the model probs `pmono`/`pquarter`/`phalf` (NA
at untested sites) — for per-stage TASSEL-parity checks.
