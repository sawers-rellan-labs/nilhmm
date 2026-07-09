# FSFHap caller driver: per-family x chromosome parent-calling + imputation

The `caller = "fsfhap"` path for
[`call_states()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_states.md).
Unlike the per-sample callers, FSFHap **pools each family**: for every
`family x chr` it builds the taxa x sites genotype matrix, runs BC
parent-calling
([`.fsfhap_call_parents_bc()`](https://sawers-rellan-labs.github.io/nilhmm/reference/dot-fsfhap_call_parents_bc.md),
stage 1) then the 5-state EM imputation + gap-fill
([`.fsfhap_impute()`](https://sawers-rellan-labs.github.io/nilhmm/reference/dot-fsfhap_impute.md),
stage 3), and emits per-marker states in the common REF/HET/ALT frame
(major allele = the recurrent parent in a backcross → `0` REF, `1` HET,
`2` ALT/donor).

## Usage

``` r
.fsfhap_states(
  data,
  phet,
  source = "nilHMM",
  donor = NA_character_,
  has_donor = FALSE,
  route = "bc",
  min_gametes = 200L,
  max_missing = 1,
  max_iter = 50L,
  threads = 1L
)
```

## Arguments

- data:

  Long table with `name, chr, pos, g` (`g` in `{0,1,2,3}`) and a
  `family` grouping column; optional `donor`.

- phet:

  Design-derived expected heterozygosity
  ([`.fsfhap_phet()`](https://sawers-rellan-labs.github.io/nilhmm/reference/dot-fsfhap_phet.md)).

- source, donor:

  Output `source`/`donor` labels (donor per-`name` if `has_donor`).

- has_donor:

  Whether `data` carries a per-row `donor` column.

- route:

  `"bc"` (backcross, BC1) or `"biparental"` (BiparentalHaplotypeFinder,
  non-BC1) — the design-selected parent-calling route.

- min_gametes, max_missing, max_iter:

  Stage-1/3 knobs (TASSEL 200 / 1.0 / 50).

- threads:

  Fan-out width for the family x chromosome units. Each unit is fully
  independent (chromosome-separate processing, family-pooled), so it
  maps cleanly onto
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  (unix only; serial [`lapply()`](https://rdrr.io/r/base/lapply.html)
  elsewhere). This is FSFHap's coarse-grained parallelism — the
  single-threaded per-unit stages 1+3 are untouched, so results are
  identical to `threads = 1L`.

## Value

`data.frame(source, donor, name, chr, pos, state)`, `state` in
`{0,1,2}`; uncalled markers (dropped by the filters, or still missing
after gap-fill) are simply absent. Feed to
[`to_segments()`](https://sawers-rellan-labs.github.io/nilhmm/reference/to_segments.md).
