# Top-level ancestry-calling API

Resolves a named caller (and optional design preset) into emission +
duration specs, runs
[`fit()`](https://sawers-rellan-labs.github.io/nilhmm/reference/fit.md)
/
[`decode()`](https://sawers-rellan-labs.github.io/nilhmm/reference/decode.md)
per (sample, chromosome), and returns common-schema segment calls.
Data-agnostic: pass the observations in; this function never reads
paths.

## Usage

``` r
call_ancestry(
  data,
  caller = c("nnil", "rtiger", "binhmm", "atlas"),
  design = NULL,
  r = 0.01,
  err = 0.01,
  conc = 20,
  fit_means = FALSE,
  p_switch = 0.01,
  f_1 = NULL,
  f_2 = NULL,
  source = "nilHMM",
  donor = NA_character_,
  parallel = FALSE,
  threads = 1L,
  seed = 1L,
  postprocess = TRUE,
  emission = NULL,
  bin_size = 1e+06,
  cluster_method = c("gauss", "gmm", "kmeans", "rebmix"),
  joint_clust = FALSE,
  obs_weights = FALSE,
  atlas_thresh = 0.95,
  atlas_het = 0.25,
  atlas_min_reads = 5L,
  germ = 0.05,
  gert = 0.1,
  p = 0.5,
  mr = 0.1,
  nir = 0.01
)
```

## Arguments

- data:

  Long observation table with columns `name, chr, pos` and either
  `n_ref, n_alt` read counts (from
  [`read_counts()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_counts.md);
  for the count/rtiger/binhmm/ atlas paths) or a pre-called
  hard-genotype column `g` in `{0,1,2,3}` (from
  [`read_vcf_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_vcf_gt.md);
  auto-selects `caller = "nnil"`, `emission = "gt"`). Optionally
  `donor`.

- caller:

  One of `"nnil"`, `"rtiger"`, `"binhmm"`, `"atlas"`.

- design:

  Breeding-design key for priors (e.g. `"BC2S2"`, `"BC2S3"`). Required
  unless `f_1`/`f_2` are supplied.

- r, err, conc, fit_means, p_switch:

  Caller parameters forwarded to
  [`caller_spec()`](https://sawers-rellan-labs.github.io/nilhmm/reference/caller_spec.md).
  Explicit formals (not `...`) so that, e.g., `r` is never
  partial-matched to another argument.

- f_1, f_2:

  Single-locus priors, used when `design` is `NULL`.

- source:

  Value for the output `source` column.

- donor:

  Donor/taxon label when `data` has no `donor` column.

- parallel:

  If `TRUE`, decode samples across cores via RcppParallel (only the
  fixed-means batched path; control thread count with
  [`RcppParallel::setThreadOptions()`](https://rdrr.io/pkg/RcppParallel/man/setThreadOptions.html)).
  Identical results to serial.

- threads, seed:

  RTIGER caller only: E-step threads and the seed for its randomized
  init.

- postprocess:

  RTIGER caller only: apply the border re-placement (default TRUE).

- emission:

  Optional emission override (`"count"`, `"gt"`) for the `nnil` caller;
  `NULL` uses the caller's default.

- bin_size, cluster_method:

  `binhmm` caller only: genomic bin width in bp (default 1 Mb) and the
  per-bin genotyping backend — `"gauss"` (default: the anchored 3-state
  Gaussian-emission HMM; fixes the HET over-call and high-coverage
  fragmentation of the K=3 clustering), or the original
  cluster-then-smooth route `"gmm"` / `"kmeans"` / `"rebmix"` (the last
  needs the suggested rebmix package, for bit-exact rpubs reproduction).
  `joint_clust`/`obs_weights` apply only to the cluster backends.

- joint_clust:

  `binhmm` caller only: if `TRUE`, pool all samples' bins and learn one
  shared set of REF/HET/ALT clusters on raw alt-freq (borrows strength
  across the cohort; cf. `get_joint_ancestry_calls.R`); if `FALSE`
  (default), cluster each sample independently on its
  informative-count-weighted alt-freq. This is joint *clustering* only;
  the HMM stays per-sample.

- obs_weights:

  `binhmm` caller only, `joint_clust = TRUE` only: if `TRUE`, weight the
  (gmm) clustering fit by each bin's informative-variant count (weights
  the influence, not the alt-freq value).

- atlas_thresh, atlas_het, atlas_min_reads:

  `atlas` caller only: GOOGA genotype-call thresholds on the donor read
  fraction — homozygous call when a parent's fraction \>= `atlas_thresh`
  (0.95), HET when both parents \>= `atlas_het` (0.25), and a minimum of
  `atlas_min_reads` (5) informative reads per gene (else missing). For
  `atlas`, `n_ref`/`n_alt` are the recurrent/donor competitive-alignment
  read counts (ambiguous excluded upstream).

- germ, gert, p, mr, nir:

  Genotype-error rates for the `gt` (categorical) emission (Holland's
  nNIL model): `germ` error on true homozygotes, `gert` on true
  heterozygotes, `p` fraction of hom errors called het, `mr` missing
  rate, `nir` non-informative-marker rate. Raising `germ`/`gert` toward
  the platform's real error rate is the native cure for
  over-fragmentation (isolated miscalled markers are absorbed as errors
  rather than opening 1-marker segments); calibrate to a clean control,
  don't crank blindly. Used by the gt emission (`emission = "gt"`, a `g`
  genotype input, or the `atlas` caller).

## Value

data.frame in the common schema
(`source, donor, name, chr, start_bp, end_bp, state`).

## Examples

``` r
# Toy single-NIL count table: a REF stretch then a donor (ALT) block.
set.seed(1)
toy <- data.frame(
  name  = "NIL1", chr = 1L, pos = seq_len(40L) * 1e5L,
  n_ref = c(rpois(20, 8), rpois(20, 4)),
  n_alt = c(rpois(20, 0), rpois(20, 4)))
call_ancestry(toy, caller = "nnil", design = "BC2S2", r = 1e-4, err = 0.01)
#>   source donor name chr start_bp  end_bp state
#> 1 nilHMM  <NA> NIL1   1   100000 2000000     0
#> 2 nilHMM  <NA> NIL1   1  2100000 4000000     1
```
