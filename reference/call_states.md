# Coordinate-free ancestry decode (per-observation states)

Resolves a named caller (and optional design preset) into emission +
duration specs and runs
[`fit()`](https://sawers-rellan-labs.github.io/nilhmm/reference/fit.md)
/
[`decode()`](https://sawers-rellan-labs.github.io/nilhmm/reference/decode.md)
per (sample, chromosome), returning one ancestry **state per observation
unit** — coordinate-free: no genomic intervals. Feed the result to
[`to_segments()`](https://sawers-rellan-labs.github.io/nilhmm/reference/to_segments.md)
to reinstate coordinates and collapse into common-schema segments;
[`call_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md)
chains the two. Data-agnostic: pass the observations in; this function
never reads paths.

## Usage

``` r
call_states(
  data,
  caller = c("nnil", "rtiger", "binhmm", "atlas", "lbimpute", "fsfhap"),
  family = NULL,
  phet = NULL,
  design = NULL,
  rrate = 0.01,
  rigidity = NULL,
  err = 0.01,
  conc = 20,
  fit_means = FALSE,
  xrate = 0.01,
  f_1 = NULL,
  f_2 = NULL,
  source = "nilHMM",
  donor = NA_character_,
  parallel = FALSE,
  threads = 1L,
  seed = 1L,
  postprocess = TRUE,
  min_cov = 1L,
  emission = NULL,
  bin_size = 1e+06,
  cluster_method = c("gauss", "gmm", "kmeans", "rebmix"),
  joint_clust = FALSE,
  obs_weights = FALSE,
  atlas_thresh = 0.95,
  atlas_het = 0.25,
  atlas_min_reads = 5L,
  genotypeerr = 0.05,
  recombdist = NULL,
  drp = FALSE,
  unit = c("bp", "cm"),
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

  One of `"nnil"`, `"rtiger"`, `"binhmm"`, `"atlas"`, `"lbimpute"`,
  `"fsfhap"`.

- family:

  `fsfhap` caller only: per-row family grouping (a vector length
  `nrow(data)`, or supply a `family` column on `data`). FSFHap pools
  each family.

- phet:

  `fsfhap` caller only: expected heterozygosity. If `NULL`, derived from
  `design` as `(1 - F)/2` (BC1S4 → 0.03125); the `fsfhap` port supports
  **BC1 designs** (`"BC1S{m}"`) — other designs route to unported TASSEL
  paths.

- design:

  Breeding-design key for priors (e.g. `"BC2S2"`, `"BC2S3"`). Required
  unless `f_1`/`f_2` are supplied.

- rrate:

  Count/geometric callers (`nnil`, `atlas`): expected per-marker
  **recombination rate** for the geometric duration (self-stay =
  `1 - rrate`). A resolution hyperparameter, not an MLE. Holland's nNIL
  sets it to `2 * total_cM / (100 * n_markers)` (~`30 / n_markers` for a
  1500 cM maize map; the factor 2 is the RIL-like doubling for the
  backcross + self meioses, Haldane & Waddington), and it is optimizable
  from data by KS-vs-sim
  ([`calibrate_r()`](https://sawers-rellan-labs.github.io/nilhmm/reference/calibrate_r.md)).

- rigidity:

  `rtiger` caller only: integer minimum run length (e.g. `5`).

- err, conc, fit_means:

  Count-emission parameters forwarded to
  [`caller_spec()`](https://sawers-rellan-labs.github.io/nilhmm/reference/caller_spec.md).

- xrate:

  Exit rate of nilHMM's **rigidity duration**
  ([`duration_rigidity()`](https://sawers-rellan-labs.github.io/nilhmm/reference/duration_rigidity.md)):
  the per-marker switch probability at the free (post-minimum-run) state
  — the geometric tail beyond the enforced run. A nilHMM construct,
  **not** a RTIGER parameter; the faithful `caller = "rtiger"` port uses
  only `rigidity` (plus `threads`/`seed`/`postprocess`) and ignores
  `xrate`.

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

- threads:

  Fan-out width for coarse-grained parallelism. `rtiger`: E-step
  threads. `fsfhap`: family x chromosome units decoded in parallel
  ([`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html),
  unix only; identical results to `threads = 1L`). `lbimpute`:
  per-(name, chr) decode fan-out. Other callers ignore it (use
  `parallel` for the batched count path).

- seed:

  RTIGER caller only: the seed for its randomized init.

- postprocess:

  RTIGER caller only: apply the border re-placement (default TRUE).

- min_cov:

  Drop no-coverage units before decoding (default `1L`; `0L` keeps
  everything, the old behaviour). "No coverage" is caller-specific but
  the intent is uniform — never make a confident call from no data, and
  keep all callers on the same support for comparability:

  - `nnil` (count) and `rtiger`: markers with `n_ref + n_alt < min_cov`;

  - `nnil` gt path: missing genotypes (`g == 3`);

  - `binhmm`: bins with fewer than `min_cov` informative markers
    (`ninf`). Zero-coverage units carry no emission signal — they only
    slow decoding, marginally inflate `nnil` fragmentation, and dilute
    `rtiger`'s rigidity run. (The `atlas` caller has its own
    `atlas_min_reads` gate and is unaffected.)

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

- genotypeerr, recombdist, drp, unit:

  `lbimpute` caller only (Fragoso et al. 2014): `genotypeerr` is the
  coverage-independent genotyping-error rate (LB-Impute `genotypeerr`,
  default 0.05) that bounds every emission to
  `[genotypeerr, 1 - genotypeerr]`; `err` doubles as LB-Impute's
  per-read `readerr` (its default is 0.05, higher than the shared
  `err = 0.01`). `unit` chooses the coordinate the transition decays
  over: `"bp"` (default, the faithful LB-Impute model – a uniform
  genome-wide recombination rate over physical distance) or `"cm"`
  (map-aware – uses a `cm` column of genetic positions per marker, so
  the local recombination rate, e.g. maize centromeric suppression, is
  captured). Output coordinates are always bp (`pos`); cM only feeds the
  transition. `recombdist` is the coordinate distance over which the
  recombination probability equalizes (LB-Impute `recombdist`); it
  shares units with the coordinate, so its default is unit-aware (`1e7`
  bp = ~50 cM in maize, or `50` cM) and larger values mean stiffer paths
  – a validation layer warns on a bp/cM `recombdist` mismatch. `drp`
  (LB-Impute `-dr`): when `TRUE`, a homozygous-\>homozygous switch is
  priced as a single recombination rather than a double event (use for
  inbred / RIL populations). For `lbimpute`, `design`/`f_1,f_2` only
  seed the start distribution (flat if absent) and zero-coverage markers
  are kept so the transition sees true marker spacing.

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

data.frame of per-unit states (`source, donor, name, chr, pos, state`;
the `binhmm` caller additionally carries per-bin `start_bp, end_bp`).
Pass to
[`to_segments()`](https://sawers-rellan-labs.github.io/nilhmm/reference/to_segments.md)
for intervals.

## Examples

``` r
# Toy single-NIL count table: a REF stretch then a donor (ALT) block.
set.seed(1)
toy <- data.frame(
  name  = "NIL1", chr = 1L, pos = seq_len(40L) * 1e5L,
  n_ref = c(rpois(20, 8), rpois(20, 4)),
  n_alt = c(rpois(20, 0), rpois(20, 4)))
st <- call_states(toy, caller = "nnil", design = "BC2S2", rrate = 1e-4, err = 0.01)
to_segments(st)                       # or, in one step:
#>   source donor name chr start_bp  end_bp state
#> 1 nilHMM  <NA> NIL1   1   100000 2000000     0
#> 2 nilHMM  <NA> NIL1   1  2100000 4000000     1
call_ancestry(toy, caller = "nnil", design = "BC2S2", rrate = 1e-4, err = 0.01)
#>   source donor name chr start_bp  end_bp state
#> 1 nilHMM  <NA> NIL1   1   100000 2000000     0
#> 2 nilHMM  <NA> NIL1   1  2100000 4000000     1
```
