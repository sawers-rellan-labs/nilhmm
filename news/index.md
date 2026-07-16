# Changelog

## nilHMM 0.3.0

### Caller family renamed to match the paper grid (breaking)

The named callers are now the explicit coordinates of the engine’s
(emission × duration) grid, matching the paper’s Fig. 1 / caller table:

- **`nnil` is now the categorical `gt` caller** (gt + geometric),
  Holland’s original on hard genotype calls. The old count-emission
  behaviour moved to the new **`bbnil`** caller (count/BetaBinomial +
  geometric) — the low-coverage count extension. **Migration:**
  `call_ancestry(counts, caller = "nnil")` → `caller = "bbnil"`;
  hard-call/`GT` inputs stay on `caller = "nnil"`.
- New **`catiger`** caller: categorical `gt` emission + rigidity
  duration (the `gt`-side counterpart of `rtiger`).
- Removed the `emission = c("count", "gt")` override on
  [`call_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md)/
  [`call_states()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_states.md).
  Each grid caller now pins its own emission, so the override is
  redundant — choose the caller instead.
- The no-HMM per-site **genotype** baseline (the paper’s het-excess
  “control”) is
  [`call_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_gt.md),
  a *genotype* caller — **not** an ancestry caller. It is deliberately
  not dispatchable through
  [`call_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md),
  keeping the terminology wall between the ancestry mosaic and per-site
  genotypes: `call_gt(prior = "flat")` is the maximum-likelihood call,
  `call_gt(prior = "hwe")` the HWE MAP (het-excess). A consumer that
  wants it in an ancestry comparison converts it
  ([`to_segments()`](https://sawers-rellan-labs.github.io/nilhmm/reference/to_segments.md))
  itself; the package provides no genotype→mosaic shortcut.
- [`caller_sweep()`](https://sawers-rellan-labs.github.io/nilhmm/reference/caller_sweep.md)
  renamed its count option `"nnil"` → `"bbnil"` (the swept `rrate` grid
  is the count/geometric caller).
- **GOOGA transcript callers split into `googa` and `atlas`.** The GOOGA
  source and the Flagel 2019 / Veltsos 2024 methods use a
  recombination-fraction (geometric) F2 HMM with no rigidity, so the
  faithful reproduction is now **`googa`** (gt + geometric). The
  **`atlas`** name is retained for this work’s transcript caller — the
  same GOOGA competitive-alignment thresholding decoded with the
  **rigidity** duration. **Migration:** the old `caller = "atlas"` (gt +
  geometric) is now `caller = "googa"`; `caller = "atlas"` now means
  gt + rigidity.

### Pedigree-aware ancestry calling

- New `call_ancestry(caller = "pedigree", pedigree = ...)`: a
  family-coupled caller that couples relatives by loopy belief
  propagation over the (pedigree × genome) grid (kernel
  `src/pedigree_bp.cpp`). It shares the engine’s emission axis — read
  counts give a depth-aware count/BetaBinomial *de novo* call (missing =
  zero depth → flat), while a hard-call `state`/`g` column gives the
  categorical gt path. Requires `design`; groups by `family` from the
  pedigree and dispatches per family (like `fsfhap`). New args:
  `pedigree`, `ped_format`, `ped_maxiter/ped_tol/ped_lambda`. Benchmark:
  on a ZEAL-structured simulated cohort (82 founders, ~17 lines/founder,
  skim depth) the caller reaches **parity** with the single-chain
  callers and fills the fully-uncovered marker grid, but does **not**
  beat `rtiger` — the per-marker transmission factor pools siblings only
  marginally, not over inherited blocks (`design/PEDIGREE_HMM.md` §2). A
  block-coherent accuracy gain awaits the phased copy-switch kernel
  (§18); this is correct/principled infrastructure, not an accuracy
  improvement.
- [`refine_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/refine_ancestry.md)
  is now a thin wrapper over the same shared `.pedigree_states()` kernel
  for the hard-call refinement use — unchanged behaviour, return shape,
  and `emission = c("gt", "count")` modes.

### Genotype calling — clean break (breaking)

- Removed the `call_gl` alias entirely. Use `call_gt`. (Consumers must
  migrate; there is no longer a soft-deprecated fallback.)
- Removed the `prior = "breeding"` + `f` deprecation shim and the `f`
  formal. `prior` now accepts only `"flat"`, `"hwe"` (+ `af`), or a
  length-3 numeric vector `c(f_REF, f_HET, f_ALT)`. Build design priors
  with
  [`design_prior()`](https://sawers-rellan-labs.github.io/nilhmm/reference/design_prior.md).

## nilHMM 0.2.0

### Genotype calling

- `call_gt` replaces `call_gl` as the per-site genotype caller. The old
  name was tied to the GL emission stage shared with every caller, not
  to what it does (call a per-site genotype); `call_gl` remains as a
  soft-deprecated alias.
- `call_gt(prior = )` is now polymorphic — `"flat"` (argmax-GL / ML),
  `"hwe"` (+ `af`), or a fixed `c(f_REF, f_HET, f_ALT)` vector
  (renormalized, covering both design and custom priors). The old
  `prior = "breeding"` + `f` API folds into the numeric-vector path via
  a deprecation shim (still warns); the `f` formal is removed.
- `design_prior(design)` returns the design’s single-locus frequencies
  as the `c(f_REF, f_HET, f_ALT)` vector
  [`call_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_gt.md)
  consumes as its `prior`, so the design prior is derived rather than
  typed (e.g. `call_gt(n_ref, n_alt, prior = design_prior("BC2S3"))`).

### Documentation

- `R CMD check` is clean (Status OK): documented
  [`call_states()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_states.md)’s
  `rtiger_fit` argument, fixed “Lost braces” in the `fsfhap_*` /
  `read_hapmap` man pages, and re-attached the `.fsfhap_biparental_call`
  roxygen block (was landing on the `.FSFHAP_BHF_WINDOW` constant).

## nilHMM 0.1.0

First development release: a unified R + Rcpp reimplementation of the
nilHMM ancestry-calling family for Near-Isogenic Lines and related
backcross / full-sib populations.

### Engine

- A single duration-aware 3-state (REF / HET / ALT) HMM with two
  swappable axes — emission (`emission_count`, `emission_gt`) and
  duration (`duration_geometric`, `duration_rigidity`, `duration_hsmm`)
  — plus breeding-design priors (`design_priors`).
- One-verb API
  [`call_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md)
  (and the coordinate-free
  [`call_states()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_states.md)
  /
  [`to_segments()`](https://sawers-rellan-labs.github.io/nilhmm/reference/to_segments.md)
  pair) returning a common segment schema
  (`source, donor, name, chr, start_bp, end_bp, state`).
- Data-agnostic throughout: functions take `(data, params)`; paths and
  sample lists stay in pipeline scripts.

### Callers

- `nnil` — Holland’s nNIL count/genotype caller (BetaBinomial or
  categorical genotype-error emission).
- `rtiger` — a Julia-free port of the RTIGER rigidity HMM (EM +
  Viterbi + border re-placement).
- `binhmm` — per-bin calling with an anchored 3-state Gaussian-emission
  HMM (GMM / k-means / rebmix backends available).
- `atlas` — a GOOGA-style caller for competitive-alignment RNA-seq
  counts.
- `lbimpute` — a native port of LB-Impute (Fragoso et al. 2014) for very
  low-coverage data: coverage-aware emission + distance-dependent
  transition, bp or cM units;
  [`caller_sweep()`](https://sawers-rellan-labs.github.io/nilhmm/reference/caller_sweep.md)
  calibrates `recombdist` (exact per value).
- `fsfhap` — a native port of TASSEL’s FSFHap (Swarts et al. 2014) for
  full-sib families. Both parent-calling routes (backcross and
  `BiparentalHaplotypeFinder`) feed a 5-state Viterbi-training EM
  imputation + gap-fill. Bit-exact vs TASSEL on backcross and RIL/inbred
  populations. Cross-family/chromosome parallelism via `threads=`
  (bit-identical to serial); on real TeoNAM it runs 4.9x (single family)
  to 7.6x (five families) faster than TASSEL, and faster even
  single-threaded. See `design/FSFHAP_PORT.md`.

### Input / output

- Readers that normalise formats to the observation table: `read_counts`
  (TSV / GATK table / VCF `AD`), `read_vcf_gt`, `read_hapmap` (TASSEL
  HapMap), `read_plink` (`.bed`/`.bim`/`.fam`), `read_pedigree` (FSFHap
  TSV or PLINK `.fam`).
- Writers: `write_common_schema` and `write_vcf_impute` (LB-Impute-style
  imputed VCF).

### Genotype calling & downstream utilities

- `call_gl` — a linkage-free per-site GATK-style genotype-likelihood
  caller.
- `interpolate_genotype` — genotype densification onto a target marker
  grid.
- `pairwise_distance` + `select_independent` — LD-based marker thinning
  (FastIndep port).

### Documentation

- Four vignettes: *Getting started with nilHMM*, *The callers*, *The
  engine*, and *Full-sib families with the fsfhap caller*.
- A pkgdown site at <https://sawers-rellan-labs.github.io/nilhmm/>.
