# CLAUDE.md

Guidance for Claude Code (claude.ai/code) working in this repository.

## What this is

**nilHMM** (`Package: nilHMM`) is an **R + Rcpp** package that calls ancestry
(donor introgressions) in Near-Isogenic Lines from sequencing data. It is a single
**duration-aware 3-state (REF / HET / ALT) HMM engine** with two swappable axes,
whose combinations express a family of named callers.

> IMPORTANT: This package is **R + Rcpp**, not Python. It began as a Python package
> (Jim Holland's [nNIL](https://github.com/ncsumaize/nNIL) methodology); that
> original is **frozen at `legacy/python/`** (tag `python-final`) and is not the
> live code. Ignore any Python-era description — the authoritative overview is
> `README.md` and the design-of-record docs under `design/`.

## The one-verb API

```r
call_ancestry(data, caller = ..., design = ..., r = ..., err = ...)
```

Everything else is a building block. The package is **data-agnostic**: every
function takes `(data, params)` and returns calls — **no hardcoded paths, sample
lists, or mounts**. Pipeline scripts (in the companion `zealhmm` repo, not here)
own file paths, sample lists, and output locations.

Output is the **common segment schema**:
`source, donor, name, chr, start_bp, end_bp, state` with `state ∈ {REF, HET, ALT}`.

## Architecture: one engine, two axes

Numeric states: **0 = REF** (recurrent hom, e.g. B73), **1 = HET**, **2 = ALT**
(donor hom, e.g. teosinte).

- **Emission** — `emission_count()` (BetaBinomial over ref/alt read depths;
  `err`, `conc`, optional `fit_means` EM) or `emission_gt()` (categorical genotype
  model with an explicit genotyping-error model: `germ, gert, p, mr, nir`).
- **Duration** — `duration_geometric()`, `duration_rigidity()` (minimum-run-length
  prior), or `duration_hsmm()`.
- **Design priors** — `design_priors()` sets the state frequencies (`f_1`, `f_2`)
  from a breeding design (e.g. `"BC2S2"`).
- **Engine** — `fit()` (EM / parameter fitting) then `decode()` (Viterbi /
  posteriors) → `to_segments()` (RLE to the common schema).

### The callers (all share the REF/HET/ALT chain + design priors)

| caller | emission | duration | typical input | lineage |
|---|---|---|---|---|
| `nnil` | `count` or `gt` | geometric | allelic read counts (skim/BrB) or called `GT` | Holland nNIL |
| `rtiger` | `count` | rigidity | allelic read counts | RTIGER (Julia-free port) |
| `binhmm` | anchored Gaussian on binned alt-freq | per-bin HMM | allelic read counts | "ancestry by bins" |
| `atlas` | `gt` (GOOGA thresholds) | geometric | competitive-alignment RNA-seq counts | GOOGA |

See `README.md` for per-caller detail and `?call_ancestry`.

### Downstream utilities (not callers)

Data-agnostic helpers for the companion `zealhmm` QTL-mapping pipeline; they take
`(data, params)` like everything else and hold no paths:

- `interpolate_genotype()` — densify a complete genotype block onto a target
  marker grid by flanking-marker interpolation in genetic distance (cM); modes
  `continuous` (Tian 2011), `step` (Chen/TeoNAM), `round`. Deterministic hard-call
  densification, **not** ancestry inference.
- `pairwise_distance()` + `select_independent()` — LD-based marker thinning per
  chromosome: build a relatedness matrix (`r2` / `mi` / `vi`; similarity vs
  distance carried on `attr(, "kind")`) and select a maximal independent set
  (FastIndep port; the greedy set is bit-identical to the FastIndep CLI). A
  `max_markers` guard (default `7000`, ceiling via option
  `nilHMM.marker_hard_cap`) refuses oversized O(n²) matrices before allocation.

## Repository layout

- `R/` — R layer. Entry `engine.R` (`call_ancestry`, `fit`, `decode`,
  `to_segments`); `callers.R` (`caller_spec` — per-caller definitions);
  `emissions.R` (`emission_count`/`emission_gt`); `duration.R` (`duration_*`);
  `presets_design.R` (`design_priors`, `cm_to_mb`) / `presets_regime.R`
  (`select_emission`); `io.R` (`read_counts`, `read_vcf_gt`); `calibrate.R`
  (`calibrate_r`); `sweep.R` (`caller_sweep`); `rtiger.R`, `binhmm.R`, `atlas.R`;
  `map.R` (`load_map`, stub); `interpolate_genotype.R` (genotype densification);
  `pairwise_distance.R` / `select_independent.R` (LD marker thinning);
  `plot.R`, `nilHMM-package.R`, `RcppExports.R`.
- `src/` — Rcpp engine: `emission_count.cpp`, `forward_backward.cpp`, `viterbi.cpp`,
  `viterbi_batch.cpp`, `viterbi_batch_par.cpp` (RcppParallel/TBB), `segments.cpp`
  (RLE), `rtiger.cpp`, `interpolate_genotype.cpp`, `fast_indep.cpp` (FastIndep port),
  `pairwise_distance.cpp` (r2/MI/VI), plus generated `RcppExports.cpp` and `Makevars`.
- `design/` — **design of record**: `architecture.md`, `Implementation.md`,
  `REFACTOR_R_PACKAGE.md`, `RTIGER_PORT.md`, `VALIDATION.md`, `package_structure.md`,
  `Zv_RTIGER_divergence.md`, `BRB_run_findings.md`. Read these before non-trivial changes.
- `tests/testthat/` + `tests/fixtures/baseline_pre_refactor/` — regression baselines
  (SHA256-summed CSVs from the pre-refactor code; callers must reproduce them).
- `calibration/`, `data-raw/`, `inst/`, `man/` (roxygen-generated), `legacy/` (frozen Python).

## C++ / build workflow

- Rcpp bindings are generated — **never hand-edit** `src/RcppExports.cpp` or
  `R/RcppExports.R`. After adding/changing a `// [[Rcpp::export]]` function, run
  `Rcpp::compileAttributes()` then `devtools::document()`.
- C++ style (see `src/segments.cpp`): `#include <Rcpp.h>`, `using namespace Rcpp;`,
  `//'` roxygen with `@param`/`@return`/`@keywords internal`, `// [[Rcpp::export]]`,
  `*_cpp` naming. Process **one chromosome at a time**; markers must be sorted.
- `Makevars` links the TBB runtime via `RcppParallel::RcppParallelLibs()`
  (`SystemRequirements: GNU make`). `nilHMM-package.R` imports RcppParallel first so
  the bundled TBB loads before nilHMM's DLL (otherwise `dlopen` fails on clean install).
- Typical dev loop: `devtools::load_all()` → `devtools::test()` →
  `devtools::document()` → `R CMD check`. `Rcpp::compileAttributes()` after C++ changes.

## Design constraints

- **Chromosome-separate processing** (linkage assumptions); markers **sorted by
  chr, position**.
- **State coding** 0/1/2 = REF/HET/ALT; missing observation handled by the emission.
- **Data-agnostic** — keep paths/sample lists out of package functions.
- `rtiger` is a **Julia-free** port (no Julia needed); `rebmix` is an optional
  `Suggests` (bit-exact `binhmm` clustering reproduction only).
- `load_map()` is currently a **stub** (Task 4) — the bundled B73 v5 consensus map
  is not yet wired in.

## Related

Companion analysis repo `zealhmm` installs nilHMM as a dependency and hosts the
caller comparison and methods paper. Cite the nNIL population paper (Zhong et al.
2025, bioRxiv 10.1101/2025.01.29.635337) per `CITATION.cff`.
