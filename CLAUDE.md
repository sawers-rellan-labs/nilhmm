# CLAUDE.md

Guidance for Claude Code (claude.ai/code) working in this repository.

## What this is

**nilHMM** (`Package: nilHMM`) is an **R + Rcpp** package that calls ancestry
(donor introgressions) in Near-Isogenic Lines from sequencing data. It is a single
**duration-aware 3-state (REF / HET / ALT) HMM engine** with two swappable axes,
whose combinations express a family of named callers.

> IMPORTANT: This package is **R + Rcpp**, not Python. It began as the maintainer's
> Python optimization of Jim Holland's [nNIL](https://github.com/ncsumaize/nNIL)
> methodology. That Python optimization (a reorganization of Jim's code, **not**
> Jim's original) is **frozen at `legacy/python/`** (tag `python-final`) and is not
> the live code; Jim's original nNIL is at
> [`ncsumaize/nNIL`](https://github.com/ncsumaize/nNIL). Ignore any Python-era
> description — the authoritative overview is `README.md` and the design-of-record
> docs under `design/`.

## The one-verb API

```r
call_ancestry(data, caller = ..., design = ..., rrate = ..., err = ...)
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
- **Design priors** — `single_locus_expectation("BCnSm")` propagates the
  single-locus genotype-frequency vector through the backcross/selfing transition
  matrices (any design, no lookup); `breeding_prior()` wraps it as the public
  `c(REF, HET, ALT)` prior consumed by `call_gt()` and the engine's state
  frequencies. (Internal `.state_freqs()` maps a design to the engine's `f_1`/`f_2`.)
- **Engine** — `fit()` (EM / parameter fitting) then `decode()` (Viterbi /
  posteriors) → `to_segments()` (RLE to the common schema).

### The callers (all share the REF/HET/ALT chain + design priors)

The four **grid** callers are the pure (emission × duration) cells; the rest sit
off the grid.

| caller | emission | duration | typical input | lineage |
|---|---|---|---|---|
| `nnil` | `gt` (categorical) | geometric | called `GT` (hard genotypes) | Holland nNIL |
| `bbnil` | `count` (BetaBinomial) | geometric | allelic read counts (skim/BrB) | low-coverage count extension of nNIL |
| `catiger` | `gt` (categorical) | rigidity | called `GT` (hard genotypes) | categorical + rigidity |
| `rtiger` | `count` (BetaBinomial) | rigidity | allelic read counts | RTIGER (Julia-free port) |
| `binhmm` | anchored Gaussian on binned alt-freq | per-bin HMM | allelic read counts | "ancestry by bins" |
| `googa` | `gt` (GOOGA thresholds) | geometric | competitive-alignment RNA-seq counts | GOOGA (Flagel 2019 / Veltsos 2024), faithful |
| `atlas` | `gt` (GOOGA thresholds) | rigidity | competitive-alignment RNA-seq counts | this work's rigidity transcript caller |
| `lbimpute` | coverage-aware (bounded by `genotypeerr`) | distance-based (double-recomb penalty) | allelic read counts | LB-Impute (Fragoso 2014, native port) |
| `fsfhap` | genotype-error (5-state EM) | family-pooled Haldane | called `GT` for full-sib families | FSFHap (Swarts 2014, TASSEL) |
| `pedigree` | count or `gt` (input-detected) | BP over pedigree × genome | counts or hard-call `state`/`g` + a pedigree | family-coupled belief propagation |

`nnil` = categorical `gt` + geometric; `bbnil` swaps in the count/BetaBinomial
emission (same geometric duration); `catiger` swaps the rigidity duration onto the
`gt` emission; `rtiger` is count + rigidity. `googa` and `atlas` share GOOGA competitive-alignment
thresholding (`atlas_thresh`/`atlas_het`/`atlas_min_reads`): `googa` = gt +
**geometric** is the faithful reproduction (GOOGA's recombination-fraction F2 HMM
has **no** rigidity/minimum-run duration), and `atlas` = gt + **rigidity** is this
work's transcript caller. The no-HMM per-site **genotype** baseline (the paper's
"control": flat prior = maximum likelihood, HWE prior = MAP / het-excess) is
`call_gt()` — a *genotype* caller, deliberately **not** an ancestry caller
(`call_ancestry` does not dispatch it). This keeps the terminology wall between the
ancestry mosaic and genotypes (`design/TERMINOLOGY.md`); if a genotype baseline is
wanted in an ancestry comparison, the consumer converts it (`to_segments()`) itself.
`call_gt()` is also distinct from `interpolate_genotype()`, a deterministic densifier.

See `README.md` for per-caller detail and `?call_ancestry`. `lbimpute` is a native
port of LB-Impute: it keeps LB-Impute's coverage-aware emission and
distance-dependent transition (the `-dr` double-recomb penalty is `drp`) but
decodes with the engine's full-chromosome Viterbi rather than LB-Impute's
windowed forward/reverse consensus (the optimal path that window approximates).
The transition decays over a coordinate chosen by `unit`: `"bp"` (faithful,
uniform genome-wide rate) or `"cm"` (map-aware — pass a `cm` column of genetic
positions so local recombination rate is captured; the C++ `lb_viterbi_cpp` takes
a numeric `tpos` and is unit-agnostic). `recombdist` shares the coordinate's units
(unit-aware default 1e7 bp / 50 cM) and a validation layer warns on unit/`recombdist`
mismatches. Output coordinates are always bp. `write_vcf_impute()` (in `io.R`)
emits LB-Impute's imputed-VCF deliverable from the per-marker states — optional
and decoupled from the engine. `read_counts(format = "vcf_ad")` reads per-sample
allelic depths from a biallelic VCF's `AD` field (e.g. the LB-Impute example data).
`caller_sweep(caller = "lbimpute", values = <recombdist grid>)` calibrates
`recombdist`; it is **exact per value** (recombdist touches only the transition,
so the emission is computed once per run and only the Viterbi transition is
re-run over the grid, batched in C++ via `lb_viterbi_sweep_cpp`) — every swept
value equals a cold `call_ancestry(caller = "lbimpute", recombdist = v)`.

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
  `breeding_design.R` (`single_locus_expectation`, `breeding_prior`,
  `parse_design`) / `emission_regime.R`
  (`select_emission`); `io.R` (`read_counts`, `read_vcf_gt`); `calibrate.R`
  (`calibrate_r`); `sweep.R` (`caller_sweep`); `rtiger.R`, `binhmm.R`, `atlas.R`,
  `lbimpute.R`; `io.R` also holds `write_vcf_impute` (LB-Impute imputed-VCF output);
  `map.R` (`load_map`, `cm_to_mb` stub); `interpolate_genotype.R` (genotype densification);
  `pairwise_distance.R` / `select_independent.R` (LD marker thinning);
  `plot.R`, `nilHMM-package.R`, `RcppExports.R`.
- `src/` — Rcpp engine: `emission_count.cpp`, `forward_backward.cpp`, `viterbi.cpp`,
  `viterbi_batch.cpp`, `viterbi_batch_par.cpp` (RcppParallel/TBB), `segments.cpp`
  (RLE), `rtiger.cpp`, `lbimpute.cpp` (LB-Impute emission + distance-aware Viterbi),
  `interpolate_genotype.cpp`, `fast_indep.cpp` (FastIndep port),
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
- `load_map()` returns the bundled B73 v5 consensus map (`maize_map_v5`); the
  cM↔bp Marey-spline interpolators (`bp_to_cm`/`cm_to_bp`/`cm_to_mb`) default to it
  and accept any `chr/bp/cm` map.

## Related

Companion analysis repo `zealhmm` installs nilHMM as a dependency and hosts the
caller comparison and methods paper. Cite the nNIL population paper (Zhong et al.
2025, bioRxiv 10.1101/2025.01.29.635337) per `CITATION.cff`.
