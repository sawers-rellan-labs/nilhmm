# Changelog

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
