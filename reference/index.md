# Package index

## Calling ancestry

The top-level API, the named-caller registry, and parameter sweeps.

- [`call_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md)
  : Top-level ancestry-calling API (decode + segment)
- [`call_states()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_states.md)
  : Coordinate-free ancestry decode (per-observation states)
- [`caller_spec()`](https://sawers-rellan-labs.github.io/nilhmm/reference/caller_spec.md)
  : Resolve a named caller into emission + duration specs
- [`caller_sweep()`](https://sawers-rellan-labs.github.io/nilhmm/reference/caller_sweep.md)
  : Sweep a caller's segmentation parameter with an amortized fit

## Engine

Fit parameters, decode the state path, and collapse to segments.

- [`fit()`](https://sawers-rellan-labs.github.io/nilhmm/reference/fit.md)
  : Fit HMM emission/transition parameters
- [`decode()`](https://sawers-rellan-labs.github.io/nilhmm/reference/decode.md)
  : Decode the most-likely state path (Viterbi)
- [`to_segments()`](https://sawers-rellan-labs.github.io/nilhmm/reference/to_segments.md)
  : Collapse per-unit ancestry states into genomic segments

## Emissions

The swappable emission axis (count / gt) and the depth-regime selector.

- [`emission_count()`](https://sawers-rellan-labs.github.io/nilhmm/reference/emission_count.md)
  : Count (BetaBinomial) emission
- [`emission_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/emission_gt.md)
  : Genotype (categorical) emission
- [`select_emission()`](https://sawers-rellan-labs.github.io/nilhmm/reference/select_emission.md)
  : Select an emission model from sequencing depth

## Duration

The swappable duration axis (geometric / rigidity / hsmm).

- [`duration_geometric()`](https://sawers-rellan-labs.github.io/nilhmm/reference/duration_geometric.md)
  : Geometric (memoryless) duration
- [`duration_rigidity()`](https://sawers-rellan-labs.github.io/nilhmm/reference/duration_rigidity.md)
  : Rigidity duration (RTIGER reimplementation)
- [`duration_hsmm()`](https://sawers-rellan-labs.github.io/nilhmm/reference/duration_hsmm.md)
  : Explicit-duration (HSMM) sojourn – reserved

## Input / output

Readers that normalise formats to the observation table, and the
writers.

- [`read_counts()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_counts.md)
  : Read allelic counts into the engine's observation table
- [`read_vcf_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_vcf_gt.md)
  : Read hard genotype calls (VCF GT) into the engine's observation
  table
- [`read_hapmap()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_hapmap.md)
  : Read a TASSEL HapMap into the engine's observation table
- [`read_plink()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_plink.md)
  : Read a PLINK binary genotype (.bed/.bim/.fam) into the engine's
  table
- [`read_pedigree()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_pedigree.md)
  : Read a pedigree (FSFHap TSV or PLINK .fam) into a family/priors
  table
- [`write_common_schema()`](https://sawers-rellan-labs.github.io/nilhmm/reference/write_common_schema.md)
  : Write segment calls in the common schema
- [`write_vcf_impute()`](https://sawers-rellan-labs.github.io/nilhmm/reference/write_vcf_impute.md)
  : Write imputed per-marker genotypes as a VCF (LB-Impute-style
  deliverable)

## Genotype calling & downstream utilities

A linkage-free genotype caller and helpers for the QTL-mapping pipeline
(genotype densification, LD-based marker thinning).

- [`call_gl()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_gl.md)
  : Per-site genotype-likelihood (GL) genotype caller with a swappable
  prior
- [`interpolate_genotype()`](https://sawers-rellan-labs.github.io/nilhmm/reference/interpolate_genotype.md)
  : Interpolate genotypes onto a target marker grid (Tian 2011 / Chen
  2019)
- [`pairwise_distance()`](https://sawers-rellan-labs.github.io/nilhmm/reference/pairwise_distance.md)
  : Pairwise marker relatedness (r2, mutual information, or variation of
  information)
- [`select_independent()`](https://sawers-rellan-labs.github.io/nilhmm/reference/select_independent.md)
  : Select a maximal independent set of markers (LD thinning)
- [`position_distance()`](https://sawers-rellan-labs.github.io/nilhmm/reference/position_distance.md)
  : Pairwise marker distance from map/physical coordinates

## Design priors & calibration

Breeding-design priors, the fragment-size Null, and r calibration.

- [`design_priors()`](https://sawers-rellan-labs.github.io/nilhmm/reference/design_priors.md)
  : Single-locus genotype-frequency priors for a breeding design

- [`calibrate_r()`](https://sawers-rellan-labs.github.io/nilhmm/reference/calibrate_r.md)
  :

  Calibrate the duration hyperparameter `r` by KS-vs-sim

- [`expected_fragment_dist()`](https://sawers-rellan-labs.github.io/nilhmm/reference/expected_fragment_dist.md)
  : Expected fragment-size distribution (the Null / KS target)

- [`fit_design_gamma()`](https://sawers-rellan-labs.github.io/nilhmm/reference/fit_design_gamma.md)
  : Fit a design Gamma from simulated segments

- [`load_map()`](https://sawers-rellan-labs.github.io/nilhmm/reference/load_map.md)
  : Load the bundled consensus map

- [`cm_to_mb()`](https://sawers-rellan-labs.github.io/nilhmm/reference/cm_to_mb.md)
  : Convert segment sizes from cM to Mb using a map

## Plotting

- [`plot_fragment_sizes()`](https://sawers-rellan-labs.github.io/nilhmm/reference/plot_fragment_sizes.md)
  : Plot called fragment sizes against the expected Null
