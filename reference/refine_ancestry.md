# Refine per-individual ancestry calls over a pedigree

Couples relatives through the pedigree to correct per-individual
[`call_states()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_states.md)
calls: structured loopy belief propagation over the pedigree x genome
grid (see
[`pedigree_bp_cpp()`](https://sawers-rellan-labs.github.io/nilhmm/reference/pedigree_bp_cpp.md),
design/PEDIGREE_HMM.md). Two emission modes: `"gt"` (depth-blind,
[`emission_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/emission_gt.md)
over hard `state` calls) or `"count"` (depth-aware "counts-first"
pedigree calling –
[`emission_count()`](https://sawers-rellan-labs.github.io/nilhmm/reference/emission_count.md)
BetaBinomial over read depths, so a zero-coverage marker emits flat and
the pedigree fills it in weighted by how confident the relatives are).
Families are processed independently; latent ungenotyped ancestors (taxa
named as parents but absent from `mosaic`) impose chromosome continuity
across siblings.

## Usage

``` r
refine_ancestry(
  mosaic,
  pedigree,
  design = "BC2S3",
  emission = c("gt", "count"),
  err = NULL,
  gert = 0.1,
  conc = 20,
  rrate = 0.01,
  maxiter = 30L,
  tol = 1e-04,
  lambda = 0.5,
  ped_format = c("fam", "fsfhap")
)
```

## Arguments

- mosaic:

  Per-marker table keyed by `name, chr, pos` (+ optional
  `source, donor, cm`). For `emission = "gt"` needs `state` in
  `{0,1,2,3=missing}` (from
  [`call_states()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_states.md));
  for `emission = "count"` needs `n_ref, n_alt` read depths.

- pedigree:

  A pedigree: a path (read via
  [`read_pedigree()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_pedigree.md))
  or a
  [`read_pedigree()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_pedigree.md)-shaped
  data.frame (`taxon, family, parent1, parent2`). `taxon` joins to
  `mosaic$name`; a `taxon` used as a parent but absent from `mosaic` is
  a latent ancestor.

- design:

  Breeding design `"BC{n}S{m}"` -\> founder prior `pi_0` and per-node
  `meioses` (via
  [`design_priors()`](https://sawers-rellan-labs.github.io/nilhmm/reference/design_priors.md)).

- emission:

  `"gt"` (depth-blind, over hard states) or `"count"` (depth-aware
  BetaBinomial over `n_ref`/`n_alt`).

- err:

  Genotyping/read error:
  [`emission_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/emission_gt.md)
  `germ` when `emission="gt"` (default 0.05), per-read error when
  `emission="count"` (default 0.01).

- gert:

  [`emission_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/emission_gt.md)
  genotyping-error-rate-in-transmission (`"gt"` only).

- conc:

  [`emission_count()`](https://sawers-rellan-labs.github.io/nilhmm/reference/emission_count.md)
  BetaBinomial concentration (`"count"` only).

- rrate:

  Per-bp recombination fraction applied to `pos` gaps when `mosaic` has
  no `cm` column (else Haldane on `cm`).

- maxiter, tol, lambda:

  BP sweeps, convergence tolerance, damping.

- ped_format:

  Passed to
  [`read_pedigree()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_pedigree.md)
  when `pedigree` is a path.

## Value

`mosaic` with `state` replaced by the refined per-marker calls (same
columns and row order; genotyped leaves only). Feed to
[`to_segments()`](https://sawers-rellan-labs.github.io/nilhmm/reference/to_segments.md).

## See also

[`call_states()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_states.md),
[`simulate_family()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_family.md),
[`pedigree_bp_cpp()`](https://sawers-rellan-labs.github.io/nilhmm/reference/pedigree_bp_cpp.md),
[`to_segments()`](https://sawers-rellan-labs.github.io/nilhmm/reference/to_segments.md).

## Examples

``` r
if (FALSE) { # \dontrun{
sim   <- simulate_family("BC2S3", families = 10, sibs = 10, seed = 1)
obs   <- simulate_counts(sim$truth, depth = 0.5, seed = 1)
calls <- call_states(obs, caller = "rtiger", design = "BC2S3")
ref   <- refine_ancestry(calls, sim$pedigree, design = "BC2S3")
} # }
```
