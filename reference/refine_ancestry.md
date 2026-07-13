# Refine per-individual ancestry calls over a pedigree

Couples relatives through the pedigree to correct per-individual
[`call_states()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_states.md)
calls: structured loopy belief propagation over the pedigree x genome
grid (see
[`pedigree_bp_cpp()`](https://sawers-rellan-labs.github.io/nilhmm/reference/pedigree_bp_cpp.md),
design/PEDIGREE_HMM.md). Depth-blind refinement – emission is
[`emission_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/emission_gt.md)
over the hard `state` calls. Families are processed independently;
latent ungenotyped ancestors (taxa named as parents but absent from
`mosaic`) impose chromosome continuity across siblings.

## Usage

``` r
refine_ancestry(
  mosaic,
  pedigree,
  design = "BC2S3",
  err = 0.05,
  gert = 0.1,
  rrate = 0.01,
  maxiter = 30L,
  tol = 1e-04,
  lambda = 0.5,
  ped_format = c("fam", "fsfhap")
)
```

## Arguments

- mosaic:

  Per-marker state table from
  [`call_states()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_states.md):
  needs `name, chr, pos, state` (+ optional `source, donor, cm`).
  `state` in `{0 REF, 1 HET, 2 ALT, 3 missing}`.

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

- err, gert:

  [`emission_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/emission_gt.md)
  genotyping-error rates (call-level, not raw sequencing error – see the
  depth caveat in design/PEDIGREE_HMM.md).

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
