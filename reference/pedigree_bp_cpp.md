# Loopy belief propagation over one family's pedigree x genome grid

Data-agnostic BP kernel for
[`refine_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/refine_ancestry.md);
processes ONE (family, chromosome). Node fields are pre-built in R. See
design/PEDIGREE_HMM.md.

## Usage

``` r
pedigree_bp_cpp(
  M,
  parent,
  meioses,
  hasData,
  emit,
  rho,
  pimat,
  r,
  root,
  maxIters,
  tol,
  lambda
)
```

## Arguments

- M:

  Number of markers (chromosome length).

- parent:

  0-based parent index per node (`-1` for the root).

- meioses:

  Per-node accumulated meiosis count (transition block length).

- hasData:

  Per-node logical; `TRUE` for genotyped leaves.

- emit:

  Length-V list; for a `hasData` node an `M x 3` emission matrix (P(obs
  \| state)), otherwise ignored (latent nodes emit 1).

- rho:

  V x 3 marker-0 prior per node.

- pimat:

  V x 3 generation stationary per node (transition relaxes to it).

- r:

  Length `M-1` per-interval recombination fraction (base; meiosis
  compounding is applied per node).

- root:

  0-based index of the family root (founder).

- maxIters, tol, lambda:

  Message-passing iterations, convergence tolerance, damping.

## Value

Length-V list of `M x 3` posterior belief matrices; attribute `iters`
records sweeps run.
