# nilHMM

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE.md)
[![R-CMD-check](https://github.com/sawers-rellan-labs/nilhmm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sawers-rellan-labs/nilhmm/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**nilHMM** calls ancestry (donor introgressions) in Near-Isogenic Lines from
sequencing data. It is a single **duration-aware 3-state (REF / HET / ALT) HMM
engine** with two swappable axes:

- **emission** — `count` (BetaBinomial over ref/alt read depths) or `gt`
  (categorical genotype model with an explicit genotyping-error model), and
- **duration** — `geometric`, `rigidity` (a minimum-run-length prior), or `hsmm`.

Those choices, plus a per-unit calling rule, express a family of named **callers**:
`nnil`, `rtiger`, `binhmm`, `atlas`, and `lbimpute`. The package is **data-agnostic** — every
function takes `(data, params)` and returns calls; pipeline scripts own file paths
and sample lists. The whole API is one verb:

```r
call_ancestry(data, caller = ..., design = ..., r = ..., err = ...)
```

nilHMM began as a Python package on Jim Holland's
[nNIL](https://github.com/ncsumaize/nNIL) methodology; it is now an R + Rcpp
package. The Python original is frozen at [`legacy/python/`](legacy/python/)
(tag `python-final`).

## Installation

```r
# install.packages("devtools")
devtools::install_github("sawers-rellan-labs/nilhmm")
```

Building from source needs a **C++ toolchain** (the engine is Rcpp /
RcppParallel). The `rtiger` caller is a **Julia-free port** of the RTIGER HMM, so
no Julia install is required. `rebmix` is an optional `Suggests` (only for
bit-exact reproduction of the original `binhmm` clustering backend).

## Quick start

```r
library(nilHMM)

# Read allelic counts (headerless "chr pos ref n_ref alt n_alt" TSV, or a
# directory of them) into the engine's observation table, then call ancestry.
calls <- call_ancestry(
  read_counts("path/to/sample_counts.tsv"),
  caller = "nnil",     # Holland's nNIL count caller
  design = "BC2S2",    # breeding-design priors (f_1, f_2)
  r   = 1e-4,          # recombination / geometric self-transition
  err = 0.01           # baseline read error for the count emission
)

# calls is the common segment schema:
#   source, donor, name, chr, start_bp, end_bp, state   (state in {REF, HET, ALT})
head(calls)
```

For saturated-depth data with trustworthy called genotypes (e.g. target
sequencing), read the VCF `GT` field instead and let it select the categorical
emission:

```r
calls <- call_ancestry(read_vcf_gt("target.vcf.gz"),
                       caller = "nnil", design = "BC2S2")   # g-only input -> gt emission
```

## The callers

All four share the 3-state REF/HET/ALT chain and the design priors; they differ
in emission, duration, and the input they expect.

| caller | emission | duration | typical input | key parameters | lineage |
|---|---|---|---|---|---|
| **`nnil`** | `count` (BetaBinomial) or `gt` (categorical) | geometric | allelic read counts (low-cov skim / BrB), or called `GT` | `r`, `err`, `conc`, `fit_means` | Holland [nNIL](https://github.com/ncsumaize/nNIL) |
| **`rtiger`** | `count` (BetaBinomial) | rigidity (min run length) | allelic read counts | `r` (integer rigidity), `seed`, `threads` | [RTIGER](https://github.com/rfael0cm/RTIGER) (Julia-free port of [`faustovrz/RTIGER`](https://github.com/faustovrz/RTIGER)) |
| **`binhmm`** | anchored Gaussian on binned alt-freq | per-bin HMM smooth | allelic read counts | `bin_size`, `cluster_method` | "Ancestry Analysis by bins" |
| **`atlas`** | `gt` (categorical, GOOGA thresholds) | geometric | competitive-alignment recurrent/donor read counts (RNA-seq) | `atlas_thresh`, `atlas_het`, `atlas_min_reads` | GOOGA competitive alignment |
| **`lbimpute`** | coverage-aware (LB-Impute) | distance-based (double-recomb penalty) | low-coverage allelic read counts (GBS / skim, <1×) | `err`, `genotypeerr`, `recombdist`, `drp` | [LB-Impute](https://github.com/dellaporta-laboratory/LB-Impute) (Fragoso et al. 2014) |

- `nnil` — count caller pools single-read observations along a segment via a
  BetaBinomial emission (`err`, `conc`); `fit_means = TRUE` EM-fits the per-state
  alt fractions. Also runs the categorical `gt` emission (Holland's error model:
  `germ`, `gert`, `p`, `mr`, `nir`) on called genotypes.
- `rtiger` — a faithful reimplementation of the RTIGER rigidity HMM (EM + Viterbi,
  border re-placement), ported off Julia; `r` is the integer rigidity (minimum
  run length).
- `binhmm` — bins the genome (default 1 Mb), calls per-bin state with an anchored
  3-state Gaussian-emission HMM (the default `cluster_method = "gauss"` fixes HET
  over-call and high-coverage fragmentation), with the original GMM/k-means/rebmix
  clustering backends available for reproduction.
- `atlas` — a GOOGA-style caller for transcript/competitive-alignment data:
  recurrent vs donor read fractions thresholded into genotype calls, then smoothed.
- `lbimpute` — a native port of LB-Impute (Fragoso et al. 2014) for very
  low-coverage (<1×) biallelic populations: a coverage-aware emission (bounded by
  `genotypeerr` so one artifactual marker can't dominate) and a distance-dependent
  transition (recombination scales with the marker gap over `recombdist`; the
  homozygous↔homozygous switch carries a double-recombination penalty unless
  `drp = TRUE`, for RILs). The transition decays over physical bp (`unit = "bp"`,
  the faithful uniform-rate model) or, with a `cm` column of map positions,
  genetic distance (`unit = "cm"`) so the *local* recombination rate — e.g. maize
  centromeric suppression — is captured; output coordinates stay bp either way,
  and `recombdist`'s default is unit-aware (`1e7` bp / `50` cM). Decoded with the
  engine's full-chromosome Viterbi — the optimal path that LB-Impute's windowed
  forward/reverse consensus approximates. Emit an imputed VCF from the result
  with `write_vcf_impute()`.

## Related building blocks

`read_counts()`, `read_vcf_gt()` (I/O); `caller_spec()`, `emission_count()`,
`emission_gt()`, `duration_geometric()`, `duration_rigidity()`,
`duration_hsmm()`, `design_priors()`, `fit()`, `decode()` (engine internals);
`calibrate_r()`, `select_emission()` (calibration); `plot_fragment_sizes()`,
`write_common_schema()`, `write_vcf_impute()` (output). See `?call_ancestry` and the package
documentation for the full reference.

## Documentation & benchmarks

- Package site (pkgdown): <https://sawers-rellan-labs.github.io/nilhmm/>
- Caller comparison & the methods paper live in the companion
  [`zealhmm`](https://github.com/sawers-rellan-labs/zealhmm) analysis repo, which
  installs nilHMM as a dependency.

## Citation

If you use nilHMM, please cite the nNIL population paper it is built on, and (once
released) the nilHMM methods paper — see `CITATION`.

```
Zhong, T., Mullens, A., Morales, L., Swarts, K.L., Stafstrom, W.C., He, Y.,
Sermons, S.M., Yang, Q., Lopez-Zuniga, L.O., Rucker, E., Thomason, W.E.,
Nelson, R.J., Jamann, T., Balint-Kurti, P.J., & Holland, J.B. (2025).
A maize near-isogenic line population designed for gene discovery and
characterization of allelic effects. bioRxiv.
https://doi.org/10.1101/2025.01.29.635337
```

## License

MIT — see [LICENSE.md](LICENSE.md).
