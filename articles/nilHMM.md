# Getting started with nilHMM

**nilHMM** calls ancestry — donor introgressions — in Near-Isogenic
Lines and related backcross/full-sib populations from sequencing data.
Under the hood it is a single **duration-aware 3-state (REF / HET / ALT)
HMM engine** with two swappable axes (emission × duration); their
combinations, plus breeding-design priors, express a family of named
callers (the grid cells `nnil`, `bbnil`, `catiger`, `rtiger`, plus
`googa`, `atlas`, `binhmm`, `lbimpute`, `fsfhap`, `pedigree`). This
vignette covers the common workflow; see
[`vignette("callers")`](https://sawers-rellan-labs.github.io/nilhmm/articles/callers.md)
to choose a caller,
[`vignette("engine")`](https://sawers-rellan-labs.github.io/nilhmm/articles/engine.md)
for the internals, and
[`vignette("fsfhap")`](https://sawers-rellan-labs.github.io/nilhmm/articles/fsfhap.md)
for full-sib families.

``` r

library(nilHMM)
```

## Genotype calls and ancestry inference

nilHMM reconstructs ancestry with a hidden Markov model, and the way the
package is organized follows from keeping the model’s output separate
from the data it was computed on. Two different questions get asked of
the same sequencing reads.

The first is local. At a single marker, what genotype does this sample
carry? That is decided from the read counts at that one position — the
reference and alternate depths, `n_ref` and `n_alt` — and nothing else.
[`call_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_gt.md)
does this, one marker at a time, without looking at any neighbouring
marker. Its answer is a genotype (`0`, `1`, or `2`) that summarizes the
reads at that site: with a flat prior it is the maximum-likelihood call,
and with a Hardy–Weinberg prior it is the MAP call, which at low depth
leans toward heterozygous — the naive baseline the HMM is meant to
improve on.

The second question is the one the HMM answers. Given all the markers
along a chromosome, what is the most likely ancestry state — REF, HET,
or ALT — at each position? Neighbouring markers matter here, because
ancestry only changes at recombination breakpoints, so a stretch of
consecutive markers shares one state. The HMM uses that linkage to pool
the per-position observations — read counts for the count callers, or
already-called genotypes for the categorical ones — and then projects an
ancestry state back onto every input position. That projection is a
model result, not a direct reading of the observation at any single
site, and it can disagree with a per-site genotype call exactly where
linkage outweighs a noisy one.

Because a per-site summary of the reads and a linkage-based inference
across the chromosome are different things, nilHMM keeps them in
separate functions. Genotype calling is
[`call_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_gt.md),
which takes read counts. Ancestry inference is the family of callers,
run through
[`call_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md),
which return the segment schema — the mosaic. The callers that model
called genotypes (`nnil`, `catiger`) therefore expect genotypes you have
already called and will not threshold read counts on their own, so that
the hard-calling step, when you take it, stays a deliberate choice
rather than something the ancestry caller did silently.

## The one-verb API

Everything funnels through a single call:

``` r

call_ancestry(data, caller = ..., design = ..., ...)
```

`data` is a tidy long table the readers produce; `caller` picks the
model; `design` supplies the breeding-design priors. The package is
**data-agnostic** — functions take `(data, params)` and never touch file
paths, so you own where data comes from and where results go.

## A worked example

We need an observation table. In practice you read one with
[`read_counts()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_counts.md)
(a `chr pos ref n_ref alt n_alt` TSV or a directory of them),
[`read_vcf_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_vcf_gt.md),
[`read_hapmap()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_hapmap.md),
or
[`read_plink()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_plink.md).
Here we generate one with the package’s own design-driven simulator:
[`simulate_nil()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_nil.md)
builds a breeding design and runs it through
[**simcross**](https://cran.r-project.org/package=simcross) (real
meiosis and recombination on the **bundled B73 v5 consensus map**,
\[load_map()\]) to get the true donor mosaic, and
[`simulate_counts()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_counts.md)
degrades that truth to observed allelic read counts. We take a **BC2S2**
cohort — six lines on chromosomes 1–2, donor **B** crossed onto
recurrent parent **A** (so `REF` = A, `ALT` = donor B):

``` r

truth  <- simulate_nil("BC2S2", n = 6, chr = 1:2, n_markers = 300, donor = "B",
                       names = sprintf("NIL%02d", 1:6), seed = 1)   # true donor mosaic
counts <- simulate_counts(truth, depth = 6, seed = 1)[        # observed low-cov counts
  c("name", "donor", "chr", "pos", "n_ref", "n_alt")]
head(counts)
#>    name donor chr     pos n_ref n_alt
#> 1 NIL01     B   1   37410     0     4
#> 2 NIL01     B   1 1883430     0     5
#> 3 NIL01     B   1 3729450     0     6
#> 4 NIL01     B   1 5575469     0     9
#> 5 NIL01     B   1 7421489     0     4
#> 6 NIL01     B   1 9267509     0     9
```

Call ancestry with the `bbnil` count caller and BC2S2 design priors:

``` r

calls <- call_ancestry(counts, caller = "bbnil", design = "BC2S2",
                       rrate = 1e-4, err = 0.01)
head(calls)
#>   source donor  name chr  start_bp    end_bp state
#> 1 nilHMM     B NIL01   1     37410  27727705     2
#> 2 nilHMM     B NIL01   1  29573725 308322690     0
#> 3 nilHMM     B NIL01   2     98554 223047189     0
#> 4 nilHMM     B NIL01   2 224905095 232336716     2
#> 5 nilHMM     B NIL01   2 234194621 243484148     0
#> 6 nilHMM     B NIL02   1     37410   3729450     0
```

## The common segment schema

Every caller returns the same tidy **segment** table, so results are
directly comparable across callers and populations:

| column | meaning |
|----|----|
| `source` | free label for the run/method |
| `donor` | donor/taxon label (per line if the input carries one) |
| `name` | sample identifier |
| `chr` | chromosome |
| `start_bp`, `end_bp` | segment bounds (bp) |
| `state` | `REF` (0, recurrent hom), `HET` (1), or `ALT` (2, donor hom) |

``` r

table(calls$state)
#> 
#>  0  1  2 
#> 20  4  9
```

Numeric states are `0 = REF`, `1 = HET`, `2 = ALT`. The package’s
[`paint_calls()`](https://sawers-rellan-labs.github.io/nilhmm/reference/paint_calls.md)
renders those segments as a chromosome painting — one band per line,
faceted by chromosome (REF gold / HET green / ALT purple):

``` r

paint_calls(calls)
```

![](nilHMM_files/figure-html/plot-1.png)

## The callers

Four callers form the core family — the grid in panel B of the paper —
laid out over two choices. The first is the **emission**: how a marker’s
observation becomes evidence for each state, either a BetaBinomial on
the read counts or a categorical model on an already-called genotype.
The second is the **duration**: how readily the ancestry state is
allowed to switch between adjacent markers, either a plain geometric
per-marker switch probability or a rigidity model that requires a
minimum run length before a switch. The four cells are those two choices
crossed:

- `bbnil` — BetaBinomial on read counts, geometric duration. The
  read-count caller for low-coverage skim/BrB data.
- `nnil` — categorical on called genotypes, geometric duration.
  Holland’s original nNIL, for saturated genotype data.
- `rtiger` — BetaBinomial on read counts, rigidity duration. The RTIGER
  segmentation (a Julia-free port), robust at low depth.
- `catiger` — categorical on called genotypes, rigidity duration.

Read across the grid, the emission follows the data you have (read
counts vs called genotypes) and the duration sets how hard short
segments are smoothed away.

`googa` and `atlas` are the same categorical HMM applied to RNA-seq.
When reads are competitively aligned to the two parents, expression and
allele-specific expression distort the allele fraction, so a
BetaBinomial on the counts would be misled. Both callers instead
threshold the competitive counts to hard genotype calls (GOOGA’s rule)
and then run the categorical model — differing only in duration: `googa`
uses the geometric duration (a faithful reproduction of GOOGA) and
`atlas` uses the rigidity duration.

`pedigree` is the one that leaves the per-sample pattern. The callers
above decode each line on its own; `pedigree` couples relatives, passing
ancestry information across a family through the pedigree so a
thinly-covered line borrows from its better-covered relatives.

The rest are specialized or external. `binhmm` sits off the grid: it
bins the genome and calls a per-bin state from a Gaussian emission, an
option for dense data with uneven coverage. `lbimpute` and `fsfhap` are
faithful ports of established biparental imputers — LB-Impute (Fragoso
et al., a distance-dependent transition) and TASSEL’s FSFHap (per-family
haplotype imputation) — carried as comparison baselines; each keeps its
own internal model but returns the same segment schema.

[`vignette("callers")`](https://sawers-rellan-labs.github.io/nilhmm/articles/callers.md)
runs one example of each
([`vignette("fsfhap")`](https://sawers-rellan-labs.github.io/nilhmm/articles/fsfhap.md)
for the family workflow). The per-site genotype baseline of the previous
section is
[`call_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_gt.md)
—
[`?call_gt`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_gt.md).

## Next steps

- **Choosing and comparing callers:**
  [`vignette("callers")`](https://sawers-rellan-labs.github.io/nilhmm/articles/callers.md)
- **The engine (emission × duration, custom callers):**
  [`vignette("engine")`](https://sawers-rellan-labs.github.io/nilhmm/articles/engine.md)
- **Full-sib families (HapMap/pedigree → fsfhap):**
  [`vignette("fsfhap")`](https://sawers-rellan-labs.github.io/nilhmm/articles/fsfhap.md)
- Function reference:
  [`?call_ancestry`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md),
  [`?read_counts`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_counts.md),
  [`?to_segments`](https://sawers-rellan-labs.github.io/nilhmm/reference/to_segments.md).
  \`\`\`
