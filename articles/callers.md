# The callers: one example each

nilHMM draws a firm line between two different operations on the same
sequencing data. Keeping them apart is a deliberate design choice, and
it is the key to reading the rest of this vignette:

- **Genotype calling** — deciding the genotype at a marker from *that
  marker’s own read observations*, one site at a time, with **no
  linkage**. This is
  [`call_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_gt.md).
  Its output is a per-site genotype (`0/1/2`) — an *observation-level*
  call, not an inference across the genome.
- **Ancestry-mosaic inference** — reconstructing the REF/HET/ALT
  ancestry state *along the chromosome*, borrowing strength across
  neighbouring markers through the HMM’s linkage (duration) model. This
  is the family of named **callers**, run through
  [`call_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md).
  Its output is the common segment schema — the ancestry **mosaic**.

The wall between them is enforced, not incidental: an ancestry caller
never silently calls genotypes for you (`nnil`/`catiger` require
*called* genotypes; they will not threshold read counts), and
[`call_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_gt.md)
is not dispatchable through
[`call_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md).
So this vignette **starts with the genotype call**, then walks the
ancestry callers — all on the same synthetic data, all returning the
same segment schema so you can compare them directly.

``` r

library(nilHMM)
```

## Synthetic data

[`simulate_nil()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_nil.md)
builds a BC2S2 cohort (donor **B** on recurrent parent **A**) on the
bundled maize map — the **truth** — and
[`simulate_counts()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_counts.md)
degrades it to observed data carrying **both** allelic read counts
(`n_ref`/`n_alt`, for the count/BetaBinomial callers and
[`call_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_gt.md))
and a called genotype `g` in `{0,1,2,3}` (for the categorical gt callers
`nnil`/`catiger`):

``` r

truth <- simulate_nil("BC2S2", n = 8, chr = 1:2, n_markers = 300, donor = "B",
                      names = sprintf("NIL%02d", 1:8), seed = 1)
obs   <- simulate_counts(truth, depth = 6, seed = 1)
counts <- obs[c("name", "donor", "chr", "pos", "n_ref", "n_alt")]
gt     <- obs[c("name", "donor", "chr", "pos", "g")]
```

## Step 1 — the genotype call: `call_gt()`

Before any ancestry inference, you can ask the most basic question of
the data: *what is the genotype at each marker?*
[`call_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_gt.md)
answers it **per site**, deciding each `(marker, sample)` from that
marker’s own `(n_ref, n_alt)` — no linkage, no neighbours, no HMM. It is
a **genotype** caller, not an ancestry caller, so it is called directly
(never through
[`call_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md)),
and its output is a genotype `0/1/2` (with `NA` at zero depth), **not**
the ancestry mosaic.

The single-locus prior is the whole story of the two useful settings:

- `prior = "flat"` → the maximum-likelihood genotype (argmax
  genotype-likelihood). Het-blind at depth 1: one ALT read alone reads
  as ALT-hom.
- `prior = "hwe"` → the Hardy–Weinberg posterior MAP. At low depth this
  is deliberately **het-excess** (one ALT read reads as HET) — the
  classic no-HMM “control” that the ancestry callers must beat.

``` r

ml_calls <- call_gt(counts$n_ref, counts$n_alt, prior = "flat")  # maximum likelihood
hw_calls <- call_gt(counts$n_ref, counts$n_alt, prior = "hwe")   # HWE MAP (het-excess)
table(ml_calls, useNA = "ifany"); table(hw_calls, useNA = "ifany")
#> ml_calls
#>    0    1    2 <NA> 
#> 1864  214  313    9
#> hw_calls
#>    0    1    2 <NA> 
#> 1883  189  319    9
```

That is the entire genotype layer. Everything below **infers ancestry**:
it takes observations (read counts, or genotypes you have *already*
called) and reconstructs the REF/HET/ALT mosaic *along the chromosome*
using linkage. If you want a genotype baseline inside an ancestry
comparison, you convert these calls yourself
([`to_segments()`](https://sawers-rellan-labs.github.io/nilhmm/reference/to_segments.md)
on the genotype-as-state) — the package offers no genotype→mosaic
shortcut, and
[`call_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_gt.md)
is also distinct from the deterministic
[`interpolate_genotype()`](https://sawers-rellan-labs.github.io/nilhmm/reference/interpolate_genotype.md)
densifier.

## The ancestry callers

Everything from here runs through
[`call_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md)
and returns the ancestry mosaic (the segment schema). The callers share
the 3-state REF/HET/ALT chain and the breeding-design priors, and differ
in the **emission** model, the **duration** prior, and the **input**
they expect — so you can run several on the same data and compare
directly.

### The four grid callers (emission × duration)

The pure cells of the engine’s grid share the REF/HET/ALT chain and
differ only in two axes: the emission (categorical `gt` vs
count/BetaBinomial) and the duration (geometric vs rigidity).

- **`nnil`** — Holland’s nNIL: categorical `gt` emission + geometric
  duration. Takes a **called genotype** `g` column; it does *not*
  threshold read counts into genotypes (hard-calling is your explicit
  step,
  e.g. [`call_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_gt.md)).
  Error model: `germ`, `gert`, `p`, `mr`, `nir`.
- **`bbnil`** — the low-coverage count extension: count/BetaBinomial
  emission (`err`, `conc`; `fit_means = TRUE` EM-fits the per-state alt
  fractions) + the same geometric duration.
- **`catiger`** — categorical `gt` emission + rigidity duration.
- **`rtiger`** — count/BetaBinomial emission + rigidity duration
  (below).

``` r

nnil    <- call_ancestry(gt,     caller = "nnil",    design = "BC2S2")           # gt + geometric
bbnil   <- call_ancestry(counts, caller = "bbnil",   design = "BC2S2",
                         rrate = 1e-4, err = 0.01)                               # count + geometric
catiger <- call_ancestry(gt,     caller = "catiger", design = "BC2S2", rigidity = 5L)  # gt + rigidity
nrow(nnil); nrow(bbnil); nrow(catiger)
#> [1] 51
#> [1] 50
#> [1] 47
```

### `rtiger` — rigidity segmentation

A Julia-free port of the RTIGER rigidity HMM (EM + Viterbi + border
re-placement), the count-emission rigidity cell. `rigidity` is the
integer minimum run length.

``` r

rt <- call_ancestry(counts, caller = "rtiger", design = "BC2S2",
                    rigidity = 5L, seed = 1L)
head(rt, 3)
#>   source donor  name chr start_bp    end_bp state
#> 1 nilHMM     B NIL01   1    37410  27727705     2
#> 2 nilHMM     B NIL01   1 29573725 308322690     0
#> 3 nilHMM     B NIL01   2    98554 223047189     0
```

### `binhmm` — per-bin calling

Bins the genome (default 1 Mb) and calls per-bin state with an anchored
3-state Gaussian-emission HMM. Good for noisy or uneven coverage.

``` r

bh <- call_ancestry(counts, caller = "binhmm", design = "BC2S2", bin_size = 5e6)
head(bh, 3)
#>   source donor  name chr start_bp    end_bp state
#> 1 nilHMM     B NIL01   1    37410  29573725     2
#> 2 nilHMM     B NIL01   1 31419744 308322690     0
#> 3 nilHMM     B NIL01   2    98554 243484148     0
```

### `googa` / `atlas` — competitive-alignment (transcript ancestry)

For transcript / competitive-alignment data: `n_ref`/`n_alt` are the
recurrent and donor read counts, thresholded into hard genotype calls
(`atlas_thresh`, `atlas_het`, `atlas_min_reads`) then smoothed.
**`googa`** is the faithful reproduction — gt + **geometric**, matching
GOOGA’s recombination-fraction F2 HMM (Flagel 2019; Veltsos 2024), which
carries no rigidity. **`atlas`** is this work’s transcript caller: the
same thresholding with the **rigidity** duration.

``` r

gg <- call_ancestry(counts, caller = "googa", design = "BC2S2",
                    atlas_thresh = 0.95, atlas_het = 0.25, atlas_min_reads = 5L)
at <- call_ancestry(counts, caller = "atlas", design = "BC2S2", rigidity = 5L,
                    atlas_thresh = 0.95, atlas_het = 0.25, atlas_min_reads = 5L)
head(gg, 3); head(at, 3)
#>   source donor  name chr start_bp    end_bp state
#> 1 nilHMM     B NIL01   1    37410  27727705     2
#> 2 nilHMM     B NIL01   1 29573725 308322690     0
#> 3 nilHMM     B NIL01   2    98554 243484148     0
#>   source donor  name chr start_bp    end_bp state
#> 1 nilHMM     B NIL01   1    37410  27727705     2
#> 2 nilHMM     B NIL01   1 29573725 308322690     0
#> 3 nilHMM     B NIL01   2    98554 243484148     0
```

### `lbimpute` — LB-Impute port (distance-dependent transition)

A native port of LB-Impute (Fragoso et al. 2014) for biallelic
populations: a coverage-aware emission bounded by `genotypeerr`, and a
distance-dependent transition (recombination scales with the marker gap
over `recombdist`). No design priors needed.

``` r

lb <- call_ancestry(counts, caller = "lbimpute", recombdist = 1e7, genotypeerr = 0.05)
head(lb, 3)
#>   source donor  name chr start_bp    end_bp state
#> 1 nilHMM     B NIL01   1    37410  27727705     2
#> 2 nilHMM     B NIL01   1 29573725 308322690     0
#> 3 nilHMM     B NIL01   2    98554 223047189     0
```

### `fsfhap` — full-sib families

A port of TASSEL’s FSFHap. Unlike the per-line callers it **pools each
family**, so it needs a `family` grouping and a called genotype `g`. See
[`vignette("fsfhap")`](https://sawers-rellan-labs.github.io/nilhmm/articles/fsfhap.md)
for the full workflow (HapMap + pedigree input, design routing,
parallelism).

``` r

fam <- transform(gt, family = donor)          # a family grouping on the genotype table
# (a real family needs enough segregating markers; see vignette("fsfhap"))
```

## Comparing callers

Because the output schema is shared, comparison is a table join away —
e.g. how many segments and which states each caller produced on the same
simulated cohort:

``` r

summ <- function(x) c(segments = nrow(x), states = length(unique(x$state)))
rbind(nnil = summ(nnil), bbnil = summ(bbnil), catiger = summ(catiger),
      rtiger = summ(rt), binhmm = summ(bh), googa = summ(gg), atlas = summ(at),
      lbimpute = summ(lb))
#>          segments states
#> nnil           51      3
#> bbnil          50      3
#> catiger        47      3
#> rtiger         48      3
#> binhmm         32      2
#> googa          48      3
#> atlas          44      3
#> lbimpute       51      3
```

### Chromosome painting

The real payoff of the shared schema is that you can **paint the callers
on the same axes** — here the per-sample callers run above, against the
simcross **ground truth** — and eyeball whether the independent methods
recover the donor (**B**) blocks. Each row is a NIL, each column a
chromosome, and within a cell the tracks are stacked as bands (REF gold
/ HET green / ALT purple), with the noise-free truth on top.

The truth track is just
[`to_segments()`](https://sawers-rellan-labs.github.io/nilhmm/reference/to_segments.md)
on the
[`simulate_nil()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_nil.md)
table (its `state` is the noise-free donor mosaic, before
[`simulate_counts()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_counts.md)
added depth and missingness).
[`paint_calls()`](https://sawers-rellan-labs.github.io/nilhmm/reference/paint_calls.md)
then stacks it above the callers: `rbind` each track’s segments with a
column naming it, and pass that column as `track` (its levels stack
top-down, so truth goes first).

``` r

tracks <- list("nnil (gt)" = nnil, "bbnil (count)" = bbnil, catiger = catiger,
               rtiger = rt, binhmm = bh, googa = gg, atlas = at, lbimpute = lb)
comparison <- do.call(rbind, Map(function(seg, m) { seg$method <- m; seg },
                                 tracks, names(tracks)))

# simcross ground truth (the simulate_nil() table) as the top track
truth_seg <- to_segments(truth)
truth_seg$method <- "truth (simcross)"

comparison <- rbind(truth_seg, comparison)
comparison$method <- factor(comparison$method,
                            levels = c("truth (simcross)", names(tracks)))
paint_calls(comparison, track = "method", samples = sprintf("NIL%02d", 1:4))
```

![](callers_files/figure-html/paint-1.png)

The callers recover the true donor blocks — each caller’s band lines up
with the truth on top — with expected method-specific behaviour
(e.g. `binhmm` paints broader per-bin blocks). This is the same painting
used for the real coverage-sweep NILs, where the truth track is instead
an independent data source or a high-confidence call set.

For per-caller parameters and lineage, see
[`?call_ancestry`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md)
and the README caller table. \`\`\`
