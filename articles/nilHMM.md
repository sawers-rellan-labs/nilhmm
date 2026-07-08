# Getting started with nilHMM

**nilHMM** calls ancestry — donor introgressions — in Near-Isogenic
Lines and related backcross/full-sib populations from sequencing data.
Under the hood it is a single **duration-aware 3-state (REF / HET / ALT)
HMM engine** with two swappable axes (emission × duration); their
combinations, plus breeding-design priors, express a family of named
callers (`nnil`, `rtiger`, `binhmm`, `atlas`, `lbimpute`, `fsfhap`).
This vignette covers the common workflow; see
[`vignette("callers")`](https://sawers-rellan-labs.github.io/nilhmm/articles/callers.md)
to choose a caller,
[`vignette("engine")`](https://sawers-rellan-labs.github.io/nilhmm/articles/engine.md)
for the internals, and
[`vignette("fsfhap")`](https://sawers-rellan-labs.github.io/nilhmm/articles/fsfhap.md)
for full-sib families.

``` r

library(nilHMM)
```

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
Here we simulate a small cohort of allelic read counts so the vignette
is self-contained: six lines, two chromosomes, each carrying one donor
(ALT) introgression block on a REF background.

``` r

sim_counts <- function(n = 6L, m = 300L, n_chr = 2L, depth = 6, err = 0.01) {
  chr <- rep(seq_len(n_chr), each = m)
  pos <- as.integer(rep(seq_len(m), n_chr) * 1e5)
  do.call(rbind, lapply(seq_len(n), function(s) {
    truth <- integer(length(chr))                        # 0 REF / 1 HET / 2 ALT
    for (ch in seq_len(n_chr)) {
      idx <- which(chr == ch); w <- sample(20:60, 1L); s0 <- sample(idx, 1L)
      blk <- s0:min(max(idx), s0 + w); truth[blk] <- 2L; truth[blk[1]] <- 1L
    }
    d  <- rpois(length(chr), depth)
    na <- rbinom(length(chr), d, c(err, 0.5, 1 - err)[truth + 1L])
    data.frame(name = sprintf("NIL%02d", s), chr = chr, pos = pos,
               n_ref = d - na, n_alt = na)
  }))
}
counts <- sim_counts()
head(counts)
#>    name chr    pos n_ref n_alt
#> 1 NIL01   1 100000     6     1
#> 2 NIL01   1 200000     2     0
#> 3 NIL01   1 300000     4     0
#> 4 NIL01   1 400000     4     0
#> 5 NIL01   1 500000     7     0
#> 6 NIL01   1 600000     5     0
```

Call ancestry with the `nnil` count caller and BC2S2 design priors:

``` r

calls <- call_ancestry(counts, caller = "nnil", design = "BC2S2",
                       rrate = 1e-4, err = 0.01)
head(calls)
#>   source donor  name chr start_bp   end_bp state
#> 1 nilHMM  <NA> NIL01   1   100000 16600000     0
#> 2 nilHMM  <NA> NIL01   1 16700000 19000000     2
#> 3 nilHMM  <NA> NIL01   1 19100000 30000000     0
#> 4 nilHMM  <NA> NIL01   2   100000 29800000     0
#> 5 nilHMM  <NA> NIL01   2 29900000 30000000     1
#> 6 nilHMM  <NA> NIL02   1   100000  2100000     0
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
#> 22  1 11
```

Numeric states are `0 = REF`, `1 = HET`, `2 = ALT`. A quick
karyotype-style view of one line’s introgressions (base graphics, no
extra packages):

``` r

one <- calls[calls$name == "NIL01", ]
cols <- c("0" = "grey85", "1" = "goldenrod", "2" = "firebrick")
plot(NA, xlim = c(0, max(one$end_bp)), ylim = c(0.5, max(one$chr) + 0.5),
     xlab = "position (bp)", ylab = "chromosome", yaxt = "n",
     main = "NIL01 ancestry")
axis(2, at = sort(unique(one$chr)))
with(one, rect(start_bp, chr - 0.3, end_bp, chr + 0.3,
               col = cols[as.character(state)], border = NA))
legend("topright", fill = cols, legend = c("REF", "HET", "ALT"), bty = "n")
```

![](nilHMM_files/figure-html/plot-1.png)

## Which caller?

All callers share the REF/HET/ALT chain and the design priors; they
differ in the emission model, the duration prior, and the input they
expect.

| caller | input | when |
|----|----|----|
| `nnil` | allelic counts, or called `GT` | general-purpose NILs (low-cov skim or saturated GT) |
| `rtiger` | allelic counts | minimum-run-length (rigidity) segmentation |
| `binhmm` | allelic counts | per-bin calling for noisy/uneven coverage |
| `atlas` | recurrent/donor counts | competitive-alignment RNA-seq (GOOGA) |
| `lbimpute` | very low-cov counts (\<1×) | GBS/skim imputation (LB-Impute) |
| `fsfhap` | called `GT` + a `family` | full-sib families (TASSEL FSFHap) |

[`vignette("callers")`](https://sawers-rellan-labs.github.io/nilhmm/articles/callers.md)
gives a runnable example for each.

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
