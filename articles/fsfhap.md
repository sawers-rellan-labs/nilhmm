# Full-sib families with the fsfhap caller

`fsfhap` is a native port of TASSEL’s **FSFHap** (Swarts et al. 2014).
Unlike the per-line callers, it **pools each full-sib family**: the
scattered introgressions across the family’s lines collectively phase
the donor haplotype, then a 5-state Viterbi-training EM imputes and
smooths each line. It is bit-exact vs TASSEL on its intended populations
(backcross and RIL/inbred). It needs two things the per-line callers
don’t: a called genotype `g`, and a **`family` grouping**.

``` r

library(nilHMM)
```

## Input: HapMap + pedigree

FSFHap’s native inputs are a TASSEL **HapMap** and a 7-column
**pedigree**.
[`read_hapmap()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_hapmap.md)
and
[`read_pedigree()`](https://sawers-rellan-labs.github.io/nilhmm/reference/read_pedigree.md)
normalise them to the engine’s canonical tables. Here we write tiny
example files to demonstrate the readers (in practice you point them at
your files):

``` r

hmp <- tempfile(fileext = ".hmp.txt")
writeLines(c(
  "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\tL1\tL2",
  "m1\tA/C\t1\t100\t+\tNA\tNA\tNA\tNA\tNA\tNA\tA\tC",
  "m2\tA/C\t1\t200\t+\tNA\tNA\tNA\tNA\tNA\tNA\tM\tA"), hmp)
read_hapmap(hmp)          # -> name, chr, pos, g   (g in {0 A-hom, 1 het, 2 C-hom, 3 NA})
#>   name chr pos g
#> 1   L1   1 100 0
#> 2   L1   1 200 1
#> 3   L2   1 100 2
#> 4   L2   1 200 0

ped_file <- tempfile()
writeLines(c("family\ttaxon\tparent1\tparent2\tp1\tp2\tF",
             "TIL01\tL1\tW22\tTIL01_teo\t0.75\t0.25\t0.9375",
             "TIL01\tL2\tW22\tTIL01_teo\t0.75\t0.25\t0.9375"), ped_file)
read_pedigree(ped_file)   # -> taxon, family, parents, contribution1/2, F
#>   taxon family parent1   parent2 contribution1 contribution2      F
#> 1    L1  TIL01     W22 TIL01_teo          0.75          0.25 0.9375
#> 2    L2  TIL01     W22 TIL01_teo          0.75          0.25 0.9375
```

You attach the family grouping by joining on the taxon name:

``` r

data <- read_hapmap("family.hmp.txt")
ped  <- read_pedigree("family_pedigree.txt")
data$family <- ped$family[match(data$name, ped$taxon)]
stopifnot(!anyNA(data$family))            # every sample must be in the pedigree
```

## Design routing derives `phet`

Pass a `design` token (`BC{n}S{m}`); the dispatcher selects the
parent-calling route and **derives the expected heterozygosity**
`phet = (1 - F)/2` from the design’s inbreeding coefficient — no magic
constant. BC1S4 (the TeoNAM design) gives `F = 0.9375`, so
`phet = 0.03125`. `fsfhap` supports BC1 designs end-to-end; other
designs route to as-yet-unported TASSEL paths.

## Calling on a family

A real family needs enough segregating markers to phase the donor (and
to clear FSFHap’s coverage floor), so we simulate a 24-line **BC1S4**
family with [**simcross**](https://cran.r-project.org/package=simcross)
— the TeoNAM design: one backcross to the recurrent parent, then four
selfings (`F = 0.9375`, so `phet = 0.03125`). simcross gives real
recombination breakpoints; the donor blocks and residual heterozygosity
fall out of the crossing scheme rather than being placed by hand.

``` r

library(simcross)
sim_family <- function(n_lines = 24L, m = 500L, L = 100, nself = 4L, miss = 0.1) {
  rec <- create_parent(L, 1); don <- create_parent(L, 2)
  one_line <- function() {
    ind <- cross(cross(rec, don), rec)                  # F1 -> backcross to recurrent (BC1)
    for (s in seq_len(nself)) ind <- cross(ind, ind)    # four selfings -> BC1S4
    ind
  }
  map <- seq(0, L, length.out = m)
  g <- as.vector(convert2geno(lapply(seq_len(n_lines), function(i) one_line()), map)) - 1L
  g[runif(length(g)) < miss] <- 3L                      # skeleton missingness
  data.frame(name = rep(sprintf("L%02d", seq_len(n_lines)), times = m),
             family = "TIL01", chr = 1L,
             pos = rep(as.integer(map * 1e6), each = n_lines), g = as.integer(g))
}
fam <- sim_family()

calls <- call_ancestry(fam, caller = "fsfhap", design = "BC1S4")
head(calls, 4)
#>   source donor name chr start_bp   end_bp state
#> 1 nilHMM  <NA>  L01   1        0 16232464     0
#> 2 nilHMM  <NA>  L01   1 16633266 45891783     2
#> 3 nilHMM  <NA>  L01   1 46292585 98797595     0
#> 4 nilHMM  <NA>  L02   1        0  5410821     2
table(calls$state)
#> 
#>  0  1  2 
#> 32 11 25
```

## Parallelism

Each **family × chromosome** is an independent unit, so `fsfhap` fans
them out across cores with `threads` (bit-identical to serial). On a
real full-genome cohort this is a large win — the family × chromosome
units fill all cores:

``` r

identical(
  call_ancestry(fam, caller = "fsfhap", design = "BC1S4", threads = 1L),
  call_ancestry(fam, caller = "fsfhap", design = "BC1S4", threads = 2L)
)
#> [1] TRUE
```

On real TeoNAM (all five families, 51,004 markers, 50 units) the port
runs in **28.6 s on 10 threads vs TASSEL’s 218.5 s** — 7.6× faster, and
faster than TASSEL even single-threaded. See `design/FSFHAP_PORT.md` for
the full benchmark.

## When *not* to use `fsfhap`

`fsfhap` assumes homozygous founders and balanced full-sib segregation.
For NILs with **heterozygous or outbred donors** (e.g. teosinte), the
two-haplotype model breaks down — use the donor-agnostic `nnil` caller
instead. Het-heavy F2 is a documented stress case. See
[`vignette("callers")`](https://sawers-rellan-labs.github.io/nilhmm/articles/callers.md)
for the alternatives and
[`?call_ancestry`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md)
for the full parameter list. \`\`\`
