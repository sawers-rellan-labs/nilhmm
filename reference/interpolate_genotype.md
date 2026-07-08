# Interpolate genotypes onto a target marker grid (Tian 2011 / Chen 2019)

Densifies a complete genotype block onto a common target marker grid by
flanking-marker interpolation in genetic distance (cM). For a target
between called flanking markers L and R, the interpolation weight is
`w = (cM_x - cM_L) / (cM_R - cM_L)` and the continuous dosage is
`v_L + w * (v_R - v_L)` – the Tian 2011 rule (the weighted average
interpreted as the probability the SNP comes from the non-recurrent
parent). Ends are clamped to the terminal observed value (equivalent to
[`stats::approx()`](https://rdrr.io/r/stats/approxfun.html) with
`rule = 2`).

## Usage

``` r
interpolate_genotype(
  geno,
  obs,
  target,
  mode = c("continuous", "step", "round", "chen2019"),
  coord = "cm"
)
```

## Arguments

- geno:

  Numeric matrix, `nrow(geno) == nrow(obs)`, column = sample; COMPLETE
  (no `NA`). Values are the alt-allele dosage in `[0, 2]`.

- obs:

  `data.frame` with columns `chr` and the coordinate column named by
  `coord` (default `cm`), aligned row-for-row to `geno`, sorted by
  `(chr, coord)` with `coord` strictly increasing within each
  chromosome.

- target:

  `data.frame` with columns `chr` and the coordinate column named by
  `coord`, sorted by `(chr, coord)`. The coordinate may repeat (tied
  positions allowed); each row yields one output row, and markers
  sharing a coordinate get identical genotypes (see Details).

- mode:

  One of `"continuous"` (Tian 2011 dosage ramp), `"step"`, `"round"`
  (all three = distance-based JLM/GWAS densification), or `"chen2019"`
  (the composite genetic-map rule: concordant flanks fill, discordant
  flanks or chromosome ends -\> `NA`; see Details).

- coord:

  Name of the coordinate column to interpolate along; default `"cm"`
  (genetic distance). Use e.g. `"bp"` for physical-distance
  interpolation – appropriate when interpolating in cM would be
  circular, such as building a native genetic map from bp positions. The
  coordinate is just the axis to interpolate along; the arithmetic is
  unit-agnostic.

## Value

Numeric matrix `nrow(target)` x `ncol(geno)`; rownames from `target` (if
any), colnames from `geno`. Modes `continuous`/`step`/`round` never
introduce `NA`; mode `chen2019` returns `NA` at discordant-flank targets
and chromosome ends.

## Details

Genotypes are the alt/teosinte-allele dosage: `0` = REF/recurrent hom,
`1` = HET, `2` = ALT/donor hom (`continuous` returns real values in
`[0, 2]`). Three modes share one interpolation core:

- `continuous`:

  dosage ramp (Tian 2011); returns fractional dosage.

- `step`:

  nearest flanking value (`w < 0.5 ? v_L : v_R`, hard step at the cM
  midpoint, tie `w == 0.5` -\> `v_R`); the faithful Chen/TeoNAM discrete
  form – fabricates no het across a `0<->2` gap.

- `round`:

  `round(continuous)` to 0/1/2; fabricates a het band across a `0<->2`
  gap (width proportional to the cM gap) – provided to demonstrate that
  artifact.

- `chen2019`:

  the composite genetic-MAP construction rule of Chen et al. 2019
  (TeoNAM): *"the missing genotypes were imputed according to the
  flanking markers. If the flanking markers had same genotypes, the
  missing genotype was imputed as the same with flanking markers, or
  otherwise left as missing."* So concordant flanks fill; **discordant
  flanks -\> NA**; and targets past the terminal observed marker
  (chromosome ends) -\> **NA** (no rule=2 clamping). Order-based: the
  coordinate only locates the flanks, no distance weight enters the
  decision. Use this for MAP construction; use
  `continuous`/`step`/`round` for the Tian 2011 JLM/GWAS densification,
  which imputes by the genetic distance of the flanking markers. Unlike
  the other three modes, `chen2019` **can emit `NA`**.

The function densifies **one dense block** to a target grid. For
block-sparse per-family data (each family genotyped on its own marker
subset), loop over families – interpolating each family's dense block
onto the shared grid – and `cbind` the results; that loop lives in the
consumer pipeline, not here.

Tied target positions: `obs$cm` must be strictly increasing (source
flanks must be distinct), but **`target$cm` may repeat** – the target
grid does NOT need unique cM. Every target row gets its own output row,
so pass the full marker set to densify all of it; do not pre-collapse to
unique cM (that silently drops markers). Two target markers at the same
cM are, by the map, at one genetic location (zero modelled recombination
between them, i.e. perfect LD), so they receive an **identical genotype
vector** – deterministic, and correct for the interpolation model, not a
bug. This is common in recombination-cold centromeric regions (a whole
Mb mapping to one cM). If a downstream step needs non-duplicate columns,
thin the twins there
([`select_independent()`](https://sawers-rellan-labs.github.io/nilhmm/reference/select_independent.md)
on an
`r2`/[`position_distance()`](https://sawers-rellan-labs.github.io/nilhmm/reference/position_distance.md)
matrix); interpolation cannot manufacture sub-cM resolution the map
lacks.

`geno` columns are assumed COMPLETE (no `NA`) – true for imputed truth
and for HMM-caller output. There is no `NA`-aware flanking search.

## Examples

``` r
obs <- data.frame(chr = 1L, cm = c(0, 1))
geno <- matrix(c(0, 2), nrow = 2, dimnames = list(NULL, "S1"))
target <- data.frame(chr = 1L, cm = c(0, 0.5, 1))
interpolate_genotype(geno, obs, target, "continuous")  # 0, 1, 2
#>      S1
#> [1,]  0
#> [2,]  1
#> [3,]  2

# Interpolate along physical distance (bp) instead of cM.
obs_bp    <- data.frame(chr = 1L, bp = c(1e6, 3e6))
target_bp <- data.frame(chr = 1L, bp = c(1e6, 2e6, 3e6))
interpolate_genotype(geno, obs_bp, target_bp, "continuous", coord = "bp")  # 0, 1, 2
#>      S1
#> [1,]  0
#> [2,]  1
#> [3,]  2

# Chen 2019 composite-map rule: concordant flanks fill, discordant/ends -> NA.
obs2    <- data.frame(chr = 1L, cm = c(0, 2))
geno2   <- matrix(c(0, 0,  0, 2), nrow = 2,
                  dimnames = list(NULL, c("concord", "discord")))
target2 <- data.frame(chr = 1L, cm = c(0, 1, 2, 3))
interpolate_genotype(geno2, obs2, target2, "chen2019")
#>      concord discord
#> [1,]       0       0
#> [2,]       0      NA
#> [3,]       0       2
#> [4,]      NA      NA
#            concord discord
# cm0 (obs)        0       0
# cm1 (mid)        0      NA   # discordant flanks -> NA
# cm2 (obs)        0       2
# cm3 (end)       NA      NA   # chromosome end -> NA
```
