# Interpolate a complete genotype block onto a target cM grid (single chr)

Linear flanking-marker interpolation in genetic distance for one
chromosome. For a target at cM `t` between called flanking markers L (at
`obs_cm[jj]`) and R (at `obs_cm[jj+1]`), the weight is
`w = (t - cM_L) / (cM_R - cM_L)` and the continuous dosage is
`vL + w*(vR - vL)`. Ends are clamped to the terminal observed value
([`stats::approx`](https://rdrr.io/r/stats/approxfun.html) rule = 2);
tied flanking positions (`denom == 0`) collapse to `w = 0`.

## Usage

``` r
interp_geno_cpp(obs_cm, G, target_cm, mode)
```

## Arguments

- obs_cm:

  Numeric cM of the observed markers, ascending, length k.

- G:

  Numeric k x n genotype matrix (row = observed marker, col = sample),
  COMPLETE (no NA); values are the alt/teosinte-allele dosage in 0 to 2.

- target_cm:

  Numeric cM of the target grid, ascending, length M.

- mode:

  Interpolation mode: 0 = continuous dosage ramp (Tian 2011), 1 = step /
  nearest flanking value (`w < 0.5 ? vL : vR`, tie `w == 0.5` -\> vR;
  Chen/TeoNAM densification), 2 = round(continuous) to 0/1/2, 3 = Chen
  2019 composite-map rule (concordant flanks fill, discordant flanks or
  chromosome ends -\> NA; distance-independent).

## Value

Numeric M x n matrix of interpolated genotypes. Mode 3 may contain NA
(discordant flanks / chromosome ends); modes 0-2 never introduce NA.
