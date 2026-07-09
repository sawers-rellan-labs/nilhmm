# Design prior as a length-3 genotype-frequency vector

Thin accessor over
[`design_priors()`](https://sawers-rellan-labs.github.io/nilhmm/reference/design_priors.md)
that returns the single-locus genotype frequencies as the
`c(f_REF, f_HET, f_ALT)` vector consumed directly by
[`call_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_gt.md)
as its `prior`. Keeps the design prior **derived, never typed**.

## Usage

``` r
design_prior(design)
```

## Arguments

- design:

  Design key (e.g. `"BC2S2"`, `"BC2S3"`), as for
  [`design_priors()`](https://sawers-rellan-labs.github.io/nilhmm/reference/design_priors.md).

## Value

Named numeric length-3 vector `c(REF, HET, ALT)` summing to 1.

## See also

[`design_priors()`](https://sawers-rellan-labs.github.io/nilhmm/reference/design_priors.md),
[`call_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_gt.md)

## Examples

``` r
design_prior("BC2S3")                          # c(REF = .8594, HET = .0312, ALT = .1094)
#>    REF    HET    ALT 
#> 0.8594 0.0312 0.1094 
call_gt(0, 1, prior = design_prior("BC2S3"))   # design prior resists the het flip -> 2
#> [1] 2
```
