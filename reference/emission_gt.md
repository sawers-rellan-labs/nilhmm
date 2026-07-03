# Genotype (categorical) emission

Categorical over `{0, 1, 2, missing}` with a genotype-error matrix. The
cheap equivalent of `count` at saturated depth (Holland's native GT
path; the MolBreeding regime). Hard call, so no `(n, k)` cost blow-up.

## Usage

``` r
emission_gt(germ = 0.05, gert = 0.1, p = 0.5, mr = 0.1, nir = 0.01)
```

## Arguments

- germ:

  Error rate on true homozygotes.

- gert:

  Error rate on true heterozygotes.

- p:

  Fraction of homozygous errors that call as heterozygous.

- mr:

  Missing-genotype rate.

- nir:

  Non-informative-marker rate.

## Value

An emission spec for
[`fit()`](https://sawers-rellan-labs.github.io/nilhmm/reference/fit.md).

## Examples

``` r
emission_gt(germ = 0.05, gert = 0.10)
#> $type
#> [1] "gt"
#> 
#> $germ
#> [1] 0.05
#> 
#> $gert
#> [1] 0.1
#> 
#> $p
#> [1] 0.5
#> 
#> $mr
#> [1] 0.1
#> 
#> $nir
#> [1] 0.01
#> 
#> attr(,"class")
#> [1] "nilHMM_emission_gt" "nilHMM_emission"   
```
