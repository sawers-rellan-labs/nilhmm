# Genotype-frequency prior implied by a breeding design

The public design-prior contract: the `c(REF, HET, ALT)` vector of
single-locus genotype frequencies the breeding scheme implies, consumed
directly by
[`call_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_gt.md)
as its `prior` and by the engine as its state frequencies. A thin
wrapper over
[`single_locus_expectation()`](https://sawers-rellan-labs.github.io/nilhmm/reference/single_locus_expectation.md)
– callers should depend on this stable name rather than on the genetics
primitive behind it.

## Usage

``` r
breeding_prior(design, f1 = c(0, 1, 0))
```

## Arguments

- design:

  Design key of the form `"BC<n>S<m>"` (e.g. `"BC2S2"`, `"BC2S3"`).

- f1:

  Starting F1 genotype-frequency vector, passed through to
  [`single_locus_expectation()`](https://sawers-rellan-labs.github.io/nilhmm/reference/single_locus_expectation.md);
  defaults to inbred parents `c(0, 1, 0)`.

## Value

Named numeric length-3 vector `c(REF, HET, ALT)` summing to 1.

## See also

[`single_locus_expectation()`](https://sawers-rellan-labs.github.io/nilhmm/reference/single_locus_expectation.md),
[`call_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_gt.md)

## Examples

``` r
breeding_prior("BC2S3")                          # c(REF = .8594, HET = .0312, ALT = .1094)
#>      REF      HET      ALT 
#> 0.859375 0.031250 0.109375 
call_gt(0, 1, prior = breeding_prior("BC2S3"))   # design prior resists the het flip -> 2
#> [1] 2
```
