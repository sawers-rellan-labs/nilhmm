# Single-locus (Mendelian) genotype expectation for a breeding design

Propagates the single-locus genotype-frequency vector through the
backcross and selfing transition matrices for the design named by
`design` (`"BC<n>S<m>"`): `n` backcrosses to the recurrent parent, then
`m` generations of selfing. The genotype basis is `(AA, Aa, aa)` with
`AA` the recurrent (backcross-target) homozygote and `aa` the donor
homozygote, so the recurrent allele is `A` (frequency `p_A`) and the
donor allele is `a`. The donor allele frequency after the backcrosses is
`p_a = (aa + Aa/2) * 0.5^n` in terms of the starting F1 composition
`f1 = (AA, Aa, aa)` – `0.5^(n + 1)` for the default fully-het F1 – and
selfing leaves it invariant. The backcross matrix `B` mates each
generation to `AA` (`aa -> Aa`, `Aa -> half AA / half Aa`) and the
selfing matrix `S` splits a quarter of the remaining `Aa` into equal
parts `AA`/`aa` each generation.

## Usage

``` r
single_locus_expectation(design, f1 = c(0, 1, 0))
```

## Arguments

- design:

  Design key of the form `"BC<n>S<m>"` (e.g. `"BC2S2"`, `"BC1S4"`).

- f1:

  Starting genotype-frequency vector in `(AA, Aa, aa)` order (recurrent
  hom, het, donor hom). Defaults to `c(0, 1, 0)`, the fully-heterozygous
  F1 of inbred parents.

## Value

Named numeric vector `c(REF, HET, ALT)` summing to 1, mapping the
`(AA, Aa, aa)` genotypes to the REF/HET/ALT ancestry states.

## Details

The starting composition `f1` is a free argument, so the expectation is
exact for a non-inbred F1, not only the fully-heterozygous `c(0, 1, 0)`
case that inbred parents give. At the inbred default, BC2S2 lands at REF
0.844 / HET 0.0625 / ALT 0.0938.

## See also

[`parse_design()`](https://sawers-rellan-labs.github.io/nilhmm/reference/parse_design.md),
[`breeding_prior()`](https://sawers-rellan-labs.github.io/nilhmm/reference/breeding_prior.md)

## Examples

``` r
single_locus_expectation("BC2S2")                       # inbred parents
#>     REF     HET     ALT 
#> 0.84375 0.06250 0.09375 
single_locus_expectation("BC2S2", f1 = c(0.5, 0.5, 0))  # non-inbred F1
#>      REF      HET      ALT 
#> 0.921875 0.031250 0.046875 
```
