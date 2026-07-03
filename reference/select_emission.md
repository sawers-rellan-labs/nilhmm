# Select an emission model from sequencing depth

Selector rule (S5): saturated (\>=~20x) -\> `gt`; intermediate (~1-20x)
-\> `count`. Cost basis: BetaBinomial cost is proportional to the number
of distinct `(n, k)` pairs, i.e. to coverage; above ~20-30x a count
emission is effectively a hard call, so `gt` is equivalent and cheaper
(memory `rtiger-betabinomial-cost`).

## Usage

``` r
select_emission(depth)
```

## Arguments

- depth:

  Mean per-marker sequencing depth.

## Value

An emission spec
([`emission_count()`](https://sawers-rellan-labs.github.io/nilhmm/reference/emission_count.md)
/
[`emission_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/emission_gt.md)).

## Examples

``` r
class(select_emission(2))    # ~2x  -> count emission
#> [1] "nilHMM_emission_count" "nilHMM_emission"      
class(select_emission(30))   # ~30x -> gt (hard-call) emission
#> [1] "nilHMM_emission_gt" "nilHMM_emission"   
```
