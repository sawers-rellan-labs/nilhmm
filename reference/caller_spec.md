# Resolve a named caller into emission + duration specs

The named cells of the shared engine's grid (emission x duration; the gt
row also varies by hard-call thresholding lineage):

- `nnil` : categorical `gt` emission + geometric duration – Holland's
  original NIL caller on hard genotype calls (1/3-2/3 dosage cutoffs).

- `bbnil` : count/BetaBinomial emission + geometric duration – the
  low-coverage read-count extension of nNIL; the self-transition is the
  smoother.

- `catiger` : categorical `gt` emission + rigidity duration (S7) – the
  Holland-threshold categorical + minimum-run-length caller.

- `rtiger` : count/BetaBinomial emission + rigidity duration (S7).

- `googa` : categorical `gt` emission + geometric duration, with GOOGA
  competitive-alignment thresholds – the faithful GOOGA/Veltsos
  reproduction.

- `atlas` : categorical `gt` emission + rigidity duration, with GOOGA
  thresholds – this work's transcript caller (rigidity variant of
  `googa`). The GOOGA thresholding is applied in
  [`call_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md).

## Usage

``` r
caller_spec(
  caller = c("nnil", "bbnil", "catiger", "rtiger", "googa", "atlas"),
  rrate = 0.01,
  rigidity = NULL,
  err = 0.01,
  conc = 20,
  fit_means = FALSE,
  xrate = 0.01,
  germ = 0.05,
  gert = 0.1,
  p = 0.5,
  mr = 0.1,
  nir = 0.01,
  ...
)
```

## Arguments

- caller:

  One of `"nnil"`, `"bbnil"`, `"catiger"`, `"rtiger"`, `"googa"`,
  `"atlas"`.

- rrate:

  Geometric callers (`nnil`/`bbnil`/`googa`): expected per-marker
  recombination rate (self-stay = `1 - rrate`). Holland's nNIL sets it
  to `2 * total_cM / (100 * n_markers)`.

- rigidity:

  Rigidity callers (`catiger`/`rtiger`/`atlas`): integer minimum run
  length (e.g. `5`).

- err:

  Count-emission baseline error.

- conc:

  Count-emission BetaBinomial concentration.

- fit_means:

  EM-fit emission means (count emission; S10).

- xrate:

  Exit rate of the **rigidity duration**
  ([`duration_rigidity()`](https://sawers-rellan-labs.github.io/nilhmm/reference/duration_rigidity.md)):
  free-state (post-minimum-run) switch probability. A nilHMM construct,
  not a RTIGER parameter.

- germ, gert, p, mr, nir:

  Genotype-error rates for the gt emission (the gt callers
  `nnil`/`catiger`/`googa`/`atlas`; Holland's nNIL error model): hom
  error, het error, hom-error-\>het fraction, missing rate,
  non-informative-marker rate.

- ...:

  Ignored extra args (e.g. `f_1`/`f_2` consumed by
  [`call_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md)).

## Value

`list(emission, duration)`.

## Examples

``` r
caller_spec("nnil", rrate = 1e-4)      # gt emission + geometric duration
#> $emission
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
#> 
#> $duration
#> $type
#> [1] "geometric"
#> 
#> $r
#> [1] 1e-04
#> 
#> attr(,"class")
#> [1] "nilHMM_duration_geometric" "nilHMM_duration"          
#> 
caller_spec("bbnil", rrate = 1e-4)     # count emission + geometric duration
#> $emission
#> $type
#> [1] "count"
#> 
#> $err
#> [1] 0.01
#> 
#> $conc
#> [1] 20
#> 
#> $fit_means
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "nilHMM_emission_count" "nilHMM_emission"      
#> 
#> $duration
#> $type
#> [1] "geometric"
#> 
#> $r
#> [1] 1e-04
#> 
#> attr(,"class")
#> [1] "nilHMM_duration_geometric" "nilHMM_duration"          
#> 
caller_spec("catiger", rigidity = 5)   # gt emission + rigidity duration
#> $emission
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
#> 
#> $duration
#> $type
#> [1] "rigidity"
#> 
#> $r
#> [1] 5
#> 
#> $p_switch
#> [1] 0.01
#> 
#> attr(,"class")
#> [1] "nilHMM_duration_rigidity" "nilHMM_duration"         
#> 
caller_spec("rtiger", rigidity = 5)    # count emission + rigidity duration
#> $emission
#> $type
#> [1] "count"
#> 
#> $err
#> [1] 0.01
#> 
#> $conc
#> [1] 20
#> 
#> $fit_means
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "nilHMM_emission_count" "nilHMM_emission"      
#> 
#> $duration
#> $type
#> [1] "rigidity"
#> 
#> $r
#> [1] 5
#> 
#> $p_switch
#> [1] 0.01
#> 
#> attr(,"class")
#> [1] "nilHMM_duration_rigidity" "nilHMM_duration"         
#> 
str(caller_spec("googa"))              # gt emission + geometric (faithful GOOGA)
#> List of 2
#>  $ emission:List of 6
#>   ..$ type: chr "gt"
#>   ..$ germ: num 0.05
#>   ..$ gert: num 0.1
#>   ..$ p   : num 0.5
#>   ..$ mr  : num 0.1
#>   ..$ nir : num 0.01
#>   ..- attr(*, "class")= chr [1:2] "nilHMM_emission_gt" "nilHMM_emission"
#>  $ duration:List of 2
#>   ..$ type: chr "geometric"
#>   ..$ r   : num 0.01
#>   ..- attr(*, "class")= chr [1:2] "nilHMM_duration_geometric" "nilHMM_duration"
str(caller_spec("atlas", rigidity = 5))# gt emission + rigidity (this-work transcript)
#> List of 2
#>  $ emission:List of 6
#>   ..$ type: chr "gt"
#>   ..$ germ: num 0.05
#>   ..$ gert: num 0.1
#>   ..$ p   : num 0.5
#>   ..$ mr  : num 0.1
#>   ..$ nir : num 0.01
#>   ..- attr(*, "class")= chr [1:2] "nilHMM_emission_gt" "nilHMM_emission"
#>  $ duration:List of 3
#>   ..$ type    : chr "rigidity"
#>   ..$ r       : int 5
#>   ..$ p_switch: num 0.01
#>   ..- attr(*, "class")= chr [1:2] "nilHMM_duration_rigidity" "nilHMM_duration"
```
