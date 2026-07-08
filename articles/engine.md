# The engine: one chain, two swappable axes

[`call_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md)
is a thin convenience layer. Underneath, every count/geometric caller is
the **same 3-state (REF/HET/ALT) HMM** parameterised along two
independent axes plus a set of design priors:

    caller  =  emission model  ×  duration prior  +  design priors

- **Emission** — how an observation maps to a state likelihood:
  [`emission_count()`](https://sawers-rellan-labs.github.io/nilhmm/reference/emission_count.md)
  (BetaBinomial over ref/alt depths) or
  [`emission_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/emission_gt.md)
  (categorical genotype-error model).
- **Duration** — the transition prior on run lengths:
  [`duration_geometric()`](https://sawers-rellan-labs.github.io/nilhmm/reference/duration_geometric.md),
  [`duration_rigidity()`](https://sawers-rellan-labs.github.io/nilhmm/reference/duration_rigidity.md)
  (minimum run length), or
  [`duration_hsmm()`](https://sawers-rellan-labs.github.io/nilhmm/reference/duration_hsmm.md).
- **Design priors** —
  [`design_priors()`](https://sawers-rellan-labs.github.io/nilhmm/reference/design_priors.md)
  turns a breeding design into the state frequencies (`f_1`, `f_2`).

This vignette drives the engine directly with
[`fit()`](https://sawers-rellan-labs.github.io/nilhmm/reference/fit.md)
→
[`decode()`](https://sawers-rellan-labs.github.io/nilhmm/reference/decode.md)
→
[`to_segments()`](https://sawers-rellan-labs.github.io/nilhmm/reference/to_segments.md),
which is exactly what
[`call_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md)
chains for you.

``` r

library(nilHMM)
```

## The building blocks

Each axis is a small constructor returning a spec:

``` r

design_priors("BC2S2")     # state frequencies from the breeding design
#> $g
#> [1] "BC2S2"
#> 
#> $f_1
#> [1] 0.0625
#> 
#> $f_2
#> [1] 0.0938
emission_count(err = 0.01, conc = 20)
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
duration_geometric(rrate = 1e-4)
#> $type
#> [1] "geometric"
#> 
#> $r
#> [1] 1e-04
#> 
#> attr(,"class")
#> [1] "nilHMM_duration_geometric" "nilHMM_duration"
```

## Fit and decode by hand

[`fit()`](https://sawers-rellan-labs.github.io/nilhmm/reference/fit.md)
takes a per-(sample, chromosome) observation frame — `n` = depth, `a` =
alt-read count — an emission spec, a duration spec, and priors, and
returns a fitted model.
[`decode()`](https://sawers-rellan-labs.github.io/nilhmm/reference/decode.md)
runs Viterbi and returns one **state per observation** (coordinate-free:
`0/1/2`).

We simulate one line with
[**simcross**](https://cran.r-project.org/package=simcross) (a BC2S2
backcross of a donor **B** onto a recurrent parent **A**, on an 80 cM
chromosome) and layer read counts, so the observations come from a real
recombination track:

``` r

library(simcross)
map <- seq(0, 80, length.out = 24)
rec <- create_parent(80, 1); don <- create_parent(80, 2)
bc2s2 <- function() {
  ind <- cross(rec, don)                       # F1
  for (b in 1:2) ind <- cross(ind, rec)        # two backcrosses to recurrent
  for (s in 1:2) ind <- cross(ind, ind)        # two selfings
  ind
}
G <- convert2geno(lapply(1:12, function(i) bc2s2()), map) - 1L   # lines x markers, 0/1/2
g <- G[which.max(rowMeans(G == 2L)), ]                           # the line with the most donor
depth <- rpois(length(g), 8)
obs   <- data.frame(n = depth, a = rbinom(length(g), depth, c(0.01, 0.5, 0.99)[g + 1L]))

model <- fit(obs, emission_count(err = 0.01),
             duration_geometric(1e-4), priors = design_priors("BC2S2"))
rbind(simulated = g, decoded = decode(model, obs))               # decode recovers the track
#>           [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#> simulated    2    2    2    2    2    0    0    0    0     0     0     0     0
#> decoded      2    2    2    2    2    0    0    0    0     0     0     0     0
#>           [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24]
#> simulated     0     0     2     2     2     2     0     0     0     0     0
#> decoded       0     0     2     2     2     2     0     0     0     0     0
```

## Coordinates back on: `to_segments()`

[`decode()`](https://sawers-rellan-labs.github.io/nilhmm/reference/decode.md)
is coordinate-free. Attach `name/chr/pos` and run
[`to_segments()`](https://sawers-rellan-labs.github.io/nilhmm/reference/to_segments.md)
to collapse the per-marker states into the common segment schema (this
is the second half of
[`call_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md)):

``` r

states <- data.frame(name = "S1", donor = "B", chr = 1L,
                     pos = as.integer(map * 1e6),      # cM map -> bp coordinates
                     state = decode(model, obs))
to_segments(states)
#>   source donor name chr start_bp   end_bp state
#> 1 nilHMM     B   S1   1        0 13913043     2
#> 2 nilHMM     B   S1   1 17391304 48695652     0
#> 3 nilHMM     B   S1   1 52173913 62608695     2
#> 4 nilHMM     B   S1   1 66086956 80000000     0
```

## Swapping an axis

Because the axes are independent, changing the caller’s behaviour is a
matter of swapping one spec. Keep the same count emission but replace
the geometric duration with a **rigidity** prior (minimum run length 4)
— the RTIGER-style segmentation. Rigidity enforces a minimum run length,
so it resists short state changes; on this clean, well-covered line the
call is essentially unchanged, but on noisy low-coverage data it
suppresses the single-marker flicker that a geometric prior would
otherwise emit as spurious one-marker segments.

``` r

rigid <- fit(obs, emission_count(err = 0.01),
             duration_rigidity(rigidity = 4L), priors = design_priors("BC2S2"))
decode(rigid, obs)
#>  [1] 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0 2 2 2 2 0 0 0 0 0
```

Or keep the duration and swap the **emission** to the categorical
genotype model
([`emission_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/emission_gt.md))
when you have hard calls rather than read counts — the same axis
machinery, a different likelihood.

## Where the named callers fit

[`caller_spec()`](https://sawers-rellan-labs.github.io/nilhmm/reference/caller_spec.md)
bundles these axis choices into the named presets, and
[`call_states()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_states.md)
/
[`call_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md)
resolve a caller to its `(emission, duration, priors)` and run the
`fit`/`decode`/`to_segments` loop per (sample, chromosome) for you:

``` r

spec <- caller_spec("nnil", design = "BC2S2")
names(spec)
#> [1] "emission" "duration"
```

So `nnil` is `emission_count × duration_geometric`, `rtiger` is
`emission_count × duration_rigidity`, and the genotype callers swap in
`emission_gt`. The `binhmm`, `lbimpute`, and `fsfhap` callers carry
extra structure (binning, a coverage-aware emission, family pooling) and
have their own drivers, but they emit the same segment schema. See
[`vignette("callers")`](https://sawers-rellan-labs.github.io/nilhmm/articles/callers.md)
for a runnable example of each and
[`?caller_spec`](https://sawers-rellan-labs.github.io/nilhmm/reference/caller_spec.md)
for the axis defaults. \`\`\`
