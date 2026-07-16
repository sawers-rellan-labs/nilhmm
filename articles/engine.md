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

We use
[`simulate_nil()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_nil.md)
(a BC2S2 backcross of a donor **B** onto a recurrent parent **A**, on
chromosome 1 of the bundled map) and
[`simulate_counts()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_counts.md),
so the observations come from a real recombination track:

``` r

# a BC2S2 cohort on chromosome 1; take the line with the most donor as our example
truth <- simulate_nil("BC2S2", n = 12, chr = 1, n_markers = 24, donor = "B", seed = 7)
line  <- truth[truth$name == names(which.max(tapply(truth$state == 2L, truth$name, mean))), ]
obs   <- simulate_counts(line, depth = 8, seed = 7)
o     <- data.frame(n = obs$n_ref + obs$n_alt, a = obs$n_alt)    # fit()/decode(): n = depth, a = alt

model <- fit(o, emission_count(err = 0.01),
             duration_geometric(1e-4), priors = design_priors("BC2S2"))
rbind(simulated = line$state, decoded = decode(model, o))        # decode recovers the track
#>           [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#> simulated    0    0    2    2    2    2    0    0    0     0     0     0     0
#> decoded      0    0    2    2    2    2    0    0    0     0     0     0     0
#>           [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24]
#> simulated     0     2     2     0     0     0     1     0     0     0     0
#> decoded       0     2     2     0     0     0     0     0     0     0     0
```

## Coordinates back on: `to_segments()`

[`decode()`](https://sawers-rellan-labs.github.io/nilhmm/reference/decode.md)
is coordinate-free. Attach `name/chr/pos` and run
[`to_segments()`](https://sawers-rellan-labs.github.io/nilhmm/reference/to_segments.md)
to collapse the per-marker states into the common segment schema (this
is the second half of
[`call_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md)):

``` r

states <- data.frame(name = "S1", donor = "B", chr = line$chr, pos = line$pos,
                     state = decode(model, o))
to_segments(states)
#>   source donor name chr  start_bp    end_bp state
#> 1 nilHMM     B   S1   1     37410  13441118     0
#> 2 nilHMM     B   S1   1  26844826  67055950     2
#> 3 nilHMM     B   S1   1  80459657 174285612     0
#> 4 nilHMM     B   S1   1 187689320 201093028     2
#> 5 nilHMM     B   S1   1 214496736 308322690     0
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

rigid <- fit(o, emission_count(err = 0.01),
             duration_rigidity(rigidity = 4L), priors = design_priors("BC2S2"))
decode(rigid, o)
#>  [1] 1 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0
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

So the four grid cells are `nnil` = `emission_gt × duration_geometric`,
`bbnil` = `emission_count × duration_geometric`, `catiger` =
`emission_gt × duration_rigidity`, and `rtiger` =
`emission_count × duration_rigidity`. The `binhmm`, `lbimpute`, and
`fsfhap` callers carry extra structure (binning, a coverage-aware
emission, family pooling) and have their own drivers, but they emit the
same segment schema. See
[`vignette("callers")`](https://sawers-rellan-labs.github.io/nilhmm/articles/callers.md)
for a runnable example of each and
[`?caller_spec`](https://sawers-rellan-labs.github.io/nilhmm/reference/caller_spec.md)
for the axis defaults. \`\`\`
