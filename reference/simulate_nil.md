# Simulate the true donor mosaic of a cross via simcross

Design-driven truth generator: build the pedigree for a `"BC{n}S{m}"`
design (F1, then `n` backcrosses to the recurrent parent, then `m`
selfings), simulate it with simcross (real meiosis / crossover
interference on a genetic map), and read off each line's **true** donor
dosage at the markers. The result is a per-marker table with `state` in
`{0 REF, 1 HET, 2 ALT}` (ALT = donor) that feeds
[`to_segments()`](https://sawers-rellan-labs.github.io/nilhmm/reference/to_segments.md)
(truth segments),
[`paint_calls()`](https://sawers-rellan-labs.github.io/nilhmm/reference/paint_calls.md)
(a truth track), or
[`simulate_counts()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_counts.md)
(degrade to observed read counts). Recurrent parent = REF = A; donor =
ALT.

## Usage

``` r
simulate_nil(
  design = "BC2S2",
  n = 1L,
  map = NULL,
  n_markers = 2000L,
  chr = NULL,
  m = 10L,
  p = 0,
  seed = NULL,
  donor = "B",
  names = NULL,
  pedigree = NULL,
  nil_id = NULL
)
```

## Arguments

- design:

  Breeding design `"BC{n}S{m}"` (e.g. `"BC2S2"`, `"BC1S4"`). Ignored
  when `pedigree` is supplied.

- n:

  Number of lines to simulate (independent meioses).

- map:

  Reference genetic map with `chr`, `cm`, `bp`. `NULL` (default) uses
  the bundled
  [`load_map()`](https://sawers-rellan-labs.github.io/nilhmm/reference/load_map.md)
  (`maize_map_v5`).

- n_markers:

  Approximate marker count to sample from `map`
  ([`build_marker_grid()`](https://sawers-rellan-labs.github.io/nilhmm/reference/build_marker_grid.md)).

- chr:

  Optional integer vector restricting to a subset of chromosomes.

- m, p:

  simcross Stahl interference (`m = 10`, `p = 0` ~ maize-like).

- seed:

  Optional RNG seed for reproducibility.

- donor:

  Donor/ALT label for the output `donor` column.

- names:

  Optional length-`n` sample names (default `sim0001` ...).

- pedigree, nil_id:

  Escape hatch: a ready simcross pedigree (`id, mom, dad, sex, gen`) and
  the id to sample; bypasses `design`.

## Value

data.frame `source, donor, name, chr, pos, cm, state` (`state` 0/1/2),
one row per (line, marker), sorted by `name, chr, pos`.

## Details

Data-agnostic: defaults to the bundled B73 v5 consensus map
([`load_map()`](https://sawers-rellan-labs.github.io/nilhmm/reference/load_map.md)),
sampled to `n_markers` by
[`build_marker_grid()`](https://sawers-rellan-labs.github.io/nilhmm/reference/build_marker_grid.md);
pass your own `map` to override. Hardcodes no design parameters (the
design token drives the pedigree). For designs beyond backcross-self
(RIL, AIL, MAGIC, ...), pass a ready simcross `pedigree` + `nil_id` and
`design` is used only to label output.

## See also

[`load_map()`](https://sawers-rellan-labs.github.io/nilhmm/reference/load_map.md),
[`build_marker_grid()`](https://sawers-rellan-labs.github.io/nilhmm/reference/build_marker_grid.md),
[`simulate_counts()`](https://sawers-rellan-labs.github.io/nilhmm/reference/simulate_counts.md),
[`to_segments()`](https://sawers-rellan-labs.github.io/nilhmm/reference/to_segments.md).

## Examples

``` r
if (requireNamespace("simcross", quietly = TRUE)) {
  truth <- simulate_nil("BC2S2", n = 2, chr = 1:2, n_markers = 200, seed = 1)
  head(truth)
  to_segments(truth)          # true donor segments
}
#>    source donor    name chr  start_bp    end_bp state
#> 1     sim     B sim0001   1     37410  27810859     2
#> 2     sim     B sim0001   1  30588204 308322690     0
#> 3     sim     B sim0001   2     98554 221103863     0
#> 4     sim     B sim0001   2 223901399 232294005     2
#> 5     sim     B sim0001   2 235091541 243484148     0
#> 6     sim     B sim0002   1     37410   2814755     0
#> 7     sim     B sim0002   1   5592100  41697583     1
#> 8     sim     B sim0002   1  44474928  50029618     2
#> 9     sim     B sim0002   1  52806963 308322690     0
#> 10    sim     B sim0002   2     98554 243484148     0
```
