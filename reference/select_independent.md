# Select a maximal independent set of markers (LD thinning)

Choose a large set of markers in which no two are more related than
`threshold` – a maximal independent set on the relatedness graph, via a
faithful port of FastIndep (deterministic greedy heuristic + optional
stochastic runs). Intended for pruning a densified genotype matrix per
chromosome before joint-linkage mapping.

## Usage

``` r
select_independent(
  x,
  threshold,
  n_runs = 1L,
  seed = 1L,
  method = c("r2", "mi", "vi"),
  sense = c("auto", "similarity", "distance"),
  max_markers = 7000L,
  ...
)
```

## Arguments

- x:

  Genotype matrix (markers x samples) or a precomputed square
  relatedness matrix (see Details).

- threshold:

  Relatedness cutoff defining an edge (sense per Details).

- n_runs:

  Total runs: `1` = greedy only; `> 1` adds `n_runs - 1` stochastic runs
  and records the distinct sets found.

- seed:

  Seed for the stochastic runs.

- method:

  Measure for the genotype-matrix path and for the edge sense when `x`
  has no `attr(, "kind")`: one of `"r2"` (default), `"mi"`, `"vi"`.

- sense:

  Edge/threshold sense: `"auto"` (default) infers it from
  `attr(x, "kind")` then `method` (current behaviour); `"similarity"`
  (edge if `>= threshold`) or `"distance"` (edge if `<= threshold`)
  force it, so any precomputed matrix – e.g. a cM distance from
  [`position_distance()`](https://sawers-rellan-labs.github.io/nilhmm/reference/position_distance.md)
  – can be pruned with the correct sense regardless of `method` or
  attributes.

- max_markers:

  Safety guard on the marker count (default `7000L`, the validated
  scale). Selection builds/consumes an O(n^2) markers x markers matrix,
  so more than `max_markers` markers errors cleanly before any large
  allocation (never a `malloc` failure); large-but-allowed sizes warn.
  Raise it (up to `getOption("nilHMM.marker_hard_cap")`, default
  `30000`) only with the memory to spare – the intended workflow prunes
  per chromosome.

- ...:

  Passed to
  [`pairwise_distance()`](https://sawers-rellan-labs.github.io/nilhmm/reference/pairwise_distance.md)
  on the genotype path (e.g. `base`).

## Value

Character vector of the selected marker names (the largest independent
set), or 1-based integer indices when `x` has no names. Attributes:
`attr(, "sets")` (list of all distinct sets, as names/indices, greedy
first), `attr(, "size_dist")` (named integer vector, set size -\>
count), `attr(, "kind")` the sense used.

## Details

`x` is EITHER a genotype matrix (rows = markers, cols = samples) OR a
precomputed square relatedness matrix. It is treated as precomputed when
it carries an `attr(, "kind")` (as returned by
[`pairwise_distance()`](https://sawers-rellan-labs.github.io/nilhmm/reference/pairwise_distance.md))
or is a square symmetric matrix; otherwise it is treated as genotypes
and `pairwise_distance(x, method)` is computed. To force one
interpretation, pass the matrix you want.

**Edge / threshold sense** (a set is independent iff no pair is an
edge):

- `r2` (similarity):

  edge if `r2 >= threshold`. `t ~ 0.2` = strict GWAS-QC pruning;
  **`t ~ 0.5` is recommended for JLM** (keeps resolution, roughly VIF 2
  / PLINK `--indep`).

- `mi` (similarity):

  edge if `MI >= threshold` (nats or bits). For Gaussian-equivalent
  reasoning in r2 units, `t_MI = -0.5 * ln(1 - r2)` nats, so
  `r2 = 0.5 <=> ~0.347 nats` – given as a note; the code does no
  conversion.

- `vi` (distance, a metric):

  edge if `VI <= threshold` (small VI = related).

The sense is taken from `attr(x, "kind")` when present, else from
`method` (`vi` is a distance; `r2`/`mi` are similarities). Pass `sense`
to override this inference explicitly – e.g. `sense = "distance"` on a
precomputed coordinate/cM matrix (see
[`position_distance()`](https://sawers-rellan-labs.github.io/nilhmm/reference/position_distance.md))
prunes it as a distance without abusing `method = "vi"`.

## RNG fidelity

The greedy set (`n_runs = 1`, or the first reported set) is
deterministic and bit-identical to the FastIndep CLI. The stochastic
runs (`n_runs > 1`) are algorithmically faithful and reproducible from
`seed` (self-contained PRNG), but individual random sets are not
bit-identical to the CLI's Mersenne-Twister stream.

## See also

[`pairwise_distance()`](https://sawers-rellan-labs.github.io/nilhmm/reference/pairwise_distance.md),
[`fast_indep_cpp()`](https://sawers-rellan-labs.github.io/nilhmm/reference/fast_indep_cpp.md)

## Examples

``` r
set.seed(1)
geno <- matrix(sample(0:2, 8 * 30, replace = TRUE), nrow = 8,
               dimnames = list(paste0("m", 1:8), NULL))
select_independent(geno, threshold = 0.5, method = "r2")
#> [1] "m1" "m2" "m3" "m4" "m5" "m6" "m7" "m8"
#> attr(,"sets")
#> attr(,"sets")[[1]]
#> [1] "m1" "m2" "m3" "m4" "m5" "m6" "m7" "m8"
#> 
#> attr(,"size_dist")
#> 8 
#> 1 
#> attr(,"kind")
#> [1] "similarity"
```
