# Find a large independent set in a thresholded similarity/distance matrix

Faithful port of FastIndep (Abraham 2013): the deterministic greedy
heuristic plus `n_runs - 1` stochastic runs. Two nodes are "related" (an
edge) when the matrix entry crosses `threshold` in the sense given by
`distance`; an independent set has no edges among its members. The
matrix diagonal is ignored.

## Usage

``` r
fast_indep_cpp(sim, threshold, n_runs, seed, distance = FALSE)
```

## Arguments

- sim:

  Symmetric numeric matrix. A similarity (higher = more related,
  `distance = false`) or a distance (lower = more related,
  `distance = true`).

- threshold:

  Edge cutoff. Similarity: edge if `sim[i,j] >= threshold`. Distance:
  edge if `sim[i,j] <= threshold`.

- n_runs:

  Total runs: 1 = greedy only; `> 1` adds `n_runs - 1` stochastic runs
  and returns the distinct sets found.

- seed:

  Seed for the stochastic runs (self-contained splitmix64 PRNG;
  independent of R's RNG). The greedy run is deterministic and
  seed-free.

- distance:

  If `true`, `sim` is a distance and the edge sense is inverted (edge if
  `<= threshold`). Default `false` (similarity).

## Value

A list: `best` (1-based indices of the largest independent set found,
the greedy set when it is largest), `sets` (list of distinct sets,
greedy first then by decreasing size), `size_dist` (named integer
vector, set size -\> number of distinct sets of that size).
