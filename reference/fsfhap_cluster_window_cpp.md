# FSFHap stage 2a: cluster a window of parent-called haplotypes

Faithful port of `HaplotypeClusterer.makeClusters` + `HaplotypeCluster`
consensus, on the parent-origin frame (`g` in 0 A-hom, 1 het, 2 C-hom, 3
missing). Clusters group taxa whose window haplotypes are 0-distance
(identical modulo missing); a haplotype 0-distance to members of several
clusters joins all of them with fractional score `1/count`. Clusters are
returned sorted by score (desc) then size (desc) —
`HaplotypeCluster.compareTo`.

## Usage

``` r
fsfhap_cluster_window_cpp(
  Gw,
  maxdiff = 0L,
  merge = FALSE,
  move_biggest = FALSE,
  max_het = -1L
)
```

## Arguments

- Gw:

  Integer matrix, taxa x window-sites, canonical `g` in 0,1,2,3.

- maxdiff:

  Distance threshold for `merge`/`move_biggest` (TASSEL
  `maxDifferenceScore`, 0 on the BC/finder path).

- merge:

  Apply `mergeClusters(maxdiff)` after `makeClusters` (clusterWindow
  does this only when `maxdiff > 0`).

- move_biggest:

  Apply `moveAllHaplotypesToBiggestCluster(maxdiff)`.

- max_het:

  If `>= 0`, drop clusters with more than `max_het` heterozygous sites
  (`removeHeterozygousClusters`; the finder passes `maxdiff + 5`).

## Value

List: `size`, `score` (per cluster); `majority`, `unanimous` (clusters x
sites consensus, `3` = N); `members` (list of 1-based taxon indices per
cluster). For 0-distance clusters `majority == unanimous`; they diverge
only after merges.
