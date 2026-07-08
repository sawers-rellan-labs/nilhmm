# LB-Impute distance-aware full-chromosome Viterbi

Decodes the most-likely REF/HET/ALT path under LB-Impute's
distance-dependent transition (`FindPath2`). Between markers a physical
distance `d` bp apart, the stay probability is \\p_s = 0.5(1 +
e^{-d/recombdist})\\ and the single recombination probability is \\p_r =
0.5(1 - e^{-d/recombdist})\\. Homozygous \<-\> heterozygous transitions
cost one recombination (\\p_r\\); homozygous -\> the OTHER homozygous
state costs two (\\p_r^2\\) unless `drp = TRUE` (LB-Impute `-dr`), which
prices it as a single event (for inbred / RIL populations). Transition
weights are LB-Impute's exact (un-normalized) model; the Viterbi argmax
over full-chromosome paths supersedes the original's windowed
best-path + forward/reverse consensus.

## Usage

``` r
lb_viterbi_cpp(log_init, log_emit, tpos, recombdist, drp)
```

## Arguments

- log_init:

  Length-3 vector of log initial-state probabilities (REF/HET/ALT).

- log_emit:

  T x 3 matrix of log emissions (from
  [`lb_emission_loglik_cpp()`](https://sawers-rellan-labs.github.io/nilhmm/reference/lb_emission_loglik_cpp.md)).

- tpos:

  Numeric length-T transition coordinate per marker, non-decreasing. The
  genetic (cM) or physical (bp) coordinate the transition decays over;
  the arithmetic is unit-agnostic, so `tpos` and `recombdist` must share
  units.

- recombdist:

  Coordinate distance (same units as `tpos`) over which the
  recombination probability equalizes (LB-Impute `recombdist`; 1e7 bp,
  or ~50 cM).

- drp:

  If `TRUE`, a homozygous-\>homozygous switch is priced as a single
  recombination rather than a double event.

## Value

Integer length-T most-likely state path (0 = REF, 1 = HET, 2 = ALT).
