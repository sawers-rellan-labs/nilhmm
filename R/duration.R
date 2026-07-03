# Duration layer -- the transition machinery (S4). Orthogonal to emission.
# geometric = memoryless self-transition; rigidity = hard minimum run length
# via state expansion (this IS RTIGER's rigidity, reimplemented per S7);
# hsmm = explicit-duration sojourn (reserved, NOT for baking in the fitted
# Gamma -- that would be the S7 circularity trap).

#' Geometric (memoryless) duration
#'
#' Per-marker recombination rate `rrate` (self-stay = `1 - rrate`); the nilHMM
#' transition. `rrate` is a resolution hyperparameter, **not** an MLE --
#' calibrate by KS-vs-sim (memory `rigidity-not-mle`, [calibrate_r()]).
#'
#' @param rrate Per-marker recombination / switch rate between adjacent markers.
#' @return A duration spec for [fit()].
#' @examples
#' duration_geometric(rrate = 1e-4)
#' @export
duration_geometric <- function(rrate = 0.01) {
  structure(list(type = "geometric", r = rrate),
            class = c("nilHMM_duration_geometric", "nilHMM_duration"))
}

#' Rigidity duration (RTIGER reimplementation)
#'
#' Enforces a hard minimum run length of `rigidity` markers by expanding each
#' state into a chain of `rigidity` sub-states (phase-type / Erlang). Reimplements
#' the subset of RTIGER we use (BetaBinomial emission + rigidity + Viterbi +
#' KS-calibrated `r`); `optimize_R` is intentionally NOT ported (it
#' over-rigidifies, memory `optimize-R-overrigidifies-nil-sim`). Validate
#' against `tests/fixtures/baseline_pre_refactor/rtiger_rigidity_ref/`.
#'
#' @param rigidity Minimum run length in markers (integer >= 1).
#' @param xrate Per-marker switch / exit probability at the free (post-minimum)
#'   state -- the geometric tail beyond the enforced minimum run. `rigidity` of 1
#'   reduces the expansion to a plain geometric transition with this exit rate.
#' @return A duration spec for [fit()].
#' @examples
#' duration_rigidity(rigidity = 5, xrate = 2e-3)
#' @export
duration_rigidity <- function(rigidity = 5L, xrate = 0.01) {
  r <- as.integer(rigidity)
  if (r < 1L) stop("duration_rigidity(): rigidity must be >= 1")
  structure(list(type = "rigidity", r = r, p_switch = xrate),
            class = c("nilHMM_duration_rigidity", "nilHMM_duration"))
}

#' Explicit-duration (HSMM) sojourn -- reserved
#'
#' Explicit sojourn distribution. Reserved for a genuine duration model; **must
#' not** be used to bake in the design's fitted Gamma (S7 circularity trap).
#'
#' @param sojourn Sojourn distribution spec (placeholder).
#' @return A duration spec for [fit()].
#' @examples
#' duration_hsmm()   # reserved placeholder spec
#' @export
duration_hsmm <- function(sojourn = NULL) {
  structure(list(type = "hsmm", sojourn = sojourn),
            class = c("nilHMM_duration_hsmm", "nilHMM_duration"))
}