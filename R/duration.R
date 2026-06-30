# Duration layer — the transition machinery (§4). Orthogonal to emission.
# geometric = memoryless self-transition; rigidity = hard minimum run length
# via state expansion (this IS RTIGER's rigidity, reimplemented per §7);
# hsmm = explicit-duration sojourn (reserved, NOT for baking in the fitted
# Gamma — that would be the §7 circularity trap).

#' Geometric (memoryless) duration
#'
#' Self-transition rate `r`; the nilHMM transition. `r` is a resolution
#' hyperparameter, **not** an MLE — calibrate by KS-vs-sim (memory
#' `rigidity-not-mle`, [calibrate_r()]).
#'
#' @param r Self-transition / recombination rate between adjacent markers.
#' @return A duration spec for [fit()].
#' @export
duration_geometric <- function(r = 0.01) {
  structure(list(type = "geometric", r = r),
            class = c("nilHMM_duration_geometric", "nilHMM_duration"))
}

#' Rigidity duration (RTIGER reimplementation)
#'
#' Enforces a hard minimum run length of `r` markers by expanding each state
#' into a chain of `r` sub-states (phase-type / Erlang machinery). Reimplements
#' the subset of RTIGER we use (BetaBinomial emission + rigidity + Viterbi +
#' KS-calibrated `r`); `optimize_R` is intentionally NOT ported (it
#' over-rigidifies, memory `optimize-R-overrigidifies-nil-sim`). Validate
#' against `tests/fixtures/baseline_pre_refactor/rtiger_rigidity_ref/`.
#'
#' @param r Minimum run length in markers (the rigidity).
#' @return A duration spec for [fit()].
#' @export
duration_rigidity <- function(r = 5L) {
  structure(list(type = "rigidity", r = as.integer(r)),
            class = c("nilHMM_duration_rigidity", "nilHMM_duration"))
}

#' Explicit-duration (HSMM) sojourn — reserved
#'
#' Explicit sojourn distribution. Reserved for a genuine duration model; **must
#' not** be used to bake in the design's fitted Gamma (§7 circularity trap).
#'
#' @param sojourn Sojourn distribution spec (placeholder).
#' @return A duration spec for [fit()].
#' @export
duration_hsmm <- function(sojourn = NULL) {
  structure(list(type = "hsmm", sojourn = sojourn),
            class = c("nilHMM_duration_hsmm", "nilHMM_duration"))
}