# Calibration (§7, §9). r is set by minimum-distance KS against simulated truth,
# NOT by maximizing a likelihood (memory `rigidity-not-mle`). The generative
# parameter is the Gamma's lambda; r is the caller knob tuned to it.
#
# BRB caveat (BRB_run_findings.md): when the caller under-calls donor, the donor
# RATE is r-invariant and KS-on-block-size rails to the grid edge. A
# rate/footprint-aware objective is an open item (§10).

#' Calibrate the duration hyperparameter `r` by KS-vs-sim
#'
#' Sweeps `r`, calls on data, and returns the `r` minimising the KS distance
#' between called donor-block sizes and the simulated-truth distribution.
#'
#' @param data Marker-level input for a taxon/group.
#' @param sim_segments Simulated-truth segments for the matching design (the KS
#'   target; see [expected_fragment_dist()] and the frozen `sim_truth/`
#'   fixtures).
#' @param r_grid Candidate `r` values to sweep.
#' @param ... Forwarded to [call_ancestry()].
#' @return data.frame of `(r, D, ...)` plus the argmin; warns if `at_grid_edge`.
#' @export
calibrate_r <- function(data, sim_segments, r_grid, ...) {
  stop("nilHMM::calibrate_r() not yet implemented (Task 4)")
}