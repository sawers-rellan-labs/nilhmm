# Presets layer 2b: breeding-design (generation) axis (§6). Orthogonal to the
# emission axis: a dataset carries BOTH a depth-regime (-> emission) and a
# breeding design (-> generation-derived quantities). e.g. skim = count x BC2S2,
# BRB = count x BC2S3.
#
# Generation determines the single-locus genotype-freq priors (f_1, f_2) and the
# expected fragment-size law (the fitted Gamma). The Gamma is CALIBRATION/NULL
# ONLY -- NEVER an engine prior (§7 circularity trap): baking it into the
# transition would make the caller reproduce it by construction and mask the
# open "callers run longer than the model" signal.

#' Single-locus genotype-frequency priors for a breeding design
#'
#' @param design Design key (e.g. `"BC1S1"`, `"BC2S2"`, `"BC2S3"`). Looked up in
#'   the bundled `breeding_designs` dataset.
#' @return `list(g, f_1, f_2)`.
#' @export
design_priors <- function(design) {
  stop("nilHMM::design_priors() not yet implemented (Task 4)")
}

#' Expected fragment-size distribution (the Null / KS target)
#'
#' Returns the fitted Gamma `(k, lambda)` for a design plus density/CDF
#' closures, for plotting the grey Null and as the KS-vs-sim calibration target.
#' Never feeds the engine transition.
#'
#' @param design Design key.
#' @return `list(k, lambda, density, cdf)`.
#' @export
expected_fragment_dist <- function(design) {
  stop("nilHMM::expected_fragment_dist() not yet implemented (Task 4)")
}

#' Fit a design Gamma from simulated segments
#'
#' Extensible entry point used by `data-raw/make_breeding_designs.R` to derive
#' the bundled `(k, lambda)` from simcross output.
#'
#' @param sim_segments data.frame of simulated segment sizes (cM).
#' @return `list(k, lambda)`.
#' @export
fit_design_gamma <- function(sim_segments) {
  stop("nilHMM::fit_design_gamma() not yet implemented (Task 4)")
}

#' Convert segment sizes from cM to Mb using a map
#'
#' cM-space Gammas are assembly-robust; anything bp/Mb is tied to B73 v5.
#'
#' @param seg Segments with cM coordinates.
#' @param map A consensus map (see [load_map()]).
#' @return `seg` with Mb coordinates added.
#' @export
cm_to_mb <- function(seg, map) {
  stop("nilHMM::cm_to_mb() not yet implemented (Task 4)")
}