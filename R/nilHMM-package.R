#' nilHMM: duration-aware HMM ancestry caller
#'
#' A unified R + Rcpp engine for calling REF/HET/ALT ancestry in NILs from
#' sequencing data. One 3-state HMM with swappable emission ([emission_count()],
#' [emission_gt()], [emission_dosage()]) and duration ([duration_geometric()],
#' [duration_rigidity()], [duration_hsmm()]) layers expresses the `nnil`,
#' `rtiger` and `skimbin` callers. See `REFACTOR_R_PACKAGE.md` for the design.
#'
#' The package is **data-agnostic**: no hardcoded paths, sample lists, or mount
#' locations. Functions take `(data, params)` and return calls; pipeline scripts
#' (not part of the package) own which files/samples and where outputs go.
#'
#' @keywords internal
#' @useDynLib nilHMM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"