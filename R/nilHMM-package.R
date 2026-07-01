#' nilHMM: duration-aware HMM ancestry caller
#'
#' A unified R + Rcpp engine for calling REF/HET/ALT ancestry in NILs from
#' sequencing data. One 3-state HMM with swappable emission ([emission_count()],
#' [emission_gt()]) and duration ([duration_geometric()], [duration_rigidity()],
#' [duration_hsmm()]) layers expresses the `nnil` and `rtiger` callers (plus the
#' `binhmm` bin/cluster caller). See `REFACTOR_R_PACKAGE.md` for the design.
#'
#' The package is **data-agnostic**: no hardcoded paths, sample lists, or mount
#' locations. Functions take `(data, params)` and return calls; pipeline scripts
#' (not part of the package) own which files/samples and where outputs go.
#'
#' @keywords internal
#' @useDynLib nilHMM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
"_PACKAGE"
# The RcppParallel import forces its namespace (and the bundled TBB runtime) to
# load before nilHMM's DLL, so viterbi_batch_par_cpp resolves @rpath/libtbb at
# load time (otherwise dlopen fails on a clean install / R CMD check).