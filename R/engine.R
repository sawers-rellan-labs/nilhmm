# Engine: the duration-aware 3-state HMM (REF/HET/ALT). Layer 1 of the
# architecture (REFACTOR_R_PACKAGE.md §4). Emission and duration are pluggable
# interfaces (see emissions.R, duration.R). Hot loops live in src/ (Rcpp).

#' Fit HMM emission/transition parameters by Baum-Welch EM
#'
#' Fits the fittable parameters of an emission model (e.g. EM-fit emission
#' means for count/ref-biased data, see §10) and transition rate where
#' applicable. Hot loop in Rcpp.
#'
#' @param data Per-chromosome observation object (see emission constructors).
#' @param emission An emission spec from [emission_count()] / [emission_gt()] /
#'   [emission_dosage()].
#' @param duration A duration spec from [duration_geometric()] /
#'   [duration_rigidity()] / [duration_hsmm()].
#' @param priors Single-locus genotype-frequency priors `list(f_1, f_2)` (from
#'   [design_priors()]); never a generation prior baked into the transition.
#' @param control EM control (max_iter, tol, init).
#' @return A fitted model object consumed by [decode()].
#' @export
fit <- function(data, emission, duration, priors, control = list()) {
  stop("nilHMM::fit() not yet implemented (Task 4)")
}

#' Decode the most-likely state path (Viterbi)
#'
#' @param model A fitted model from [fit()].
#' @param data Per-chromosome observation object.
#' @return Integer state path over markers (0 = REF, 1 = HET, 2 = ALT).
#' @export
decode <- function(model, data) {
  stop("nilHMM::decode() not yet implemented (Task 4)")
}

#' Top-level ancestry-calling API
#'
#' Convenience wrapper that resolves a named caller (and optional
#' source/design presets) into emission + duration specs, runs [fit()] /
#' [decode()] per chromosome, and returns common-schema segment calls. This is
#' the function pipeline scripts call.
#'
#' @param data Marker-level input (counts / GT / dosage) plus positions; the
#'   caller is data-agnostic about provenance.
#' @param caller One of `"nnil"`, `"rtiger"`, `"skimbin"` (see [caller_spec()]).
#' @param design Breeding-design key (e.g. `"BC2S2"`, `"BC2S3"`) selecting
#'   priors and the expected-fragment Null (see [design_priors()]).
#' @param regime Optional emission-by-depth regime override (see
#'   [select_emission()]); inferred from depth when `NULL`.
#' @param ... Caller-specific parameter overrides (`r`, `err`, `conc`, ...).
#' @return A data.frame in the common segment schema
#'   (`source, donor, name, chr, start_bp, end_bp, state`).
#' @export
call_ancestry <- function(data, caller = c("nnil", "rtiger", "skimbin"),
                          design = NULL, regime = NULL, ...) {
  caller <- match.arg(caller)
  stop("nilHMM::call_ancestry() not yet implemented (Task 4)")
}