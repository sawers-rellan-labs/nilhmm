# Engine: the duration-aware 3-state HMM (REF/HET/ALT). Layer 1 of the
# architecture (REFACTOR_R_PACKAGE.md §4). Emission and duration are pluggable
# interfaces (see emissions.R, duration.R). Hot loops live in src/ (Rcpp:
# count_emission_loglik_cpp, viterbi_log_cpp).

# --- internal: transition + init from r and state-frequency priors ----------
# Identical parameterization to the Python _build_transition: self-transition
# 1-r, off-diagonal recombination weighted by expected state frequencies. The
# generation enters ONLY through f_1/f_2 here (a prior), never a baked-in
# duration (§7).
.build_transition <- function(r, f_1, f_2) {
  f_0 <- 1 - f_1 - f_2
  if (f_0 <= 0) stop(".build_transition: f_1 + f_2 must be < 1")
  tmat <- matrix(c(
    1 - r,                       r * f_1 / (f_1 + f_2),       r * f_2 / (f_1 + f_2),
    r * f_0 / (f_0 + f_2),       1 - r,                       r * f_2 / (f_0 + f_2),
    r * f_0 / (f_0 + f_1),       r * f_1 / (f_0 + f_1),       1 - r
  ), nrow = 3, byrow = TRUE)
  list(log_start = log(c(f_0, f_1, f_2)), log_trans = log(tmat))
}

# --- internal: run-length-encode a state path into common-schema segments ----
# Boundaries where the state changes; start_bp/end_bp = first/last marker
# position of each contiguous run (same convention as the Python rle_segments
# and the RTIGER segments).
.rle_segments <- function(path, pos, chr, name, source, donor) {
  if (length(path) == 0) return(NULL)
  brk    <- which(diff(path) != 0)
  starts <- c(1L, brk + 1L)
  ends   <- c(brk, length(path))
  data.frame(
    source = source, donor = donor, name = name, chr = as.integer(chr),
    start_bp = as.integer(pos[starts]), end_bp = as.integer(pos[ends]),
    state = as.integer(path[starts]), stringsAsFactors = FALSE
  )
}

#' Fit HMM emission/transition parameters
#'
#' For a fixed-mean emission (the default, reproducing the Python count caller)
#' this resolves `theta` and the transition and is otherwise a no-op. For a
#' fittable-mean emission (`fit_means = TRUE`; required for RNA / reference-biased
#' data per §10) it Baum-Welch EM-fits the emission means. Viterbi/EM hot loops
#' are in Rcpp.
#'
#' @param obs A per-(sample, chromosome) observation table from an emission's
#'   reader (for `count`: columns `n`, `a`).
#' @param emission Emission spec ([emission_count()] / [emission_gt()] /
#'   [emission_dosage()]).
#' @param duration Duration spec ([duration_geometric()] /
#'   [duration_rigidity()] / [duration_hsmm()]).
#' @param priors Single-locus genotype-frequency priors `list(f_1, f_2)`.
#' @param control EM control (`max_iter`, `tol`).
#' @return A fitted model: `list(theta, log_start, log_trans, emission, duration)`.
#' @export
fit <- function(obs, emission, duration, priors, control = list()) {
  if (!inherits(duration, "nilHMM_duration_geometric"))
    stop("fit(): only geometric duration is implemented (Task 4); rigidity is next")
  tr <- .build_transition(duration$r, priors$f_1, priors$f_2)
  theta <- .emission_theta(emission)                 # initial / fixed means
  if (isTRUE(emission$fit_means)) {
    theta <- .em_fit_means(obs, emission, tr, theta, control)
  }
  structure(list(theta = theta, log_start = tr$log_start, log_trans = tr$log_trans,
                 emission = emission, duration = duration),
            class = "nilHMM_model")
}

#' Decode the most-likely state path (Viterbi)
#'
#' @param model A fitted model from [fit()].
#' @param obs A per-(sample, chromosome) observation table.
#' @return Integer state path over markers (0 = REF, 1 = HET, 2 = ALT).
#' @export
decode <- function(model, obs) {
  em <- .emission_loglik(model$emission, obs, model$theta)
  viterbi_log_cpp(model$log_start, model$log_trans, em)
}

#' Top-level ancestry-calling API
#'
#' Resolves a named caller (and optional design preset) into emission + duration
#' specs, runs [fit()] / [decode()] per (sample, chromosome), and returns
#' common-schema segment calls. Data-agnostic: pass the observations in; this
#' function never reads paths.
#'
#' @param data Long observation table with columns `name, chr, pos, n_ref,
#'   n_alt` (and optionally `donor`). From [read_counts()].
#' @param caller One of `"nnil"`, `"rtiger"`, `"skimbin"`.
#' @param design Breeding-design key for priors (e.g. `"BC2S2"`, `"BC2S3"`).
#'   Required unless `f_1`/`f_2` are supplied.
#' @param r,err,conc,fit_means Caller parameters forwarded to [caller_spec()].
#'   Explicit formals (not `...`) so that, e.g., `r` is never partial-matched to
#'   another argument.
#' @param f_1,f_2 Single-locus priors, used when `design` is `NULL`.
#' @param source Value for the output `source` column.
#' @param donor Donor/taxon label when `data` has no `donor` column.
#' @return data.frame in the common schema
#'   (`source, donor, name, chr, start_bp, end_bp, state`).
#' @export
call_ancestry <- function(data, caller = c("nnil", "rtiger", "skimbin"),
                          design = NULL, r = 0.01, err = 0.01, conc = 20,
                          fit_means = FALSE, f_1 = NULL, f_2 = NULL,
                          source = "nilHMM", donor = NA_character_) {
  caller <- match.arg(caller)
  spec <- caller_spec(caller, r = r, err = err, conc = conc, fit_means = fit_means)

  priors <- if (!is.null(design)) design_priors(design)
            else if (!is.null(f_1) && !is.null(f_2)) list(f_1 = f_1, f_2 = f_2)
            else stop("call_ancestry(): supply `design` or both `f_1` and `f_2`")

  req <- c("name", "chr", "pos", "n_ref", "n_alt")
  if (!all(req %in% names(data))) stop("call_ancestry(): data needs columns ", paste(req, collapse = ", "))
  has_donor <- "donor" %in% names(data)

  out <- list()
  for (nm in unique(data$name)) {
    dn <- data[data$name == nm, , drop = FALSE]
    donor_nm <- if (has_donor) dn$donor[1] else donor
    for (cc in unique(dn$chr)) {
      dc <- dn[dn$chr == cc, , drop = FALSE]
      dc <- dc[order(dc$pos), , drop = FALSE]
      obs <- list(n = dc$n_ref + dc$n_alt, a = dc$n_alt)
      model <- fit(obs, spec$emission, spec$duration, priors)
      path <- decode(model, obs)
      out[[length(out) + 1L]] <- .rle_segments(path, dc$pos, cc, nm, source, donor_nm)
    }
  }
  calls <- do.call(rbind, out)
  calls[order(calls$donor, calls$name, calls$chr, calls$start_bp), , drop = FALSE]
}
