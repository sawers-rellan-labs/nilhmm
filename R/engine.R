# Engine: the duration-aware 3-state HMM (REF/HET/ALT). Layer 1 of the
# architecture (REFACTOR_R_PACKAGE.md S4). Emission and duration are pluggable
# interfaces (see emissions.R, duration.R). Hot loops live in src/ (Rcpp:
# count_emission_loglik_cpp, viterbi_log_cpp).

# --- internal: transition + init from r and state-frequency priors ----------
# Identical parameterization to the Python _build_transition: self-transition
# 1-r, off-diagonal recombination weighted by expected state frequencies. The
# generation enters ONLY through f_1/f_2 here (a prior), never a baked-in
# duration (S7).
.build_transition <- function(r, f_1, f_2) {
  f_0 <- 1 - f_1 - f_2
  if (f_0 <= 0) stop(".build_transition: f_1 + f_2 must be < 1")
  tmat <- matrix(c(
    1 - r,                       r * f_1 / (f_1 + f_2),       r * f_2 / (f_1 + f_2),
    r * f_0 / (f_0 + f_2),       1 - r,                       r * f_2 / (f_0 + f_2),
    r * f_0 / (f_0 + f_1),       r * f_1 / (f_0 + f_1),       1 - r
  ), nrow = 3, byrow = TRUE)
  list(log_start = log(c(f_0, f_1, f_2)), log_trans = log(tmat), n_sub = 1L)
}

# --- internal: rigidity (RTIGER) transition via state expansion (S4, S7) ------
# Each macro-state (REF/HET/ALT) becomes a chain of `r` sub-states: positions
# 1..r-1 advance deterministically (forced minimum run length r), position r is
# a free state that either self-loops (stay in the macro-state) or exits to
# another macro-state's sub-state 1, with probabilities from the macro
# transition `A` (built from `p_switch` + priors). r = 1 collapses to `A`
# itself. Expanded state index = m * r + k (m, k 0-based). All sub-states of a
# macro-state share its emission (handled in decode()).
.build_transition_rigidity <- function(r, p_switch, f_1, f_2) {
  macro <- .build_transition(p_switch, f_1, f_2)
  A     <- exp(macro$log_trans)              # 3 x 3 macro transition
  start <- exp(macro$log_start)              # c(f_0, f_1, f_2)
  n <- as.integer(r); K <- 3L; S <- K * n
  idx <- function(m, k) m * n + k + 1L       # m,k 0-based -> 1-based
  Tm <- matrix(0, S, S)
  for (m in 0:2) for (k in 0:(n - 1L)) {
    if (k < n - 1L) {
      Tm[idx(m, k), idx(m, k + 1L)] <- 1                       # forced advance
    } else {
      Tm[idx(m, k), idx(m, n - 1L)] <- A[m + 1L, m + 1L]       # free state: stay
      for (mp in 0:2) if (mp != m) Tm[idx(m, k), idx(mp, 0L)] <- A[m + 1L, mp + 1L]  # exit
    }
  }
  s0 <- numeric(S); for (m in 0:2) s0[idx(m, 0L)] <- start[m + 1L]
  list(log_start = log(s0), log_trans = log(Tm), n_sub = n)
}

# Dispatch: (duration, priors) -> list(log_start, log_trans, n_sub).
.duration_transition <- function(duration, priors) {
  if (inherits(duration, "nilHMM_duration_geometric"))
    return(.build_transition(duration$r, priors$f_1, priors$f_2))
  if (inherits(duration, "nilHMM_duration_rigidity"))
    return(.build_transition_rigidity(duration$r, duration$p_switch, priors$f_1, priors$f_2))
  stop(".duration_transition(): unsupported duration '", duration$type,
       "' (hsmm is reserved, S4)")
}

# --- internal: run-length-encode a state path into common-schema segments ----
# Boundaries where the state changes; start_bp/end_bp = first/last marker
# position of each contiguous run (same convention as the Python rle_segments
# and the RTIGER segments).
.rle_segments <- function(path, pos, chr, name, source, donor) {
  if (length(path) == 0) return(NULL)
  r <- rle_segments_cpp(as.integer(path), as.integer(pos))   # boundaries in C++
  data.frame(
    source = source, donor = donor, name = name, chr = as.integer(chr),
    start_bp = r$start_bp, end_bp = r$end_bp, state = r$state,
    stringsAsFactors = FALSE
  )
}

# Batched fixed-means count calls: per chromosome, pivot the (rectangular)
# cohort to T x S count matrices, memoize the emission over distinct (n,a)
# pairs, and run one C++ viterbi_batch over all samples. Returns NULL if the
# data is not a clean rectangular grid (every sample x every position once),
# so call_ancestry can fall back to the safe per-sample loop.
.batched_calls <- function(data, emission, td, theta, source, donor, has_donor) {
  by_chr <- split(seq_len(nrow(data)), data$chr)   # row indices per chr, O(n)
  out <- vector("list", length(by_chr))
  for (i in seq_along(by_chr)) {
    ri    <- by_chr[[i]]
    nm_v  <- data$name[ri]; pos_v <- data$pos[ri]
    samples <- unique(nm_v); pos_lv <- sort(unique(pos_v))
    Tn <- length(pos_lv); S <- length(samples)
    si <- match(nm_v, samples); pii <- match(pos_v, pos_lv)
    lin <- (si - 1L) * Tn + pii                    # column-major cell id
    if (length(ri) != Tn * S || anyDuplicated(lin)) return(NULL)   # not a rectangular grid
    n_mat <- integer(Tn * S); n_mat[lin] <- as.integer(data$n_ref[ri] + data$n_alt[ri]); dim(n_mat) <- c(Tn, S)
    a_mat <- integer(Tn * S); a_mat[lin] <- as.integer(data$n_alt[ri]);                   dim(a_mat) <- c(Tn, S)
    base  <- max(n_mat) + 1
    key   <- as.double(n_mat) * base + as.double(a_mat)
    u     <- unique(key)
    em_u  <- count_emission_loglik_cpp(as.integer(u %/% base), as.integer(u %% base),
                                       theta, emission$conc)
    inv   <- matrix(match(key, u) - 1L, Tn, S)
    paths <- viterbi_batch_cpp(td$log_start, td$log_trans, em_u, inv)   # Tn x S
    seg   <- rle_segments_batch_cpp(paths, pos_lv)                      # all samples at once
    dn_by <- if (has_donor) data$donor[ri][match(samples, nm_v)] else rep(donor, S)
    out[[i]] <- data.frame(
      source = source, donor = dn_by[seg$sample], name = samples[seg$sample],
      chr = as.integer(names(by_chr)[i]),
      start_bp = seg$start_bp, end_bp = seg$end_bp, state = seg$state,
      stringsAsFactors = FALSE)
  }
  calls <- do.call(rbind, out)
  calls[order(calls$donor, calls$name, calls$chr, calls$start_bp), , drop = FALSE]
}

#' Fit HMM emission/transition parameters
#'
#' For a fixed-mean emission (the default, reproducing the Python count caller)
#' this resolves `theta` and the transition and is otherwise a no-op. For a
#' fittable-mean emission (`fit_means = TRUE`; required for RNA / reference-biased
#' data per S10) it Baum-Welch EM-fits the emission means. Viterbi/EM hot loops
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
  td <- .duration_transition(duration, priors)
  theta <- .emission_theta(emission)                 # initial / fixed means
  if (isTRUE(emission$fit_means)) {
    theta <- .em_fit_means(list(obs), emission, td, theta, control)  # single sequence
  }
  structure(list(theta = theta, log_start = td$log_start, log_trans = td$log_trans,
                 n_sub = td$n_sub, emission = emission, duration = duration),
            class = "nilHMM_model")
}

# Replicate a T x 3 macro emission across each macro-state's n_sub sub-states
# (expanded column index = m * n_sub + k). n_sub = 1 is a no-op.
.expand_emission <- function(em3, n_sub) {
  if (n_sub == 1L) return(em3)
  out <- matrix(0, nrow(em3), 3L * n_sub)
  for (m in 0:2) for (k in 0:(n_sub - 1L)) out[, m * n_sub + k + 1L] <- em3[, m + 1L]
  out
}

#' Decode the most-likely state path (Viterbi)
#'
#' @param model A fitted model from [fit()].
#' @param obs A per-(sample, chromosome) observation table.
#' @return Integer macro-state path over markers (0 = REF, 1 = HET, 2 = ALT);
#'   for rigidity the expanded sub-states are mapped back to macro-states.
#' @export
decode <- function(model, obs) {
  em <- .expand_emission(.emission_loglik(model$emission, obs, model$theta), model$n_sub)
  path <- viterbi_log_cpp(model$log_start, model$log_trans, em)
  if (model$n_sub > 1L) path %/% model$n_sub else path
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
#' @param r,err,conc,fit_means,p_switch Caller parameters forwarded to
#'   [caller_spec()]. Explicit formals (not `...`) so that, e.g., `r` is never
#'   partial-matched to another argument.
#' @param f_1,f_2 Single-locus priors, used when `design` is `NULL`.
#' @param source Value for the output `source` column.
#' @param donor Donor/taxon label when `data` has no `donor` column.
#' @return data.frame in the common schema
#'   (`source, donor, name, chr, start_bp, end_bp, state`).
#' @export
call_ancestry <- function(data, caller = c("nnil", "rtiger", "skimbin"),
                          design = NULL, r = 0.01, err = 0.01, conc = 20,
                          fit_means = FALSE, p_switch = 0.01,
                          f_1 = NULL, f_2 = NULL,
                          source = "nilHMM", donor = NA_character_) {
  caller <- match.arg(caller)
  spec <- caller_spec(caller, r = r, err = err, conc = conc,
                      fit_means = fit_means, p_switch = p_switch)

  priors <- if (!is.null(design)) design_priors(design)
            else if (!is.null(f_1) && !is.null(f_2)) list(f_1 = f_1, f_2 = f_2)
            else stop("call_ancestry(): supply `design` or both `f_1` and `f_2`")

  req <- c("name", "chr", "pos", "n_ref", "n_alt")
  if (!all(req %in% names(data))) stop("call_ancestry(): data needs columns ", paste(req, collapse = ", "))
  has_donor <- "donor" %in% names(data)

  td     <- .duration_transition(spec$duration, priors)   # same for all samples
  theta0 <- .emission_theta(spec$emission)

  # Fast path: fixed-means count caller on a rectangular cohort -> one batched
  # Viterbi per chromosome across all samples in C++ (mirrors the numpy-batched
  # Python). Falls back to the per-sample loop when fit_means / rigidity / a
  # non-rectangular grid make batching unsafe.
  if (!isTRUE(spec$emission$fit_means) && td$n_sub == 1L &&
      inherits(spec$emission, "nilHMM_emission_count")) {
    batched <- .batched_calls(data, spec$emission, td, theta0, source, donor, has_donor)
    if (!is.null(batched)) return(batched)
  }

  # Split by sample ONCE (O(n)) rather than re-scanning the frame per sample
  # (which is O(samples * rows) and dominates wall-clock on large cohorts).
  by_name <- split(data, data$name, drop = TRUE)
  out <- list()
  for (nm in names(by_name)) {
    dn <- by_name[[nm]]
    donor_nm <- if (has_donor) dn$donor[1] else donor
    # per-chromosome observation sequences for this sample (HMM decodes each chr
    # independently; emission params are shared and fit pooled across them).
    obs_list <- lapply(split(dn, dn$chr, drop = TRUE), function(dc) {
      dc <- dc[order(dc$pos), , drop = FALSE]
      list(chr = dc$chr[1], pos = dc$pos, n = dc$n_ref + dc$n_alt, a = dc$n_alt)
    })
    theta <- if (isTRUE(spec$emission$fit_means))
               .em_fit_means(obs_list, spec$emission, td, theta0) else theta0
    model <- structure(list(theta = theta, log_start = td$log_start,
                            log_trans = td$log_trans, n_sub = td$n_sub,
                            emission = spec$emission, duration = spec$duration),
                       class = "nilHMM_model")
    for (o in obs_list)
      out[[length(out) + 1L]] <- .rle_segments(decode(model, o), o$pos, o$chr, nm, source, donor_nm)
  }
  calls <- do.call(rbind, out)
  calls[order(calls$donor, calls$name, calls$chr, calls$start_bp), , drop = FALSE]
}
