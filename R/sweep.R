# caller_sweep: sweep a caller's SEGMENTATION parameter with a single shared fit.
#
# The segmentation parameter (rtiger `rigidity`, nnil `rrate`) is a prior on run
# length / switch rate, not an emission property, so the fitted emission is
# ~parameter-independent. That lets the expensive estimation run ONCE and the grid
# be swept by cheap decodes -- exact at the reference value, a close approximation
# elsewhere -- fanned across cores. Turns an N-fit sweep into 1 fit + N decodes.

#' Sweep a caller's segmentation parameter with one shared fit
#'
#' Fit once (rtiger: the joint EM; nnil: the per-sample emission means, or nothing
#' when `fit_means = FALSE`), then decode across `values`, fanning the independent
#' decodes over `threads`. Results are exact at `ref` and a close approximation for
#' other values (the emission is nearly parameter-invariant), at a fraction of a
#' cold per-value refit.
#'
#' @param data Common input: `name, chr, pos, n_ref, n_alt` (+ optional `donor`).
#' @param caller `"rtiger"` (sweeps `rigidity`) or `"nnil"` (sweeps `rrate`).
#' @param values Parameter grid to sweep.
#' @param design,f_1,f_2 Population priors (a design name, or explicit `f_1`,`f_2`).
#' @param threads Fan-out width (`parallel::mclapply` on unix; serial otherwise).
#' @param ref Reference value for the single shared fit (default `median(values)`,
#'   rounded for rtiger). The emission is fit here and reused across the grid.
#' @param min_cov Covered-marker filter before decoding (default `1L`, as in
#'   [call_states()]); `0L` keeps every marker.
#' @param err,conc,fit_means nnil count-emission parameters (fixed across the grid).
#' @param seed,postprocess rtiger fit seed and border post-processing.
#' @param source,donor Output labels.
#' @return A common-schema segment table (`source, donor, name, chr, start_bp,
#'   end_bp, state`) with an added column named for the swept parameter
#'   (`rigidity` or `rrate`) tagging each value's calls.
#' @export
caller_sweep <- function(data, caller = c("rtiger", "nnil"), values,
                         design = NULL, f_1 = NULL, f_2 = NULL, threads = 1L,
                         ref = NULL, min_cov = 1L, err = 0.01, conc = 20,
                         fit_means = FALSE, seed = 1L, postprocess = TRUE,
                         source = "nilHMM", donor = NA_character_) {
  caller <- match.arg(caller)
  if (!all(c("name", "chr", "pos") %in% names(data)))
    stop("caller_sweep(): data needs columns name, chr, pos")
  if (!length(values)) stop("caller_sweep(): `values` is empty")
  if (!all(c("n_ref", "n_alt") %in% names(data)))
    stop("caller_sweep(): needs n_ref/n_alt read counts")
  has_donor <- "donor" %in% names(data)
  fan <- function(X, FUN) {
    r <- if (threads > 1L && .Platform$OS.type == "unix")
      parallel::mclapply(X, FUN, mc.cores = threads) else lapply(X, FUN)
    if (any(vapply(r, function(x) inherits(x, "try-error") || is.null(x), logical(1))))
      stop("caller_sweep(): a decode/fit job failed")
    r
  }
  if (!is.null(min_cov) && min_cov > 0L)
    data <- data[data$n_ref + data$n_alt >= min_cov, , drop = FALSE]
  if (!nrow(data)) stop("caller_sweep(): no markers with coverage >= min_cov")

  # ---------------- rtiger: one joint EM, decode per rigidity -----------------
  if (caller == "rtiger") {
    values <- as.integer(values)
    o <- .rtiger_obs(data, has_donor, donor)
    .rtiger_check_coverage(o$obs, max(values))            # feasibility for the whole grid
    r_ref <- if (is.null(ref)) as.integer(round(stats::median(values))) else as.integer(ref)
    fit <- .rtiger_fit(o$obs, r_ref, threads = threads, seed = seed)   # ONE fit
    segs <- fan(values, function(v) {
      paths <- .rtiger_decode(o$obs, fit, v, postprocess = postprocess, threads = 1L)
      seg <- to_segments(.rtiger_assemble(o$obs, o$pos, paths, o$donor_of, source))
      seg$rigidity <- v
      seg
    })
    return(do.call(rbind, segs))
  }

  # ---------------- nnil: per-sample emission once, decode per rrate ----------
  priors <- if (!is.null(design)) design_priors(design)
            else if (!is.null(f_1) && !is.null(f_2)) list(f_1 = f_1, f_2 = f_2)
            else stop("caller_sweep(nnil): supply `design` or both `f_1`, `f_2`")
  emission <- emission_count(err, conc, fit_means)
  theta0 <- .emission_theta(emission)
  ref_rrate <- if (is.null(ref)) stats::median(values) else ref
  by_name <- split(data, data$name, drop = TRUE)

  # Phase 1 (parallel): per-sample chains + emission theta -- fixed (theta0) or
  # fit ONCE at ref_rrate. Its transition only matters for the fit_means EM.
  td_ref <- .duration_transition(
    caller_spec("nnil", rrate = ref_rrate, err = err, conc = conc, fit_means = fit_means)$duration, priors)
  samples <- fan(names(by_name), function(nm) {
    dn <- by_name[[nm]]
    obs_list <- lapply(split(dn, dn$chr, drop = TRUE), function(dc) {
      dc <- dc[order(dc$pos), , drop = FALSE]
      list(chr = dc$chr[1], pos = dc$pos, n = dc$n_ref + dc$n_alt, a = dc$n_alt)
    })
    theta <- if (isTRUE(fit_means)) .em_fit_means(obs_list, emission, td_ref, theta0) else theta0
    list(name = nm, donor = if (has_donor) dn$donor[1] else donor,
         obs_list = obs_list, theta = theta)
  })

  # Phase 2 (parallel over value x sample): build the transition per rrate, decode
  # every sample with its already-fit emission. Precompute the transitions serially.
  tds <- lapply(values, function(v) .duration_transition(
    caller_spec("nnil", rrate = v, err = err, conc = conc, fit_means = fit_means)$duration, priors))
  jobs <- expand.grid(vi = seq_along(values), si = seq_along(samples))
  rows <- fan(seq_len(nrow(jobs)), function(j) {
    s <- samples[[jobs$si[j]]]; td <- tds[[jobs$vi[j]]]
    model <- structure(list(theta = s$theta, log_start = td$log_start,
                            log_trans = td$log_trans, n_sub = td$n_sub,
                            emission = emission), class = "nilHMM_model")
    do.call(rbind, lapply(s$obs_list, function(o) data.frame(
      source = source, donor = s$donor, name = s$name, chr = as.integer(o$chr),
      pos = as.integer(o$pos), state = as.integer(decode(model, o)),
      rrate = values[jobs$vi[j]], stringsAsFactors = FALSE)))
  })
  st <- do.call(rbind, rows)
  do.call(rbind, lapply(values, function(v) {
    seg <- to_segments(st[st$rrate == v, setdiff(names(st), "rrate"), drop = FALSE])
    seg$rrate <- v
    seg
  }))
}