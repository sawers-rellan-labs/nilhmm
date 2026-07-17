# caller_sweep: sweep a caller's SEGMENTATION parameter with a single shared fit.
#
# The segmentation parameter (rtiger `rigidity`, bbnil `rrate`) is a prior on run
# length / switch rate, not an emission property, so the fitted emission is
# ~parameter-independent (with `fit_means = TRUE` the shared fit is conditioned on
# the reference parameter `ref`, then held fixed across the grid). That lets the
# expensive estimation be amortized across the grid, fanned across cores. This sweep is for CALIBRATION -- finding the best
# parameter -- so it needs segmentations that are comparable across the grid for
# scoring, not per-value fits identical to a cold refit. Two `refit` modes:
#   "none" -- fit ONCE at `ref`, reuse for every decode. Holds the emission fixed
#             so differences across the grid reflect ONLY the segmentation prior
#             (what you're calibrating); exact at `ref`, a close approximation
#             elsewhere. Cheapest (1 fit + N decodes); the recommended scan mode.
#   "cold" -- fit each value from scratch. Exact per value (== per-value
#             call_ancestry); the baseline. Refitting the emission per value can
#             add local-optimum jitter across the grid on multimodal data.
# Workflow: scan with "none" -> pick the best value -> take the exact final calls
# from a single call_ancestry() at that value (a warm-start mode was evaluated and
# dropped: on multimodal fits it reaches a different optimum than cold, so it is
# neither exact nor better than "none" for ranking).

#' Sweep a caller's segmentation parameter with an amortized fit
#'
#' Fit the emission once (rtiger: the joint EM; bbnil: per-sample means, or nothing
#' when `fit_means = FALSE`) and sweep `values`, fanning the decodes over
#' `threads`. `refit` controls the accuracy/speed trade (see details).
#'
#' @param data Common input: `name, chr, pos, n_ref, n_alt` (+ optional `donor`;
#'   `lbimpute` with `unit = "cm"` also needs a `cm` map-position column).
#' @param caller `"rtiger"` (sweeps `rigidity`), `"bbnil"` or `"nnil"` (both sweep
#'   `rrate`; `bbnil` = count/BetaBinomial emission on read counts, `nnil` =
#'   categorical genotype emission on a hard-called `g` column), or `"lbimpute"`
#'   (sweeps `recombdist`).
#' @param values Parameter grid to sweep.
#' @param refit `"none"` (fit once at `ref` and reuse -- exact at `ref`, a close
#'   approximation elsewhere; recommended for calibration, as it isolates the
#'   segmentation prior) or `"cold"` (fit each value from scratch -- exact per
#'   value, the baseline). For bbnil with `fit_means = FALSE` the emission is
#'   `rrate`-independent, so both are identical (and exact). This sweep *finds* the
#'   best value; for the exact final calls, refit once with [call_ancestry()] at
#'   the chosen value. **Ignored for `lbimpute`**: `recombdist` touches only the
#'   transition, never the emission, so every swept value is already EXACT --
#'   identical to a cold `call_ancestry(caller = "lbimpute", recombdist = v)` (the
#'   emission is computed once per run and only the Viterbi transition is re-run
#'   over the grid, batched in C++).
#' @param design,f_1,f_2 Population priors (a design name, or explicit `f_1`,`f_2`).
#'   For `lbimpute` these only seed the start distribution (flat if absent).
#' @param threads Fan-out width (`parallel::mclapply` on unix; serial otherwise).
#' @param ref Reference value for the shared fit (default `median(values)`,
#'   rounded for rtiger). Unused for `lbimpute` (no shared fit).
#' @param min_reads Minimum read depth to keep a marker before decoding (default
#'   `1L`); `0L` keeps all. **No-op for `lbimpute`**, which keeps zero-read markers
#'   (flat emission) so the distance transition sees true marker spacing.
#' @param err,conc,fit_means bbnil count-emission parameters (fixed across the grid).
#'   `err` is also LB-Impute's per-read error (`readerr`).
#' @param seed,postprocess rtiger fit seed and border post-processing.
#' @param unit,genotypeerr,drp `lbimpute` only: `unit` is the transition coordinate
#'   (`"bp"` physical / `"cm"` genetic map; output stays bp), `genotypeerr` the
#'   emission floor/ceiling, `drp` the single-vs-double homozygous-switch cost.
#'   Same validation/mismatch warnings as [call_states()].
#' @param source,donor Output labels.
#' @return A common-schema segment table with an added column named for the swept
#'   parameter (`rigidity`, `rrate`, or `recombdist`) tagging each value's calls.
#' @export
caller_sweep <- function(data, caller = c("nnil", "bbnil", "rtiger", "lbimpute"), values,
                         refit = c("none", "cold"),
                         design = NULL, f_1 = NULL, f_2 = NULL, threads = 1L,
                         ref = NULL, min_reads = 1L, err = 0.01, conc = 20,
                         fit_means = FALSE, seed = 1L, postprocess = TRUE,
                         unit = c("bp", "cm"), genotypeerr = 0.05, drp = FALSE,
                         source = "nilHMM", donor = NA_character_) {
  caller <- match.arg(caller)
  refit <- match.arg(refit)
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

  # ---------------- lbimpute: sweep `recombdist`, EXACT per value --------------
  # recombdist affects only the transition, never the emission, so each swept value
  # is identical to a cold call_ancestry(caller = "lbimpute", recombdist = v) -- no
  # refit approximation (the `refit` arg is inapplicable and ignored). Routed
  # before the min_reads filter: lbimpute keeps zero-read markers (flat emission)
  # so the distance transition sees true marker spacing (min_reads is a no-op here,
  # matching call_states for lbimpute).
  if (caller == "lbimpute") {
    unit <- match.arg(unit)
    return(.caller_sweep_lbimpute(data, values, unit, err, genotypeerr, drp,
                                  design, f_1, f_2, threads, source, donor, has_donor))
  }

  if (!is.null(min_reads) && min_reads > 0L)
    data <- data[data$n_ref + data$n_alt >= min_reads, , drop = FALSE]
  if (!nrow(data)) stop("caller_sweep(): no markers with >= min_reads reads")

  # ---------------- rtiger: EM (once / warm / cold) + decode per rigidity -----
  if (caller == "rtiger") {
    values <- as.integer(values)
    o <- .rtiger_obs(data, has_donor, donor)
    .rtiger_check_coverage(o$obs, max(values))            # feasibility for the whole grid
    r_ref <- if (is.null(ref)) as.integer(round(stats::median(values))) else as.integer(ref)
    ref_fit <- if (refit == "cold") NULL
               else .rtiger_fit(o$obs, r_ref, threads = threads, seed = seed)   # shared fit
    segs <- fan(values, function(v) {
      fit <- if (refit == "cold")
               .rtiger_fit(o$obs, v, threads = 1L, seed = seed)                 # exact, from scratch
             else ref_fit                                                        # none: reuse the ref fit
      paths <- .rtiger_decode(o$obs, fit, v, postprocess = postprocess, threads = 1L)
      seg <- to_segments(.rtiger_assemble(o$obs, o$pos, paths, o$donor_of, source))
      seg$rigidity <- v
      seg
    })
    return(do.call(rbind, segs))
  }

  # ---- nnil (categorical gt) / bbnil (count): geometric caller, decode per rrate --
  # Both are geometric-duration callers differing only in emission. nnil uses the
  # categorical genotype-confusion emission on hard-called genotypes (a `g` column
  # in {0,1,2,3}); bbnil uses the count/BetaBinomial emission on read counts.
  # caller_spec() supplies the matching emission + duration for the requested caller.
  is_gt <- caller == "nnil"
  if (is_gt && !("g" %in% names(data)))
    stop("caller_sweep(nnil): needs a hard-called `g` column (0/1/2/3). The ",
         "categorical emission does not threshold read counts -- hard-call them ",
         "first with call_gt().")
  priors <- .state_freqs(design, f_1, f_2, "caller_sweep")
  spec_of <- function(v) caller_spec(caller, rrate = v, err = err, conc = conc, fit_means = fit_means)
  ref_rrate <- if (is.null(ref)) stats::median(values) else ref
  emission <- spec_of(ref_rrate)$emission  # rrate-independent; gt emission ignores theta
  theta0 <- .emission_theta(emission)      # NA for the categorical (gt) emission
  td_of <- function(v) .duration_transition(spec_of(v)$duration, priors)
  by_name <- split(data, data$name, drop = TRUE)

  # Phase 1 (parallel): per-sample chains + (for none/warm) a reference emission
  # fit at ref_rrate. Fixed-means needs no fit; cold refits per value in phase 2.
  td_ref <- td_of(ref_rrate)
  need_ref <- isTRUE(fit_means) && refit == "none" && !is_gt   # gt has no means to fit
  samples <- fan(names(by_name), function(nm) {
    dn <- by_name[[nm]]
    obs_list <- lapply(split(dn, dn$chr, drop = TRUE), function(dc) {
      dc <- dc[order(dc$pos), , drop = FALSE]
      list(chr = dc$chr[1], pos = dc$pos, n = dc$n_ref + dc$n_alt, a = dc$n_alt,
           g = if ("g" %in% names(dc)) as.integer(dc$g) else NULL)
    })
    ref_theta <- if (need_ref) .em_fit_means(obs_list, emission, td_ref, theta0) else theta0
    list(name = nm, donor = if (has_donor) dn$donor[1] else donor,
         obs_list = obs_list, ref_theta = ref_theta)
  })

  # Phase 2 (parallel over value x sample): theta per (value, sample) per `refit`,
  # then decode. Precompute transitions serially (cheap).
  tds <- lapply(values, td_of)
  jobs <- expand.grid(vi = seq_along(values), si = seq_along(samples))
  # Segment IN-WORKER (like the rtiger path) and return compact segments, not raw
  # per-marker states. Aggregating raw markers in the master -- st <- rbind(rows)
  # over |values| x n_samples per-marker frames -- OOMs on large panels (118K x 200
  # x |values| ~ 1e8 rows, past the per-process vector cap). to_segments RLE never
  # spans a `name` change (runs break on it) and re-sorts its output by
  # (donor, name, chr, start_bp), so segmenting per (value, sample) and re-applying
  # that sort per value is byte-identical to a whole-cohort to_segments per value.
  seglist <- fan(seq_len(nrow(jobs)), function(j) {
    s <- samples[[jobs$si[j]]]; td <- tds[[jobs$vi[j]]]; v <- values[jobs$vi[j]]
    theta <- if (is_gt || !isTRUE(fit_means)) theta0   # gt / fixed-means: rrate-independent (exact)
             else if (refit == "none") s$ref_theta     # fit once at ref, reuse
             else .em_fit_means(s$obs_list, emission, td, theta0)   # cold: refit per value
    model <- structure(list(theta = theta, log_start = td$log_start,
                            log_trans = td$log_trans, n_sub = td$n_sub,
                            emission = emission), class = "nilHMM_model")
    markers <- do.call(rbind, lapply(s$obs_list, function(o) data.frame(
      source = source, donor = s$donor, name = s$name, chr = as.integer(o$chr),
      pos = as.integer(o$pos), state = as.integer(decode(model, o)),
      stringsAsFactors = FALSE)))
    seg <- to_segments(markers)
    seg$rrate <- v
    seg
  })
  do.call(rbind, lapply(seq_along(values), function(vi) {
    seg <- do.call(rbind, seglist[jobs$vi == vi])
    seg[order(seg$donor, seg$name, seg$chr, seg$start_bp), , drop = FALSE]
  }))
}