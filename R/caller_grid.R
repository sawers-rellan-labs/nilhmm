# caller_grid: sweep an arbitrary emission x transition parameter grid for the
# geometric callers (nnil / bbnil), on the unified viterbi_sweep_cpp kernel.
#
# The optimized factorial: the transition (rrate) enters only log_trans, and for
# geometric duration log_start and n_sub are rrate-independent (.build_transition),
# so one shared log_init + one 3 x 3 log_trans per rrate suffice; the emission
# enters only log_emit. So we build log_emit ONCE per (emission combo, sample) and
# batch the whole rrate grid through viterbi_sweep_cpp in C++. Each grid point is
# byte-identical to a cold call_ancestry() at that (emission, rrate) -- the
# emission is fixed-means (no EM), exactly as the geometric decode is.
#
# This is the general "sweep any parameter" driver: caller_sweep is its single-axis
# (rrate-only, default emission) special case, and it retires the ad-hoc
# emission-grid scripts (Holland's 945-config nir x germ x gert x p x r, the
# MolBreeding calibration). Rigidity (rtiger, n_sub > 1) is out of scope here.

#' Sweep an emission x transition parameter grid (geometric callers)
#'
#' Runs the full factorial of an emission-parameter grid crossed with a
#' transition (`rrate`) grid for `nnil` (categorical genotype emission) or `bbnil`
#' (count/BetaBinomial, fixed means), reusing one C++ transition sweep
#' ([viterbi_sweep_cpp()]) per (emission combo, sample). Each grid point equals a
#' cold [call_ancestry()] at that configuration.
#'
#' @param data Long table `name, chr, pos` plus the caller's emission input: a
#'   hard-called `g` column (`0/1/2/3`) for `nnil`, or `n_ref, n_alt` read counts
#'   for `bbnil` (+ optional `donor`).
#' @param caller `"nnil"` or `"bbnil"` (geometric-duration callers). Rigidity
#'   (`rtiger`) is not supported -- its transition expands the state space
#'   (`n_sub > 1`); use [caller_sweep()].
#' @param emission_grid A data.frame whose columns are emission parameters and
#'   whose rows are the combos to sweep. For `nnil`: any of `germ, gert, p, mr,
#'   nir`. For `bbnil`: `err, conc`. Absent parameters take [caller_spec()]
#'   defaults. Every column is echoed onto the output as a tag.
#' @param rrate Numeric transition grid (the segmentation switch rate), swept in
#'   C++ with the emission held fixed.
#' @param design,f_1,f_2 Population priors (a design name, or explicit `f_1,f_2`).
#' @param err,conc Count-emission parameters when not columns of `emission_grid`.
#' @param threads Fan-out width over (emission combo x sample) jobs.
#' @param min_reads Minimum depth to keep a marker (count callers; no-op for nnil).
#' @param source,donor Output labels.
#' @return A common-schema segment table (`source, donor, name, chr, start_bp,
#'   end_bp, state`) with a `rrate` column and one column per `emission_grid`
#'   parameter, tagging every segment's grid point.
#' @export
caller_grid <- function(data, caller = c("nnil", "bbnil"),
                        emission_grid, rrate,
                        design = NULL, f_1 = NULL, f_2 = NULL,
                        err = 0.01, conc = 20,
                        threads = 1L, min_reads = 1L,
                        source = "nilHMM", donor = NA_character_) {
  caller <- match.arg(caller)
  if (!all(c("name", "chr", "pos") %in% names(data)))
    stop("caller_grid(): data needs columns name, chr, pos")
  if (!length(rrate)) stop("caller_grid(): `rrate` is empty")
  eg <- as.data.frame(emission_grid, stringsAsFactors = FALSE)
  if (!nrow(eg)) stop("caller_grid(): `emission_grid` has no rows")
  is_gt <- caller == "nnil"
  has_counts <- all(c("n_ref", "n_alt") %in% names(data))
  if (is_gt) {
    if (!("g" %in% names(data)))
      stop("caller_grid(nnil): needs a hard-called `g` column (0/1/2/3) from call_gt().")
  } else if (!has_counts) {
    stop("caller_grid(bbnil): needs `n_ref`/`n_alt` read counts.")
  }
  has_donor <- "donor" %in% names(data)
  # marker-support filter, matching call_states(): the categorical (nnil) path
  # treats a missing genotype (g == 3, "no coverage") as an uncovered marker and
  # drops it when min_reads > 0 (min_reads = 0 decodes every marker); count callers
  # drop markers with < min_reads reads. Keeps the grid on the same marker support
  # a cold call_ancestry() uses.
  if (!is.null(min_reads) && min_reads > 0L) {
    data <- if (is_gt) data[data$g != 3L, , drop = FALSE]
    else data[data$n_ref + data$n_alt >= min_reads, , drop = FALSE]
  }
  if (!nrow(data)) stop("caller_grid(): no markers left after the marker-support filter")

  priors <- .state_freqs(design, f_1, f_2, "caller_grid")

  # transition grid: geometric -> log_start & n_sub are rrate-independent, so a
  # single shared log_init + one 3 x 3 log_trans per rrate.
  tds <- lapply(rrate, function(r) .duration_transition(duration_geometric(r), priors))
  if (any(vapply(tds, function(td) td$n_sub, integer(1)) != 1L))
    stop("caller_grid(): only geometric callers (n_sub = 1) are supported; ",
         "use caller_sweep() for rtiger rigidity.")
  log_start <- tds[[1]]$log_start
  trans_list <- lapply(tds, `[[`, "log_trans")
  tb <- if (is_gt) 1L else 0L # GT ties structurally -> incumbent, matching decode()

  # one emission spec per combo (fixed means; no EM)
  specs <- lapply(seq_len(nrow(eg)), function(ei) {
    args <- c(list(caller = caller, err = err, conc = conc), as.list(eg[ei, , drop = FALSE]))
    do.call(caller_spec, args)$emission
  })
  if (any(vapply(specs, function(e) isTRUE(e$fit_means), logical(1))))
    stop("caller_grid(): fit_means = TRUE is unsupported (emission must be fixed across the grid).")

  # per-sample observation lists (built once, shared across all combos)
  by_name <- split(seq_len(nrow(data)), data$name)
  samples <- lapply(names(by_name), function(nm) {
    dn <- data[by_name[[nm]], , drop = FALSE]
    obs_list <- lapply(split(seq_len(nrow(dn)), dn$chr), function(ci) {
      dc <- dn[ci, , drop = FALSE]
      dc <- dc[order(dc$pos), , drop = FALSE]
      list(chr = dc$chr[1], pos = dc$pos,
           n = if (has_counts) dc$n_ref + dc$n_alt else NULL,
           a = if (has_counts) dc$n_alt else NULL,
           g = if ("g" %in% names(dc)) as.integer(dc$g) else NULL)
    })
    list(name = nm, donor = if (has_donor) dn$donor[1] else donor, obs_list = obs_list)
  })

  fan <- function(X, FUN) {
    r <- if (threads > 1L && .Platform$OS.type == "unix")
      parallel::mclapply(X, FUN, mc.cores = threads) else lapply(X, FUN)
    if (any(vapply(r, function(x) inherits(x, "try-error") || is.null(x), logical(1))))
      stop("caller_grid(): a decode job failed")
    r
  }

  # one job per (emission combo, sample); each sweeps the whole rrate grid in C++
  jobs <- expand.grid(ei = seq_len(nrow(eg)), si = seq_along(samples))
  seglist <- fan(seq_len(nrow(jobs)), function(j) {
    ei <- jobs$ei[j]; s <- samples[[jobs$si[j]]]
    emission <- specs[[ei]]; theta <- .emission_theta(emission)
    per_chr <- lapply(s$obs_list, function(o) {
      em <- .emission_loglik(emission, o, theta) # T x 3, fixed across rrate
      paths <- viterbi_sweep_cpp(log_start, em, mode = "const",
                                 trans_list = trans_list, tie_break = tb) # T x V
      lapply(seq_along(rrate), function(vi) {
        mk <- data.frame(source = source, donor = s$donor, name = s$name,
                         chr = as.integer(o$chr), pos = as.integer(o$pos),
                         state = as.integer(paths[, vi]), stringsAsFactors = FALSE)
        seg <- to_segments(mk)
        seg$rrate <- rrate[vi]
        seg
      })
    })
    # bind chromosomes within each rrate, tag the emission-combo params
    combo <- eg[ei, , drop = FALSE]
    do.call(rbind, lapply(seq_along(rrate), function(vi) {
      seg <- do.call(rbind, lapply(per_chr, `[[`, vi))
      if (is.null(seg) || !nrow(seg)) return(NULL)
      for (col in names(combo)) seg[[col]] <- combo[[col]]
      seg
    }))
  })
  out <- do.call(rbind, seglist)
  out[order(out$donor, out$name, out$chr, out$start_bp), , drop = FALSE]
}