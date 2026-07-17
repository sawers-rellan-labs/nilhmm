# LB-Impute driver: thin R orchestration over the Rcpp kernels in
# src/lbimpute.cpp. LB-Impute (Fragoso et al. 2014, G3;
# dellaporta-laboratory/LB-Impute) is, structurally, a 3-state REF/HET/ALT HMM --
# the chain nilHMM is built around. Its distinctive pieces are the coverage-aware
# emission (lb_emission_loglik_cpp) and the distance-dependent transition with a
# double-recombination penalty (lb_viterbi_cpp), both faithful to the reference
# Java. The decoder is the engine's full-chromosome Viterbi rather than
# LB-Impute's windowed forward/reverse consensus (the "native" choice: the
# optimal path the window approximates; see src/lbimpute.cpp).
#
# Like rtiger/binhmm this is routed separately from the count engine because its
# transition is per-marker (distance-based), not the single time-homogeneous
# matrix the count/gt callers share.

# --- shared parameter resolution (used by call_states AND caller_sweep) --------

# Transition-coordinate column, with the same value-based validation call_states
# applies (checks are on VALUES, not storage type -- R literals default to double,
# so an is.integer test would misread the natural `c(1e6, 2e6)` bp input as cM).
.lbimpute_tcol <- function(data, unit, ctx = "call_states") {
  tcol <- if (unit == "cm") "cm" else "pos"
  if (unit == "cm" && !("cm" %in% names(data)))
    stop(ctx, "(): caller = 'lbimpute' with unit = 'cm' needs a `cm` column of ",
         "genetic (map) positions per marker")
  tvals <- data[[tcol]]
  if (anyNA(tvals))
    stop(ctx, "(): `", tcol, "` (lbimpute transition coordinate) contains NA")
  if (unit == "bp" && any(tvals != floor(tvals)))
    stop(ctx, "(): bp positions must be whole numbers ",
         "(fractional values look like cM -- did you mean unit = 'cm'?)")
  tcol
}

# Warn on the unit/recombdist mismatches that silently ruin the calls. `recombdist`
# may be a scalar (call_states) or a grid (caller_sweep); the checks are vectorized.
.lbimpute_warn_recombdist <- function(recombdist, unit, tvals, ctx = "call_states") {
  if (unit == "cm" && any(recombdist > 1000))
    warning(ctx, "(): unit = 'cm' but recombdist has bp-sized value(s) (> 1000); ",
            "cM `recombdist` is typically ~50. The transition will barely relax ",
            "and the path may collapse to one segment.")
  if (unit == "bp" && any(recombdist < 1000))
    warning(ctx, "(): unit = 'bp' but recombdist has cM-sized value(s) (< 1000); ",
            "bp `recombdist` is typically ~1e7. The transition will over-relax ",
            "and the path may over-fragment.")
  if (unit == "cm" && max(tvals) > 1000 && all(tvals == floor(tvals)))
    warning(ctx, "(): unit = 'cm' but the `cm` values are all whole and exceed ",
            "1000 -- they look like bp. Check the coordinate units.")
  invisible(NULL)
}

# Start distribution: design/f_1,f_2 only SEED the start (LB-Impute has no
# state-frequency prior in its transition); absent, flat.
.lbimpute_log_init <- function(design, f_1, f_2, ctx = "call_states") {
  if (is.null(design) && (is.null(f_1) || is.null(f_2))) return(log(rep(1 / 3, 3)))
  fs <- .state_freqs(design, f_1, f_2, ctx)
  log(c(1 - fs$f_1 - fs$f_2, fs$f_1, fs$f_2))
}

# --- shared decode machinery ---------------------------------------------------

# One global radix sort into decode+output order + per-(name, chr) run boundaries.
# donor is constant within a name, so (name, chr) runs stay contiguous under this
# key and the result needs no re-sort. Radix on INTEGER factor codes: a character
# multi-key order() over a cohort-sized table profiled at ~70% of runtime; radix
# on integer codes is linear. Factor levels sort in the session locale, matching
# the previous character ordering. Output coordinates are always bp (`pos`); the
# transition decays over `tcol` (bp `pos` or cM `cm`).
.lbimpute_prep <- function(data, tcol, donor, has_donor) {
  n <- nrow(data)
  name_c <- as.integer(factor(data$name))
  chr_i  <- as.integer(data$chr)
  pos_i  <- as.integer(data$pos)
  ord <- if (has_donor)
           order(as.integer(factor(data$donor)), name_c, chr_i, pos_i, method = "radix")
         else order(name_c, chr_i, pos_i, method = "radix")
  nm     <- data$name[ord]
  chr    <- chr_i[ord]
  starts <- which(c(TRUE, nm[-1] != nm[-n] | chr[-1] != chr[-n]))
  list(nm = nm, chr = chr, pos = pos_i[ord],
       nref = as.integer(data$n_ref[ord]), nalt = as.integer(data$n_alt[ord]),
       tpos = as.numeric(data[[tcol]][ord]),
       dncol = if (has_donor) as.character(data$donor[ord]) else rep(donor, n),
       starts = starts, ends = c(starts[-1] - 1L, n))
}

# Memoize the emission over DISTINCT (n_ref, n_alt) pairs, index back (read depths
# repeat heavily). Returns the T x 3 log-emission matrix for a run.
.lb_run_emission <- function(r, a, err, errg) {
  base <- max(a) + 1L                            # base > max(nalt) so the key is invertible
  key  <- as.double(r) * base + a
  u    <- unique(key)
  em_u <- lb_emission_loglik_cpp(as.integer(u %/% base), as.integer(u %% base), err, errg)
  em_u[match(key, u), , drop = FALSE]
}

.lb_fan <- function(X, FUN, threads) {
  runs <- if (threads > 1L && .Platform$OS.type == "unix")
            parallel::mclapply(X, FUN, mc.cores = threads) else lapply(X, FUN)
  if (any(vapply(runs, function(x) inherits(x, "try-error") || is.null(x), logical(1))))
    stop("call_states(): lbimpute decode failed for a (name, chr) run")
  runs
}

# Single decode: one Viterbi per (name, chr) run at a fixed recombdist. Returns the
# per-marker state table (source, donor, name, chr, pos, state).
.lbimpute_states <- function(data, err, errg, recombdist, drp, log_init,
                             source, donor, has_donor, tcol = "pos", threads = 1L) {
  if (nrow(data) == 0L)
    return(data.frame(source = character(), donor = character(), name = character(),
                      chr = integer(), pos = integer(), state = integer(),
                      stringsAsFactors = FALSE))
  p <- .lbimpute_prep(data, tcol, donor, has_donor)
  runs <- .lb_fan(seq_along(p$starts), function(g) {
    i <- p$starts[g]:p$ends[g]
    em <- .lb_run_emission(p$nref[i], p$nalt[i], err, errg)
    lb_viterbi_cpp(log_init, em, p$tpos[i], recombdist, drp)
  }, threads)
  data.frame(source = source, donor = p$dncol, name = p$nm, chr = p$chr, pos = p$pos,
             state = as.integer(unlist(runs, use.names = FALSE)), stringsAsFactors = FALSE)
}

# Batched sweep: decode every recombdist in `recombdists` reusing one emission per
# run (recombdist affects only the transition), so each value is EXACT -- identical
# to a cold .lbimpute_states at that recombdist. Returns `base` (the per-marker
# source/donor/name/chr/pos frame, in sorted order) and `states` (an n x V integer
# matrix; column k is the decode at recombdists[k]). Memory is O(n_markers x
# n_values); calibration grids are small.
#
# TODO(memory): the n_markers x n_values state matrix is held whole. Fine for
# golden-section refinement (a few values per iteration); for a very wide grid over
# the full cohort at once, chunk `recombdists` (loop caller_sweep over value blocks
# and rbind) or stream to_segments per column instead of materializing all columns.
.lbimpute_sweep <- function(data, err, errg, recombdists, drp, log_init,
                            source, donor, has_donor, tcol = "pos", threads = 1L) {
  recombdists <- as.numeric(recombdists)
  if (nrow(data) == 0L)
    return(list(base = data.frame(source = character(), donor = character(),
                                  name = character(), chr = integer(), pos = integer(),
                                  stringsAsFactors = FALSE),
                states = matrix(integer(0), 0L, length(recombdists)),
                values = recombdists))
  p <- .lbimpute_prep(data, tcol, donor, has_donor)
  runs <- .lb_fan(seq_along(p$starts), function(g) {
    i <- p$starts[g]:p$ends[g]
    em <- .lb_run_emission(p$nref[i], p$nalt[i], err, errg)
    lb_viterbi_sweep_cpp(log_init, em, p$tpos[i], recombdists, drp)   # n_i x V
  }, threads)
  list(base = data.frame(source = source, donor = p$dncol, name = p$nm,
                         chr = p$chr, pos = p$pos, stringsAsFactors = FALSE),
       states = do.call(rbind, runs), values = recombdists)
}

# caller_sweep(caller = "lbimpute"): sweep `recombdist` over `values`, returning
# the common segment schema with a `recombdist` column tagging each value. EXACT
# per value (recombdist touches only the transition), so each is identical to
# call_ancestry(caller = "lbimpute", recombdist = v). Keeps zero-coverage markers
# (flat emission) so the distance transition sees true spacing; `min_reads` is a
# no-op here, matching call_states for lbimpute.
.caller_sweep_lbimpute <- function(data, values, unit, err, errg, drp,
                                   design, f_1, f_2, threads, source, donor, has_donor) {
  values <- as.numeric(values)
  if (any(values <= 0)) stop("caller_sweep(lbimpute): all `values` (recombdist) must be > 0")
  if (nrow(data) == 0L)                              # empty -> empty common schema + recombdist
    return(data.frame(source = character(), donor = character(), name = character(),
                      chr = integer(), start_bp = integer(), end_bp = integer(),
                      state = integer(), recombdist = numeric(), stringsAsFactors = FALSE))
  tcol <- .lbimpute_tcol(data, unit, ctx = "caller_sweep")
  .lbimpute_warn_recombdist(values, unit, data[[tcol]], ctx = "caller_sweep")
  log_init <- .lbimpute_log_init(design, f_1, f_2, ctx = "caller_sweep")
  sw <- .lbimpute_sweep(data, err, errg, values, drp, log_init,
                        source, donor, has_donor, tcol, threads)
  do.call(rbind, lapply(seq_along(values), function(k) {
    seg <- to_segments(cbind(sw$base, state = sw$states[, k]))
    seg$recombdist <- values[k]
    seg
  }))
}
