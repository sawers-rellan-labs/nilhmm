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

# Decode per (sample, chromosome). Output coordinates are always bp (`pos`); the
# transition decays over `tcol` (bp `pos` or cM `cm`). Returns per-marker state
# rows (source, donor, name, chr, pos, state).
#
# Perf: sort ONCE by (donor, name, chr, pos), then decode over contiguous vector
# slices and assemble a single output data.frame -- avoiding the per-(sample,chr)
# data.frame split/subset/rbind that profiling showed dominated wall-clock (the
# same vectorized-loop win used for to_segments). The C++ decode itself was ~10%
# of the old runtime; the R marshaling was the rest.
.lbimpute_states <- function(data, err, errg, recombdist, drp, log_init,
                             source, donor, has_donor, tcol = "pos", threads = 1L) {
  n <- nrow(data)
  if (n == 0L)
    return(data.frame(source = character(), donor = character(), name = character(),
                      chr = integer(), pos = integer(), state = integer(),
                      stringsAsFactors = FALSE))

  # One global sort into decode+output order. donor is constant within a name, so
  # (name, chr) runs stay contiguous under this key -- and the result needs no
  # re-sort (it already matches the common (donor, name, chr, pos) order).
  # Sort on INTEGER factor codes with radix: a character multi-key order() over a
  # cohort-sized table is O(n log n) string-compares and profiled at ~70% of the
  # runtime; radix on integer codes is linear and drops it to noise. Factor levels
  # sort in the session locale, matching the previous character ordering.
  name_c <- as.integer(factor(data$name))
  chr_i  <- as.integer(data$chr)
  pos_i  <- as.integer(data$pos)
  ord <- if (has_donor)
           order(as.integer(factor(data$donor)), name_c, chr_i, pos_i, method = "radix")
         else order(name_c, chr_i, pos_i, method = "radix")
  nm    <- data$name[ord]
  chr   <- chr_i[ord]
  pos   <- pos_i[ord]
  nref  <- as.integer(data$n_ref[ord])
  nalt  <- as.integer(data$n_alt[ord])
  tpos  <- as.numeric(data[[tcol]][ord])
  dncol <- if (has_donor) as.character(data$donor[ord]) else rep(donor, n)

  # Run starts where (name, chr) changes -> one Viterbi sequence per run.
  starts <- which(c(TRUE, nm[-1] != nm[-n] | chr[-1] != chr[-n]))
  ends   <- c(starts[-1] - 1L, n)

  decode_run <- function(g) {
    i <- starts[g]:ends[g]
    r <- nref[i]; a <- nalt[i]
    # memoize emission over DISTINCT (n_ref, n_alt) pairs (read depths repeat).
    base <- max(a) + 1L                       # base > max(nalt) so key is invertible
    key  <- as.double(r) * base + a
    u    <- unique(key)
    em_u <- lb_emission_loglik_cpp(as.integer(u %/% base), as.integer(u %% base), err, errg)
    em   <- em_u[match(key, u), , drop = FALSE]
    lb_viterbi_cpp(log_init, em, tpos[i], recombdist, drp)   # integer states for this run
  }

  gs <- seq_along(starts)
  runs <- if (threads > 1L && .Platform$OS.type == "unix")
            parallel::mclapply(gs, decode_run, mc.cores = threads)
          else lapply(gs, decode_run)
  bad <- vapply(runs, function(x) inherits(x, "try-error") || is.null(x), logical(1))
  if (any(bad)) {
    who <- unique(paste(nm[starts[bad]], chr[starts[bad]], sep = ":"))
    stop("call_states(): lbimpute decode failed for (name:chr): ", paste(who, collapse = ", "))
  }
  # runs are in sorted-group order, so unlist rebuilds the full state vector aligned
  # to (nm, chr, pos); one data.frame, no per-group rbind.
  state <- unlist(runs, use.names = FALSE)
  data.frame(source = source, donor = dncol, name = nm, chr = chr, pos = pos,
             state = as.integer(state), stringsAsFactors = FALSE)
}
