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

# Per (sample, chromosome) decode. `data` has name, chr, pos, n_ref, n_alt
# (+ optional donor). Returns per-marker state rows (source, donor, name, chr,
# pos, state) for to_segments().
.lbimpute_states <- function(data, err, errg, recombdist, drp, log_init,
                             source, donor, has_donor, threads = 1L) {
  by_name <- split(data, data$name, drop = TRUE)

  one_sample <- function(nm) {
    dn <- by_name[[nm]]
    donor_nm <- if (has_donor) dn$donor[1] else donor
    do.call(rbind, lapply(split(dn, dn$chr, drop = TRUE), function(dc) {
      dc <- dc[order(dc$pos), , drop = FALSE]
      nref <- as.integer(dc$n_ref); nalt <- as.integer(dc$n_alt)
      # memoize the emission over DISTINCT (n_ref, n_alt) pairs, index back (the
      # count-emission trick in emissions.R): read depths repeat heavily.
      base <- max(nref + nalt) + 1L
      key  <- as.double(nref) * base + nalt
      u    <- unique(key)
      em_u <- lb_emission_loglik_cpp(as.integer(u %/% base), as.integer(u %% base), err, errg)
      em   <- em_u[match(key, u), , drop = FALSE]
      path <- lb_viterbi_cpp(log_init, em, as.integer(dc$pos), recombdist, drp)
      data.frame(source = source, donor = donor_nm, name = nm,
                 chr = as.integer(dc$chr[1]), pos = as.integer(dc$pos),
                 state = as.integer(path), stringsAsFactors = FALSE)
    }))
  }

  nms <- names(by_name)
  per <- if (threads > 1L && .Platform$OS.type == "unix")
           parallel::mclapply(nms, one_sample, mc.cores = threads)
         else lapply(nms, one_sample)
  bad <- vapply(per, function(x) inherits(x, "try-error") || is.null(x), logical(1))
  if (any(bad)) stop("call_states(): lbimpute decode failed for sample(s): ",
                     paste(nms[bad], collapse = ", "))
  states <- if (requireNamespace("data.table", quietly = TRUE))
              as.data.frame(data.table::rbindlist(per)) else do.call(rbind, per)
  states[order(states$donor, states$name, states$chr, states$pos), , drop = FALSE]
}
