# RTIGER driver: thin R orchestration over the Rcpp port in src/rtiger.cpp.
# The entire EM fit (E-step + M-steps + Brent emission + convergence loop) lives
# in C++ (rtiger_fit_cpp), faithfully reproducing the fork's fit/EM
# (rHMM_methods.jl). This R layer only builds the per-chain (k, n) observations
# and runs the Viterbi decode pass; it supersedes the earlier invented
# duration_rigidity/p_switch rtiger path.
#
# Observation convention (generateObject): per sample x chromosome, k = ref/P1
# count, n = total. States [pat,het,mat] = (ref~0.95, 0.5, ref~0.05); remapped
# to the common schema as common = rtiger_state - 1 (pat->REF, het->HET, mat->ALT).

# Flatten obs (list of samples, each a list of per-chr list(k=, n=)) into the
# flat k/n chain lists rtiger_fit_cpp consumes.
.rtiger_chains <- function(obs) {
  ks_list <- list(); ns_list <- list()
  for (sm in obs) for (ch in sm) {
    ks_list[[length(ks_list) + 1L]] <- as.integer(ch$k)
    ns_list[[length(ns_list) + 1L]] <- as.integer(ch$n)
  }
  list(ks = ks_list, ns = ns_list)
}

# Fit RTIGER parameters by the full C++ EM. Returns list(A, pi, alpha, beta, iterations).
# `threads` parallelizes the per-chain E-step (RcppParallel); result is
# deterministic for a fixed thread count (Viterbi-identical to threads=1).
# RTIGER's randomized emission init (generate_params, randomize=TRUE): jitter the
# canonical alpha=[20,20,1] beta=[1,20,20] by runif and scale by priorstrength
# ~1+-0.5. Seeded for reproducibility (RTIGER's own is unseeded). A/pi keep the
# C++ deterministic diagonal-dominant default (the emission means drive the fit;
# transition/start converge quickly).
.rtiger_init <- function(nstates = 3L, seed = 1L) {
  set.seed(seed)
  ps <- 1 + runif(1, -0.5, 0.5)
  av <- (if (nstates == 3L) c(20, 20, 1) else rep(20, nstates)) + runif(nstates, -0.5, 0.5)
  bv <- (if (nstates == 3L) c(1, 20, 20) else rep(20, nstates)) + runif(nstates, -0.5, 0.5)
  list(alpha = av * ps, beta = bv * ps)
}

# Default init is now RTIGER's seeded randomized init (faithful to the original);
# pass init_alpha/init_beta to override, or seed= to change the draw.
.rtiger_fit <- function(obs, r, nstates = 3L, eps = 0.01, max_iter = 50L, threads = 1L,
                        init_alpha = NULL, init_beta = NULL, seed = 1L) {
  if (is.null(init_alpha) || is.null(init_beta)) {
    ini <- .rtiger_init(nstates, seed); init_alpha <- ini$alpha; init_beta <- ini$beta
  }
  ch <- .rtiger_chains(obs)
  rtiger_fit_cpp(ch$ks, ch$ns, as.integer(r), as.integer(nstates), eps,
                 as.integer(max_iter), as.integer(threads),
                 as.numeric(init_alpha), as.numeric(init_beta))
}

# Total data log-likelihood of `obs` under fitted `params`: sum over chains of
# logsumexp_k(forward[k,r] + backward[k,r]) (the un-normalized gammar column,
# which is logP(O) for the chain). Used to rank multi-start fits.
.rtiger_loglik <- function(obs, params, r) {
  logPI <- log(params$pi); logA <- log(params$A)
  lse <- function(x) { m <- max(x); if (!is.finite(m)) return(-Inf); m + log(sum(exp(x - m))) }
  total <- 0
  for (sm in obs) for (ch in sm) {
    lp <- rtiger_getlogpsi_cpp(ch$k, ch$n, params$alpha, params$beta)
    LP <- rtiger_productpsi_cpp(lp, r)
    al <- rtiger_forward_cpp(logPI, LP, logA, lp, r)
    be <- rtiger_backward_cpp(LP, logA, lp, r)
    total <- total + lse(al[, r] + be[, r])
  }
  total
}

# Subsample-probe multi-start: probe M random inits on a small subsample for a
# few iters, keep the highest-likelihood one, then run the full fit warm-started
# from it. Robust to the local-optimum traps a single deterministic init can hit
# (e.g. Zv). Costs ~one full fit plus the cheap probes.
.rtiger_fit_multistart <- function(obs, r, nstates = 3L, eps = 0.01, max_iter = 50L,
                                   threads = 1L, M = 6L, sub_n = 50L, probe_iter = 10L,
                                   seed = 1L, verbose = FALSE) {
  set.seed(seed)
  sub <- obs[utils::head(names(obs), sub_n)]
  # candidate inits: the canonical one + (M-1) jittered over the het/donor mean space
  cand <- list(list(a = c(20, 20, 1), b = c(1, 20, 20)))            # canonical (= default)
  for (m in seq_len(max(0L, M - 1L))) {
    mu  <- c(runif(1, 0.90, 0.99), runif(1, 0.30, 0.70), runif(1, 0.02, 0.20))  # pat/het/donor means
    tau <- runif(3, 5, 30)
    cand[[length(cand) + 1L]] <- list(a = tau * mu, b = tau * (1 - mu))
  }
  best <- NULL; best_ll <- -Inf
  for (m in seq_along(cand)) {
    f  <- .rtiger_fit(sub, r, nstates, eps, probe_iter, threads, cand[[m]]$a, cand[[m]]$b)
    ll <- .rtiger_loglik(sub, f, r)
    if (verbose) cat(sprintf("  probe %d: loglik=%.1f means=%s\n", m, ll,
                             paste(round(f$alpha / (f$alpha + f$beta), 3), collapse = ",")))
    if (ll > best_ll) { best_ll <- ll; best <- f }
  }
  # warm-start the full fit from the winning probe's evolved emission means
  .rtiger_fit(obs, r, nstates, eps, max_iter, threads, best$alpha, best$beta)
}

# Border post-processing (port of the fork's postprocessing/findsplit/segmente/
# getVitSeg). For each pair of adjacent segments, re-place the boundary at the
# split that maximizes the per-marker emission likelihood: assign the first m-1
# of the two-segment span to the left state and the rest to the right, pick the
# best m. `psi` is the s x T per-marker log-emission (getlogpsi, NOT windowed).
# Sequential left-to-right (only the right-border branch, as in the fork).
.rtiger_postprocess <- function(path, psi) {
  Tn <- length(path)
  if (Tn < 2L) return(path)
  brk    <- which(path[-1L] != path[-Tn])
  starts <- c(1L, brk + 1L); ends <- c(brk, Tn); states <- path[starts]
  nseg   <- length(starts)
  if (nseg < 2L) return(path)
  for (i in seq_len(nseg - 1L)) {
    lo <- starts[i]; hi <- ends[i + 1L]               # span of segments i and i+1
    dec <- psi[, lo:hi, drop = FALSE]
    ls  <- cumsum(dec[states[i], ])                   # left state, prefix sums
    rs  <- rev(cumsum(rev(dec[states[i + 1L], ])))    # right state, suffix sums
    e   <- which.max(c(0, ls) + c(rs, 0))             # best split (1-based)
    ends[i] <- lo + e - 2L; starts[i + 1L] <- lo + e - 1L
  }
  out <- integer(Tn)
  for (i in seq_len(nseg)) if (starts[i] <= ends[i]) out[starts[i]:ends[i]] <- states[i]
  out
}

# Decode all chains with the fitted parameters (one Viterbi pass per chain),
# optionally applying the border post-processing (default on, as RTIGER does).
# Returns the same nested list shape as `obs`, with 1-based state paths.
.rtiger_decode <- function(obs, params, r, postprocess = TRUE) {
  logPI <- log(params$pi); logA <- log(params$A)
  lapply(obs, function(sm) lapply(sm, function(ch) {
    lp <- rtiger_getlogpsi_cpp(ch$k, ch$n, params$alpha, params$beta)
    LP <- rtiger_productpsi_cpp(lp, r)
    v  <- rtiger_viterbi_cpp(logPI, LP, lp, logA, r)  # 1-based states
    if (postprocess) .rtiger_postprocess(v, lp) else v
  }))
}

# call_ancestry(caller="rtiger") backend: build (k,n) observations from the
# common input, fit + decode the faithful RTIGER, remap states
# [pat,het,mat]=(1,2,3) -> common [REF,HET,ALT]=(0,1,2) via state-1, and RLE to
# the common segment schema. RTIGER fits across all samples in `data` (the
# consumer groups per taxon). No post-processing yet (separate item).
.call_ancestry_rtiger <- function(data, rigidity, source, donor, has_donor, threads, seed,
                                  postprocess = TRUE) {
  by_name <- split(data, data$name, drop = TRUE)
  obs <- list(); pos <- list(); donor_of <- character(0)
  for (nm in names(by_name)) {
    dn <- by_name[[nm]]
    donor_of[nm] <- if (has_donor) as.character(dn$donor[1]) else donor
    chr_l <- split(dn, dn$chr, drop = TRUE)
    obs[[nm]] <- list(); pos[[nm]] <- list()
    for (cc in names(chr_l)) {
      dc <- chr_l[[cc]]; dc <- dc[order(dc$pos), , drop = FALSE]
      obs[[nm]][[cc]] <- list(k = dc$n_ref, n = dc$n_ref + dc$n_alt)
      pos[[nm]][[cc]] <- dc$pos
    }
  }
  fit   <- .rtiger_fit(obs, rigidity, threads = threads, seed = seed)
  paths <- .rtiger_decode(obs, fit, rigidity, postprocess = postprocess)
  out <- list()
  for (nm in names(obs)) for (cc in names(obs[[nm]])) {
    macro <- paths[[nm]][[cc]] - 1L                      # 1/2/3 -> 0/1/2
    out[[length(out) + 1L]] <- .rle_segments(macro, pos[[nm]][[cc]], as.integer(cc),
                                             nm, source, donor_of[[nm]])
  }
  calls <- do.call(rbind, out)
  calls[order(calls$donor, calls$name, calls$chr, calls$start_bp), , drop = FALSE]
}
