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
.rtiger_fit <- function(obs, r, nstates = 3L, eps = 0.01, max_iter = 50L, threads = 1L) {
  ch <- .rtiger_chains(obs)
  rtiger_fit_cpp(ch$ks, ch$ns, as.integer(r), as.integer(nstates), eps,
                 as.integer(max_iter), as.integer(threads))
}

# Decode all chains with the fitted parameters (one Viterbi pass per chain).
# Returns the same nested list shape as `obs`, with 1-based state paths.
.rtiger_decode <- function(obs, params, r) {
  logPI <- log(params$pi); logA <- log(params$A)
  lapply(obs, function(sm) lapply(sm, function(ch) {
    lp <- rtiger_getlogpsi_cpp(ch$k, ch$n, params$alpha, params$beta)
    LP <- rtiger_productpsi_cpp(lp, r)
    rtiger_viterbi_cpp(logPI, LP, lp, logA, r)        # 1-based states
  }))
}
