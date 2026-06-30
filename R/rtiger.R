# RTIGER driver: the R-level orchestration of the Rcpp kernels in src/rtiger.cpp,
# faithfully reproducing the fork's EM / fit / viterbi (rHMM_methods.jl L1040-1622).
# This is the real RTIGER reproduction; it supersedes the earlier invented
# duration_rigidity/p_switch rtiger path.
#
# Observation convention (generateObject): per sample x chromosome, a 2-row
# matrix [k; n] with k = ref/P1-allele count, n = total. States [pat,het,mat] =
# (ref~0.95, 0.5, ref~0.05); remapped to the common schema [REF=0,HET=1,ALT=2]
# as common = rtiger_state - 1 (pat->REF, het->HET, mat->ALT).

# --- emission M-step: per-state BetaBinomial update (Brent on the concentration)
# Port of emissionUpdateState (rHMM_methods.jl L459-517). Brent -> stats::optimize.
.rtiger_emission_update_state <- function(i, ks, ns, ws, sumk, sumn, alpha_old, beta_old) {
  mi <- sumk / sumn
  if (is.nan(mi)) {                                   # state never supported
    mi <- alpha_old[i] / (alpha_old[i] + beta_old[i])
    return(c(a = alpha_old[i], b = beta_old[i]))
  }
  if (mi < 0.01) mi <- 0.01
  if (mi > 0.99) mi <- 0.99
  tau_i <- alpha_old[i] / mi + beta_old[i] / (1 - mi)
  if (tau_i > 100) tau_i <- 100
  Q <- function(t) {
    a <- max(t * mi, 1e-6); b <- max(t * (1 - mi), 1e-6)
    sum(ws * (lchoose(ns, ks) + lbeta(ks + a, ns - ks + b) - lbeta(a, b)))
  }
  tau <- optimize(Q, lower = max(1e-6, tau_i - 100), upper = max(tau_i + 1, 100),
                  maximum = TRUE)$maximum
  c(a = tau * mi, b = tau * (1 - mi))
}

# --- E-step for one chain: returns zeta (T x s x s) and gamma (s x T).
.rtiger_estep_chain <- function(O_k, O_n, logPI, logA, alpha, beta, r) {
  lp  <- rtiger_getlogpsi_cpp(O_k, O_n, alpha, beta)
  LP  <- rtiger_productpsi_cpp(lp, r)
  al  <- rtiger_forward_cpp(logPI, LP, logA, lp, r)
  be  <- rtiger_backward_cpp(LP, logA, lp, r)
  z   <- rtiger_zeta_cpp(al, be, logA, LP, lp, r)
  gam <- rtiger_gamma_cpp(z, al, be, r)
  list(zeta = z, gamma = gam)
}

# --- one EM iteration over all chains; returns updated (A, pi, alpha, beta).
# The per-chain E-step and the pooled suff-stat fold run in C++
# (rtiger_em_suffstats_cpp); only the M-step (cheap) stays in R. Mirrors the
# fork's EM (rHMM_methods.jl L1057-1242): transition (sumZeta), start
# (startAcc/nOb), emission (per-state distinct-(k,n) weights + sumk/sumn).
.rtiger_em <- function(ks_list, ns_list, logPI, logA, alpha, beta, r, nstates = 3L) {
  s <- nstates
  ss <- rtiger_em_suffstats_cpp(ks_list, ns_list, logPI, logA, alpha, beta,
                                as.integer(r), as.integer(nstates))

  rs <- rowSums(ss$sumZeta); rs[rs == 0] <- 1           # zero-row guard
  Anew  <- ss$sumZeta / rs
  PInew <- ss$startAcc / ss$nOb

  anew <- numeric(s); bnew <- numeric(s)
  for (st in 1:s) {
    ab <- .rtiger_emission_update_state(st, ss$kvals, ss$nvals, ss$wmat[, st],
                                        ss$sumk[st], ss$sumn[st], alpha, beta)
    anew[st] <- ab["a"]; bnew[st] <- ab["b"]
  }
  list(A = Anew, pi = PInew, alpha = anew, beta = bnew)
}

# --- fit: EM until max(|Δα|,|Δβ|) <= eps or max.iter (fit, rHMM_methods.jl L1412-1582).
# Deterministic init (randomize off, for reproducibility): transition diag-dominant,
# emission α=[20,20,1] β=[1,20,20], pi uniform (the generate_params forms).
.rtiger_fit <- function(obs, r, nstates = 3L, eps = 0.01, max_iter = 50L) {
  s <- nstates
  # flatten chains once (the C++ fold consumes flat k/n lists)
  ks_list <- list(); ns_list <- list()
  for (sm in obs) for (ch in sm) {
    ks_list[[length(ks_list) + 1L]] <- as.integer(ch$k)
    ns_list[[length(ns_list) + 1L]] <- as.integer(ch$n)
  }

  A  <- matrix(0.1, s, s) + diag(s) * 10; A <- A / rowSums(A)
  PI <- rep(1 / s, s)
  alpha <- if (s == 3L) c(20, 20, 1)  else rep(20, s)
  beta  <- if (s == 3L) c(1, 20, 20)  else rep(20, s)
  iter <- 0L
  repeat {
    nw <- .rtiger_em(ks_list, ns_list, log(PI), log(A), alpha, beta, r, nstates)
    er <- max(abs(alpha - nw$alpha), abs(beta - nw$beta))
    A <- nw$A; PI <- nw$pi; alpha <- nw$alpha; beta <- nw$beta
    iter <- iter + 1L
    if (round(er, 6) <= eps || iter >= max_iter) break
  }
  list(A = A, pi = PI, alpha = alpha, beta = beta, iterations = iter)
}

# --- decode all chains with the fitted parameters (viterbi); returns 1-based paths.
.rtiger_decode <- function(obs, params, r) {
  logPI <- log(params$pi); logA <- log(params$A)
  lapply(obs, function(sm) lapply(sm, function(ch) {
    lp <- rtiger_getlogpsi_cpp(ch$k, ch$n, params$alpha, params$beta)
    LP <- rtiger_productpsi_cpp(lp, r)
    rtiger_viterbi_cpp(logPI, LP, lp, logA, r)        # 1-based states
  }))
}
