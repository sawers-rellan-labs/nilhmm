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
# Pooled streaming sufficient statistics, mirroring EM (rHMM_methods.jl L1057-1242):
# transitionMultiple (sumZeta), startMultiple (startAcc/nOb), emission (per-state
# distinct-(k,n) weights W + weighted totals sumk/sumn).
.rtiger_em <- function(obs, logPI, logA, alpha, beta, r, nstates = 3L) {
  s <- nstates
  sumZeta <- matrix(0, s, s)
  startAcc <- numeric(s)
  nOb <- 0L
  W    <- replicate(s, new.env(parent = emptyenv()), simplify = FALSE)  # (k,n) -> Σγ
  sumk <- numeric(s); sumn <- numeric(s)

  for (sm in obs) for (ch in sm) {
    O_k <- ch$k; O_n <- ch$n
    Tc <- length(O_k)
    e  <- .rtiger_estep_chain(O_k, O_n, logPI, logA, alpha, beta, r)
    z  <- e$zeta; gam <- e$gamma                       # z: T x s x s ; gam: s x T

    # transition: Σ over band (r+1):(T-r+1) of zeta[t,,]
    band <- (r + 1):(Tc - r + 1)
    sumZeta <- sumZeta + apply(z[band, , , drop = FALSE], c(2, 3), sum)
    # start: Σ gamma[,1]
    startAcc <- startAcc + gam[, 1]
    nOb <- nOb + 1L
    # emission: per state, group γ by distinct (k,n); accumulate weighted totals
    key <- paste(O_k, O_n, sep = "_")
    for (st in 1:s) {
      agg <- tapply(gam[st, ], key, sum)
      for (nm in names(agg)) {
        cur <- W[[st]][[nm]]
        W[[st]][[nm]] <- (if (is.null(cur)) 0 else cur) + agg[[nm]]
      }
      sumk[st] <- sumk[st] + sum(O_k * gam[st, ])
      sumn[st] <- sumn[st] + sum(O_n * gam[st, ])
    }
  }

  # M-step: transition (row-normalize, zero-row guard), start, emission.
  rs <- rowSums(sumZeta); rs[rs == 0] <- 1
  Anew <- sumZeta / rs
  PInew <- startAcc / nOb
  anew <- numeric(s); bnew <- numeric(s)
  for (st in 1:s) {
    keys <- ls(W[[st]])
    kn <- do.call(rbind, strsplit(keys, "_", fixed = TRUE))
    ks <- as.integer(kn[, 1]); ns <- as.integer(kn[, 2])
    ws <- vapply(keys, function(nm) W[[st]][[nm]], numeric(1))
    ab <- .rtiger_emission_update_state(st, ks, ns, ws, sumk[st], sumn[st], alpha, beta)
    anew[st] <- ab["a"]; bnew[st] <- ab["b"]
  }
  list(A = Anew, pi = PInew, alpha = anew, beta = bnew)
}

# --- fit: EM until max(|Δα|,|Δβ|) <= eps or max.iter (fit, rHMM_methods.jl L1412-1582).
# Deterministic init (randomize off, for reproducibility): transition diag-dominant,
# emission α=[20,20,1] β=[1,20,20], pi uniform (the generate_params forms).
.rtiger_fit <- function(obs, r, nstates = 3L, eps = 0.01, max_iter = 50L) {
  s <- nstates
  A  <- matrix(0.1, s, s) + diag(s) * 10; A <- A / rowSums(A)
  PI <- rep(1 / s, s)
  alpha <- if (s == 3L) c(20, 20, 1)  else rep(20, s)
  beta  <- if (s == 3L) c(1, 20, 20)  else rep(20, s)
  iter <- 0L
  repeat {
    nw <- .rtiger_em(obs, log(PI), log(A), alpha, beta, r, nstates)
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
