# Emission models -- the axis that distinguishes callers (S4, S5).
# Each constructor returns a spec: state-conditional emission log-likelihoods
# plus which parameters are EM-fittable. Depth drives the choice (S5):
# saturated (>=~20x) -> gt; intermediate (~1-20x) -> count; imputed -> dosage.

#' Count (BetaBinomial) emission
#'
#' Emission on `(n_ref, n_alt)` read counts. State means `theta`; depth-0 markers
#' emit flat. **Means are FITTABLE** (EM-fit / reference-bias-corrected), not the
#' fixed `c(err, 0.5, 1 - err)` -- the BRB run showed fixed `theta_ALT = 1 - err`
#' collapses ALT->HET under reference bias (S10, BRB_run_findings.md).
#'
#' @param err Baseline genotyping/sequencing error (initialises `theta`).
#' @param conc BetaBinomial concentration (overdispersion); near-no-op on BRB
#'   (S10) but retained for the regime axis.
#' @param fit_means If `TRUE`, EM-fit the state means; if `FALSE` (default) use
#'   the fixed `c(err, 0.5, 1 - err)` that reproduces the Python baseline.
#'   Reference-biased data (RNA / BRB) needs `TRUE` (S10).
#' @return An emission spec for [fit()].
#' @export
emission_count <- function(err = 0.01, conc = 20, fit_means = FALSE) {
  structure(list(type = "count", err = err, conc = conc, fit_means = fit_means),
            class = c("nilHMM_emission_count", "nilHMM_emission"))
}

#' Genotype (categorical) emission
#'
#' Categorical over `{0, 1, 2, missing}` with a genotype-error matrix. The cheap
#' equivalent of `count` at saturated depth (Holland's native GT path; the
#' MolBreeding regime). Hard call, so no `(n, k)` cost blow-up.
#'
#' @param germ Error rate on true homozygotes.
#' @param gert Error rate on true heterozygotes.
#' @param p Fraction of homozygous errors that call as heterozygous.
#' @param mr Missing-genotype rate.
#' @param nir Non-informative-marker rate.
#' @return An emission spec for [fit()].
#' @export
emission_gt <- function(germ = 0.05, gert = 0.10, p = 0.5, mr = 0.10, nir = 0.01) {
  structure(list(type = "gt", germ = germ, gert = gert, p = p, mr = mr, nir = nir),
            class = c("nilHMM_emission_gt", "nilHMM_emission"))
}

#' Dosage (Gaussian/Beta) emission
#'
#' Continuous emission centred at `0 / 1 / 2` with variance set by imputation
#' uncertainty (Skim-BIN style; imputed-dosage sources).
#'
#' @param sd_dosage Per-state dosage standard deviation (imputation uncertainty).
#' @return An emission spec for [fit()].
#' @export
emission_dosage <- function(sd_dosage = 0.25) {
  structure(list(type = "dosage", sd_dosage = sd_dosage),
            class = c("nilHMM_emission_dosage", "nilHMM_emission"))
}

# --- internal emission interface (consumed by fit()/decode()) ----------------

# Initial / fixed REF/HET/ALT expected alt-fractions for an emission.
.emission_theta <- function(emission) {
  if (inherits(emission, "nilHMM_emission_count"))
    return(c(emission$err, 0.5, 1 - emission$err))
  stop(".emission_theta(): only the count emission is implemented (Task 4)")
}

# T x K log-emission matrix for a per-(sample, chromosome) observation table.
.emission_loglik <- function(emission, obs, theta) {
  if (inherits(emission, "nilHMM_emission_count"))
    return(count_emission_loglik_cpp(as.integer(obs$n), as.integer(obs$a),
                                     theta, emission$conc))
  stop(".emission_loglik(): only the count emission is implemented (Task 4)")
}

# Baum-Welch EM for the count-emission state means (S10 fix for reference-biased
# data). E-step: forward-backward posteriors per chromosome sequence. M-step:
# theta_s = sum(gamma_s * a) / sum(gamma_s * n) (the posterior-weighted alt
# fraction; exact for a Binomial, the mean update for BetaBinomial at fixed
# conc). Pooled across `obs_list` (a sample's chromosomes) so theta_ALT is
# estimated from all donor evidence, not one chromosome. A state that is never
# visited (denominator ~ 0) keeps its current mean. Transitions are held fixed
# (only means are fit); rigidity sub-state posteriors are summed back to macro.
.em_fit_means <- function(obs_list, emission, td, theta, control = list()) {
  max_iter <- if (!is.null(control$max_iter)) control$max_iter else 100L
  tol      <- if (!is.null(control$tol)) control$tol else 1e-5
  eps      <- 1e-4
  ns <- td$n_sub
  for (iter in seq_len(max_iter)) {
    num <- numeric(3); den <- numeric(3)
    for (o in obs_list) {
      em <- .expand_emission(count_emission_loglik_cpp(as.integer(o$n), as.integer(o$a),
                                                       theta, emission$conc), ns)
      g <- forward_backward_cpp(td$log_start, td$log_trans, em)   # T x (3*ns)
      gm <- if (ns > 1L)
              sapply(0:2, function(m) rowSums(g[, (m * ns + 1L):(m * ns + ns), drop = FALSE]))
            else g
      num <- num + colSums(gm * o$a)
      den <- den + colSums(gm * o$n)
    }
    new_theta <- theta
    upd <- den > eps
    new_theta[upd] <- pmin(pmax(num[upd] / den[upd], eps), 1 - eps)
    delta <- max(abs(new_theta - theta))
    theta <- new_theta
    if (delta < tol) break
  }
  theta
}