# Emission models — the axis that distinguishes callers (§4, §5).
# Each constructor returns a spec: state-conditional emission log-likelihoods
# plus which parameters are EM-fittable. Depth drives the choice (§5):
# saturated (>=~20x) -> gt; intermediate (~1-20x) -> count; imputed -> dosage.

#' Count (BetaBinomial) emission
#'
#' Emission on `(n_ref, n_alt)` read counts. State means `theta`; depth-0 markers
#' emit flat. **Means are FITTABLE** (EM-fit / reference-bias-corrected), not the
#' fixed `c(err, 0.5, 1 - err)` — the BRB run showed fixed `theta_ALT = 1 - err`
#' collapses ALT->HET under reference bias (§10, BRB_run_findings.md).
#'
#' @param err Baseline genotyping/sequencing error (initialises `theta`).
#' @param conc BetaBinomial concentration (overdispersion); near-no-op on BRB
#'   (§10) but retained for the regime axis.
#' @param fit_means If `TRUE`, EM-fit the state means; if `FALSE` (default) use
#'   the fixed `c(err, 0.5, 1 - err)` that reproduces the Python baseline.
#'   Reference-biased data (RNA / BRB) needs `TRUE` (§10).
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

# Baum-Welch EM for the emission means (§10). Implemented in the next Task-4 step.
.em_fit_means <- function(obs, emission, tr, theta, control) {
  stop(".em_fit_means(): fittable emission means not yet implemented; ",
       "use fit_means = FALSE (the baseline-reproducing fixed means) for now")
}