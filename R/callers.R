# Callers = named (emission x duration) combinations (S2, S8). "nilHMM" is the
# PACKAGE; the callers are explicit methods inside it. This kills the
# package-vs-caller ambiguity.

#' Resolve a named caller into emission + duration specs
#'
#' - `nnil`    : Holland's nNIL -- count/gt emission + geometric duration.
#' - `rtiger`  : rigidity mode -- count emission + rigidity duration (S7).
#' - `atlas`   : GOOGA competitive-alignment caller -- gt (categorical) emission
#'   + geometric duration; the per-unit call uses GOOGA fraction thresholds
#'   (applied in [call_ancestry()]).
#'
#' @param caller One of `"nnil"`, `"rtiger"`, `"atlas"`.
#' @param r Duration hyperparameter: geometric self-transition rate for
#'   `nnil`/`atlas`; integer rigidity (minimum run length) for `rtiger`.
#' @param err Count-emission baseline error.
#' @param conc Count-emission BetaBinomial concentration.
#' @param fit_means EM-fit emission means (count emission; S10).
#' @param p_switch Free-state switch probability for the `rtiger` rigidity tail.
#' @param germ,gert,p,mr,nir Genotype-error rates for the gt emission (the `atlas`
#'   caller; Holland's nNIL error model): hom error, het error, hom-error->het
#'   fraction, missing rate, non-informative-marker rate.
#' @param ... Ignored extra args (e.g. `f_1`/`f_2` consumed by [call_ancestry()]).
#' @return `list(emission, duration)`.
#' @examples
#' caller_spec("nnil", r = 1e-4)          # count emission + geometric duration
#' caller_spec("rtiger", r = 5)           # count emission + rigidity duration
#' str(caller_spec("atlas"))              # gt emission + geometric duration
#' @export
caller_spec <- function(caller = c("nnil", "rtiger", "atlas"),
                        r = 0.01, err = 0.01, conc = 20, fit_means = FALSE,
                        p_switch = 0.01, germ = 0.05, gert = 0.10, p = 0.5,
                        mr = 0.10, nir = 0.01, ...) {
  caller <- match.arg(caller)
  switch(caller,
    nnil    = list(emission = emission_count(err, conc, fit_means),
                   duration = duration_geometric(r)),
    rtiger  = list(emission = emission_count(err, conc, fit_means),
                   duration = duration_rigidity(r, p_switch)),
    atlas   = list(emission = emission_gt(germ, gert, p, mr, nir),
                   duration = duration_geometric(r))
  )
}
