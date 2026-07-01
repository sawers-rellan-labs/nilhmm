# Callers = named (emission x duration) combinations (S2, S8). "nilHMM" is the
# PACKAGE; the callers are explicit methods inside it. This kills the
# package-vs-caller ambiguity.

#' Resolve a named caller into emission + duration specs
#'
#' - `nnil`    : Holland's nNIL -- count/gt emission + geometric duration.
#' - `rtiger`  : rigidity mode -- count emission + rigidity duration (S7).
#'
#' @param caller One of `"nnil"`, `"rtiger"`.
#' @param r Duration hyperparameter: geometric self-transition rate for
#'   `nnil`; integer rigidity (minimum run length) for `rtiger`.
#' @param err Count-emission baseline error.
#' @param conc Count-emission BetaBinomial concentration.
#' @param fit_means EM-fit emission means (count emission; S10).
#' @param p_switch Free-state switch probability for the `rtiger` rigidity tail.
#' @param ... Ignored extra args (e.g. `f_1`/`f_2` consumed by [call_ancestry()]).
#' @return `list(emission, duration)`.
#' @export
caller_spec <- function(caller = c("nnil", "rtiger"),
                        r = 0.01, err = 0.01, conc = 20, fit_means = FALSE,
                        p_switch = 0.01, ...) {
  caller <- match.arg(caller)
  switch(caller,
    nnil    = list(emission = emission_count(err, conc, fit_means),
                   duration = duration_geometric(r)),
    rtiger  = list(emission = emission_count(err, conc, fit_means),
                   duration = duration_rigidity(r, p_switch))
  )
}
