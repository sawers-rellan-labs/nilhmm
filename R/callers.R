# Callers = named (emission x duration) combinations (§2, §8). "nilHMM" is the
# PACKAGE; the callers are explicit methods inside it. This kills the
# package-vs-caller ambiguity.

#' Resolve a named caller into emission + duration specs
#'
#' - `nnil`    : Holland's nNIL — count/gt emission + geometric duration.
#' - `rtiger`  : rigidity mode — count emission + rigidity duration (§7).
#' - `skimbin` : Skim-BIN — dosage/count emission + geometric duration.
#'
#' @param caller One of `"nnil"`, `"rtiger"`, `"skimbin"`.
#' @param r Duration hyperparameter (self-transition rate / rigidity).
#' @param err Count-emission baseline error.
#' @param conc Count-emission BetaBinomial concentration.
#' @param fit_means EM-fit emission means (count emission; §10).
#' @param ... Ignored extra args (e.g. `f_1`/`f_2` consumed by [call_ancestry()]).
#' @return `list(emission, duration)`.
#' @export
caller_spec <- function(caller = c("nnil", "rtiger", "skimbin"),
                        r = 0.01, err = 0.01, conc = 20, fit_means = FALSE, ...) {
  caller <- match.arg(caller)
  switch(caller,
    nnil    = list(emission = emission_count(err, conc, fit_means),
                   duration = duration_geometric(r)),
    rtiger  = list(emission = emission_count(err, conc, fit_means),
                   duration = duration_rigidity(r)),
    skimbin = list(emission = emission_dosage(),
                   duration = duration_geometric(r))
  )
}
