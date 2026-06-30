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
#' @param ... Parameter overrides forwarded to the emission/duration
#'   constructors (`r`, `err`, `conc`, ...).
#' @return `list(emission, duration)`.
#' @export
caller_spec <- function(caller = c("nnil", "rtiger", "skimbin"), ...) {
  caller <- match.arg(caller)
  stop("nilHMM::caller_spec() not yet implemented (Task 4)")
}
