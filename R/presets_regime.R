# Presets layer 2a: emission-by-depth regime (S5). Convenience/defaults on top
# of the engine, NOT the engine. See docs/emission_by_depth_regime.

#' Select an emission model from sequencing depth
#'
#' Selector rule (S5): saturated (>=~20x) -> `gt`; intermediate (~1-20x) ->
#' `count`. Cost basis: BetaBinomial cost is proportional to the number of
#' distinct `(n, k)` pairs, i.e. to coverage; above ~20-30x a count emission is
#' effectively a hard call, so `gt` is equivalent and cheaper (memory
#' `rtiger-betabinomial-cost`).
#'
#' @param depth Mean per-marker sequencing depth.
#' @return An emission spec ([emission_count()] / [emission_gt()]).
#' @examples
#' class(select_emission(2))    # ~2x  -> count emission
#' class(select_emission(30))   # ~30x -> gt (hard-call) emission
#' @export
select_emission <- function(depth) {
  if (depth >= 20) return(emission_gt())    # saturated -> hard genotype call
  emission_count()                          # intermediate -> counts
}
