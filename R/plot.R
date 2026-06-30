# Plotting (S8). Fragment-size distributions against the grey Null (the design's
# expected Gamma from expected_fragment_dist() -- calibration/Null only, S7),
# and chromosome-painting of called segments.

#' Plot called fragment sizes against the expected Null
#'
#' @param calls Segment calls (common schema).
#' @param design Design key for the grey Null (see [expected_fragment_dist()]).
#' @param space One of `"cM"`, `"Mb"`.
#' @return A ggplot object.
#' @export
plot_fragment_sizes <- function(calls, design, space = c("Mb", "cM")) {
  space <- match.arg(space)
  stop("nilHMM::plot_fragment_sizes() not yet implemented (Task 4)")
}