# Plotting (S8). Fragment-size distributions against the grey Null (the design's
# expected Gamma from expected_fragment_dist() -- calibration/Null only, S7),
# and chromosome-painting of called segments.

#' Plot called fragment sizes against the expected Null
#'
#' @param calls Segment calls (common schema).
#' @param design Design key for the grey Null (see [expected_fragment_dist()]).
#' @param space One of `"cM"`, `"Mb"`.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' # Planned (Task 4): called fragment sizes vs the expected Null.
#' calls <- call_ancestry(read_counts("counts/"), caller = "nnil", design = "BC2S2")
#' plot_fragment_sizes(calls, design = "BC2S2", space = "Mb")
#' }
#' @export
plot_fragment_sizes <- function(calls, design, space = c("Mb", "cM")) {
  space <- match.arg(space)
  stop("nilHMM::plot_fragment_sizes() not yet implemented (Task 4)")
}