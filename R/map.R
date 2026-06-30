# Map utilities (S6). The consensus map is version-namespaced (maize_map_v5)
# and overridable; warn on assembly mismatch. Powers position files (cM<->Mb)
# and optional map-aware (position-dependent) transitions (open item, S10).

#' Load the bundled consensus map
#'
#' @param version Map version key (default `"v5"`); version-namespaced and
#'   overridable.
#' @return A map data.frame (marker -> chr, cM, bp).
#' @export
load_map <- function(version = "v5") {
  stop("nilHMM::load_map() not yet implemented (Task 4)")
}