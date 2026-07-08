# Map utilities (S6). The consensus map is version-namespaced (maize_map_v5)
# and overridable; warn on assembly mismatch. Powers position files (cM<->Mb)
# and optional map-aware (position-dependent) transitions (open item, S10).

#' Load the bundled consensus map
#'
#' Returns the version-namespaced bundled consensus map (marker -> chr, cM, bp).
#' Only `"v5"` (B73 v5, [maize_map_v5]) is bundled; the map is overridable
#' everywhere it is used (e.g. pass your own `map` to [simulate_nil()]).
#'
#' @param version Map version key (default `"v5"`).
#' @return The map data.frame `locus, chr, cm, bp` (with an `"assembly"` attr).
#' @seealso [maize_map_v5], [simulate_nil()], [build_marker_grid()]
#' @examples
#' map <- load_map("v5")
#' head(map)
#' @export
load_map <- function(version = "v5") {
  if (!identical(version, "v5"))
    stop("load_map(): only the bundled 'v5' consensus map is available; got '", version, "'")
  e <- new.env()
  utils::data("maize_map_v5", package = "nilHMM", envir = e)
  e[["maize_map_v5"]]
}