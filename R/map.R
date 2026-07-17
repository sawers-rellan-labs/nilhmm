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

# Per-chromosome monotone Hyman spline mapping column `from` -> column `to` of a
# consensus map (a Marey map), clamped to each chromosome's observed `from`
# range; duplicate `from` values collapsed by mean. Returns a vectorized
# function(chr, x) -> interpolated `to`. The one place the bp<->cM interpolation
# lives; bp_to_cm/cm_to_bp/cm_to_mb and build_marker_grid all route through it.
.marey <- function(map, from, to) {
  map <- as.data.frame(map, stringsAsFactors = FALSE)
  if (!all(c("chr", from, to) %in% names(map)))
    stop(".marey(): map needs columns chr, ", from, ", ", to)
  fns <- lapply(split(map, map$chr), function(d) {
    y <- tapply(d[[to]], d[[from]], mean)          # collapse duplicate `from`
    x <- as.numeric(names(y)); o <- order(x)
    s <- stats::splinefun(x[o], as.numeric(y)[o], method = "hyman")
    rng <- range(x)
    function(v) s(pmin(pmax(v, rng[1]), rng[2]))   # clamp to observed range
  })
  function(chr, x) {
    x <- as.numeric(x); chr <- rep(as.character(chr), length.out = length(x))
    out <- numeric(length(x))
    for (ch in unique(chr)) {
      if (is.null(fns[[ch]])) stop("map interpolation: chromosome not in map: ", ch)
      out[chr == ch] <- fns[[ch]](x[chr == ch])
    }
    out
  }
}

#' Map interpolators: physical <-> genetic position
#'
#' Build a per-chromosome monotone Hyman spline (a Marey map) through a consensus
#' map's `(bp, cm)` pairs, clamped to each chromosome's observed range, and return
#' a vectorized `function(chr, x)`. Duplicate coordinates are collapsed by mean.
#' The default map is the bundled B73 v5 consensus map; pass any data frame with
#' `chr`, `bp`, `cm` columns to use your own (e.g. a population's native map).
#'
#' @param map A consensus map with columns `chr`, `bp`, `cm`. Defaults to
#'   [load_map()] (bundled B73 v5).
#' @return `bp_to_cm`: `function(chr, bp) -> cm`; `cm_to_bp`: `function(chr, cm) -> bp`.
#' @examples
#' to_cm <- bp_to_cm()          # bundled map
#' to_cm(1L, 1e6)
#' @seealso [cm_to_mb()], [load_map()]
#' @name map_interpolators
#' @export
bp_to_cm <- function(map = load_map()) .marey(map, "bp", "cm")

#' @rdname map_interpolators
#' @export
cm_to_bp <- function(map = load_map()) .marey(map, "cm", "bp")

#' Project segment cM coordinates to physical Mb
#'
#' Converts segment coordinates in genetic (cM) space to physical Mb via the
#' inverse Marey spline [cm_to_bp()]. cM coordinates are assembly-robust; the
#' bp/Mb they map to are tied to the map's assembly (bundled = B73 v5).
#'
#' @param seg Segments with columns `chr`, `start_cm`, `end_cm`.
#' @param map A consensus map (`chr`, `bp`, `cm`); defaults to [load_map()].
#' @return `seg` with `start_mb` and `end_mb` added.
#' @examples
#' cm_to_mb(data.frame(chr = 1L, start_cm = 0, end_cm = 10))
#' @seealso [cm_to_bp()], [bp_to_cm()]
#' @export
cm_to_mb <- function(seg, map = load_map()) {
  f <- cm_to_bp(map)
  seg$start_mb <- f(seg$chr, seg$start_cm) / 1e6
  seg$end_mb   <- f(seg$chr, seg$end_cm) / 1e6
  seg
}