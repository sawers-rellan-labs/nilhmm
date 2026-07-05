# Physical/genetic distance matrix builder for FastIndep marker thinning. A
# sibling of pairwise_distance(): instead of a genotype-derived relatedness, it
# builds the markers x markers matrix of |pos_i - pos_j| from a coordinate vector
# (cM or bp). Feeding it to select_independent() with sense = "distance" and a
# threshold t prunes to a maximal set of markers pairwise > t apart -- Chen's cM
# thinning, expressed on the same maximal-independent-set machinery as r2/MI/VI.
# Data-agnostic: the caller supplies the coordinates and owns the per-chromosome
# split (deriving cM for a specific dataset -- liftover, consensus-map spline --
# stays in the consumer repo).

#' Pairwise marker distance from map/physical coordinates
#'
#' Build a symmetric markers x markers **distance** matrix of the absolute
#' coordinate difference `|pos_i - pos_j|` from a vector of marker positions
#' (genetic cM or physical bp). The result is handed to [select_independent()]
#' with `sense = "distance"`: thresholding at `t` links markers `<= t` apart, so
#' the selected independent set is a maximal set of markers pairwise `> t` apart
#' -- distance-based thinning on the same machinery as r2/MI/VI relatedness.
#'
#' @details
#' Distance sense: an edge exists when `|pos_i - pos_j| <= threshold`, so the
#' independent set has all pairwise gaps `> threshold`. For a cM map, `threshold`
#' is in cM (e.g. `0.1` ~ Chen's 0.1-cM thinning). On a one-dimensional
#' coordinate the maximal independent set coincides with a greedy min-spacing
#' sweep (greedy is optimal on a line); the matrix path exists to give a uniform
#' interface across r2/MI/VI/cM rather than to beat the sweep.
#'
#' If `chr` is supplied, cross-chromosome pairs are set to `Inf` (never within
#' any finite `threshold`, hence never linked), so a genome-wide matrix still
#' prunes correctly. Per-chromosome use is preferred and cheaper (the matrix is
#' O(n^2)); the `chr` argument is a convenience for a single pooled call.
#'
#' @param pos Numeric vector of marker coordinates (cM or bp). `names(pos)`, if
#'   present, are carried onto both dimensions of the result.
#' @param chr Optional chromosome label per marker (same length as `pos`; any
#'   type comparable with `!=`). When given, cross-chromosome pairs are `Inf`.
#' @param method Label recorded in `attr(, "method")` (default `"cm"`); use
#'   `"pos"` or `"bp"` for physical coordinates. Does not change the arithmetic.
#' @param max_markers Safety guard on the marker count (default `7000L`, the
#'   validated scale). The output is an O(n^2) markers x markers matrix, so a
#'   request above `max_markers` errors cleanly (before any large allocation)
#'   rather than risking a `malloc` failure; large-but-allowed sizes warn. Raise
#'   it (up to the ceiling in `getOption("nilHMM.marker_hard_cap")`, default
#'   `30000`) only if the memory is available; the intended workflow prunes per
#'   chromosome, well under the default.
#' @return Symmetric numeric matrix (markers x markers) of `|pos_i - pos_j|`,
#'   zero diagonal, with `attr(, "kind") = "distance"` and `attr(, "method")`
#'   the label. Dimnames from `names(pos)` when present. Cross-chromosome pairs
#'   are `Inf` when `chr` is supplied.
#' @seealso [select_independent()], [pairwise_distance()]
#' @examples
#' cm <- c(m1 = 0, m2 = 0.05, m3 = 0.2, m4 = 0.9)
#' d <- position_distance(cm)
#' attr(d, "kind")   # "distance"
#' # Prune to markers pairwise > 0.1 cM apart:
#' select_independent(d, threshold = 0.1, sense = "distance")
#' @export
position_distance <- function(pos, chr = NULL, method = "cm",
                              max_markers = 7000L) {
  if (!is.numeric(pos))
    stop("position_distance(): `pos` must be a numeric vector.")
  n <- length(pos)
  if (n < 1L) stop("position_distance(): `pos` has no markers.")
  if (!is.null(chr) && length(chr) != n)
    stop("position_distance(): `chr` must have the same length as `pos`.")
  if (anyNA(pos)) stop("position_distance(): `pos` contains NA.")
  .check_marker_budget(n, max_markers, "position_distance")

  out <- abs(outer(as.numeric(pos), as.numeric(pos), "-"))
  if (!is.null(chr)) {
    cross <- outer(chr, chr, "!=")
    out[cross] <- Inf
  }
  diag(out) <- 0                      # |pos_i - pos_i| = 0 (also on-chromosome)

  nm <- names(pos)
  if (!is.null(nm)) dimnames(out) <- list(nm, nm)
  attr(out, "kind")   <- "distance"
  attr(out, "method") <- method
  out
}
