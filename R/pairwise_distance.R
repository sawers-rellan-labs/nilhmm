# Pairwise marker-relatedness matrix builder for LD-based marker thinning. Feeds
# select_independent(): compute a markers x markers matrix under a pluggable
# measure (r2 / MI / VI), then prune to a maximal independent set. Data-agnostic:
# the caller supplies the genotype matrix (e.g. the densified TeoNAM block) and
# owns the per-chromosome split. See select_independent() for the selector.

# Absolute safety ceiling on the marker count, independent of the caller's
# `max_markers`. Both the r2/MI/VI matrix (8*n^2 bytes) and FastIndep's separation
# sets (up to ~4*n^2 bytes) are O(n^2), so an unbounded n would malloc-fail. This
# ceiling (~7 GB for the n x n double matrix alone) refuses such a request in R,
# before any large allocation, rather than risking a crash. Per-chromosome pruning
# keeps n well under this. Overridable per session via option `nilHMM.marker_hard_cap`.
.marker_hard_cap <- function() {
  as.integer(getOption("nilHMM.marker_hard_cap", 30000L))
}

# Bytes an n x n double matrix would occupy (the dominant O(n^2) allocation).
.nn_bytes <- function(n) 8 * (as.numeric(n))^2

# Guard the marker count before any O(n^2) allocation. `max_markers` is the
# caller-facing, configurable limit (default 7000, the validated scale); it may be
# raised up to the hard cap when the memory is available. Errors cleanly (never a
# malloc failure) when exceeded, and warns for large-but-allowed sizes.
.check_marker_budget <- function(n, max_markers, fn = "select_independent") {
  max_markers <- as.integer(max_markers)
  hard <- .marker_hard_cap()
  gb <- function(x) .nn_bytes(x) / 2^30
  if (is.na(max_markers) || max_markers < 1L)
    stop(sprintf("%s(): `max_markers` must be a positive integer.", fn))
  if (max_markers > hard)
    stop(sprintf(paste0("%s(): `max_markers` (%d) exceeds the safety ceiling of %d ",
                        "markers (~%.1f GB for the n x n matrix). Raise ",
                        "getOption('nilHMM.marker_hard_cap') only if you truly have the memory."),
                 fn, max_markers, hard, gb(hard)))
  if (n > max_markers)
    stop(sprintf(paste0("%s(): %d markers exceeds `max_markers` (%d). The relatedness ",
                        "matrix is O(n^2) (~%.1f GB at n=%d). Prune per chromosome ",
                        "(intended usage), or raise `max_markers` (up to %d) if the memory is available."),
                 fn, n, max_markers, gb(n), n, hard))
  if (.nn_bytes(n) > 2^30)   # > ~1 GB: large but allowed -> warn
    warning(sprintf(paste0("%s(): building an %d x %d matrix (~%.1f GB, O(n^2) time and ",
                           "memory). Consider per-chromosome pruning."),
                    fn, n, n, gb(n)), call. = FALSE)
  invisible(TRUE)
}

#' Pairwise marker relatedness (r2, mutual information, or variation of information)
#'
#' Compute a symmetric markers x markers relatedness matrix from a
#' markers x samples dosage matrix, under a pluggable measure. The measure is
#' computed over the samples where both markers are non-missing (any `NA`/`NaN`).
#' The result is handed to [select_independent()] for maximal-independent-set
#' marker thinning; the measure is thresholded in its own units (no r2<->MI
#' conversion).
#'
#' @details
#' \describe{
#'   \item{`r2` (default, a similarity)}{`cor(g_i, g_j)^2`, the squared Pearson
#'     correlation of the dosage rows -- the PLINK `--indep-pairwise` convention.
#'     High = related. A constant (zero-variance) marker contributes `0`.}
#'   \item{`mi` (a similarity)}{plug-in mutual information from the joint genotype
#'     contingency table, with the Miller-Madow bias correction
#'     `I_MM = I_plugin + (m_X + m_Y - m_XY - 1) / (2N)`. High = related.}
#'   \item{`vi` (a distance, a true metric)}{variation of information
#'     `VI(X,Y) = H(X) + H(Y) - 2*I(X,Y)`, on MM-corrected entropies. LOW =
#'     related; `VI` is a metric (obeys the triangle inequality), unlike r2/MI.}
#' }
#' For near-linearly-dependent genotypes (e.g. RIL ancestry dosages) `r2` and `mi`
#' rank pairs almost identically; `mi`/`vi` earn their keep on non-monotone or
#' multiallelic data and for the metric/embedding direction.
#'
#' @param geno Numeric matrix, rows = markers, cols = samples; dosages (typically
#'   `0`/`1`/`2`). Missing entries are `NA`/`NaN`. Row names (marker ids) are
#'   carried onto both dimensions of the result.
#' @param method One of `"r2"` (default), `"mi"`, `"vi"`.
#' @param base For `mi`/`vi`: information unit, `"nats"` (default) or `"bits"`.
#'   Ignored for `r2`.
#' @param max_markers Safety guard on the marker count (default `7000L`, the
#'   validated scale). The output is an O(n^2) markers x markers matrix, so a
#'   request above `max_markers` errors cleanly (before any large allocation)
#'   rather than risking a `malloc` failure; large-but-allowed sizes warn. Raise
#'   it (up to the ceiling in `getOption("nilHMM.marker_hard_cap")`, default
#'   `30000`) only if the memory is available; the intended workflow prunes per
#'   chromosome, well under the default.
#' @return Symmetric numeric matrix (markers x markers) with `attr(, "kind")`
#'   `"similarity"` (r2, mi) or `"distance"` (vi) and `attr(, "method")` the
#'   measure. Dimnames from `rownames(geno)`.
#' @seealso [select_independent()]
#' @examples
#' set.seed(1)
#' geno <- matrix(sample(0:2, 5 * 20, replace = TRUE), nrow = 5,
#'                dimnames = list(paste0("m", 1:5), NULL))
#' d <- pairwise_distance(geno, "r2")
#' attr(d, "kind")   # "similarity"
#' @export
pairwise_distance <- function(geno, method = c("r2", "mi", "vi"),
                              base = c("nats", "bits"), max_markers = 7000L) {
  method <- match.arg(method)
  base   <- match.arg(base)
  method_int <- match(method, c("r2", "mi", "vi")) - 1L      # 0/1/2
  base_int   <- match(base, c("nats", "bits")) - 1L          # 0/1

  geno <- as.matrix(geno)
  if (!is.numeric(geno)) storage.mode(geno) <- "double"
  if (nrow(geno) < 1L) stop("pairwise_distance(): `geno` has no markers (rows).")
  .check_marker_budget(nrow(geno), max_markers, "pairwise_distance")

  out <- pairwise_distance_cpp(geno, method_int, base_int)
  rn <- rownames(geno)
  if (!is.null(rn)) dimnames(out) <- list(rn, rn)
  attr(out, "kind")   <- if (method == "vi") "distance" else "similarity"
  attr(out, "method") <- method
  out
}
