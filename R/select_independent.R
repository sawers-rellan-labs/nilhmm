# Maximal-independent-set marker thinning (FastIndep port). Selects a large set
# of markers no two of which are more related than a threshold -- LD-based
# pruning for joint-linkage QTL mapping. Data-agnostic and per-chromosome: the
# caller supplies genotypes (or a precomputed relatedness matrix) for one
# chromosome at a time. The C++ core is fast_indep_cpp(); the matrix builder is
# pairwise_distance().

#' Select a maximal independent set of markers (LD thinning)
#'
#' Choose a large set of markers in which no two are more related than
#' `threshold` -- a maximal independent set on the relatedness graph, via a
#' faithful port of FastIndep (deterministic greedy heuristic + optional
#' stochastic runs). Intended for pruning a densified genotype matrix per
#' chromosome before joint-linkage mapping.
#'
#' @details
#' `x` is EITHER a genotype matrix (rows = markers, cols = samples) OR a
#' precomputed square relatedness matrix. It is treated as precomputed when it
#' carries an `attr(, "kind")` (as returned by [pairwise_distance()]) or is a
#' square symmetric matrix; otherwise it is treated as genotypes and
#' `pairwise_distance(x, method)` is computed. To force one interpretation, pass
#' the matrix you want.
#'
#' **Edge / threshold sense** (a set is independent iff no pair is an edge):
#' \describe{
#'   \item{`r2` (similarity)}{edge if `r2 >= threshold`. `t ~ 0.2` = strict
#'     GWAS-QC pruning; **`t ~ 0.5` is recommended for JLM** (keeps resolution,
#'     roughly VIF 2 / PLINK `--indep`).}
#'   \item{`mi` (similarity)}{edge if `MI >= threshold` (nats or bits). For
#'     Gaussian-equivalent reasoning in r2 units, `t_MI = -0.5 * ln(1 - r2)`
#'     nats, so `r2 = 0.5 <=> ~0.347 nats` -- given as a note; the code does no
#'     conversion.}
#'   \item{`vi` (distance, a metric)}{edge if `VI <= threshold` (small VI =
#'     related).}
#' }
#' The sense is taken from `attr(x, "kind")` when present, else from `method`
#' (`vi` is a distance; `r2`/`mi` are similarities). Pass `sense` to override
#' this inference explicitly -- e.g. `sense = "distance"` on a precomputed
#' coordinate/cM matrix (see [position_distance()]) prunes it as a distance
#' without abusing `method = "vi"`.
#'
#' @section RNG fidelity:
#' The greedy set (`n_runs = 1`, or the first reported set) is deterministic and
#' bit-identical to the FastIndep CLI. The stochastic runs (`n_runs > 1`) are
#' algorithmically faithful and reproducible from `seed` (self-contained PRNG),
#' but individual random sets are not bit-identical to the CLI's Mersenne-Twister
#' stream.
#'
#' @param x Genotype matrix (markers x samples) or a precomputed square
#'   relatedness matrix (see Details).
#' @param threshold Relatedness cutoff defining an edge (sense per Details).
#' @param n_runs Total runs: `1` = greedy only; `> 1` adds `n_runs - 1` stochastic
#'   runs and records the distinct sets found.
#' @param seed Seed for the stochastic runs.
#' @param method Measure for the genotype-matrix path and for the edge sense when
#'   `x` has no `attr(, "kind")`: one of `"r2"` (default), `"mi"`, `"vi"`.
#' @param sense Edge/threshold sense: `"auto"` (default) infers it from
#'   `attr(x, "kind")` then `method` (current behaviour); `"similarity"` (edge if
#'   `>= threshold`) or `"distance"` (edge if `<= threshold`) force it, so any
#'   precomputed matrix -- e.g. a cM distance from [position_distance()] -- can be
#'   pruned with the correct sense regardless of `method` or attributes.
#' @param max_markers Safety guard on the marker count (default `7000L`, the
#'   validated scale). Selection builds/consumes an O(n^2) markers x markers
#'   matrix, so more than `max_markers` markers errors cleanly before any large
#'   allocation (never a `malloc` failure); large-but-allowed sizes warn. Raise it
#'   (up to `getOption("nilHMM.marker_hard_cap")`, default `30000`) only with the
#'   memory to spare -- the intended workflow prunes per chromosome.
#' @param ... Passed to [pairwise_distance()] on the genotype path (e.g. `base`).
#' @return Character vector of the selected marker names (the largest independent
#'   set), or 1-based integer indices when `x` has no names. Attributes:
#'   `attr(, "sets")` (list of all distinct sets, as names/indices, greedy first),
#'   `attr(, "size_dist")` (named integer vector, set size -> count),
#'   `attr(, "kind")` the sense used.
#' @seealso [pairwise_distance()], [fast_indep_cpp()]
#' @examples
#' set.seed(1)
#' geno <- matrix(sample(0:2, 8 * 30, replace = TRUE), nrow = 8,
#'                dimnames = list(paste0("m", 1:8), NULL))
#' select_independent(geno, threshold = 0.5, method = "r2")
#' @export
select_independent <- function(x, threshold, n_runs = 1L, seed = 1L,
                               method = c("r2", "mi", "vi"),
                               sense = c("auto", "similarity", "distance"),
                               max_markers = 7000L, ...) {
  method <- match.arg(method)
  sense  <- match.arg(sense)
  if (missing(threshold) || !is.numeric(threshold) || length(threshold) != 1L)
    stop("select_independent(): `threshold` must be a single number.")

  x <- as.matrix(x)
  if (!is.numeric(x)) storage.mode(x) <- "double"

  # Decide whether `x` is a precomputed relatedness matrix or raw genotypes.
  # attr("kind") (from pairwise_distance) is authoritative; otherwise a square
  # symmetric numeric matrix is taken as precomputed, else it is genotypes.
  kind <- attr(x, "kind", exact = TRUE)
  is_precomputed <- !is.null(kind) ||
    (nrow(x) == ncol(x) && nrow(x) > 0L && isSymmetric(unname(x)))

  # Guard the marker count before any O(n^2) allocation. For genotypes the count
  # is the rows; for a precomputed matrix it is the (square) dimension. The
  # genotype path re-checks inside pairwise_distance() with the same limit.
  n_markers <- nrow(x)
  .check_marker_budget(n_markers, max_markers, "select_independent")

  if (is_precomputed) {
    sim <- x
    if (is.null(kind))
      kind <- if (method == "vi") "distance" else "similarity"
  } else {
    sim  <- pairwise_distance(x, method = method, max_markers = max_markers, ...)
    kind <- attr(sim, "kind", exact = TRUE)
  }

  if (nrow(sim) != ncol(sim))
    stop("select_independent(): relatedness matrix must be square.")

  # `sense` overrides the inferred `kind` when set; "auto" keeps the inference.
  if (sense != "auto") kind <- sense
  distance <- identical(kind, "distance")
  res <- fast_indep_cpp(sim, as.numeric(threshold), as.integer(n_runs),
                        as.integer(seed), distance)

  # Map 1-based indices to marker names when the matrix is named.
  nm <- rownames(sim)
  to_ids <- function(idx) if (is.null(nm)) idx else nm[idx]

  best <- to_ids(res$best)
  attr(best, "sets")      <- lapply(res$sets, to_ids)
  attr(best, "size_dist") <- res$size_dist
  attr(best, "kind")      <- kind
  best
}
