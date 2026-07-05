# position_distance(): the |pos_i - pos_j| distance-matrix builder, and its use
# through select_independent(..., sense = "distance") for coordinate/cM thinning.
# On a 1-D coordinate the maximal independent set (all pairs > t apart) has the
# same cardinality as a greedy min-spacing sweep, which is provably optimal on a
# line -- so we cross-check size and validity against that independent reference.

# Greedy left-to-right min-spacing sweep: keep a marker iff it is > t beyond the
# last kept one. Optimal-cardinality set of points pairwise > t apart on a line.
greedy_minspace <- function(pos, t) {
  o <- order(pos)
  p <- unname(pos)[o]
  keep <- logical(length(p))
  last <- -Inf
  for (i in seq_along(p)) {
    if (p[i] - last > t) { keep[i] <- TRUE; last <- p[i] }
  }
  nm <- names(pos)
  if (is.null(nm)) o[keep] else nm[o][keep]
}

test_that("position_distance returns |pos_i - pos_j| exactly", {
  pos <- c(0, 0.05, 0.2, 0.9, 1.5)
  d <- position_distance(pos)
  ref <- abs(outer(pos, pos, "-"))
  expect_equal(d, ref, ignore_attr = TRUE, tolerance = 1e-12)
  expect_identical(attr(d, "kind"), "distance")
  expect_identical(attr(d, "method"), "cm")
  expect_true(isSymmetric(unname(d)))
  expect_equal(diag(d), rep(0, length(pos)))
})

test_that("names(pos) are carried onto both dimensions", {
  pos <- c(m1 = 0, m2 = 0.3, m3 = 1.1)
  d <- position_distance(pos)
  expect_identical(rownames(d), c("m1", "m2", "m3"))
  expect_identical(colnames(d), c("m1", "m2", "m3"))
  expect_equal(d["m1", "m3"], 1.1, tolerance = 1e-12)
})

test_that("cross-chromosome pairs are Inf, within-chromosome finite", {
  pos <- c(0, 1, 5, 0, 2)
  chr <- c(1, 1, 1, 2, 2)
  d <- position_distance(pos, chr = chr)
  # within chr 1
  expect_equal(d[1, 3], 5)
  # within chr 2
  expect_equal(d[4, 5], 2)
  # across chromosomes -> Inf
  expect_true(is.infinite(d[1, 4]))
  expect_true(is.infinite(d[3, 5]))
  # diagonal still 0 (same chromosome as itself)
  expect_equal(diag(d), rep(0, 5))
  # symmetry preserved (Inf both ways)
  expect_true(isSymmetric(unname(d)))
})

test_that("method label is settable and does not change the arithmetic", {
  pos <- c(100, 250, 900)
  d <- position_distance(pos, method = "bp")
  expect_identical(attr(d, "method"), "bp")
  expect_equal(d[1, 3], 800)
})

test_that("input validation: NA, non-numeric, chr length mismatch, chr NA", {
  expect_error(position_distance(c(1, NA, 3)), "NA")
  expect_error(position_distance(letters[1:3]), "numeric")
  expect_error(position_distance(c(1, 2, 3), chr = c(1, 1)), "same length")
  expect_error(position_distance(c(1, 2, 3), chr = c(1, NA, 2)), "`chr` contains NA")
})

test_that("cM pruning: selected set is pairwise > t and matches greedy sweep size", {
  pos <- c(m1 = 0, m2 = 0.05, m3 = 0.2, m4 = 0.55, m5 = 0.6,
           m6 = 0.9, m7 = 1.5, m8 = 1.52, m9 = 2.4, m10 = 2.9)
  t <- 0.1
  d <- position_distance(pos)
  sel <- select_independent(d, threshold = t, n_runs = 200L, sense = "distance")

  # All selected pairs are strictly more than t cM apart.
  sub <- d[sel, sel, drop = FALSE]
  offdiag <- sub[upper.tri(sub)]
  expect_true(all(offdiag > t))

  # Optimal cardinality: equals the greedy min-spacing reference (optimal on a line).
  ref <- greedy_minspace(pos, t)
  expect_equal(length(sel), length(ref))
  expect_identical(attr(sel, "kind"), "distance")
})

test_that("sense = 'distance' override works on a raw (attr-less) matrix", {
  pos <- c(0, 0.05, 0.2, 0.9, 1.0, 1.5, 2.4)
  t <- 0.1
  d <- position_distance(pos)
  raw <- unname(d)
  attr(raw, "kind") <- NULL
  attr(raw, "method") <- NULL          # a plain square symmetric numeric matrix

  sel <- select_independent(raw, threshold = t, sense = "distance")
  sub <- raw[sel, sel, drop = FALSE]
  expect_true(all(sub[upper.tri(sub)] > t))   # pruned as a distance
  expect_identical(attr(sel, "kind"), "distance")

  # Without the override, an attr-less square matrix defaults to the similarity
  # sense -- a different (wrong-for-distance) interpretation, proving the override.
  sel_auto <- select_independent(raw, threshold = t, sense = "auto")
  expect_identical(attr(sel_auto, "kind"), "similarity")
  expect_false(identical(sort(as.character(sel)), sort(as.character(sel_auto))))
})

test_that("sense = 'auto' preserves current inference (r2 regression)", {
  set.seed(6)
  geno <- matrix(sample(0:2, 10 * 40, replace = TRUE), nrow = 10,
                 dimnames = list(paste0("snp", 1:10), NULL))
  base_sel <- select_independent(geno, threshold = 0.5, n_runs = 5L, method = "r2")
  auto_sel <- select_independent(geno, threshold = 0.5, n_runs = 5L, method = "r2",
                                 sense = "auto")
  expect_identical(auto_sel, base_sel)
  expect_identical(attr(auto_sel, "kind"), "similarity")
})

test_that("sense override can force a vi distance matrix to be read as similarity", {
  # Not a recommended use, but proves the override is authoritative over attr(kind).
  set.seed(7)
  geno <- matrix(sample(0:2, 12 * 40, replace = TRUE), nrow = 12,
                 dimnames = list(paste0("m", 1:12), NULL))
  d_vi <- pairwise_distance(geno, "vi")             # attr(kind) = "distance"
  sel <- select_independent(d_vi, threshold = 0.5, sense = "similarity")
  expect_identical(attr(sel, "kind"), "similarity")
  esim <- d_vi >= 0.5; diag(esim) <- FALSE
  expect_false(any(esim[sel, sel, drop = FALSE]))   # independent in the forced sense
})
