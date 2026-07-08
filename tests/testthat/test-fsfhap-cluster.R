# FSFHap port — stage 2a: HaplotypeClusterer.makeClusters + consensus (window).
# Faithful-port checks of the 0-distance clustering, multi-membership fractional
# scoring, the g-code distance metric, consensus, and score/size sort order.

int_mat <- function(x) { storage.mode(x) <- "integer"; x }

test_that("identical haplotypes collapse to one cluster (score == size == n)", {
  G <- int_mat(matrix(0L, 5, 6))                 # 5 taxa, all A-hom
  r <- fsfhap_cluster_window_cpp(G)
  expect_equal(length(r$size), 1L)
  expect_equal(r$size, 5L)
  expect_equal(r$score, 5)
  expect_equal(as.integer(r$majority[1, ]), rep(0L, 6))
  expect_equal(r$members[[1]], 1:5)
})

test_that("two distinct blocks -> two clusters, sorted by score desc", {
  G <- int_mat(rbind(matrix(0L, 4, 6), matrix(2L, 2, 6)))   # 4 A-hom, 2 C-hom
  r <- fsfhap_cluster_window_cpp(G)
  expect_equal(length(r$size), 2L)
  expect_equal(r$size, c(4L, 2L))                # larger/higher-score first
  expect_equal(r$score, c(4, 2))
  expect_equal(as.integer(r$majority[1, ]), rep(0L, 6))
  expect_equal(as.integer(r$majority[2, ]), rep(2L, 6))
})

test_that("an all-missing haplotype joins every cluster with fractional score", {
  G <- int_mat(rbind(matrix(0L, 3, 6), matrix(2L, 2, 6), matrix(3L, 1, 6)))
  r <- fsfhap_cluster_window_cpp(G)
  expect_equal(r$size, c(4L, 3L))                # taxon 6 in both clusters
  expect_equal(r$score, c(3.5, 2.5))             # 3 + .5  and  2 + .5
  expect_true(6L %in% r$members[[1]] && 6L %in% r$members[[2]])
})

test_that("het-vs-hom is distance 1 -> a lone het splits into its own cluster", {
  # 4 taxa A-hom; taxon 5 identical except a het at site 1 (distance 1 > 0)
  G <- rbind(matrix(0L, 4, 6), c(1L, rep(0L, 5))); G <- int_mat(G)
  r <- fsfhap_cluster_window_cpp(G)
  expect_equal(length(r$size), 2L)
  expect_equal(r$size, c(4L, 1L))
  expect_equal(r$members[[2]], 5L)
})

test_that("missing sites are tolerated within a cluster (0 distance)", {
  # taxon 2 identical to taxon 1 except N at two sites -> still 0 distance
  G <- int_mat(rbind(rep(2L, 6), c(2L, 3L, 2L, 3L, 2L, 2L), rep(2L, 6)))
  r <- fsfhap_cluster_window_cpp(G)
  expect_equal(length(r$size), 1L)
  expect_equal(r$size, 3L)
  expect_equal(as.integer(r$majority[1, ]), rep(2L, 6))   # N filled by the others
})

test_that("majority == unanimous for makeClusters (0-distance) output", {
  G <- int_mat(rbind(matrix(0L, 3, 5), matrix(2L, 3, 5)))
  r <- fsfhap_cluster_window_cpp(G)
  expect_equal(r$majority, r$unanimous)
})

# ---- clusterer completions: mergeClusters / moveToBiggest / removeHet --------

test_that("mergeClusters merges clusters within maxdiff, not beyond", {
  # block A (all 0) and A' (site1 = 2) differ by one hom site -> distance 2
  G <- int_mat(rbind(matrix(0L, 3, 5), cbind(2L, matrix(0L, 3, 4))))
  expect_equal(length(fsfhap_cluster_window_cpp(G, maxdiff = 0, merge = TRUE)$size), 2L)
  m2 <- fsfhap_cluster_window_cpp(G, maxdiff = 2, merge = TRUE)
  expect_equal(length(m2$size), 1L)
  expect_equal(m2$size, 6L)
})

test_that("moveAllHaplotypesToBiggestCluster absorbs near haplotypes at maxdiff", {
  G <- int_mat(rbind(matrix(0L, 3, 5), cbind(2L, matrix(0L, 3, 4))))   # dist 2 apart
  expect_equal(length(fsfhap_cluster_window_cpp(G, maxdiff = 0, move_biggest = TRUE)$size), 2L)
  expect_equal(length(fsfhap_cluster_window_cpp(G, maxdiff = 2, move_biggest = TRUE)$size), 1L)
})

test_that("removeHeterozygousClusters drops het-heavy clusters", {
  G <- int_mat(rbind(matrix(1L, 4, 20), matrix(0L, 4, 20)))   # 4 all-het taxa + 4 A-hom
  r <- fsfhap_cluster_window_cpp(G, max_het = 5L)             # 20 het sites > 5 -> drop
  expect_equal(length(r$size), 1L)
  expect_equal(as.integer(r$majority[1, ]), rep(0L, 20))      # the surviving cluster is A-hom
})

test_that("removeHet keeps clusters below the het threshold", {
  G <- int_mat(rbind(matrix(0L, 4, 20), matrix(2L, 4, 20)))   # two homozygous blocks, 0 het sites
  expect_equal(length(fsfhap_cluster_window_cpp(G, max_het = 5L)$size), 2L)
})
