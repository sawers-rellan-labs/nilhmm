test_that("load_map returns the bundled B73 v5 consensus map", {
  m <- load_map()
  expect_setequal(names(m), c("locus", "chr", "cm", "bp"))
  expect_equal(sort(unique(m$chr)), 1:10)
  expect_identical(attr(m, "assembly"), "B73v5")
  expect_error(load_map("v99"), "only the bundled")
})

test_that("build_marker_grid samples the map to ~n_markers, cM interpolated", {
  g <- build_marker_grid(load_map(), n_markers = 500)
  expect_setequal(names(g), c("chr", "pos", "cm"))
  expect_equal(sort(unique(g$chr)), 1:10)
  expect_lt(abs(nrow(g) - 500), 60)                     # approximately the target
  expect_false(is.unsorted(g$pos[g$chr == 1L]))         # sorted within a chromosome
  expect_error(build_marker_grid(load_map()[, c("chr", "cm")]), "needs columns")
})

test_that("simulate_nil generates design-consistent truth on the bundled map", {
  skip_if_not_installed("simcross")
  tr <- simulate_nil("BC1S4", n = 4, chr = 1:2, n_markers = 200, seed = 1)
  expect_setequal(names(tr), c("source", "donor", "name", "chr", "pos", "cm", "state"))
  expect_setequal(unique(tr$chr), 1:2)                  # honoured the chr subset
  expect_equal(length(unique(tr$name)), 4L)
  expect_true(all(tr$state %in% 0:2))
  expect_identical(unique(tr$donor), "B")

  fr <- tabulate(tr$state + 1L, 3) / nrow(tr)           # BC1S4: recurrent-dominated, low het
  expect_gt(fr[1], 0.5)
  expect_lt(fr[2], 0.15)
  expect_gt(fr[3], 0)
  # deterministic under a pinned seed
  expect_identical(tr, simulate_nil("BC1S4", n = 4, chr = 1:2, n_markers = 200, seed = 1))
  # truth is directly segmentable
  expect_true(all(c("start_bp", "end_bp", "state") %in% names(to_segments(tr))))
})

test_that("simulate_counts degrades truth to counts + a masked genotype", {
  skip_if_not_installed("simcross")
  tr <- simulate_nil("BC2S2", n = 2, chr = 1, n_markers = 150, seed = 1)
  obs <- simulate_counts(tr, depth = 5, p_missing = 0.1, seed = 1)
  expect_setequal(names(obs), c("name", "donor", "chr", "pos", "n_ref", "n_alt", "g"))
  expect_true(all(obs$g %in% 0:3))
  expect_true(any(obs$g == 3L))                          # some missing
  expect_true(all(obs$n_ref >= 0L & obs$n_alt >= 0L))
  # missing <=> zero total depth
  expect_true(all((obs$n_ref + obs$n_alt == 0L) == (obs$g == 3L)))
  expect_error(simulate_counts(tr[, c("name", "chr")]), "needs columns")
})
