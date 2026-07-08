# FSFHap stage 1a/1b — frozen TASSEL parity gate (CI, no JVM).
#
# agent/fsfhap_tassel_parity.R runs TASSEL FSFHap once and freezes the fixture
# (genotype matrix + positions) alongside TASSEL's logged cardinalities. This
# test replays our port on that exact input and asserts our counts equal TASSEL's
# — so a regression in the segregation binomial or the same-tag R^2 filter fails
# CI without needing TASSEL installed. Regenerate the fixture with that harness.

test_that("stage 1a/1b counts match frozen TASSEL cardinalities", {
  fixf <- testthat::test_path("..", "fixtures", "fsfhap_tassel", "parity_fixture.rds")
  skip_if_not(file.exists(fixf), "TASSEL parity fixture not present (run agent/fsfhap_tassel_parity.R)")
  fx <- readRDS(fixf)

  res <- .fsfhap_call_parents_bc(fx$G, fx$pos, max_missing = 1.0, min_rsq = 0.8,
                                 min_gametes = 200L, min_r = 0.0)

  # polybits = segregating sites (stage 1a)
  expect_equal(res$n_seg, fx$tassel$polybits)
  # filteredBits = polybits AND same-tag (stage 1b); == called snps at minR=0
  expect_equal(res$n_kept_sites, fx$tassel$filteredBits)
  expect_equal(res$n_kept_sites, fx$tassel$called)
  # the fixture must actually exercise the same-tag filter, else 1b is untested
  expect_lt(fx$tassel$filteredBits, fx$tassel$polybits)
})
