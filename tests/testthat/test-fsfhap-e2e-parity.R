# FSFHap port — frozen END-TO-END TASSEL parity gate (stages 1 + 3, no JVM).
#
# agent/fsfhap_e2e_parity.R runs TASSEL FSFHap once and freezes the raw input +
# TASSEL's imputed parent calls. This replays our full pipeline (call->g -> stage 1
# BC parent-calling -> stage 3 5-state EM + fillgaps) and asserts PER-CELL identity
# with TASSEL — so any regression in the segregation test, same-tag filter, recode,
# EM Viterbi, transition scaling, or gap-fill fails CI without needing TASSEL.

test_that("full pipeline reproduces TASSEL imputed parent calls per cell", {
  fixf <- testthat::test_path("..", "fixtures", "fsfhap_tassel", "e2e_fixture.rds")
  skip_if_not(file.exists(fixf), "e2e parity fixture absent (run agent/fsfhap_e2e_parity.R)")
  fx <- readRDS(fixf)

  s1 <- .fsfhap_call_parents_bc(fx$G0, fx$pos, max_missing = 1.0, min_gametes = 200L, min_r = 0.0)
  imp <- .fsfhap_impute(s1$G, s1$pos, phet = fx$phet, fill_gaps = TRUE)
  dimnames(imp) <- list(fx$taxa[s1$keep_taxa], fx$markers[s1$keep_sites])

  # align on common markers/taxa; compare where both non-missing (3 = N)
  Tt <- fx$tassel_parents
  cm <- intersect(colnames(imp), rownames(Tt))
  ct <- intersect(rownames(imp), colnames(Tt))
  expect_gt(length(cm), 100L)                    # kept-site set must be substantial
  O  <- imp[ct, cm, drop = FALSE]
  Tm <- t(Tt[cm, ct, drop = FALSE])              # -> taxa x markers
  both <- O != 3L & Tm != 3L
  expect_gt(sum(both), 10000L)
  expect_equal(mean(O[both] == Tm[both]), 1)     # bit-exact parent calls
})
