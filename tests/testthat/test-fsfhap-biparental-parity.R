# FSFHap port — frozen TASSEL parity gate for the BIPARENTAL route (stage 2b, no JVM).
#
# agent/fsfhap_biparental_parity.R runs TASSEL FSFHap (BiparentalHaplotypeFinder, via
# an F2 pedigree at contribution 0.5) once and freezes the fixture + TASSEL's imputed
# parent calls. This replays our biparental pipeline (preFilterSites -> reconstruct
# parental haplotypes -> recode -> 5-state EM impute) and asserts near-exact per-cell
# concordance -- so a regression in the seed/extend/assign/recode/impute chain fails
# CI without TASSEL. (LD filter off on both sides for this fixture, matching -minR 0.)

test_that("biparental route reproduces TASSEL imputed parent calls (per cell)", {
  fixf <- testthat::test_path("..", "fixtures", "fsfhap_tassel", "biparental_fixture.rds")
  skip_if_not(file.exists(fixf), "biparental parity fixture absent (run agent/fsfhap_biparental_parity.R)")
  fx <- readRDS(fixf)

  bip <- .fsfhap_biparental_call(fx$G, fx$pos, min_r2 = 0.0, min_gametes = 200L)
  expect_gt(length(bip$keep_sites), 100L)                       # seeded + reconstructed
  imp <- .fsfhap_impute(bip$G, bip$pos, phet = 0.5, fill_gaps = TRUE)
  dimnames(imp) <- list(fx$taxa[bip$keep_taxa], fx$markers[bip$keep_sites])

  Tt <- fx$tassel_parents
  # site-set match: our kept markers == TASSEL's (incl. the seed window)
  expect_setequal(fx$markers[bip$keep_sites], rownames(Tt))
  cm <- intersect(colnames(imp), rownames(Tt)); ct <- intersect(rownames(imp), colnames(Tt))
  O <- imp[ct, cm, drop = FALSE]; Tm <- t(Tt[cm, ct, drop = FALSE])
  both <- O != 3L & Tm != 3L
  expect_gt(sum(both), 10000L)
  expect_gt(mean(O[both] == Tm[both]), 0.999)                   # near-exact (a lone Viterbi tie is tolerated)
})

test_that("biparental route is BIT-EXACT vs TASSEL on a RIL (its target population)", {
  # RIL (inbred, ~all homozygous) is the population BiparentalHaplotypeFinder targets;
  # F2 (~50% het) is a documented stress case. On RIL the port is exact (task #5).
  fixf <- testthat::test_path("..", "fixtures", "fsfhap_tassel", "ril_fixture.rds")
  skip_if_not(file.exists(fixf), "RIL parity fixture absent (run agent/fsfhap_ril_parity.R)")
  fx <- readRDS(fixf)
  bip <- .fsfhap_biparental_call(fx$G, fx$pos, min_r2 = 0.0, min_gametes = 200L)
  imp <- .fsfhap_impute(bip$G, bip$pos, phet = 0.0, fill_gaps = TRUE)
  dimnames(imp) <- list(fx$taxa[bip$keep_taxa], fx$markers[bip$keep_sites])
  Tt <- fx$tassel_parents
  expect_setequal(fx$markers[bip$keep_sites], rownames(Tt))     # exact site-set match
  cm <- intersect(colnames(imp), rownames(Tt)); ct <- intersect(rownames(imp), colnames(Tt))
  O <- imp[ct, cm, drop = FALSE]; Tm <- t(Tt[cm, ct, drop = FALSE]); both <- O != 3L & Tm != 3L
  expect_gt(sum(both), 30000L)
  expect_equal(mean(O[both] == Tm[both]), 1)                    # bit-exact
})
