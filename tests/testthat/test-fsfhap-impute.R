# FSFHap port — stage 3: 5-state EM Viterbi imputation (imputeUsingViterbiFiveState)
# + fillGapsInAlignment. Deterministic-behavior checks; end-to-end TASSEL parity
# lives in the parity harness (imputed HapMap is exportable).

int_mat <- function(x) { storage.mode(x) <- "integer"; x }

# a backcross-ish family: each taxon is A-hom then C-hom at a random breakpoint
mk_family <- function(ntaxa = 20L, nsites = 30L, seed = 1L) {
  set.seed(seed)
  G <- matrix(0L, ntaxa, nsites)
  for (t in seq_len(ntaxa)) {
    bp <- sample(5:(nsites - 5), 1L)
    if (runif(1) < 0.5) G[t, bp:nsites] <- 2L
  }
  int_mat(G)
}
pos30 <- as.integer(seq_len(30L) * 1e5L)

test_that("clean family: imputation recovers the true genotypes, no missing", {
  G <- mk_family()
  r <- fsfhap_impute_five_state_cpp(G, pos30, phet = 0.03125, max_iter = 50L)
  expect_true(all(r$imputed == G))
  expect_false(any(r$imputed == 3L))
  expect_lte(r$iters, 50L)
})

test_that("5-state HMM smooths single-site genotyping errors (het in a hom block)", {
  G <- mk_family(); truth <- G
  G[3, 10] <- 1L; G[7, 20] <- 1L                 # spurious hets inside hom runs
  r <- fsfhap_impute_five_state_cpp(G, pos30, phet = 0.03125, max_iter = 50L)
  expect_equal(r$imputed[3, 10], truth[3, 10])   # corrected back to hom
  expect_equal(r$imputed[7, 20], truth[7, 20])
})

test_that("the 5-state decodes only observed sites; missing stay N until fillgaps", {
  G <- mk_family(); truth <- G
  G[5, 12] <- 3L; G[5, 13] <- 3L                 # missing inside a run
  r <- fsfhap_impute_five_state_cpp(G, pos30, phet = 0.03125, max_iter = 50L)
  # missing sites are NOT in the observation sequence -> not decoded -> stay N(3)
  expect_equal(r$imputed[5, 12], 3L)
  expect_equal(r$imputed[5, 13], 3L)
  # observed sites still recover truth
  expect_true(all(r$imputed[5, G[5, ] != 3L] == truth[5, G[5, ] != 3L]))
})

test_that("phet is a required (design-derived) argument on .fsfhap_impute", {
  G <- mk_family()
  expect_error(.fsfhap_impute(G, pos30), "argument \"phet\"")
})

test_that(".fsfhap_phet derives (1-F)/2, falls back outside [0,1]", {
  expect_equal(.fsfhap_phet(0.9375), 0.03125)    # BC1S4
  expect_equal(.fsfhap_phet(0.75), 0.125)        # BC2S2
  expect_equal(.fsfhap_phet(NA), 0.07)           # fallback
  expect_equal(.fsfhap_phet(2), 0.07)            # out of range
})

# ---- fillGapsInAlignment ----------------------------------------------------

test_that("fillgaps fills missing runs between EQUAL flanks, not different ones", {
  G <- int_mat(matrix(c(0, 3, 3, 0,  2, 3, 0,  1, 3, 1), nrow = 1))
  f <- fsfhap_fill_gaps_cpp(G)[1, ]
  expect_equal(f[1:4], c(0L, 0L, 0L, 0L))        # 0 _ _ 0 -> filled
  expect_equal(f[5:7], c(2L, 3L, 0L))            # 2 _ 0 -> different flanks, unfilled
  expect_equal(f[8:10], c(1L, 1L, 1L))           # 1 _ 1 -> filled
})

test_that("fillgaps leaves leading/trailing missing untouched", {
  G <- int_mat(matrix(c(3, 3, 0, 0, 3, 3), nrow = 1))
  expect_equal(fsfhap_fill_gaps_cpp(G)[1, ], c(3L, 3L, 0L, 0L, 3L, 3L))
})

test_that(".fsfhap_impute chains EM + fillgaps and carries EM attributes", {
  G <- mk_family()
  out <- .fsfhap_impute(G, pos30, phet = 0.03125)
  expect_equal(dim(out), dim(G))
  expect_true(is.numeric(attr(out, "iters")))
  expect_equal(dim(attr(out, "emission")), c(5L, 3L))
})
