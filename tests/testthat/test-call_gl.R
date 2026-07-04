# call_gl(): closed-form per-site GL genotype caller with a swappable prior. The
# het-excess control for the coverage-degradation study -- assert the decision rule
# exactly (worked-out log-posterior math), not by eyeball.

test_that("single ALT read: flat is het-blind (hom-ALT), HWE calls HET (het excess)", {
  # (n_ref=0, n_alt=1): flat -> argmax GL = hom-ALT (L2 > L1 > L0).
  expect_equal(call_gl(0, 1, prior = "flat"), 2L)
  # HWE with moderate p=0.3: P1 = log(.5)+log(2*.3*.7) beats P2 = log(.99)+2log(.3).
  expect_equal(call_gl(0, 1, prior = "hwe", af = 0.3), 1L)
  # even at the low teosinte AF of a TeoNAM family (p~0.06) a single ALT read is
  # still HET -- the hom-ALT prior p^2 is smaller still. This IS the het excess.
  expect_equal(call_gl(0, 1, prior = "hwe", af = 0.06), 1L)
})

test_that("single REF read: both priors call REF (hom-REF)", {
  # (n_ref=1, n_alt=0): flat -> 0 (L0 > L1 > L2).
  expect_equal(call_gl(1, 0, prior = "flat"), 0L)
  # HWE p=0.3: P0 = log(.99)+2log(.7) is the max -> 0.
  expect_equal(call_gl(1, 0, prior = "hwe", af = 0.3), 0L)
})

test_that("high depth: the data dominates any prior", {
  # 10 ALT / 0 REF -> hom-ALT regardless of prior (even a hostile low-ALT prior).
  for (pr in list(list(prior = "flat"),
                  list(prior = "hwe", af = 0.01),
                  list(prior = "breeding", f = c(0.98, 0.01, 0.01)))) {
    expect_equal(do.call(call_gl, c(list(0, 10), pr)), 2L, info = pr$prior)
    expect_equal(do.call(call_gl, c(list(10, 0), pr)), 0L, info = pr$prior)
  }
  # balanced high depth -> HET under any reasonable prior.
  expect_equal(call_gl(8, 8, prior = "hwe", af = 0.3), 1L)
})

test_that("zero depth -> NA (missing), preserved in matrix shape", {
  expect_true(is.na(call_gl(0, 0, prior = "flat")))
  expect_true(is.na(call_gl(0, 0, prior = "hwe", af = 0.2)))
  R <- matrix(c(0, 5, 0, 1), 2, 2); A <- matrix(c(0, 0, 3, 0), 2, 2)
  out <- call_gl(R, A, prior = "flat")
  expect_equal(dim(out), c(2L, 2L))
  expect_true(is.na(out[1, 1]))          # (0,0) -> NA
  expect_equal(out[2, 1], 0L)            # (5 ref, 0 alt) -> REF
  expect_equal(out[1, 2], 2L)            # (0 ref, 3 alt) -> ALT
  expect_equal(out[2, 2], 0L)            # (1 ref, 0 alt) -> REF
})

test_that("vectorized matrix input returns the right shape; vector in -> vector out", {
  set.seed(1)
  M <- 50; N <- 4
  R <- matrix(rpois(M * N, 3), M, N); A <- matrix(rpois(M * N, 3), M, N)
  af <- runif(M, 0.05, 0.5)
  out <- call_gl(R, A, prior = "hwe", af = af)
  expect_equal(dim(out), c(M, N))
  expect_true(all(out %in% c(0L, 1L, 2L) | is.na(out)))
  # a plain vector is treated as one sample and returns a plain vector
  vout <- call_gl(R[, 1], A[, 1], prior = "hwe", af = af)
  expect_null(dim(vout))
  expect_equal(length(vout), M)
  expect_equal(unname(vout), unname(out[, 1]))   # same call, column vs vector
})

test_that("af length mismatch and bad error are rejected", {
  expect_error(call_gl(c(0, 1), c(1, 0), prior = "hwe", af = 0.3), "length")
  expect_error(call_gl(1, 0, prior = "flat", error = 0.6), "error")
  expect_error(call_gl(1, 0, prior = "breeding", f = c(0.5, 0.5)), "breeding")
  expect_error(call_gl(matrix(0, 2, 2), matrix(0, 2, 3), prior = "flat"), "shape")
})

test_that("af = NULL estimates p per marker from the reads", {
  # marker 1: all-ALT reads across samples -> p_hat ~ 1 -> hom-ALT where covered.
  # marker 2: all-REF -> p_hat ~ 0 -> hom-REF.
  R <- rbind(c(0, 0, 0), c(5, 4, 6)); A <- rbind(c(5, 4, 6), c(0, 0, 0))
  out <- call_gl(R, A, prior = "hwe", af = NULL)
  expect_true(all(out[1, ] == 2L))
  expect_true(all(out[2, ] == 0L))
})

test_that("return = 'post' gives normalized posteriors summing to 1; NA at zero depth", {
  post <- call_gl(c(0, 0, 2), c(1, 0, 8), prior = "hwe", af = c(0.3, 0.3, 0.3),
                  return = "post")
  expect_equal(dim(post), c(3L, 1L, 3L))
  # row 1 (single ALT) -> HET is the modal posterior
  expect_equal(unname(which.max(post[1, 1, ])), 2L)
  expect_equal(sum(post[1, 1, ]), 1, tolerance = 1e-12)
  # row 2 (zero depth) -> NA
  expect_true(all(is.na(post[2, 1, ])))
  # 'dosage' return is the hard call as double
  d <- call_gl(c(0, 1), c(1, 0), prior = "hwe", af = c(0.3, 0.3), return = "dosage")
  expect_type(d, "double")
  expect_equal(unname(d), c(1, 0))   # (0 ref,1 alt)->HET=1 ; (1 ref,0 alt)->REF=0
})
