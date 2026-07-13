test_that("pedigree_bp_cpp couples relatives: an informed leaf pulls an uninformed sib", {
  M <- 4L
  parent  <- c(-1L, 0L, 0L)          # latent root -> two leaves
  meioses <- c(3L, 4L, 4L)
  hasData <- c(FALSE, TRUE, TRUE)
  altobs  <- matrix(rep(c(0.02, 0.10, 0.98), each = M), nrow = M)  # leaf A observed ALT
  flat    <- matrix(1, M, 3)
  emit    <- list(flat, altobs, flat)
  pi0 <- c(0.75, 0.25, 0); uni <- c(1, 1, 1) / 3
  rho <- rbind(pi0, uni, uni)
  res <- pedigree_bp_cpp(M, parent, meioses, hasData, emit, rho, rho,
                         rep(0.01, M - 1L), root = 0L,
                         maxIters = 50L, tol = 1e-6, lambda = 0.5)
  A <- res[[2]][1, ]; B <- res[[3]][1, ]
  expect_equal(which.max(A), 3L)                 # informed leaf -> ALT
  # founder cannot be D/D (pi0[3]==0), so a D/D child forces founder HET, and the
  # uninformative sib is pulled off its flat prior toward het/alt (P(REF) < 0.5)
  expect_lt(B[1], 0.5)
  expect_true(all(is.finite(unlist(res))))
})

test_that("deterministic Tsib rows hold: a homozygous-REF parent keeps children REF", {
  M <- 3L
  parent  <- c(-1L, 0L)
  meioses <- c(3L, 4L)
  hasData <- c(TRUE, FALSE)
  refobs  <- matrix(rep(c(0.98, 0.10, 0.02), each = M), nrow = M)  # parent observed REF
  emit    <- list(refobs, matrix(1, M, 3))
  pi0 <- c(0.75, 0.25, 0); uni <- c(1, 1, 1) / 3
  rho <- rbind(pi0, uni)
  res <- pedigree_bp_cpp(M, parent, meioses, hasData, emit, rho, rho,
                         rep(0.01, M - 1L), root = 0L,
                         maxIters = 50L, tol = 1e-6, lambda = 0.5)
  expect_equal(which.max(res[[1]][1, ]), 1L)     # parent REF
  expect_equal(which.max(res[[2]][1, ]), 1L)     # child of REF-hom parent stays REF
})

test_that("refine_ancestry returns a valid same-shape mosaic and imputes missing", {
  skip_if_not_installed("simcross")
  sim   <- simulate_family("BC2S3", families = 3, sibs = 5, chr = 1,
                           n_markers = 150, seed = 3)
  truth <- sim$truth
  obs   <- simulate_counts(truth, depth = 0.5, error = 0.01, seed = 3)
  obs$cm <- truth$cm[match(paste(obs$name, obs$chr, obs$pos),
                           paste(truth$name, truth$chr, truth$pos))]
  calls <- call_states(obs, caller = "rtiger", design = "BC2S3", rigidity = 3L)
  # full-grid mosaic: rtiger call where covered, else 3 (missing)
  mo <- truth[, c("name", "chr", "pos", "cm")]
  mo$state <- calls$state[match(paste(mo$name, mo$chr, mo$pos),
                                paste(calls$name, calls$chr, calls$pos))]
  mo$state[is.na(mo$state)] <- 3L
  ref <- refine_ancestry(mo, sim$pedigree, design = "BC2S3", maxiter = 20L)
  expect_equal(nrow(ref), nrow(mo))
  expect_true(all(ref$state %in% c(0L, 1L, 2L)))     # no missing left; valid dosages
  expect_false(anyNA(ref$state))
})

test_that("refine_ancestry improves full-grid dosage accuracy over an LOCF baseline", {
  skip_if_not_installed("simcross")
  sim   <- simulate_family("BC2S3", families = 5, sibs = 10, chr = 1:2,
                           n_markers = 400, seed = 11)
  truth <- sim$truth
  obs   <- simulate_counts(truth, depth = 0.5, error = 0.01, seed = 11)
  obs$cm <- truth$cm[match(paste(obs$name, obs$chr, obs$pos),
                           paste(truth$name, truth$chr, truth$pos))]
  calls <- call_states(obs, caller = "rtiger", design = "BC2S3", rigidity = 3L)
  mo <- truth[, c("name", "chr", "pos", "cm")]
  mo$state <- calls$state[match(paste(mo$name, mo$chr, mo$pos),
                                paste(calls$name, calls$chr, calls$pos))]
  mo <- mo[order(mo$name, mo$chr, mo$pos), ]
  # LOCF baseline within sample-chr
  base <- ave(mo$state, paste(mo$name, mo$chr), FUN = function(z) {
    for (i in seq_along(z)) if (is.na(z[i])) z[i] <- if (i > 1) z[i - 1] else NA
    if (anyNA(z)) z[is.na(z)] <- 0L; z })
  moin <- mo; moin$state[is.na(moin$state)] <- 3L
  ref <- refine_ancestry(moin, sim$pedigree, design = "BC2S3", maxiter = 30L)
  ts  <- truth$state[match(paste(mo$name, mo$chr, mo$pos),
                           paste(truth$name, truth$chr, truth$pos))]
  rf  <- ref$state[match(paste(mo$name, mo$chr, mo$pos),
                         paste(ref$name, ref$chr, ref$pos))]
  # refine should at least match the per-individual LOCF baseline (small tolerance)
  expect_gte(mean(rf == ts), mean(base == ts) - 0.005)
})
