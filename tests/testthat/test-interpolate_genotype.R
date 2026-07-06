# interpolate_genotype(): flanking-marker genotype interpolation in cM
# (Tian 2011 continuous / Chen-TeoNAM step / round). Deterministic hard-call
# densification onto a target grid -- NOT ancestry inference.

test_that("concordant fill: a single observed state fills every target (all modes)", {
  obs <- data.frame(chr = 1L, cm = c(0, 1, 2, 3))
  target <- data.frame(chr = 1L, cm = c(-1, 0.4, 1.7, 3, 5))
  for (state in c(0, 2)) {
    geno <- matrix(state, nrow = 4, ncol = 2, dimnames = list(NULL, c("A", "B")))
    for (mode in c("continuous", "step", "round")) {
      out <- interpolate_genotype(geno, obs, target, mode)
      expect_equal(unname(out), matrix(state, nrow = nrow(target), ncol = 2),
                   info = paste(mode, state))
    }
  }
})

test_that("empty target grid returns a 0-row matrix with geno's columns", {
  obs <- data.frame(chr = 1L, cm = c(0, 1))
  geno <- matrix(c(0, 2, 1, 1), nrow = 2, dimnames = list(NULL, c("S1", "S2")))
  target <- data.frame(chr = integer(0), cm = numeric(0))
  out <- interpolate_genotype(geno, obs, target, "continuous")
  expect_equal(dim(out), c(0L, 2L))
  expect_identical(colnames(out), c("S1", "S2"))
})

test_that("continuous ramp matches hand-computed values; step/round at midpoint", {
  obs <- data.frame(chr = 1L, cm = c(0, 1))
  geno <- matrix(c(0, 2), nrow = 2, dimnames = list(NULL, "S1"))
  target <- data.frame(chr = 1L, cm = c(0, 0.25, 0.5, 0.75, 1))

  cont <- interpolate_genotype(geno, obs, target, "continuous")
  expect_equal(as.numeric(cont), c(0, 0.5, 1.0, 1.5, 2.0))

  step <- interpolate_genotype(geno, obs, target, "step")
  # w = 0, .25, .5, .75, 1 -> vL,vL,vR(tie),vR,vR = 0,0,2,2,2
  expect_equal(as.numeric(step), c(0, 0, 2, 2, 2))

  rnd <- interpolate_genotype(geno, obs, target, "round")
  # C++ std::round is round-half-away-from-zero: round(0,.5,1,1.5,2) -> 0,1,1,2,2
  expect_equal(as.numeric(rnd), c(0, 1, 1, 2, 2))
})

test_that("step fabricates no het across a 0<->2 gap; round does", {
  obs <- data.frame(chr = 1L, cm = c(0, 10))
  geno <- matrix(c(0, 2), nrow = 2, dimnames = list(NULL, "S1"))
  target <- data.frame(chr = 1L, cm = seq(0, 10, by = 0.5))

  step <- as.numeric(interpolate_genotype(geno, obs, target, "step"))
  expect_true(all(step %in% c(0, 2)))
  expect_false(any(step == 1))

  rnd <- as.numeric(interpolate_genotype(geno, obs, target, "round"))
  expect_true(all(rnd %in% c(0, 1, 2)))
  expect_true(any(rnd == 1))   # fabricated het band
})

test_that("ends are clamped to the terminal observed value (all modes)", {
  obs <- data.frame(chr = 1L, cm = c(2, 4, 6))
  geno <- matrix(c(0, 1, 2), nrow = 3, dimnames = list(NULL, "S1"))
  target <- data.frame(chr = 1L, cm = c(-5, 0, 10, 100))  # all outside [2, 6]
  for (mode in c("continuous", "step", "round")) {
    out <- as.numeric(interpolate_genotype(geno, obs, target, mode))
    expect_equal(out, c(0, 0, 2, 2), info = mode)  # left->first, right->last
  }
})

test_that("a target exactly on an observed cM returns the observed value", {
  obs <- data.frame(chr = 1L, cm = c(0, 1, 2))
  geno <- matrix(c(0, 1, 2), nrow = 3, dimnames = list(NULL, "S1"))
  target <- data.frame(chr = 1L, cm = c(0, 1, 2))
  for (mode in c("continuous", "step", "round")) {
    out <- as.numeric(interpolate_genotype(geno, obs, target, mode))
    expect_equal(out, c(0, 1, 2), info = mode)
  }
})

test_that("tied target cM is allowed: every marker imputed, twins share genotypes", {
  # Centromeric plateau: two target markers at the SAME cM (m2a, m2b). The target
  # grid need NOT be unique in cM -- each marker still gets its own output row, and
  # the twins receive an identical genotype vector (perfect LD in the map).
  obs  <- data.frame(chr = 1L, cm = c(0, 1, 2))
  geno <- matrix(c(0, 0, 2, 0, 0, 2), ncol = 2, dimnames = list(NULL, c("S1", "S2")))
  target <- data.frame(chr = 1L, cm = c(0.0, 1.0, 1.0, 2.0),
                       row.names = c("m1", "m2a", "m2b", "m3"))
  for (mode in c("continuous", "step", "round")) {
    out <- interpolate_genotype(geno, obs, target, mode)
    expect_equal(nrow(out), nrow(target), info = mode)          # nothing dropped
    expect_equal(unname(out["m2a", ]), unname(out["m2b", ]), info = mode)  # twins identical
  }
})

test_that("no interpolation bleeds across a chromosome boundary", {
  obs <- data.frame(chr = c(1L, 1L, 2L, 2L), cm = c(0, 1, 0, 1))
  # chr1: 0 -> 0 ; chr2: 2 -> 2. A midpoint must stay within its own chr.
  geno <- matrix(c(0, 0, 2, 2), nrow = 4, dimnames = list(NULL, "S1"))
  target <- data.frame(chr = c(1L, 2L), cm = c(0.5, 0.5))
  out <- as.numeric(interpolate_genotype(geno, obs, target, "continuous"))
  expect_equal(out, c(0, 2))  # if it bled, chr1's midpoint would drift toward 2
})

test_that("continuous mode equals stats::approx(rule = 2) per sample", {
  set.seed(1)
  obs <- data.frame(chr = 1L, cm = sort(cumsum(runif(40, 0.1, 2))))
  geno <- matrix(runif(40 * 5, 0, 2), nrow = 40,
                 dimnames = list(NULL, paste0("S", 1:5)))
  target <- data.frame(chr = 1L,
                       cm = sort(runif(120, min(obs$cm) - 5, max(obs$cm) + 5)))
  out <- interpolate_genotype(geno, obs, target, "continuous")
  for (j in 1:5) {
    ref <- stats::approx(obs$cm, geno[, j], xout = target$cm, rule = 2)$y
    expect_equal(out[, j], ref, info = paste("sample", j))
  }
})

test_that("row/col names are carried through", {
  obs <- data.frame(chr = 1L, cm = c(0, 1))
  geno <- matrix(c(0, 2), nrow = 2, dimnames = list(NULL, "SAMP"))
  target <- data.frame(chr = 1L, cm = c(0, 0.5, 1))
  rownames(target) <- c("t1", "t2", "t3")
  out <- interpolate_genotype(geno, obs, target, "continuous")
  expect_equal(colnames(out), "SAMP")
  expect_equal(rownames(out), c("t1", "t2", "t3"))
})

test_that("validation: NA, mismatched nrow, and unsorted obs are rejected", {
  obs <- data.frame(chr = 1L, cm = c(0, 1, 2))
  geno <- matrix(c(0, 1, 2), nrow = 3)
  target <- data.frame(chr = 1L, cm = c(0, 1))

  g_na <- geno; g_na[2] <- NA
  expect_error(interpolate_genotype(g_na, obs, target, "step"), "complete")
  expect_error(interpolate_genotype(matrix(0, nrow = 2), obs, target, "step"),
               "nrow")
  obs_bad <- data.frame(chr = 1L, cm = c(0, 0, 1))  # tied -> not strictly increasing
  expect_error(
    interpolate_genotype(geno, obs_bad, target, "step"),
    "sorted|strictly increasing")
})

test_that("coord = 'bp' is identical to the same data with the column named cm", {
  set.seed(1)
  # a two-chromosome block with a mix of REF/HET/ALT so all modes are exercised
  obs_bp <- data.frame(chr = c(1L, 1L, 1L, 2L, 2L), bp = c(1e6, 3e6, 7e6, 2e6, 9e6))
  geno   <- matrix(sample(0:2, 5 * 3, replace = TRUE), nrow = 5,
                   dimnames = list(NULL, c("A", "B", "C")))
  target_bp <- data.frame(chr = c(1L, 1L, 1L, 1L, 2L, 2L, 2L),
                          bp = c(1e6, 2e6, 5e6, 7e6, 2e6, 2e6, 9e6))  # incl a tied target
  # rename the bp column to cm -> default coord must give the identical result
  obs_cm    <- setNames(obs_bp,    sub("^bp$", "cm", names(obs_bp)))
  target_cm <- setNames(target_bp, sub("^bp$", "cm", names(target_bp)))
  for (mode in c("continuous", "step", "round")) {
    a <- interpolate_genotype(geno, obs_bp, target_bp, mode, coord = "bp")
    b <- interpolate_genotype(geno, obs_cm, target_cm, mode)  # default coord = "cm"
    expect_equal(a, b, info = mode)
  }
})

test_that("coord errors: missing coord column and non-increasing coord within a chr", {
  geno <- matrix(c(0, 1, 2), nrow = 3, dimnames = list(NULL, "S1"))
  obs  <- data.frame(chr = 1L, bp = c(1e6, 2e6, 3e6))
  target <- data.frame(chr = 1L, bp = c(1e6, 2e6))
  # obs/target lack the requested coord column
  expect_error(interpolate_genotype(geno, obs, target, "step", coord = "cm"),
               "columns `chr` and `cm`")
  # non-strictly-increasing coord within a chromosome
  obs_bad <- data.frame(chr = 1L, bp = c(1e6, 1e6, 3e6))  # tied bp
  expect_error(interpolate_genotype(geno, obs_bad, target, "step", coord = "bp"),
               "sorted|strictly increasing")
  # coord must be a single column name
  expect_error(interpolate_genotype(geno, obs, target, "step", coord = c("bp", "cm")),
               "single column name")
})
