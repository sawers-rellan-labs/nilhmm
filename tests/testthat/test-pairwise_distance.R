# Pairwise relatedness builder (r2 / MI / VI). r2 is cross-checked against
# stats::cor^2; MI/VI are cross-checked against an independent R reimplementation
# of the same Miller-Madow estimator, plus the metric properties of VI.

# Independent R reference for the MM-corrected estimator over complete pairs.
mm_ref <- function(x, y, base = c("nats", "bits")) {
  base <- match.arg(base)
  ok <- is.finite(x) & is.finite(y)
  x <- round(x[ok]); y <- round(y[ok]); N <- length(x)
  px <- table(x) / N; py <- table(y) / N
  pxy <- table(x, y) / N; pxy <- pxy[pxy > 0]
  H <- function(p) -sum(p * log(p))
  Hx <- H(px) + (length(px) - 1) / (2 * N)
  Hy <- H(py) + (length(py) - 1) / (2 * N)
  Hxy <- H(as.numeric(pxy)) + (sum(pxy > 0) - 1) / (2 * N)
  I <- Hx + Hy - Hxy
  vi <- Hx + Hy - 2 * I
  div <- if (base == "bits") log(2) else 1
  list(mi = I / div, vi = vi / div, Hx = Hx / div)
}

test_that("r2 equals stats::cor(t(geno))^2", {
  set.seed(1)
  geno <- matrix(sample(0:2, 15 * 50, replace = TRUE), nrow = 15,
                 dimnames = list(paste0("m", 1:15), NULL))
  got <- pairwise_distance(geno, "r2")
  ref <- cor(t(geno))^2
  expect_equal(got, ref, ignore_attr = TRUE, tolerance = 1e-10)
  expect_identical(attr(got, "kind"), "similarity")
  expect_identical(attr(got, "method"), "r2")
  expect_equal(diag(got), setNames(rep(1, 15), rownames(geno)))
})

test_that("r2 is 0 for a constant (zero-variance) marker", {
  geno <- rbind(a = rep(1, 20), b = c(0, 1, 2, rep(1, 17)))
  got <- pairwise_distance(geno, "r2")
  expect_equal(got["a", "b"], 0)
})

test_that("mi matches an independent MM reference and MI(X,X) = H(X)", {
  set.seed(2)
  geno <- matrix(sample(0:2, 8 * 40, replace = TRUE), nrow = 8,
                 dimnames = list(paste0("m", 1:8), NULL))
  got <- pairwise_distance(geno, "mi")
  for (i in 1:8) for (j in 1:8) {
    ref <- mm_ref(geno[i, ], geno[j, ])$mi
    expect_equal(got[i, j], ref, tolerance = 1e-10, info = paste(i, j))
  }
  # MI(X,X) = H(X) (both MM-corrected)
  for (i in 1:8)
    expect_equal(got[i, i], mm_ref(geno[i, ], geno[i, ])$Hx, tolerance = 1e-10)
  expect_true(isSymmetric(unname(got)))
})

test_that("mi hand-computed 3x3 joint table", {
  # X and Y with a known 3x3 joint distribution (each cell occupied).
  x <- c(0,0,0, 1,1,1, 2,2,2)
  y <- c(0,1,2, 0,1,2, 0,1,2)          # X, Y independent by construction
  geno <- rbind(X = x, Y = y)
  got <- pairwise_distance(geno, "mi")
  ref <- mm_ref(x, y)
  expect_equal(got["X", "Y"], ref$mi, tolerance = 1e-12)
  # Plug-in MI of independent uniform marginals is exactly 0; MM adds
  # (mX+mY-mXY-1)/(2N) = (3+3-9-1)/18 = -4/18 -> a small negative estimate.
  expect_equal(got["X", "Y"], -4 / 18, tolerance = 1e-12)
})

test_that("mi in bits equals mi in nats divided by log(2)", {
  set.seed(3)
  geno <- matrix(sample(0:2, 6 * 30, replace = TRUE), nrow = 6)
  nats <- pairwise_distance(geno, "mi", base = "nats")
  bits <- pairwise_distance(geno, "mi", base = "bits")
  expect_equal(bits, nats / log(2), tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("vi is a metric: VI(X,X)=0, symmetric, triangle inequality", {
  set.seed(4)
  geno <- matrix(sample(0:2, 10 * 60, replace = TRUE), nrow = 10,
                 dimnames = list(paste0("m", 1:10), NULL))
  got <- pairwise_distance(geno, "vi")
  expect_identical(attr(got, "kind"), "distance")
  expect_equal(diag(got), setNames(rep(0, 10), rownames(geno)), tolerance = 1e-12)
  expect_true(isSymmetric(unname(got)))
  expect_true(all(got >= -1e-12))
  # triangle inequality on all triples
  ok <- TRUE
  for (a in 1:10) for (b in 1:10) for (c in 1:10)
    if (got[a, c] > got[a, b] + got[b, c] + 1e-9) ok <- FALSE
  expect_true(ok)
  # VI = H(X) + H(Y) - 2 I  (consistency with the mi builder)
  mi <- pairwise_distance(geno, "mi")
  expect_equal(got[1, 2], mi[1, 1] + mi[2, 2] - 2 * mi[1, 2], tolerance = 1e-10)
})

test_that("missing (NA) entries use complete pairs only", {
  x <- c(0, 1, 2, 0, 1, 2, NA, NA)
  y <- c(0, 1, 2, 0, 1, 2, 0, 1)
  geno <- rbind(X = x, Y = y)
  got <- pairwise_distance(geno, "r2")
  ref <- cor(x, y, use = "complete.obs")^2
  expect_equal(got["X", "Y"], ref, tolerance = 1e-10)
})

test_that("sense wiring: r2 and vi each yield an independent set in their sense", {
  set.seed(5)
  geno <- matrix(sample(0:2, 14 * 50, replace = TRUE), nrow = 14,
                 dimnames = list(paste0("m", 1:14), NULL))
  d_r2 <- pairwise_distance(geno, "r2")
  d_vi <- pairwise_distance(geno, "vi")

  sel_r2 <- select_independent(geno, threshold = 0.5, method = "r2")
  sel_vi <- select_independent(geno, threshold = 0.5, method = "vi")

  er2 <- d_r2 >= 0.5; diag(er2) <- FALSE
  evi <- d_vi <= 0.5; diag(evi) <- FALSE
  expect_false(any(er2[sel_r2, sel_r2, drop = FALSE]))
  expect_false(any(evi[sel_vi, sel_vi, drop = FALSE]))
})

test_that("precomputed matrix path honours attr(kind) for the sense", {
  set.seed(6)
  geno <- matrix(sample(0:2, 12 * 40, replace = TRUE), nrow = 12,
                 dimnames = list(paste0("m", 1:12), NULL))
  d_vi <- pairwise_distance(geno, "vi")            # carries attr(kind)="distance"
  sel <- select_independent(d_vi, threshold = 0.5) # sense taken from attr
  expect_identical(attr(sel, "kind"), "distance")
  evi <- d_vi <= 0.5; diag(evi) <- FALSE
  expect_false(any(evi[sel, sel, drop = FALSE]))
})
