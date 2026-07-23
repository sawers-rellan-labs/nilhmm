# Parity gate for the unified sweep kernel (viterbi_sweep_cpp), Scope A.
# Before lb_viterbi_sweep_cpp / lb_viterbi_cpp can be retired, the generalized
# kernel must reproduce them byte-for-byte, and its new "const" mode must match
# the generic viterbi_log_cpp per transition matrix.

test_that("distance mode is byte-identical to lb_viterbi_sweep_cpp", {
  set.seed(1)
  T <- 300L
  log_emit <- matrix(log(runif(T * 3)), T, 3) # continuous -> no emission ties
  log_init <- log(rep(1 / 3, 3))
  tpos <- as.numeric(cumsum(sample(1:1000, T, replace = TRUE))) # strictly increasing
  recombdists <- c(1e3, 5e3, 1e4, 5e4, 1e5)
  for (drp in c(FALSE, TRUE)) {
    old <- lb_viterbi_sweep_cpp(log_init, log_emit, tpos, recombdists, drp)
    new <- viterbi_sweep_cpp(log_init, log_emit,
      mode = "distance",
      tpos = tpos, recombdists = recombdists, drp = drp
    )
    expect_identical(new, old)
  }
})

test_that("distance grid-of-1 equals lb_viterbi_cpp single decode", {
  set.seed(2)
  T <- 250L
  log_emit <- matrix(log(runif(T * 3)), T, 3)
  log_init <- log(rep(1 / 3, 3))
  tpos <- as.numeric(cumsum(sample(1:1000, T, replace = TRUE)))
  for (drp in c(FALSE, TRUE)) {
    for (rd in c(1e3, 1e5)) {
      single <- lb_viterbi_cpp(log_init, log_emit, tpos, rd, drp)
      swept <- viterbi_sweep_cpp(log_init, log_emit,
        mode = "distance",
        tpos = tpos, recombdists = rd, drp = drp
      )[, 1]
      expect_identical(as.integer(swept), as.integer(single))
    }
  }
})

test_that("const mode matches viterbi_log_cpp per transition matrix (any S)", {
  set.seed(3)
  for (S in c(3L, 4L)) {
    T <- 220L
    log_emit <- matrix(log(runif(T * S)), T, S)
    log_init <- log(rep(1 / S, S))
    mk_trans <- function(sw) { # symmetric switch-rate matrix, rows sum to 1
      M <- matrix(log(sw / (S - 1)), S, S)
      diag(M) <- log(1 - sw)
      M
    }
    trans_list <- lapply(c(0.001, 0.01, 0.1), mk_trans)
    swept <- viterbi_sweep_cpp(log_init, log_emit, mode = "const", trans_list = trans_list)
    for (vi in seq_along(trans_list)) {
      ref <- viterbi_log_cpp(log_init, trans_list[[vi]], log_emit) # tie_break = 0 default
      expect_identical(as.integer(swept[, vi]), as.integer(ref))
    }
  }
})

test_that("const mode honours tie_break, matching viterbi_log_cpp on tying emissions", {
  set.seed(4)
  S <- 3L
  T <- 240L
  log_emit <- matrix(log(runif(T * S)), T, S)
  miss <- sample(T, 40) # GT-like: missing -> all states equal
  log_emit[miss, ] <- log(1 / 3)
  het <- sample(setdiff(seq_len(T), miss), 40) # het -> REF == ALT
  log_emit[het, 1] <- log_emit[het, 3]
  log_init <- log(rep(1 / S, S))
  mk <- function(sw) {
    M <- matrix(log(sw / (S - 1)), S, S)
    diag(M) <- log(1 - sw)
    M
  }
  trans_list <- lapply(c(0.01, 0.1), mk)
  for (tb in c(0L, 1L)) {
    swept <- viterbi_sweep_cpp(log_init, log_emit,
      mode = "const", trans_list = trans_list, tie_break = tb
    )
    for (vi in seq_along(trans_list)) {
      ref <- viterbi_log_cpp(log_init, trans_list[[vi]], log_emit, tie_break = tb)
      expect_identical(as.integer(swept[, vi]), as.integer(ref))
    }
  }
})

test_that("empty / single-marker inputs are handled", {
  li <- log(rep(1 / 3, 3))
  e0 <- matrix(numeric(0), 0, 3)
  expect_equal(dim(viterbi_sweep_cpp(li, e0, mode = "distance", tpos = numeric(0), recombdists = 1e4)), c(0L, 1L))
  e1 <- matrix(log(runif(3)), 1, 3)
  d1 <- viterbi_sweep_cpp(li, e1, mode = "distance", tpos = 0, recombdists = c(1e3, 1e4))
  expect_equal(dim(d1), c(1L, 2L))
})
