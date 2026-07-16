# Scaffold-level checks: the package loads, the API surface exists, and the
# Rcpp Viterbi workhorse decodes a trivial chain correctly. Caller/engine
# behaviour is exercised by the §9.4 regression tests once Task 4 lands.

test_that("public API functions are exported", {
  for (fn in c("call_ancestry", "fit", "decode", "caller_spec",
               "emission_count", "emission_gt",
               "duration_geometric", "duration_rigidity", "duration_hsmm",
               "calibrate_r", "read_counts", "write_common_schema")) {
    expect_true(is.function(get(fn, envir = asNamespace("nilHMM"))), info = fn)
  }
})

test_that("emission/duration constructors return tagged specs", {
  expect_s3_class(emission_count(), "nilHMM_emission")
  expect_s3_class(emission_gt(), "nilHMM_emission")
  expect_s3_class(duration_geometric(), "nilHMM_duration")
  expect_s3_class(duration_rigidity(5L), "nilHMM_duration")
  expect_identical(duration_rigidity(5)$r, 5L)  # coerced to integer
})

test_that("the count and gt emissions run via call_ancestry", {
  # count emission (bbnil) takes read counts; gt emission (nnil) takes CALLED
  # genotypes -- it does not threshold counts, so feed it a `g` column.
  cnt <- data.frame(name = "s", chr = 1L, pos = 1:6, n_ref = c(2L,0L,1L,2L,0L,1L),
                    n_alt = c(0L,2L,1L,0L,2L,0L), donor = "Zx")
  gt  <- data.frame(name = "s", chr = 1L, pos = 1:6, g = c(0L,2L,1L,0L,2L,0L), donor = "Zx")
  bb <- call_ancestry(cnt, caller = "bbnil", design = "BC2S2", rrate = 0.01)
  nn <- call_ancestry(gt,  caller = "nnil",  design = "BC2S2", rrate = 0.01)
  expect_true(all(bb$state %in% 0:2))
  expect_true(all(nn$state %in% 0:2))
})

test_that("viterbi_log_cpp decodes a trivial 2-state chain", {
  # State 0 strongly favoured everywhere -> all-zero path.
  K <- 2; T <- 5
  log_init  <- log(c(0.9, 0.1))
  log_trans <- log(matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE))
  log_emit  <- matrix(log(c(0.95, 0.05)), nrow = T, ncol = K, byrow = TRUE)
  path <- nilHMM:::viterbi_log_cpp(log_init, log_trans, log_emit)
  expect_equal(path, rep(0L, T))

  # Emission flips to favour state 1 -> all-one path.
  log_emit2 <- matrix(log(c(0.05, 0.95)), nrow = T, ncol = K, byrow = TRUE)
  path2 <- nilHMM:::viterbi_log_cpp(log_init, log_trans, log_emit2)
  expect_equal(path2, rep(1L, T))
})