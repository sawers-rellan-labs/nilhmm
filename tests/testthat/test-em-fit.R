# EM-fit emission means (fit_means = TRUE, S10). The mechanism is validated by
# recovering known means from synthetic data (it is a deliberate divergence from
# the fixed-means baseline, so there is no baseline to match). The fixed-means
# default is covered by the strict golden-slice regression.

test_that("forward_backward posteriors are normalized", {
  td <- nilHMM:::.build_transition(1e-3, 0.0625, 0.0938)
  em <- nilHMM:::count_emission_loglik_cpp(c(1L, 2L, 0L, 1L), c(0L, 2L, 0L, 1L),
                                           c(0.01, 0.5, 0.99), 20)
  g <- nilHMM:::forward_backward_cpp(td$log_start, td$log_trans, em)
  expect_equal(unname(rowSums(g)), rep(1, 4), tolerance = 1e-8)
  expect_true(all(g >= 0 & g <= 1))
})

test_that("EM recovers known emission means from a perturbed start", {
  set.seed(42)
  true_theta <- c(0.02, 0.5, 0.95)
  blocks <- rep(c(0, 1, 2, 0, 2, 1, 0), each = 400)
  n <- rpois(length(blocks), 2)
  a <- rbinom(length(blocks), n, true_theta[blocks + 1])
  obs <- list(list(chr = 1, pos = seq_along(blocks), n = n, a = a))
  td  <- nilHMM:::.build_transition(1e-3, 0.0625, 0.0938)
  emc <- emission_count(err = 0.01, conc = 20, fit_means = TRUE)

  fit_theta <- nilHMM:::.em_fit_means(obs, emc, td, c(0.10, 0.40, 0.70))
  expect_equal(fit_theta, true_theta, tolerance = 0.05)
  expect_true(all(fit_theta > 0 & fit_theta < 1))
})

test_that("a never-visited state keeps its initial mean (no NaN)", {
  # all-REF data: HET/ALT never supported -> their means must stay finite
  obs <- list(list(chr = 1, pos = 1:200, n = rep(1L, 200), a = rep(0L, 200)))
  td  <- nilHMM:::.build_transition(1e-3, 0.0625, 0.0938)
  emc <- emission_count(err = 0.01, conc = 20, fit_means = TRUE)
  th  <- nilHMM:::.em_fit_means(obs, emc, td, c(0.01, 0.5, 0.99))
  expect_false(anyNA(th))
  expect_true(all(th > 0 & th < 1))
})

test_that("fit_means = TRUE runs end-to-end and stays well-defined on skim", {
  raw <- read_golden_counts(golden_slice_path("skim", "counts", "PN14_SID1259.chr1.tsv"))
  counts <- data.frame(name = "PN14_SID1259", chr = 1L, pos = raw$pos,
                       n_ref = raw$n_ref, n_alt = raw$n_alt, donor = "Zd")
  got <- call_ancestry(counts, caller = "nnil", design = "BC2S2", rrate = 7e-6, fit_means = TRUE)
  expect_true(all(c("source", "donor", "name", "chr", "start_bp", "end_bp", "state") %in% names(got)))
  expect_true(all(got$state %in% 0:2))
  expect_true(all(got$end_bp >= got$start_bp))
})
