# binhmm caller (R/binhmm.R): bin -> K=3 cluster -> HMM-smooth, the rpubs
# "Ancestry Analysis by bins" pipeline wired in as caller = "binhmm". Distinct
# from the engine callers (it clusters per-bin ALT fraction, not reads).

# (a) cluster relabeling — regressions for the two degenerate-clustering bugs ---
test_that("bimodal signal never collapses a clear-ALT group into REF", {
  # a GMM component can win nothing on clean bimodal data; the survivors must
  # still span the extremes so the high group is ALT, not silently REF.
  af <- c(rep(0, 20), rep(0.02, 15), rep(0.95, 15))
  st <- nilHMM:::.binhmm_cluster(af, "gmm")
  expect_true(all(st[af > 0.5] == 2L))       # 0.95 bins -> ALT
  expect_true(all(st[af == 0]  == 0L))       # zero-ALT bins -> REF
})

test_that("clustering does not crash when non-zero bins have < K distinct values", {
  expect_silent(s1 <- nilHMM:::.binhmm_cluster(c(rep(0, 20), rep(0.5, 5)), "kmeans"))
  expect_true(all(s1 %in% 0:2))
  s2 <- nilHMM:::.binhmm_cluster(c(rep(0, 20), rep(0.02, 5), rep(0.9, 5)), "kmeans")
  expect_true(all(s2[c(26:30)] == 2L))       # the high, distinct group -> ALT
})

# (b) end-to-end via the top-level API -----------------------------------------
test_that("call_ancestry(caller='binhmm') recovers a donor block in the common schema", {
  set.seed(11)
  pos <- seq(5000L, 200e6, by = 5000L)               # ~200 1-Mb bins on chr1
  inblock <- pos >= 40e6 & pos <= 70e6
  n_alt <- ifelse(inblock, rbinom(length(pos), 12, 0.45), rbinom(length(pos), 12, 0.01))
  data <- data.frame(name = "S1", donor = "Zx", chr = 1L, pos = as.integer(pos),
                     n_ref = 12L - n_alt, n_alt = as.integer(n_alt))

  calls <- call_ancestry(data, caller = "binhmm", design = "BC2S3", source = "test")

  expect_named(calls, c("source", "donor", "name", "chr", "start_bp", "end_bp", "state"))
  expect_true(all(calls$state %in% 0:2))
  alt <- calls[calls$state == 2L, ]
  expect_gt(nrow(alt), 0L)                            # the block is detected as ALT
  # the ALT call overlaps the simulated 40-70 Mb donor block
  expect_true(any(alt$start_bp <= 70e6 & alt$end_bp >= 40e6))
})

test_that("the rebmix backend runs and returns valid calls (when rebmix is installed)", {
  skip_if_not_installed("rebmix")
  set.seed(11)
  pos <- seq(5000L, 200e6, by = 5000L)
  inblock <- pos >= 40e6 & pos <= 70e6
  n_alt <- ifelse(inblock, rbinom(length(pos), 12, 0.45), rbinom(length(pos), 12, 0.01))
  data <- data.frame(name = "S1", chr = 1L, pos = as.integer(pos),
                     n_ref = 12L - n_alt, n_alt = as.integer(n_alt))
  calls <- call_ancestry(data, caller = "binhmm", design = "BC2S3",
                         cluster_method = "rebmix", donor = "Zx")
  expect_named(calls, c("source", "donor", "name", "chr", "start_bp", "end_bp", "state"))
  expect_true(all(calls$state %in% 0:2))
  expect_true(any(calls$state[calls$start_bp <= 70e6 & calls$end_bp >= 40e6] == 2L))
})

test_that("the rebmix backend errors clearly when rebmix is absent", {
  skip_if(requireNamespace("rebmix", quietly = TRUE), "rebmix is installed")
  d <- data.frame(name = "s", chr = 1L, pos = as.integer(seq(5000L, 50e6, 5000L)),
                  n_ref = 5L, n_alt = 0L)
  expect_error(call_ancestry(d, caller = "binhmm", design = "BC2S3", cluster_method = "rebmix"),
               "rebmix")
})

test_that("binhmm is deterministic for the kmeans backend", {
  set.seed(2)
  pos <- seq(5000L, 100e6, by = 5000L)
  n_alt <- ifelse(pos >= 20e6 & pos <= 35e6, rbinom(length(pos), 10, 0.4),
                  rbinom(length(pos), 10, 0.01))
  data <- data.frame(name = "S1", chr = 1L, pos = as.integer(pos),
                     n_ref = 10L - n_alt, n_alt = as.integer(n_alt))
  a <- call_ancestry(data, caller = "binhmm", f_1 = 0.0312, f_2 = 0.1094,
                     cluster_method = "kmeans", donor = "Zx")
  b <- call_ancestry(data, caller = "binhmm", f_1 = 0.0312, f_2 = 0.1094,
                     cluster_method = "kmeans", donor = "Zx")
  expect_identical(a, b)
})
