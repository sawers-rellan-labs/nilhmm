# GOOGA / ATLAS callers (R/atlas.R): competitive-alignment ancestry. Hard fraction
# call (recurrent = n_ref, donor = n_alt; ambiguous excluded upstream) -> gt
# (categorical confusion) emission HMM. `googa` = gt + geometric (faithful GOOGA),
# `atlas` = gt + rigidity (this work); both share the GOOGA thresholding.

test_that("GOOGA hard call maps donor fraction + depth to REF/HET/ALT/missing", {
  a <- c(0L, 10L, 5L, 1L, 2L, 8L)          # donor reads
  n <- c(10L, 10L, 10L, 2L, 10L, 10L)      # total informative reads
  g <- nilHMM:::.googa_gt_call(a, n)       # thresh 0.95, het 0.25, min 5
  expect_equal(g[1], 0L)   # f=0.00 -> REF
  expect_equal(g[2], 2L)   # f=1.00 -> ALT
  expect_equal(g[3], 1L)   # f=0.50 -> HET (both parents >= 0.25)
  expect_equal(g[4], 3L)   # depth 2 < 5 -> missing
  expect_equal(g[5], 3L)   # f=0.20 -> ambiguous (donor<0.25, recur<0.95) -> missing
  expect_equal(g[6], 3L)   # f=0.80 -> ambiguous (recur<0.25, donor<0.95) -> missing
})

test_that("googa/atlas recover a donor block; control is all-REF", {
  set.seed(1)
  pos <- seq(1e6, 200e6, by = 1e6)                     # 200 genes on chr1
  inblock <- pos >= 40e6 & pos <= 70e6
  n_alt <- ifelse(inblock, rbinom(length(pos), 20, 0.97), rbinom(length(pos), 20, 0.03))
  data <- data.frame(name = "S1", donor = "TIL18", chr = 1L, pos = as.integer(pos),
                     n_ref = 20L - n_alt, n_alt = as.integer(n_alt))
  ctrl <- data.frame(name = "B73", chr = 1L, pos = as.integer(pos), n_ref = 20L, n_alt = 0L)
  for (cl in c("googa", "atlas")) {
    calls <- call_ancestry(data, caller = cl, design = "BC2S3", source = cl)
    expect_named(calls, c("source", "donor", "name", "chr", "start_bp", "end_bp", "state"))
    expect_true(all(calls$state %in% 0:2), info = cl)
    alt <- calls[calls$state == 2L, ]
    expect_gt(nrow(alt), 0L)
    expect_true(any(alt$start_bp <= 70e6 & alt$end_bp >= 40e6), info = cl)  # donor block -> ALT

    cc <- call_ancestry(ctrl, caller = cl, design = "BC2S3")
    expect_gt(nrow(cc), 0L)
    expect_true(all(cc$state == 0L), info = cl)                             # pure recurrent -> all REF
  }
})

test_that("googa/atlas gate on min informative reads (low-depth genes -> no donor)", {
  pos <- seq(1e6, 100e6, by = 1e6)
  # every gene has only 3 informative reads (< min 5), even if donor-leaning
  data <- data.frame(name = "L", chr = 1L, pos = as.integer(pos), n_ref = 1L, n_alt = 2L)
  for (cl in c("googa", "atlas")) {
    calls <- call_ancestry(data, caller = cl, design = "BC2S3")
    expect_gt(nrow(calls), 0L)
    expect_true(all(calls$state == 0L), info = cl)                          # all missing -> smoothed to REF
  }
})
