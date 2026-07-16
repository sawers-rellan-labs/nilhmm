# gt (categorical genotype) emission.

test_that("gt emission matrix is a valid categorical distribution (Holland's emimat)", {
  M <- nilHMM:::.gt_emimat(emission_gt(germ = 0.05, gert = 0.10, p = 0.5, mr = 0.10, nir = 0.01))
  expect_equal(dim(M), c(3L, 4L))                 # 3 states x {0,1,2,missing}
  expect_equal(unname(rowSums(M)), rep(1, 3), tolerance = 1e-12)   # each row a distribution
  expect_true(all(M >= 0))
})

test_that("the gt caller runs end-to-end and returns valid calls", {
  raw <- read_golden_counts(golden_slice_path("skim", "counts", "PN14_SID1259.chr1.tsv"))
  # nnil is the categorical gt caller: it consumes CALLED genotypes, not read
  # counts. Hard-call the counts explicitly here (the user's decision) to a `g`
  # column, then run nnil.
  g <- with(raw, ifelse(n_ref + n_alt == 0L, 3L,
              ifelse(n_alt == 0L, 0L, ifelse(n_ref == 0L, 2L, 1L))))
  gtin <- data.frame(name = "PN14_SID1259", chr = as.integer(sub("^chr", "", raw$chr)),
                     pos = raw$pos, g = as.integer(g), donor = "Zd")
  gt <- call_ancestry(gtin, "nnil", design = "BC2S2", rrate = 7e-6)
  expect_true(all(gt$state %in% 0:2))
  expect_true(all(gt$end_bp >= gt$start_bp))
  expect_identical(unique(gt$chr), 1L)
})

test_that("select_emission picks emission by depth regime", {
  expect_s3_class(select_emission(0.4),  "nilHMM_emission_count")
  expect_s3_class(select_emission(30),   "nilHMM_emission_gt")
})
