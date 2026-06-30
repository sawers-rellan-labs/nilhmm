# gt (categorical genotype) and dosage (Gaussian) emissions.

test_that("gt emission matrix is a valid categorical distribution (Holland's emimat)", {
  M <- nilHMM:::.gt_emimat(emission_gt(germ = 0.05, gert = 0.10, p = 0.5, mr = 0.10, nir = 0.01))
  expect_equal(dim(M), c(3L, 4L))                 # 3 states x {0,1,2,missing}
  expect_equal(unname(rowSums(M)), rep(1, 3), tolerance = 1e-12)   # each row a distribution
  expect_true(all(M >= 0))
})

test_that("dosage emission peaks at the matching state; missing is flat", {
  em <- emission_dosage(sd_dosage = 0.25)
  ll <- nilHMM:::.emission_loglik(em, list(d = c(0, 1, 2, NA)), nilHMM:::.emission_theta(em))
  expect_equal(which.max(ll[1, ]), 1L)            # d=0 -> REF
  expect_equal(which.max(ll[2, ]), 2L)            # d=1 -> HET
  expect_equal(which.max(ll[3, ]), 3L)            # d=2 -> ALT
  expect_equal(unname(ll[4, ]), c(0, 0, 0))       # missing -> flat
})

test_that("gt and dosage callers run end-to-end and return valid calls", {
  raw <- read_golden_counts(golden_slice_path("skim", "counts", "PN14_SID1259.chr1.tsv"))
  counts <- data.frame(name = "PN14_SID1259", chr = as.integer(sub("^chr", "", raw$chr)),
                       pos = raw$pos, n_ref = raw$n_ref, n_alt = raw$n_alt, donor = "Zd")
  gt  <- call_ancestry(counts, "nnil", design = "BC2S2", r = 7e-6, emission = "gt")
  dos <- call_ancestry(counts, "skimbin", design = "BC2S2", r = 7e-6)   # dosage emission
  for (g in list(gt, dos)) {
    expect_true(all(g$state %in% 0:2))
    expect_true(all(g$end_bp >= g$start_bp))
    expect_identical(unique(g$chr), 1L)
  }
})

test_that("select_emission picks emission by depth regime", {
  expect_s3_class(select_emission(0.4),  "nilHMM_emission_count")
  expect_s3_class(select_emission(30),   "nilHMM_emission_gt")
  expect_s3_class(select_emission(0.4, imputed = TRUE), "nilHMM_emission_dosage")
})
