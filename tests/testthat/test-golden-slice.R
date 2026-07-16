# Golden-slice fixture well-formedness. The decode-vs-expected regression
# (run the engine on these counts, compare to expected_calls_chr1.csv) lands in
# Task 4 once call_ancestry() is live -- see the skipped test at the bottom,
# which is the intended shape of that check.

test_that("golden-slice files are present", {
  for (rel in c("brb/counts/PN1_SID25.chr1.tsv", "brb/counts/PN1_SID32.chr1.tsv",
                "brb/counts/PN1_SID41.chr1.tsv", "brb/counts/PN1_SID11.chr1.tsv",
                "brb/expected_calls_chr1.csv",
                "skim/counts/PN14_SID1259.chr1.tsv", "skim/counts/PN3_SID235.chr1.tsv",
                "skim/counts/PN17_SID1630.chr1.tsv", "skim/expected_calls_chr1.csv")) {
    expect_true(file.exists(golden_slice_path(rel)), info = rel)
  }
})

test_that("count slices are chr1-only with the expected 6-column schema", {
  for (f in c(golden_slice_path("brb", "counts", "PN1_SID25.chr1.tsv"),
              golden_slice_path("skim", "counts", "PN14_SID1259.chr1.tsv"))) {
    d <- read_golden_counts(f)
    expect_named(d, c("chr", "pos", "ref", "n_ref", "alt", "n_alt"))
    expect_identical(unique(d$chr), "chr1")
    expect_true(all(d$n_ref >= 0) && all(d$n_alt >= 0))
    expect_gt(nrow(d), 1000L)
  }
})

test_that("expected calls use the common schema and only chr1", {
  for (f in c(golden_slice_path("brb", "expected_calls_chr1.csv"),
              golden_slice_path("skim", "expected_calls_chr1.csv"))) {
    calls <- read_golden_expected(f)
    expect_named(calls, c("source", "donor", "name", "chr", "start_bp", "end_bp", "state"))
    expect_identical(unique(calls$chr), 1L)
    expect_true(all(calls$state %in% c(0L, 1L, 2L)))
    expect_true(all(calls$end_bp >= calls$start_bp))
  }
})

test_that("brb slice carries donor (non-REF) signal to exercise HET/ALT", {
  calls <- read_golden_expected(golden_slice_path("brb", "expected_calls_chr1.csv"))
  expect_gt(sum(calls$state != 0L), 0L)
})

# Strict regression: the R count caller (fixed means) must be bit-identical to the
# frozen Python calls on chr1. The Python baseline decoded EVERY panel marker,
# including zero-coverage ones, so reproduction requires min_reads = 0L (the
# documented "decode every marker / old behaviour" flag); the live default of
# min_reads = 1L drops no-read markers and shifts segment borders on purpose
# (commit "min_cov: uniform no-coverage filtering across all callers"). Both skim
# and BRB used fixed emission means, so both reproduce exactly. (fit_means=TRUE,
# the §10 fix for BRB's ALT collapse, is a deliberate FUTURE divergence and is
# tested separately once implemented.)
seg_cols <- c("chr", "start_bp", "end_bp", "state")

reproduces <- function(src, sample, donor, design, r) {
  raw <- read_golden_counts(golden_slice_path(src, "counts", paste0(sample, ".chr1.tsv")))
  counts <- data.frame(name = sample, chr = as.integer(sub("^chr", "", raw$chr)),
                       pos = raw$pos, n_ref = raw$n_ref, n_alt = raw$n_alt, donor = donor)
  got <- call_ancestry(counts, caller = "bbnil", design = design, rrate = r, min_reads = 0L)
  exp <- read_golden_expected(golden_slice_path(src, "expected_calls_chr1.csv"))
  exp <- exp[exp$name == sample, ]
  expect_equal(got[, seg_cols], exp[, seg_cols], ignore_attr = TRUE, info = paste(src, sample))
}

test_that("count caller reproduces the skim baseline (BC2S2, fixed means)", {
  reproduces("skim", "PN14_SID1259", "Zd", "BC2S2", r = 7e-6)
  reproduces("skim", "PN3_SID235",   "Zx", "BC2S2", r = 7e-6)
  reproduces("skim", "PN17_SID1630", "Zx", "BC2S2", r = 7e-6)
})

test_that("count caller reproduces the BRB baseline (BC2S3, fixed means)", {
  reproduces("brb", "PN1_SID25", "Zd", "BC2S3", r = 1e-8)
  reproduces("brb", "PN1_SID32", "Zd", "BC2S3", r = 1e-8)
  reproduces("brb", "PN1_SID11", "Zv", "BC2S3", r = 1e-8)
})

test_that("threaded decode (parallel=TRUE) is identical to serial", {
  RcppParallel::setThreadOptions(numThreads = 2)
  on.exit(RcppParallel::setThreadOptions(numThreads = 1), add = TRUE)
  raw <- read_golden_counts(golden_slice_path("skim", "counts", "PN14_SID1259.chr1.tsv"))
  counts <- data.frame(name = "PN14_SID1259", chr = 1L, pos = raw$pos,
                       n_ref = raw$n_ref, n_alt = raw$n_alt, donor = "Zd")
  ser <- call_ancestry(counts, "bbnil", design = "BC2S2", rrate = 7e-6, parallel = FALSE)
  par <- call_ancestry(counts, "bbnil", design = "BC2S2", rrate = 7e-6, parallel = TRUE)
  expect_equal(ser, par)
})

test_that("call_ancestry rejects missing priors and bad columns", {
  d <- data.frame(name = "s", chr = 1L, pos = 1L, n_ref = 0L, n_alt = 0L)
  expect_error(call_ancestry(d, caller = "bbnil"), "supply")               # count caller, no design
  expect_error(call_ancestry(data.frame(x = 1), caller = "bbnil", design = "BC2S2"), "columns")
})
