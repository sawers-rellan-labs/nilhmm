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

test_that("[Task 4] engine reproduces the skim golden slice", {
  skip("call_ancestry() not yet implemented (Task 4); strict concordance target")
  # counts <- read_golden_counts(golden_slice_path("skim","counts","PN14_SID1259.chr1.tsv"))
  # got    <- call_ancestry(counts, caller = "nnil", design = "BC2S2")
  # exp    <- read_golden_expected(golden_slice_path("skim","expected_calls_chr1.csv"))
  # expect concordance(got, exp[exp$name=="PN14_SID1259",]) above tolerance
})
