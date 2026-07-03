# Smoke test: each of the four callers (nnil, rtiger, binhmm, atlas) returns a
# valid common-schema segment table on the bundled golden-slice fixture. This is
# the "all four callers run end-to-end" gate (plan A4); the strict per-marker
# regressions live in the caller-specific test files.

schema_cols <- c("source", "donor", "name", "chr", "start_bp", "end_bp", "state")

# One real chr1 skim sample (biallelic AD counts, ~7k markers) as (name, chr,
# pos, n_ref, n_alt). Enough markers/positions for binning + EM to run.
smoke_counts <- function() {
  raw <- read_golden_counts(golden_slice_path("skim", "counts", "PN14_SID1259.chr1.tsv"))
  data.frame(name = "PN14_SID1259", chr = as.integer(sub("^chr", "", raw$chr)),
             pos = raw$pos, n_ref = raw$n_ref, n_alt = raw$n_alt, donor = "Zd",
             stringsAsFactors = FALSE)
}

expect_valid_calls <- function(calls, info) {
  expect_s3_class(calls, "data.frame")
  expect_true(all(schema_cols %in% names(calls)), info = info)
  expect_gt(nrow(calls), 0L)
  expect_true(all(calls$state %in% c(0L, 1L, 2L)), info = info)
  expect_true(all(calls$end_bp >= calls$start_bp), info = info)
  expect_true(all(calls$chr == 1L), info = info)
}

test_that("all four callers return valid common-schema calls", {
  counts <- smoke_counts()
  expect_valid_calls(
    call_ancestry(counts, caller = "nnil", design = "BC2S2", r = 7e-6),
    "nnil")
  expect_valid_calls(
    call_ancestry(counts, caller = "rtiger", design = "BC2S2", r = 5L, seed = 1L),
    "rtiger")
  expect_valid_calls(
    call_ancestry(counts, caller = "binhmm", design = "BC2S2"),
    "binhmm")
  expect_valid_calls(
    call_ancestry(counts, caller = "atlas", design = "BC2S2"),
    "atlas")
})

test_that("nnil also runs on a hard-genotype (g-only) input via the gt emission", {
  counts <- smoke_counts()
  g <- ifelse(counts$n_alt == 0L, 0L,
              ifelse(counts$n_ref == 0L, 2L, 1L))
  gt <- data.frame(name = counts$name, chr = counts$chr, pos = counts$pos, g = g)
  expect_valid_calls(
    call_ancestry(gt, caller = "nnil", design = "BC2S2"),
    "nnil-gt")
})
