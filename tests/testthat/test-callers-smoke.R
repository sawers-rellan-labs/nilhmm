# Smoke test: each grid caller (nnil, bbnil, catiger, rtiger) and the off-axis/
# wrapper callers (binhmm, googa, atlas) returns a valid common-schema segment
# table on the bundled golden-slice fixture. This is the "all callers run
# end-to-end" gate (plan A4); the strict per-marker regressions live in the
# caller-specific test files. (The no-HMM genotype baseline is call_gt(), a
# genotype caller -- not dispatchable here; see test-call_gt.R.)

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

test_that("the count-input grid and wrapper callers return valid calls", {
  counts <- smoke_counts()
  # grid cells: nnil (gt+geom, g derived from counts), bbnil (count+geom),
  # catiger (gt+rigidity), rtiger (count+rigidity)
  expect_valid_calls(
    call_ancestry(counts, caller = "nnil", design = "BC2S2", rrate = 7e-6), "nnil")
  expect_valid_calls(
    call_ancestry(counts, caller = "bbnil", design = "BC2S2", rrate = 7e-6), "bbnil")
  expect_valid_calls(
    call_ancestry(counts, caller = "catiger", design = "BC2S2", rigidity = 5L), "catiger")
  expect_valid_calls(
    call_ancestry(counts, caller = "rtiger", design = "BC2S2", rigidity = 5L, seed = 1L), "rtiger")
  # off-axis / wrapper; googa (gt + geometric) and atlas (gt + rigidity) share the
  # GOOGA competitive-alignment thresholding
  expect_valid_calls(
    call_ancestry(counts, caller = "binhmm", design = "BC2S2"), "binhmm")
  expect_valid_calls(
    call_ancestry(counts, caller = "googa", design = "BC2S2"), "googa")
  expect_valid_calls(
    call_ancestry(counts, caller = "atlas", design = "BC2S2", rigidity = 5L), "atlas")
})

test_that("the no-HMM per-site genotype baselines are call_gt(), not call_ancestry callers", {
  # ml/hwemap are GENOTYPE callers (call_gt), deliberately NOT ancestry callers --
  # the terminology wall between the ancestry mosaic and per-site genotypes.
  counts <- smoke_counts()
  expect_error(call_ancestry(counts, caller = "ml", design = "BC2S2"), "arg")
  expect_error(call_ancestry(counts, caller = "hwemap", design = "BC2S2"), "arg")
})

test_that("the gt callers (nnil, catiger) also run on a hard-genotype (g-only) input", {
  counts <- smoke_counts()
  g <- ifelse(counts$n_alt == 0L, 0L,
              ifelse(counts$n_ref == 0L, 2L, 1L))
  gt <- data.frame(name = counts$name, chr = counts$chr, pos = counts$pos, g = g)
  expect_valid_calls(
    call_ancestry(gt, caller = "nnil", design = "BC2S2"), "nnil-gt")
  expect_valid_calls(
    call_ancestry(gt, caller = "catiger", design = "BC2S2", rigidity = 5L), "catiger-gt")
})

# The het-blind (flat) vs het-excess (HWE) genotype-calling behaviour lives in
# test-call_gt.R -- it is a genotype-caller property, not an ancestry-caller one.
