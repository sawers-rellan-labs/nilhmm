# Smoke test: each grid caller (nnil, bbnil, catiger, rtiger), the off-axis/
# wrapper callers (binhmm, atlas), and the no-HMM baselines (ml, hwemap) returns
# a valid common-schema segment table on the bundled golden-slice fixture. This is
# the "all callers run end-to-end" gate (plan A4); the strict per-marker
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

test_that("the count-input grid, wrapper, and baseline callers return valid calls", {
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
  # no-HMM per-site baselines
  expect_valid_calls(
    call_ancestry(counts, caller = "ml", design = "BC2S2"), "ml")
  expect_valid_calls(
    call_ancestry(counts, caller = "hwemap", design = "BC2S2"), "hwemap")
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

test_that("ml is het-blind and hwemap is het-excess at low depth", {
  # A 2-sample cohort, one REF read and one ALT read per site -> pooled allele
  # frequency 0.5. At depth 1, ml (flat prior) makes a HOMOZYGOUS call from the
  # lone read (het-blind); hwemap (HWE prior, p=0.5) flips every single-read site
  # to HET (the het-excess the HMM callers must beat).
  d <- rbind(
    data.frame(name = "altread", chr = 1L, pos = 1:4, n_ref = 0L, n_alt = 1L, donor = "Zx"),
    data.frame(name = "refread", chr = 1L, pos = 1:4, n_ref = 1L, n_alt = 0L, donor = "Zx"))
  ml <- call_ancestry(d, caller = "ml", design = "BC2S2")
  hw <- call_ancestry(d, caller = "hwemap", design = "BC2S2")
  expect_false(any(ml$state == 1L))                 # het-blind: no HET calls
  expect_setequal(unique(ml$state), c(0L, 2L))      # lone read -> homozygous
  expect_true(all(hw$state == 1L))                  # het-excess: every site -> HET
})
