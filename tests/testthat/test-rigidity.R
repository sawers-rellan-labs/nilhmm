# Two distinct things are tested here:
#  (a) the engine's generic `duration_rigidity` option (min-run via state
#      expansion) — a standalone HMM duration, NOT RTIGER;
#  (b) the `rtiger` caller, which is the faithful RTIGER port (src/rtiger.cpp,
#      R/rtiger.R) — its own EM/Viterbi, validated distributionally against
#      RTIGER's calls_taxa_r5 elsewhere (agent/rtiger_fulltaxa.R). Here we only
#      smoke-test that the caller runs and is deterministic.

rig_counts <- function(sample = "PN14_SID1259", donor = "Zd") {
  raw <- read_golden_counts(golden_slice_path("skim", "counts", paste0(sample, ".chr1.tsv")))
  data.frame(name = sample, chr = as.integer(sub("^chr", "", raw$chr)),
             pos = raw$pos, n_ref = raw$n_ref, n_alt = raw$n_alt, donor = donor)
}

# (a) engine duration_rigidity --------------------------------------------------
test_that("duration_rigidity expands each macro-state into r sub-states", {
  model <- nilHMM:::fit(list(n = 1L, a = 0L), emission_count(),
                        duration_rigidity(5L, 0.01), design_priors("BC2S2"))
  expect_equal(model$n_sub, 5L)
  expect_equal(length(model$log_start), 3L * 5L)   # 3 macro-states x r
})

test_that("duration_rigidity enforces the minimum run length (interior segments)", {
  counts <- rig_counts()
  obs   <- list(n = counts$n_ref + counts$n_alt, a = counts$n_alt)
  model <- nilHMM:::fit(obs, emission_count(), duration_rigidity(8L, 0.01), design_priors("BC2S2"))
  path  <- nilHMM:::decode(model, obs)
  expect_true(all(path %in% 0:2))
  runs <- rle(path)$lengths
  interior <- runs[-c(1, length(runs))]            # boundary runs may be truncated
  if (length(interior)) expect_gte(min(interior), 8L)
})

# (b) faithful RTIGER caller ----------------------------------------------------
test_that("rtiger caller runs end-to-end and returns valid common-schema calls", {
  counts <- rig_counts()
  got <- call_ancestry(counts, caller = "rtiger", r = 5L)
  expect_named(got, c("source", "donor", "name", "chr", "start_bp", "end_bp", "state"))
  expect_true(all(got$state %in% 0:2))             # pat/het/mat -> REF/HET/ALT
  expect_true(all(got$end_bp >= got$start_bp))
  expect_identical(unique(got$chr), 1L)
})

test_that("rtiger caller is deterministic for a fixed seed", {
  counts <- rig_counts()
  a <- call_ancestry(counts, caller = "rtiger", r = 5L, seed = 1L)
  b <- call_ancestry(counts, caller = "rtiger", r = 5L, seed = 1L)
  expect_equal(a, b)
})
