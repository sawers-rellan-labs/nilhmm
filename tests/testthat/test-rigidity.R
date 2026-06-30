# Rigidity duration (RTIGER mode, §7) — structural property tests. Exact
# distributional concordance with RTIGER's calls_taxa_r5.csv is the §9.4
# (Task 5) job; RTIGER EM-fits emissions in Julia, so this engine is validated
# against it distributionally, not bit-identically.

# A skim sample's chr1 counts, as a call_ancestry input frame.
rig_counts <- function(sample = "PN14_SID1259", donor = "Zd") {
  raw <- read_golden_counts(golden_slice_path("skim", "counts", paste0(sample, ".chr1.tsv")))
  data.frame(name = sample, chr = as.integer(sub("^chr", "", raw$chr)),
             pos = raw$pos, n_ref = raw$n_ref, n_alt = raw$n_alt, donor = donor)
}

test_that("rigidity expands each macro-state into r sub-states", {
  spec  <- caller_spec("rtiger", r = 5L, p_switch = 0.01)
  model <- nilHMM:::fit(list(n = 1L, a = 0L), spec$emission, spec$duration,
                        design_priors("BC2S2"))
  expect_equal(model$n_sub, 5L)
  expect_equal(length(model$log_start), 3L * 5L)   # 3 macro-states x r
})

test_that("r = 1 rigidity equals a plain geometric transition", {
  counts <- rig_counts()
  geo <- call_ancestry(counts, "nnil",   design = "BC2S2", r = 0.01)
  rig <- call_ancestry(counts, "rtiger", design = "BC2S2", r = 1L, p_switch = 0.01)
  expect_equal(geo[, c("start_bp", "end_bp", "state")],
               rig[, c("start_bp", "end_bp", "state")], ignore_attr = TRUE)
})

test_that("rigidity enforces the minimum run length (interior segments)", {
  counts <- rig_counts()
  spec  <- caller_spec("rtiger", r = 8L, p_switch = 0.01)
  obs   <- list(n = counts$n_ref + counts$n_alt, a = counts$n_alt)
  model <- nilHMM:::fit(obs, spec$emission, spec$duration, design_priors("BC2S2"))
  path  <- nilHMM:::decode(model, obs)
  expect_true(all(path %in% 0:2))                  # mapped back to macro-states
  runs <- rle(path)$lengths
  interior <- runs[-c(1, length(runs))]            # boundary runs may be truncated
  if (length(interior)) expect_gte(min(interior), 8L)
})

test_that("higher rigidity yields no more segments than lower", {
  counts <- rig_counts()
  n_seg <- function(r) nrow(call_ancestry(counts, "rtiger", design = "BC2S2",
                                          r = r, p_switch = 0.01))
  expect_lte(n_seg(20L), n_seg(2L))
})
