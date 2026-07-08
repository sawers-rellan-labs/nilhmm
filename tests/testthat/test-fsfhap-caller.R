# FSFHap port — caller wiring: call_ancestry(caller = "fsfhap") end to end.
# Checks the family x chr driver, the design-routing dispatcher (design -> route +
# phet), and the common REF/HET/ALT segment schema. Per-cell correctness is the
# e2e parity test; here we check integration, schema, and guards.

# self-contained BC1 family. NOTE: needs enough markers that taxa survive the
# TASSEL-faithful min_gametes=200 coverage filter (>100 kept sites/taxon), so use
# a realistic marker count. Each taxon = A-hom with a random donor(C) block.
mk_long <- function(ntaxa = 30L, nsites = 300L, seed = 3L) {
  set.seed(seed)
  taxa <- sprintf("T%03d", seq_len(ntaxa)); pos <- as.integer(seq_len(nsites) * 1e5L)
  g <- matrix(0L, ntaxa, nsites)
  for (t in seq_len(ntaxa)) { s0 <- sample(1:(nsites - 60L), 1L); g[t, s0:(s0 + 50L)] <- 2L }
  g[cbind(sample(ntaxa, 20), sample(nsites, 20))] <- 3L        # sprinkle missing
  data.frame(name = rep(taxa, times = nsites), family = "F1", chr = 1L,
             pos = rep(pos, each = ntaxa), g = as.integer(g), stringsAsFactors = FALSE)
}

test_that("call_ancestry(caller='fsfhap') returns the common segment schema", {
  seg <- call_ancestry(mk_long(), caller = "fsfhap", design = "BC1S4")
  expect_named(seg, c("source", "donor", "name", "chr", "start_bp", "end_bp", "state"))
  expect_true(all(seg$state %in% c(0L, 1L, 2L)))
  expect_true(nrow(seg) > 0)
  expect_true(all(seg$end_bp >= seg$start_bp))
})

test_that("phet is design-derived (BC1S4 -> 0.03125); explicit phet overrides", {
  d <- mk_long()
  # both run without error; design derives phet, explicit phet is accepted
  expect_silent(s1 <- call_ancestry(d, caller = "fsfhap", design = "BC1S4"))
  expect_silent(s2 <- call_ancestry(d, caller = "fsfhap", phet = 0.03125))
  expect_identical(s1, s2)                                     # same phet -> same calls
})

test_that("dispatcher routes BC1 -> bc, non-BC1 -> biparental; malformed design errors", {
  d <- mk_long()
  expect_silent(call_ancestry(d, caller = "fsfhap", design = "BC1S4"))   # backcross route
  expect_silent(call_ancestry(d, caller = "fsfhap", design = "BC2S2"))   # biparental route (no error)
  expect_error(call_ancestry(d, caller = "fsfhap", design = "junk"), "BC\\{n\\}S\\{m\\}")
})

test_that("biparental route (BC2) returns the common schema on a recombinant family", {
  set.seed(5); nt <- 40L; ns <- 400L; g <- matrix(0L, nt, ns)
  for (t in seq_len(nt)) {                              # RIL-like: alternating hom-A/hom-B blocks
    st <- 1L; cur <- sample(c(0L, 2L), 1L); bps <- sort(sample(20:380, sample(2:5, 1L)))
    for (bp in c(bps, ns)) { g[t, st:bp] <- cur; cur <- ifelse(cur == 0L, 2L, 0L); st <- bp + 1L }
  }
  dat <- data.frame(name = rep(sprintf("T%03d", seq_len(nt)), times = ns), family = "F",
                    chr = 1L, pos = rep(as.integer(seq_len(ns) * 1e5L), each = nt), g = as.integer(g))
  seg <- call_ancestry(dat, caller = "fsfhap", design = "BC2S2")   # -> biparental route
  expect_named(seg, c("source", "donor", "name", "chr", "start_bp", "end_bp", "state"))
  expect_true(nrow(seg) > 0)
  expect_true(all(seg$state %in% c(0L, 1L, 2L)))
})

test_that("fsfhap needs a g column, a family, and a design/phet", {
  d <- mk_long()
  expect_error(call_ancestry(d[, c("name", "family", "chr", "pos")], caller = "fsfhap", design = "BC1S4"),
               "needs a `g`|read counts")
  expect_error(call_ancestry(d, caller = "fsfhap"), "needs `design`")
  expect_error(call_ancestry(d[, setdiff(names(d), "family")], caller = "fsfhap", design = "BC1S4"),
               "family")
})

test_that("family can be supplied as a vector argument", {
  d <- mk_long(); d$family <- NULL
  fam <- rep("F1", nrow(d))
  expect_silent(seg <- call_ancestry(d, caller = "fsfhap", design = "BC1S4", family = fam))
  expect_true(nrow(seg) > 0)
})
