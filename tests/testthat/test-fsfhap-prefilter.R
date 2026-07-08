# FSFHap port — stage 2b preFilterSites step 1: filterSnpsByTag + computeRForMissingness.
# Same-tag thinning by presence-correlation (<0.7 within 64bp) + MAF/missing/het gates.

int_mat <- function(x) { storage.mode(x) <- "integer"; x }

test_that("filterSnpsByTag drops same-tag SNPs, keeps distant, gates on quality", {
  n <- 20L
  s1 <- c(rep(0L, 10), rep(2L, 10))     # MAF 0.5
  s2 <- s1                              # identical presence, 40bp from s1 -> same tag, drop
  s3 <- c(rep(2L, 8), rep(0L, 12))      # MAF 0.4, 300bp away -> kept
  s4 <- rep(0L, 20)                     # monomorphic (MAF 0) -> quality drop
  G <- int_mat(cbind(s1, s2, s3, s4)); pos <- as.integer(c(100, 140, 300, 360))
  keep <- fsfhap_filter_snps_by_tag_cpp(G, pos, min_maf = 0.05, max_missing = 0.8, max_het = 1.0)
  expect_equal(keep, c(TRUE, FALSE, TRUE, FALSE))
})

test_that("within-64bp but LOW presence-correlation SNPs are kept", {
  # each site biallelic (MAF 0.5) but present in DISJOINT taxon halves ->
  # anti-correlated missingness (r ~ -1 < 0.7) -> not same tag -> both kept
  s1 <- rep(3L, 40); s1[1:10]  <- 0L; s1[11:20] <- 2L   # present in taxa 1-20
  s2 <- rep(3L, 40); s2[21:30] <- 0L; s2[31:40] <- 2L   # present in taxa 21-40
  G <- int_mat(cbind(s1, s2)); pos <- as.integer(c(100, 150))   # 50bp apart (<64)
  keep <- fsfhap_filter_snps_by_tag_cpp(G, pos, min_maf = 0.05, max_missing = 0.9, max_het = 1.0)
  expect_true(keep[1])                                # head
  expect_true(keep[2])                                # anti-correlated presence -> not same tag
})

test_that("filterSnpsByTag gates out high-missing and high-het sites", {
  n <- 20L
  ok   <- c(rep(0L, 10), rep(2L, 10))                 # good
  miss <- c(rep(0L, 3), rep(2L, 3), rep(3L, 14))      # 70% missing
  het  <- c(rep(1L, 18), 0L, 2L)                      # ~90% het
  G <- int_mat(cbind(ok, miss, het)); pos <- as.integer(c(100, 3000, 6000))  # all far apart
  keep <- fsfhap_filter_snps_by_tag_cpp(G, pos, min_maf = 0.05, max_missing = 0.5, max_het = 0.5)
  expect_equal(keep, c(TRUE, FALSE, FALSE))
})

test_that("filterSnpsByTag validates pos length", {
  G <- int_mat(matrix(0L, 5, 3))
  expect_error(fsfhap_filter_snps_by_tag_cpp(G, 1:2, 0.05, 0.8, 1.0), "length\\(pos\\)")
})

# ---- full preFilterSites (filterSnpsByTag -> het-dev -> biallelic -> LD) -----

mk_pref <- function(n = 40L, nm = 60L, seed = 1L) {
  set.seed(seed)
  true <- matrix(0L, nm, n)
  for (j in seq_len(n)) { s0 <- sample(1:(nm - 15L), 1L); true[s0:(s0 + 14L), j] <- 2L }
  G <- int_mat(t(true)); attr(G, "pos") <- as.integer(seq_len(nm) * 1e5L)  # 100kb -> same-tag keeps all
  G
}

test_that("preFilterSites drops monomorphic and low-LD sites; LD-off keeps a superset", {
  G <- mk_pref(); pos <- attr(G, "pos")
  G[, 10] <- 0L                                       # monomorphic -> biallelic drop
  keep_ld  <- fsfhap_prefilter_sites_cpp(G, pos, 0.05, 0.2, 5, 0.2)
  keep_off <- fsfhap_prefilter_sites_cpp(G, pos, 0.05, 0.2, 5, 0.0)
  expect_false(keep_ld[10]); expect_false(keep_off[10])   # monomorphic dropped either way
  expect_gte(sum(keep_off), sum(keep_ld))                 # LD off keeps a superset
  expect_true(all(which(keep_ld) %in% which(keep_off)))
})

test_that("preFilterSites het-deviation drops a hugely-het site", {
  G <- mk_pref(); pos <- attr(G, "pos")
  G[, 25] <- 1L                                       # all-het site: pHet=1, far above mean
  keep <- fsfhap_prefilter_sites_cpp(G, pos, 0.05, 0.2, 5, 0.0)  # LD off to isolate het-dev
  expect_false(keep[25])
})

test_that("preFilterSites validates pos length", {
  G <- int_mat(matrix(0L, 5, 3))
  expect_error(fsfhap_prefilter_sites_cpp(G, 1:2, 0.05, 0.2, 5, 0.2), "length\\(pos\\)")
})
