# FastIndep port: maximal-independent-set selection. The greedy heuristic is
# deterministic and RNG-free, so it must reproduce the FastIndep CLI exactly
# (cross-checked when the reference binary is present). The stochastic runs use a
# self-contained PRNG (not the CLI's Mersenne Twister); we test their invariants
# (validity, determinism) rather than CLI bit-identity -- see the package docs.

# Locate the reference FastIndep binary, if installed (for the CLI cross-check).
# Discovery is env/PATH-based so the test stays portable; set FASTINDEP_BIN (e.g.
# in ~/.Renviron) to point at a local build, otherwise the cross-check skips.
fastindep_bin <- function() {
  cand <- c(Sys.getenv("FASTINDEP_BIN"), unname(Sys.which("fastindep")))
  cand <- cand[nzchar(cand)]
  hit <- cand[file.exists(cand)]
  if (length(hit)) hit[1] else NA_character_
}

# Write a symmetric similarity matrix in FastIndep's text format: line 1 = names,
# each row = name + values, diagonal forced to 10000 (the ignored sentinel).
write_fastindep <- function(sim, path) {
  nm <- rownames(sim)
  con <- file(path, "w")
  on.exit(close(con))
  writeLines(paste(nm, collapse = " "), con)
  for (i in seq_len(nrow(sim))) {
    row <- sim[i, ]
    row[i] <- 10000
    writeLines(paste(c(nm[i], formatC(row, format = "g", digits = 8)), collapse = " "), con)
  }
}

# Parse the first "Set Size =" row of a FastIndep output file (the deterministic
# greedy set, always printed first), returning the member names.
parse_greedy <- function(out_path, names_pool) {
  lines <- readLines(out_path)
  hit <- grep("Set Size =", lines, value = TRUE)
  if (!length(hit)) return(character(0))
  toks <- strsplit(trimws(hit[1]), "[[:space:]]+")[[1]]
  intersect(toks, names_pool)
}

test_that("greedy set is bit-identical to the FastIndep CLI", {
  bin <- fastindep_bin()
  skip_if(is.na(bin), "FastIndep reference binary not found")

  set.seed(11)
  n <- 18
  nm <- sprintf("m%02d", seq_len(n))
  A <- matrix(runif(n * n), n, n)
  sim <- (A + t(A)) / 2                       # symmetric, values in (0,1)
  dimnames(sim) <- list(nm, nm)

  thr <- 0.5
  dir <- withr::local_tempdir()
  inp <- file.path(dir, "in.txt"); outp <- file.path(dir, "out.txt")
  write_fastindep(sim, inp)
  # -n 5 so the CLI prints set membership (n=1 prints only the size).
  status <- system2(bin, c("-t", thr, "-n", "5", "-s", "1", "-i", inp, "-o", outp),
                    stdout = FALSE, stderr = FALSE)
  skip_if_not(file.exists(outp), "FastIndep produced no output")

  cli_greedy <- sort(parse_greedy(outp, nm))
  got <- select_independent(sim, threshold = thr, n_runs = 1L)  # greedy only
  expect_gt(length(cli_greedy), 0L)
  expect_identical(sort(got), cli_greedy)
})

test_that("every returned set is independent (no pair crosses the threshold)", {
  set.seed(2)
  n <- 25
  A <- matrix(runif(n * n), n, n); sim <- (A + t(A)) / 2
  thr <- 0.6
  res <- fast_indep_cpp(sim, thr, 30L, 7L, FALSE)
  edge <- sim >= thr; diag(edge) <- FALSE
  for (s in res$sets) {
    expect_false(any(edge[s, s, drop = FALSE]),
                 info = paste("set", paste(s, collapse = ",")))
  }
})

test_that("distance sense (vi) inverts the edge rule correctly", {
  # A distance matrix: small = related, edge if <= threshold.
  n <- 12
  set.seed(3)
  A <- matrix(runif(n * n), n, n); d <- (A + t(A)) / 2; diag(d) <- 0
  thr <- 0.4
  res <- fast_indep_cpp(d, thr, 1L, 1L, TRUE)   # distance = TRUE
  edge <- d <= thr; diag(edge) <- FALSE
  expect_false(any(edge[res$best, res$best, drop = FALSE]))
})

test_that("known graph: maximum independent set is found", {
  # 6-cycle 1-2-3-4-5-6-1: the two maximum independent sets are {1,3,5},{2,4,6}.
  n <- 6
  sim <- matrix(0, n, n)
  edges <- rbind(c(1,2), c(2,3), c(3,4), c(4,5), c(5,6), c(6,1))
  for (k in seq_len(nrow(edges))) { i <- edges[k,1]; j <- edges[k,2]; sim[i,j] <- sim[j,i] <- 1 }
  res <- fast_indep_cpp(sim, 0.5, 1L, 1L, FALSE)
  expect_length(res$best, 3L)
  # An alternating triple of the cycle
  expect_true(setequal(res$best, c(1,3,5)) || setequal(res$best, c(2,4,6)))
})

test_that("sets are distinct: greedy is not duplicated when also found stochastically", {
  # 6-cycle: the only two maximum independent sets are {1,3,5} and {2,4,6}; with
  # many runs the greedy set is also discovered stochastically, exercising the
  # dedup path that must not list it twice.
  n <- 6
  sim <- matrix(0, n, n)
  edges <- rbind(c(1,2), c(2,3), c(3,4), c(4,5), c(5,6), c(6,1))
  for (k in seq_len(nrow(edges))) { i <- edges[k,1]; j <- edges[k,2]; sim[i,j] <- sim[j,i] <- 1 }
  res <- fast_indep_cpp(sim, 0.5, 100L, 1L, FALSE)
  keys <- vapply(res$sets, function(s) paste(sort(s), collapse = ","), character(1))
  expect_identical(anyDuplicated(keys), 0L)                 # no set listed twice
  expect_equal(length(res$sets), sum(res$size_dist))        # sets count == distinct total
})

test_that("singletons (unrelated to all) always join the set", {
  # Node 5 relates to nobody (a singleton) -> must be in every set. Nodes 1-4
  # form a clique (all related). Best set = {5} + one clique member.
  n <- 5
  sim <- matrix(0, n, n)
  for (i in 1:4) for (j in 1:4) if (i != j) sim[i, j] <- 1  # 1-4 mutually related
  res <- fast_indep_cpp(sim, 0.5, 1L, 1L, FALSE)
  expect_true(5 %in% res$best)
  expect_length(res$best, 2L)
})

test_that("determinism: same seed reproduces, different seed may differ", {
  set.seed(4)
  n <- 40
  A <- matrix(runif(n * n), n, n); sim <- (A + t(A)) / 2
  a <- fast_indep_cpp(sim, 0.5, 50L, 123L, FALSE)
  b <- fast_indep_cpp(sim, 0.5, 50L, 123L, FALSE)
  expect_identical(a$sets, b$sets)
  expect_identical(a$best, b$best)
})

test_that("greedy set is seed-independent (RNG-free)", {
  set.seed(5)
  n <- 30
  A <- matrix(runif(n * n), n, n); sim <- (A + t(A)) / 2
  g1 <- fast_indep_cpp(sim, 0.5, 1L, 1L, FALSE)$best
  g2 <- fast_indep_cpp(sim, 0.5, 1L, 99999L, FALSE)$best
  expect_identical(g1, g2)
})

test_that("select_independent returns marker names and carries attributes", {
  set.seed(6)
  geno <- matrix(sample(0:2, 10 * 40, replace = TRUE), nrow = 10,
                 dimnames = list(paste0("snp", 1:10), NULL))
  sel <- select_independent(geno, threshold = 0.5, n_runs = 5L, method = "r2")
  expect_type(sel, "character")
  expect_true(all(sel %in% rownames(geno)))
  expect_true(!is.null(attr(sel, "size_dist")))
  expect_true(is.list(attr(sel, "sets")))
  expect_identical(attr(sel, "kind"), "similarity")
})

test_that("marker-count guard errors cleanly before any large allocation", {
  # Direct checks on the budget helper (no matrix is allocated).
  expect_error(.check_marker_budget(8000, 7000), "exceeds `max_markers`")
  expect_error(.check_marker_budget(100, 40000), "safety ceiling")
  expect_error(.check_marker_budget(100, 0), "positive integer")
  expect_silent(.check_marker_budget(7000, 7000))            # at the default: fine
  # Large but allowed (max raised): warns about the O(n^2) size, still no alloc.
  expect_warning(.check_marker_budget(12000, 20000), "O\\(n\\^2\\)")
})

test_that("select_independent enforces max_markers (both paths)", {
  # Precomputed square matrix path.
  sim <- matrix(0.1, 4, 4); diag(sim) <- 1
  expect_error(select_independent(sim, threshold = 0.5, max_markers = 3),
               "exceeds `max_markers`")
  # Genotype path.
  geno <- matrix(sample(0:2, 6 * 10, replace = TRUE), nrow = 6)
  expect_error(select_independent(geno, threshold = 0.5, max_markers = 4),
               "exceeds `max_markers`")
  # Default (7000) admits a small input.
  expect_silent(sel <- select_independent(geno, threshold = 0.5))
})

test_that("hard-cap ceiling is honoured and option-overridable", {
  sim <- matrix(0.1, 5, 5); diag(sim) <- 1
  withr::local_options(nilHMM.marker_hard_cap = 3L)
  # Default max_markers (7000) now exceeds the lowered ceiling -> refuse.
  expect_error(select_independent(sim, threshold = 0.5), "safety ceiling")
  # Even a permissible count is refused because max_markers > ceiling.
  expect_error(select_independent(sim, threshold = 0.5, max_markers = 5),
               "safety ceiling")
})

test_that("scales to a large matrix in reasonable time", {
  set.seed(7)
  n <- 2000
  A <- matrix(runif(n * n), n, n); sim <- (A + t(A)) / 2
  t0 <- Sys.time()
  res <- fast_indep_cpp(sim, 0.5, 5L, 1L, FALSE)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  expect_gt(length(res$best), 0L)
  expect_lt(elapsed, 20)   # generous ceiling; typically well under a second
})
