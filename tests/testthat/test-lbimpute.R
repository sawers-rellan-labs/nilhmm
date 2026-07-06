# LB-Impute caller (native port of Fragoso et al. 2014): coverage-aware emission
# + distance-dependent transition (double-recomb penalty), full-chromosome Viterbi.

test_that("lb_emission is bounded to [errg, 1-errg] and flat at zero coverage", {
  em <- lb_emission_loglik_cpp(nref = c(10L, 0L, 0L), nalt = c(0L, 10L, 0L),
                               err = 0.05, errg = 0.05)
  p <- exp(em)
  expect_true(all(p >= 0.05 - 1e-9 & p <= 0.95 + 1e-9))
  # strong REF evidence -> REF is the max; strong ALT evidence -> ALT is the max
  expect_equal(which.max(em[1, ]), 1L)
  expect_equal(which.max(em[2, ]), 3L)
  # zero coverage -> all three states equal (flat emission)
  expect_equal(em[3, 1], em[3, 2]); expect_equal(em[3, 2], em[3, 3])
})

test_that("lb_emission favours HET at balanced depth", {
  em <- lb_emission_loglik_cpp(nref = 20L, nalt = 20L, err = 0.05, errg = 0.05)
  expect_equal(which.max(em[1, ]), 2L)   # HET
})

test_that("lb_viterbi recovers a REF -> ALT -> REF block", {
  pos <- seq_len(60L) * 1e5L
  nref <- c(rep(10L, 20), rep(0L, 20), rep(10L, 20))
  nalt <- c(rep(0L, 20), rep(10L, 20), rep(0L, 20))
  em <- lb_emission_loglik_cpp(nref, nalt, err = 0.05, errg = 0.05)
  path <- lb_viterbi_cpp(log(rep(1/3, 3)), em, as.integer(pos),
                         recombdist = 1e7, drp = FALSE)
  expect_equal(path[10], 0L); expect_equal(path[30], 2L); expect_equal(path[50], 0L)
  expect_equal(sort(unique(path)), c(0L, 2L))   # no spurious HET
})

test_that("lb_viterbi errors on unsorted positions and bad dims", {
  em <- lb_emission_loglik_cpp(c(5L, 5L), c(0L, 0L), 0.05, 0.05)
  expect_error(lb_viterbi_cpp(log(rep(1/3, 3)), em, c(2, 1), 1e7, FALSE), "non-decreasing")
  expect_error(lb_viterbi_cpp(log(rep(1/3, 3)), em, 1, 1e7, FALSE), "length")
})

test_that("call_ancestry(caller='lbimpute') returns the common segment schema", {
  toy <- data.frame(
    name = "NIL1", chr = 1L, pos = seq_len(40L) * 1e5L,
    n_ref = c(rep(10L, 20), rep(0L, 20)),
    n_alt = c(rep(0L, 20), rep(10L, 20)))
  seg <- call_ancestry(toy, caller = "lbimpute", err = 0.05,
                       genotypeerr = 0.05, recombdist = 1e7)
  expect_named(seg, c("source", "donor", "name", "chr", "start_bp", "end_bp", "state"))
  # a clean REF half then ALT half -> two segments
  expect_equal(nrow(seg), 2L)
  expect_equal(seg$state, c(0L, 2L))
})

test_that("lbimpute cM (map-aware) transition: bp-equivalent map reproduces bp calls", {
  # A perfectly uniform 5 cM/Mb map makes cM a linear rescale of bp, so with a
  # matched recombdist (1e7 bp <-> 50 cM) the cM path must equal the bp path.
  pos <- seq_len(40L) * 1e5L
  toy <- data.frame(name = "NIL1", chr = 1L, pos = pos, cm = pos * 5e-6,
                    n_ref = c(rep(10L, 20), rep(0L, 20)),
                    n_alt = c(rep(0L, 20), rep(10L, 20)))
  bp <- call_ancestry(toy, caller = "lbimpute", recombdist = 1e7, unit = "bp")
  cm <- call_ancestry(toy, caller = "lbimpute", recombdist = 50, unit = "cm")
  expect_equal(bp$state, cm$state)
  expect_equal(bp$start_bp, cm$start_bp)   # output stays in bp
})

test_that("lbimpute cM path differs where the map is non-uniform (centromere)", {
  # Two markers physically far apart but at the SAME cM (recombination-suppressed)
  # should NOT be allowed to switch under the cM transition (d_cm = 0 forbids it),
  # whereas the bp transition sees a large gap and permits the switch.
  pos <- as.integer(c(1e6, 5e7, 9e7))               # wide bp spacing
  cm  <- c(10, 10, 10)                              # flat cM: no recombination here
  toy <- data.frame(name = "NIL1", chr = 1L, pos = pos, cm = cm,
                    n_ref = c(10L, 0L, 10L), n_alt = c(0L, 10L, 0L))
  cmcall <- call_ancestry(toy, caller = "lbimpute", unit = "cm", recombdist = 50)
  # zero cM gap -> transition forbidden -> single homozygous block, no ALT dip
  expect_equal(nrow(cmcall), 1L)
  expect_false(2L %in% cmcall$state)
})

test_that("lbimpute unit validation errors and warns", {
  base <- data.frame(name = "NIL1", chr = 1L, pos = c(1e5, 2e5),
                     n_ref = c(5L, 5L), n_alt = c(0L, 0L))
  # unit = 'cm' without a cm column
  expect_error(call_ancestry(base, caller = "lbimpute", unit = "cm"), "needs a `cm`")
  # fractional bp -> looks like cM
  frac <- transform(base, pos = c(1.5, 2.5))
  expect_error(call_ancestry(frac, caller = "lbimpute", unit = "bp"), "whole numbers")
  # bp coords with a cM-sized recombdist -> over-fragment warning
  expect_warning(call_ancestry(base, caller = "lbimpute", unit = "bp", recombdist = 50),
                 "over-relax")
  # cM coords with a bp-sized recombdist -> collapse warning
  cmdat <- data.frame(name = "NIL1", chr = 1L, pos = c(1e5, 2e5), cm = c(0.5, 1.5),
                      n_ref = c(5L, 5L), n_alt = c(0L, 0L))
  expect_warning(call_ancestry(cmdat, caller = "lbimpute", unit = "cm", recombdist = 1e7),
                 "bp-sized")
})

test_that("lbimpute default recombdist is unit-aware", {
  # cM default (50) must NOT trip the bp-sized warning; bp default (1e7) must not warn.
  cmdat <- data.frame(name = "NIL1", chr = 1L, pos = c(1e5, 2e5), cm = c(0.5, 1.5),
                      n_ref = c(5L, 5L), n_alt = c(0L, 0L))
  expect_no_warning(call_ancestry(cmdat, caller = "lbimpute", unit = "cm"))
  bpdat <- data.frame(name = "NIL1", chr = 1L, pos = c(1e5, 2e5),
                      n_ref = c(5L, 5L), n_alt = c(0L, 0L))
  expect_no_warning(call_ancestry(bpdat, caller = "lbimpute", unit = "bp"))
})

test_that("lbimpute works without a design (flat start) and honours donor column", {
  toy <- data.frame(
    name = c(rep("A", 10), rep("B", 10)),
    donor = c(rep("teo1", 10), rep("teo2", 10)),
    chr = 1L, pos = rep(seq_len(10L) * 1e5L, 2),
    n_ref = c(rep(8L, 10), rep(0L, 10)),
    n_alt = c(rep(0L, 10), rep(8L, 10)))
  seg <- call_ancestry(toy, caller = "lbimpute")   # no design supplied
  expect_true(all(c("teo1", "teo2") %in% seg$donor))
  expect_equal(seg$state[seg$name == "A"], 0L)
  expect_equal(seg$state[seg$name == "B"], 2L)
})

test_that("drp=TRUE eases homozygous->homozygous switching vs drp=FALSE", {
  # a single-marker ALT dip inside a REF run at wide spacing: the double-recomb
  # penalty (drp=FALSE) should be at least as reluctant to switch as drp=TRUE.
  pos <- seq_len(11L) * 5e6L
  nref <- rep(6L, 11); nalt <- rep(0L, 11)
  nref[6] <- 0L; nalt[6] <- 6L
  em <- lb_emission_loglik_cpp(nref, nalt, 0.05, 0.05)
  p_pen  <- lb_viterbi_cpp(log(rep(1/3, 3)), em, as.integer(pos), 1e7, FALSE)
  p_free <- lb_viterbi_cpp(log(rep(1/3, 3)), em, as.integer(pos), 1e7, TRUE)
  expect_true(sum(p_free != 0L) >= sum(p_pen != 0L))
})

test_that("write_vcf_impute emits a valid biallelic GT VCF", {
  st <- data.frame(name = c("NIL1", "NIL1", "NIL2", "NIL2"),
                   chr = 1L, pos = c(1e5, 2e5, 1e5, 2e5),
                   state = c(0L, 2L, 1L, 0L))
  f <- tempfile(fileext = ".vcf")
  write_vcf_impute(st, f, ref = "A", alt = "T")
  ln <- readLines(f)
  hdr <- ln[startsWith(ln, "#CHROM")]
  expect_match(hdr, "NIL1\tNIL2")
  body <- ln[!startsWith(ln, "#")]
  expect_equal(length(body), 2L)                       # two markers
  expect_match(body[1], "GT\t0/0\t0/1")               # marker 1: NIL1 REF, NIL2 het
  expect_match(body[2], "GT\t1/1\t0/0")               # marker 2: NIL1 ALT, NIL2 REF
})

test_that("vectorized decode matches a naive per-sequence reference (bp and cM)", {
  set.seed(11)
  n_samp <- 12L; pos <- seq_len(60L) * 1e5L
  cohort <- do.call(rbind, lapply(seq_len(n_samp), function(s) {
    st <- rep(sample(0:2, 3, replace = TRUE), times = c(20, 20, 20))
    nr <- ifelse(st == 2, 0L, ifelse(st == 1, 4L, 8L)) + rpois(60, 0.3)
    na <- ifelse(st == 0, 0L, ifelse(st == 1, 4L, 8L)) + rpois(60, 0.3)
    data.frame(name = paste0("S", s), donor = paste0("d", s %% 3), chr = 1L,
               pos = pos, cm = pos * 5e-6, n_ref = nr, n_alt = na)
  }))
  li <- log(rep(1/3, 3))
  # naive reference: independently decode each (name, chr) sequence.
  naive <- function(tc, rd) {
    parts <- lapply(split(cohort, cohort$name), function(dn) {
      do.call(rbind, lapply(split(dn, dn$chr), function(dc) {
        dc <- dc[order(dc$pos), ]
        em <- lb_emission_loglik_cpp(dc$n_ref, dc$n_alt, 0.01, 0.05)
        st <- lb_viterbi_cpp(li, em, as.numeric(dc[[tc]]), rd, FALSE)
        data.frame(donor = dc$donor, name = dc$name, chr = dc$chr, pos = dc$pos,
                   state = as.integer(st), stringsAsFactors = FALSE)
      }))
    })
    ref <- do.call(rbind, parts)
    ref[order(ref$donor, ref$name, ref$chr, ref$pos), c("name", "chr", "pos", "state")]
  }
  for (tc in c("pos", "cm")) {
    rd <- if (tc == "cm") 50 else 1e7
    got <- nilHMM:::.lbimpute_states(cohort, 0.01, 0.05, rd, FALSE, li,
                                     "nilHMM", NA_character_, has_donor = TRUE, tcol = tc)
    ref <- naive(tc, rd)
    rownames(got) <- NULL; rownames(ref) <- NULL
    expect_equal(got[, c("name", "chr", "pos", "state")], ref, info = tc)
  }
})

test_that("ragged input (per-sample marker subsets) decodes correctly", {
  # sample B is missing marker 2 -- a non-rectangular cohort. A all REF, B all ALT.
  ragged <- data.frame(
    name = c("A", "A", "A", "B", "B"),
    chr = 1L, pos = c(1e5, 2e5, 3e5, 1e5, 3e5),
    n_ref = c(8L, 8L, 8L, 0L, 0L), n_alt = c(0L, 0L, 0L, 8L, 8L))
  seg <- call_ancestry(ragged, caller = "lbimpute")
  expect_equal(seg$state[seg$name == "A"], 0L)        # A: single REF segment
  expect_equal(seg$state[seg$name == "B"], 2L)        # B: single ALT segment
})

test_that("lb_viterbi_sweep_cpp column k == lb_viterbi_cpp at recombdists[k]", {
  n <- 40L
  nref <- c(rep(8L, 20), rep(0L, 20)); nalt <- c(rep(0L, 20), rep(8L, 20))
  em  <- lb_emission_loglik_cpp(nref, nalt, 0.01, 0.05)
  pos <- as.numeric(seq_len(n) * 1e5L)
  li  <- log(rep(1/3, 3))
  grid <- c(5e5, 1e6, 5e6, 1e7, 5e7)
  sw <- lb_viterbi_sweep_cpp(li, em, pos, grid, FALSE)
  expect_equal(dim(sw), c(n, length(grid)))
  for (k in seq_along(grid))
    expect_equal(sw[, k], lb_viterbi_cpp(li, em, pos, grid[k], FALSE), info = as.character(grid[k]))
})

test_that("lb_viterbi_sweep_cpp guards: non-decreasing tpos, empty grid", {
  em <- lb_emission_loglik_cpp(c(5L, 5L), c(0L, 0L), 0.05, 0.05)
  expect_error(lb_viterbi_sweep_cpp(log(rep(1/3, 3)), em, c(2, 1), 1e7, FALSE), "non-decreasing")
  z <- lb_viterbi_sweep_cpp(log(rep(1/3, 3)), em, c(1, 2), numeric(0), FALSE)
  expect_equal(dim(z), c(2L, 0L))
})

test_that("caller_sweep(lbimpute) == cold call_ancestry per value (bp and cm)", {
  set.seed(3)
  pos <- seq_len(50L) * 1e5L
  cohort <- do.call(rbind, lapply(1:6, function(s) {
    st <- rep(sample(0:2, 3, replace = TRUE), c(20, 15, 15))
    nr <- ifelse(st == 2, 0L, ifelse(st == 1, 4L, 8L)) + rpois(50, 0.3)
    na <- ifelse(st == 0, 0L, ifelse(st == 1, 4L, 8L)) + rpois(50, 0.3)
    data.frame(name = paste0("S", s), donor = paste0("d", s %% 2), chr = 1L,
               pos = pos, cm = pos * 5e-6, n_ref = nr, n_alt = na)
  }))
  for (u in c("bp", "cm")) {
    grid <- if (u == "cm") c(20, 50, 100) else c(5e6, 1e7, 2e7)
    sw <- caller_sweep(cohort, caller = "lbimpute", values = grid, unit = u)
    expect_true("recombdist" %in% names(sw))
    for (v in grid) {
      got  <- sw[sw$recombdist == v, setdiff(names(sw), "recombdist"), drop = FALSE]
      cold <- call_ancestry(cohort, caller = "lbimpute", recombdist = v, unit = u)
      rownames(got) <- NULL; rownames(cold) <- NULL
      expect_equal(got, cold, info = paste(u, v))     # exact per value
    }
  }
})

test_that("caller_sweep(lbimpute) edge cases: empty, single marker, all-zero coverage", {
  li_grid <- c(5e6, 1e7)
  # empty data -> empty segment table (graceful, like call_ancestry)
  empty <- data.frame(name = character(), chr = integer(), pos = integer(),
                      n_ref = integer(), n_alt = integer())
  e <- caller_sweep(empty, caller = "lbimpute", values = li_grid)
  expect_equal(nrow(e), 0L)
  # single-marker sample + an all-zero-coverage sample
  d <- data.frame(name = c("one", "zed", "zed"), chr = 1L,
                  pos = c(1e5, 1e5, 2e5), n_ref = c(8L, 0L, 0L), n_alt = c(0L, 0L, 0L))
  sw <- caller_sweep(d, caller = "lbimpute", values = li_grid)
  expect_equal(sort(unique(sw$recombdist)), sort(li_grid))
  expect_true(all(c("one", "zed") %in% sw$name))
  for (v in li_grid) {                                # still exact vs cold on odd shapes
    got  <- sw[sw$recombdist == v, setdiff(names(sw), "recombdist"), drop = FALSE]
    cold <- call_ancestry(d, caller = "lbimpute", recombdist = v)
    rownames(got) <- NULL; rownames(cold) <- NULL
    expect_equal(got, cold, info = as.character(v))
  }
})

test_that("caller_sweep(lbimpute) rejects non-positive recombdist and honours warnings", {
  d <- data.frame(name = "A", chr = 1L, pos = c(1e5, 2e5),
                  n_ref = c(5L, 5L), n_alt = c(0L, 0L))
  expect_error(caller_sweep(d, caller = "lbimpute", values = c(1e7, -1)), "must be > 0")
  expect_warning(caller_sweep(d, caller = "lbimpute", values = c(50, 1e7), unit = "bp"),
                 "over-relax")     # 50 is cM-sized for a bp sweep
})

test_that("write_vcf_impute rejects duplicate (name, chr, pos) rows", {
  st <- data.frame(name = "NIL1", chr = 1L, pos = c(1e5, 1e5),
                   state = c(0L, 2L))
  expect_error(write_vcf_impute(st, tempfile(fileext = ".vcf")), "duplicate")
})

test_that("lbimpute rejects f_1 + f_2 >= 1", {
  toy <- data.frame(name = "NIL1", chr = 1L, pos = c(1e5, 2e5),
                    n_ref = c(5L, 5L), n_alt = c(0L, 0L))
  expect_error(call_ancestry(toy, caller = "lbimpute", f_1 = 0.6, f_2 = 0.6), "< 1")
})

test_that("write_vcf_impute round-trips through read_vcf_gt", {
  st <- data.frame(name = "NIL1", chr = 1L, pos = c(1e5, 2e5, 3e5),
                   state = c(0L, 1L, 2L))
  f <- tempfile(fileext = ".vcf")
  write_vcf_impute(st, f, ref = "A", alt = "T")
  g <- read_vcf_gt(f)
  expect_equal(g$g, c(0L, 1L, 2L))   # state 0/1/2 -> dosage 0/1/2
})
