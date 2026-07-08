#!/usr/bin/env Rscript
# Tier-0 speed benchmark for the FSFHap port (design/FSFHAP_PORT.md).
#
# Head-to-head wall-clock: stock TASSEL FSFHapImputationPlugin  vs  the native
# nilHMM port, on the SAME input, across a small->large size ladder. JVM startup
# is measured separately (a fixed cost TASSEL pays per invocation that the port
# avoids). Reports seconds, throughput (taxa*markers/s), peak RSS, and speedup.
#
# Runnable NOW for the TASSEL baseline; the port side is GUARDED and self-skips
# until call_ancestry(caller = "fsfhap") exists, so the same harness times both
# without edits once the port lands.
#
# The synthetic generator emits a genuinely-segregating backcross full-sib family
# (recurrent A + donor-C introgression blocks + phet hets + skeleton missingness)
# so FSFHap runs its full pipeline (seed window -> phase -> impute) rather than
# short-circuiting. It also records the true breakpoints -> reuse as the seed for
# the Tier-3 simulated-family fixture.
#
# Usage:
#   Rscript agent/benchmark_fsfhap.R [rungs] [reps] [design]
#     rungs  : comma list of ladder keys (default "xs,s,m"); "full" is slow
#     reps   : timing repetitions per cell (default 1; TASSEL is slow)
#     design : breeding design token "BC{nbc}S{nself}" (default "BC1S4")
#   Override the data source with real input instead of the simulator:
#     FSFHAP_HAPMAP=/path/skeleton.hmp.txt FSFHAP_PED=/path/pedigree.txt Rscript ...
#
# NOTHING about the design is hardcoded: parent contributions (p1/p2) and the
# inbreeding coefficient F are DERIVED from the BC{n}S{m} token, and the donor
# introgression amount simulated + the pedigree written both follow from it.
#   p2 (donor contribution) = 0.5^(nbc+1)   [= design_priors donor dosage]
#   p1 (recurrent)          = 1 - p2
#   F  (inbreeding)         = 1 - 0.5^nself
#   phet (expected het, the TASSEL -phet flag) defaults to 0.5^(nself+1)
#        (BC1S4 -> 0.03125, matching the TeoNAM/Chen run); override with FSFHAP_PHET.
#
# Env:
#   TASSEL_PIPELINE  path to run_pipeline.pl (default /Applications/TASSEL 5/...)
#   TASSEL_XMX       -Xmx heap (default 10g)
#   FSFHAP_DESIGN    design token (overrides the 3rd positional arg)
#   FSFHAP_PHET      expected heterozygosity (overrides the derived default)
#   FSFHAP_LADDER    "n_chr,mpc,n_lines" -> defines a "custom" rung you can run

suppressMessages({
  ok_pkg <- requireNamespace("data.table", quietly = TRUE)
  # load_all so the port becomes callable the moment it exists; tolerate a repo
  # that doesn't build yet (TASSEL-only baseline still runs).
  try(devtools::load_all(".", quiet = TRUE), silent = TRUE)
})
if (!ok_pkg) stop("need data.table")
library(data.table)

a      <- commandArgs(trailingOnly = TRUE)
RUNGS  <- if (length(a) >= 1) strsplit(a[1], ",")[[1]] else c("xs", "s", "m")
REPS   <- if (length(a) >= 2) as.integer(a[2]) else 1L
TASSEL <- Sys.getenv("TASSEL_PIPELINE", "/Applications/TASSEL 5/run_pipeline.pl")
XMX    <- Sys.getenv("TASSEL_XMX", "10g")
# Port fan-out width for the family x chromosome units (matches TASSEL's core use
# so the comparison is like-for-like). Default: all cores, capped so it never
# exceeds the number of independent units on a rung.
THREADS <- { e <- Sys.getenv("FSFHAP_THREADS")
  if (nzchar(e)) as.integer(e) else max(1L, parallel::detectCores()) }
OUTDIR <- file.path("benchmarks", "benchmark_fsfhap_out"); dir.create(OUTDIR, FALSE, TRUE)
IS_MAC <- Sys.info()[["sysname"]] == "Darwin"

# --- breeding design: parse BC{nbc}S{nself} -> contributions + F (no hardcode) -
parse_design <- function(s) {
  m <- regmatches(s, regexec("^BC([0-9]+)S([0-9]+)$", s))[[1]]
  if (length(m) != 3L) stop("design must look like 'BC1S4'; got: ", s)
  nbc <- as.integer(m[2]); nself <- as.integer(m[3]); p2 <- 0.5^(nbc + 1)
  phet_env <- Sys.getenv("FSFHAP_PHET")
  list(design = s, nbc = nbc, nself = nself, p1 = 1 - p2, p2 = p2,
       F = 1 - 0.5^nself,
       phet = if (nzchar(phet_env)) as.numeric(phet_env) else 0.5^(nself + 1))
}
DESIGN_TOKEN <- Sys.getenv("FSFHAP_DESIGN",
                           if (length(a) >= 3) a[3] else "BC1S4")
DES <- parse_design(DESIGN_TOKEN)

# --- size ladder: small -> full-genotype/full-family -------------------------
# (n_chr, markers-per-chr, n_lines). Tune to your data; "full" ~ TeoNAM scale.
LADDER <- list(
  xs   = list(n_chr = 1L,  mpc = 500L,   n_lines = 20L),   # CI smoke
  s    = list(n_chr = 1L,  mpc = 5000L,  n_lines = 50L),   # 1 chr full
  m    = list(n_chr = 5L,  mpc = 5000L,  n_lines = 100L),  # multi-chr
  full = list(n_chr = 10L, mpc = 5100L,  n_lines = 100L)   # ~51k markers, full family
)
# custom rung from FSFHAP_LADDER="n_chr,mpc,n_lines" (run it with rung "custom")
.lad <- Sys.getenv("FSFHAP_LADDER")
if (nzchar(.lad)) { v <- as.integer(strsplit(.lad, ",")[[1]])
  LADDER$custom <- list(n_chr = v[1], mpc = v[2], n_lines = v[3]) }

# --- synthetic backcross full-sib family, donor dosage + het from the design --
# donor_frac (= design p2) sets the expected genome fraction that is donor(C);
# phet (= design phet) sets residual heterozygosity. Returns a marker table + an
# A/M/C/N call matrix (markers x lines) + true state (0/1/2) matrix (Tier-3 seed).
sim_family <- function(n_chr, mpc, n_lines, donor_frac, phet, miss = 0.15, seed = 1L) {
  set.seed(seed)
  chr <- rep(seq_len(n_chr), each = mpc)
  pos <- as.integer(rep(seq_len(mpc), n_chr) * 1e5L)
  n_mark <- length(chr)
  markers <- sprintf("m%06d", seq_len(n_mark))
  true <- matrix(0L, n_mark, n_lines)               # 0=A(recurrent), 2=C(donor), 1=het
  seg_w <- 0.12                                     # mean segment width (frac of chr)
  for (j in seq_len(n_lines)) for (ch in seq_len(n_chr)) {
    idx <- which(chr == ch)
    # segments so E[donor fraction] ~= donor_frac: n = donor_frac / seg_w
    nseg <- rpois(1, donor_frac / seg_w)
    if (nseg > 0) for (b in seq_len(nseg)) {
      w <- max(2L, rpois(1, length(idx) * seg_w))
      s0 <- sample(idx, 1L); rng <- s0:min(max(idx), s0 + w - 1L)
      true[rng, j] <- 2L
    }
  }
  # het band at donor/recurrent boundaries + random phet
  het <- (true == 2L) & (rbind(FALSE, head(true, -1) == 0L) | rbind(tail(true, -1) == 0L, FALSE))
  true[het] <- 1L
  true[matrix(runif(n_mark * n_lines) < phet, n_mark, n_lines) & true == 0L] <- 1L
  call <- matrix(c("A", "M", "C")[true + 1L], n_mark, n_lines)
  call[matrix(runif(n_mark * n_lines) < miss, n_mark, n_lines)] <- "N"   # skeleton gaps
  colnames(call) <- sprintf("L%03d", seq_len(n_lines))
  list(markers = markers, chr = chr, pos = pos, call = call, true = true,
       n_mark = n_mark, n_lines = n_lines)
}

# --- write TASSEL Hapmap + 7-col FSFHap pedigree (teonam_fsfhap.R format) -----
write_tassel_inputs <- function(fam, tag) {
  hd <- data.table(`rs#` = fam$markers, alleles = "A/C", chrom = as.integer(fam$chr),
                   pos = as.integer(fam$pos), strand = "+", `assembly#` = NA,
                   center = NA, protLSID = NA, assayLSID = NA, panelLSID = NA, QCcode = NA)
  hmp <- cbind(hd, as.data.table(fam$call))
  hmpf <- file.path(OUTDIR, sprintf("skeleton_%s.hmp.txt", tag))
  fwrite(hmp, hmpf, sep = "\t", quote = FALSE, na = "N")
  ped <- data.table(family = "SIM", taxon = colnames(fam$call), parent1 = "A",
                    parent2 = "C", p1 = DES$p1, p2 = DES$p2, F = DES$F)   # from design
  pedf <- file.path(OUTDIR, sprintf("pedigree_%s.txt", tag))
  fwrite(ped, pedf, sep = "\t", quote = FALSE, col.names = TRUE)
  list(hmp = hmpf, ped = pedf)
}

# --- timed subprocess with best-effort peak RSS (macOS /usr/bin/time -l) ------
timed_run <- function(cmd) {
  if (IS_MAC) {
    out <- suppressWarnings(system2("/usr/bin/time", c("-l", "sh", "-c", shQuote(cmd)),
                                    stdout = TRUE, stderr = TRUE))
    el  <- system.time(NULL)  # placeholder; we parse wall from our own clock below
    rss <- as.numeric(sub("\\s*maximum resident set size.*", "",
                          grep("maximum resident set size", out, value = TRUE)[1]))
    attr(out, "rss_bytes") <- rss
    out
  } else {
    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE); NULL
  }
}

run_tassel <- function(inp, tag) {
  out <- file.path(OUTDIR, sprintf("imputed_%s", tag))
  log <- file.path(OUTDIR, sprintf("fsfhap_%s.log", tag))
  cmd <- sprintf(
    '"%s" -Xmx%s -h "%s" -FSFHapImputationPlugin -pedigrees "%s" -bc true -phet %g -fillgaps true -maxMissing 1.0 -minMaf 0.0 -minR 0.0 -logfile "%s" -endPlugin -export "%s" -exportType Hapmap',
    TASSEL, XMX, inp$hmp, inp$ped, DES$phet, log, out)
  t <- system.time(res <- timed_run(cmd))[["elapsed"]]
  list(sec = t, rss = if (!is.null(res)) attr(res, "rss_bytes") else NA_real_)
}

# JVM+classload+IO baseline: import a tiny hapmap and re-export, no FSFHap plugin.
tassel_startup <- function(inp) {
  out <- file.path(OUTDIR, "startup_export")
  cmd <- sprintf('"%s" -Xmx%s -h "%s" -export "%s" -exportType Hapmap',
                 TASSEL, XMX, inp$hmp, out)
  system.time(timed_run(cmd))[["elapsed"]]
}

# --- port side (GUARDED: self-skips until the caller exists) ------------------
port_available <- function() {
  isTRUE(tryCatch("fsfhap" %in% eval(formals(nilHMM::call_states)$caller),
                  error = function(e) FALSE))
}
# synthetic fam -> canonical (name, chr, pos, g) long table + family column
fam_to_dt <- function(fam) {
  g <- fam$true; g[fam$call == "N"] <- 3L               # g in {0,1,2,3=missing}
  data.table(name = rep(colnames(fam$call), each = fam$n_mark), family = "SIM",
             chr = fam$chr, pos = fam$pos, g = as.integer(g))
}
# time the port on a prebuilt canonical table (real or synthetic); `design` routes
# + derives phet, `threads` fans out the family x chromosome units.
run_port <- function(dt, design, threads = 1L) {
  gc(reset = TRUE)
  t <- system.time(
    nilHMM::call_ancestry(dt, caller = "fsfhap", design = design, threads = threads)
  )[["elapsed"]]
  list(sec = t, rss = sum(gc()[, 6]) / 1024)           # peak Mb -> GB (approx)
}

# --- run the ladder -----------------------------------------------------------
have_tassel <- file.exists(TASSEL)
have_port   <- port_available()
override_hmp <- Sys.getenv("FSFHAP_HAPMAP"); override_ped <- Sys.getenv("FSFHAP_PED")
cat(sprintf("TASSEL: %s | port caller: %s | rungs: %s | reps: %d\n",
            if (have_tassel) "found" else "MISSING (baseline skipped)",
            if (have_port) "available" else "not yet implemented (skipped)",
            paste(RUNGS, collapse = ","), REPS))
cat(sprintf("design %s: nbc=%d nself=%d -> p1=%.4f p2=%.4f F=%.4f phet=%g\n",
            DES$design, DES$nbc, DES$nself, DES$p1, DES$p2, DES$F, DES$phet))

# --- REAL-DATA head-to-head (FSFHAP_HAPMAP [+ FSFHAP_PED]) --------------------
# When FSFHAP_HAPMAP is set, benchmark the port on the REAL genotypes (loaded
# through the package's own read_hapmap/read_pedigree adapters) instead of the
# synthetic ladder. Family grouping comes from either:
#   (a) FSFHAP_PED  -> read_pedigree(); design DERIVED from contribution/F
#       (no hardcode: p2 = 0.5^(nbc+1), F = 1 - 0.5^nself); taxa filtered to it.
#   (b) no pedigree -> taxon-name prefix via FSFHAP_FAMILY_REGEX (capture group 1;
#       default "^([A-Za-z]+[0-9]+)" -> "TIL01" from "TIL01A001"), ALL hapmap taxa;
#       design = DESIGN_TOKEN. TASSEL gets a pedigree synthesized from that grouping
#       + the design contributions (written to OUTDIR, a benchmark artifact).
# Either way TASSEL runs on the same HapMap + pedigree, so it is a true head-to-head.
if (nzchar(override_hmp)) {
  if (!have_port) stop("real-data mode needs the port caller (devtools::load_all failed?)")
  fam_regex <- Sys.getenv("FSFHAP_FAMILY_REGEX", "^([A-Za-z]+[0-9]+)")
  if (nzchar(override_ped)) {                              # (a) pedigree-driven
    ped <- nilHMM::read_pedigree(override_ped)
    nbc   <- max(0L, round(log2(1 / ped$contribution2[1])) - 1L)
    nself <- round(log2(1 / (1 - ped$F[1])))
    design <- if (all(is.finite(c(nbc, nself)))) sprintf("BC%dS%d", nbc, nself) else DESIGN_TOKEN
    ped_path <- override_ped; grouping <- sprintf("pedigree %s", basename(override_ped))
    samples <- ped$taxon
  } else {                                                 # (b) taxon-prefix grouping
    hdr <- strsplit(readLines(override_hmp, n = 1L), "\t", fixed = TRUE)[[1]]
    taxa <- hdr[-(1:11)]                                   # 11 fixed HapMap metadata columns
    fam_of <- vapply(regmatches(taxa, regexec(fam_regex, taxa)),
                     function(m) if (length(m) >= 2L) m[2] else NA_character_, character(1))
    if (anyNA(fam_of)) stop("FSFHAP_FAMILY_REGEX matched no family for some taxa; adjust the regex")
    ped <- data.frame(taxon = taxa, family = fam_of, stringsAsFactors = FALSE)
    design <- DESIGN_TOKEN; samples <- taxa
    grouping <- sprintf("taxon-prefix /%s/", fam_regex)
    # synthesize a TASSEL pedigree from the grouping + design contributions
    ped_path <- file.path(OUTDIR, "pedigree_real_prefix.txt")
    fwrite(data.table(family = fam_of, taxon = taxa, parent1 = "P1",
                      parent2 = paste0(fam_of, "_donor"), p1 = DES$p1, p2 = DES$p2, F = DES$F),
           ped_path, sep = "\t", quote = FALSE, col.names = TRUE)
  }
  fams <- unique(ped$family); n_chr_real <- NA_integer_
  cat(sprintf("\n[REAL] hapmap=%s\n       grouping=%s -> %d taxa, %d famil%s (%s)\n       design=%s\n",
              override_hmp, grouping, length(samples), length(fams),
              if (length(fams) == 1) "y" else "ies", paste(head(fams, 6), collapse = ","), design))

  t_io <- system.time(data <- nilHMM::read_hapmap(override_hmp, samples = samples))[["elapsed"]]
  data$family <- ped$family[match(data$name, ped$taxon)]
  if (anyNA(data$family)) stop("real-data mode: HapMap taxa not covered by the grouping")
  n_taxa <- length(unique(data$name)); n_chr_real <- length(unique(data$chr))
  n_mark <- nrow(data) / n_taxa; cells <- as.numeric(n_taxa) * n_mark
  n_units <- length(fams) * n_chr_real
  eff_threads <- min(THREADS, n_units)
  cat(sprintf("       read_hapmap %.1fs -> %d taxa x %d markers x %d chr = %s cells (%d units)\n",
              t_io, n_taxa, n_mark, n_chr_real, format(cells, big.mark = ","), n_units))

  start_sec <- tass_sec <- tass_rss <- NA_real_
  if (have_tassel) {
    inp <- list(hmp = override_hmp, ped = ped_path)
    start_sec <- tassel_startup(inp)
    tr <- replicate(REPS, unlist(run_tassel(inp, "real")), simplify = FALSE)
    tass_sec <- median(vapply(tr, `[[`, 0, "sec")); tass_rss <- median(vapply(tr, `[[`, 0, "rss"))
    cat(sprintf("  TASSEL: total=%.1fs  startup~%.1fs  compute~%.1fs  peakRSS=%s\n",
                tass_sec, start_sec, max(0, tass_sec - start_sec),
                if (is.na(tass_rss)) "NA" else sprintf("%.1f GB", tass_rss / 2^30)))
  }
  port_sec <- median(replicate(REPS, run_port(data, design, threads = 1L)$sec))
  cat(sprintf("  port (1 thread):    %.2fs\n", port_sec))
  port_par_sec <- NA_real_
  if (eff_threads > 1L) {
    port_par_sec <- median(replicate(REPS, run_port(data, design, threads = eff_threads)$sec))
    cat(sprintf("  port (%d threads):  %.2fs  (%.1fx over serial)\n",
                eff_threads, port_par_sec, port_sec / port_par_sec))
  }
  best_port <- if (!is.na(port_par_sec)) port_par_sec else port_sec
  compute <- if (is.na(tass_sec)) NA_real_ else max(0, tass_sec - start_sec)
  res <- data.table(
    rung = "real", n_taxa = n_taxa, n_markers = n_mark, n_families = length(fams),
    n_chr = n_chr_real, n_units = n_units, cells = cells, io_read_s = round(t_io, 2),
    tassel_total_s = round(tass_sec, 2), tassel_compute_s = round(compute, 2),
    tassel_cells_per_s = round(cells / compute),
    port_s = round(port_sec, 2), port_threads = eff_threads, port_par_s = round(port_par_sec, 2),
    port_cells_per_s = round(cells / best_port),
    speedup_vs_compute = round(compute / best_port, 1))
  cat("\n============= Tier-0 REAL-DATA benchmark: FSFHap port vs TASSEL =============\n")
  print(res)
  csv <- file.path(OUTDIR, "benchmark_fsfhap_real.csv"); fwrite(res, csv)
  cat(sprintf("\nwrote %s\n", csv))
  quit(save = "no")
}

rows <- list()
for (key in RUNGS) {
  cfg <- LADDER[[key]]; if (is.null(cfg)) { cat("unknown rung:", key, "\n"); next }
  fam <- sim_family(cfg$n_chr, cfg$mpc, cfg$n_lines, donor_frac = DES$p2, phet = DES$phet)
  cells <- as.numeric(fam$n_mark) * fam$n_lines
  cat(sprintf("\n[%s] %d lines x %d markers (%d chr) = %s cells\n", key,
              fam$n_lines, fam$n_mark, cfg$n_chr, format(cells, big.mark = ",")))

  inp <- if (nzchar(override_hmp) && nzchar(override_ped))
    list(hmp = override_hmp, ped = override_ped) else write_tassel_inputs(fam, key)

  tass_sec <- tass_rss <- start_sec <- NA_real_
  if (have_tassel) {
    start_sec <- tassel_startup(inp)
    tr <- replicate(REPS, unlist(run_tassel(inp, key)), simplify = FALSE)
    tass_sec <- median(vapply(tr, `[[`, 0, "sec"))
    tass_rss <- median(vapply(tr, `[[`, 0, "rss"))
    cat(sprintf("  TASSEL: total=%.1fs  startup~%.1fs  compute~%.1fs  peakRSS=%s\n",
                tass_sec, start_sec, max(0, tass_sec - start_sec),
                if (is.na(tass_rss)) "NA" else sprintf("%.1f GB", tass_rss / 2^30)))
  }

  port_sec <- port_par_sec <- NA_real_
  # effective threads: never more than the independent (family x chr) units on
  # this rung (1 family here, so = n_chr), so idle cores aren't reported as work.
  eff_threads <- min(THREADS, cfg$n_chr)
  if (have_port) {
    dt <- fam_to_dt(fam)
    pr <- replicate(REPS, run_port(dt, DESIGN_TOKEN, threads = 1L)$sec)
    port_sec <- median(pr)
    cat(sprintf("  port (1 thread):   %.2fs\n", port_sec))
    if (eff_threads > 1L) {
      pp <- replicate(REPS, run_port(dt, DESIGN_TOKEN, threads = eff_threads)$sec)
      port_par_sec <- median(pp)
      cat(sprintf("  port (%d threads):  %.2fs  (%.1fx over serial)\n",
                  eff_threads, port_par_sec, port_sec / port_par_sec))
    }
  } else {
    cat("  port:   TODO (caller not implemented) — TASSEL-only baseline\n")
  }

  # the port's best wall-clock (parallel if it ran, else serial) drives the speedup
  best_port <- if (!is.na(port_par_sec)) port_par_sec else port_sec
  compute <- if (is.na(tass_sec)) NA_real_ else max(0, tass_sec - start_sec)
  rows[[key]] <- data.table(
    rung = key, n_lines = fam$n_lines, n_markers = fam$n_mark, cells = cells,
    tassel_total_s = round(tass_sec, 2), tassel_startup_s = round(start_sec, 2),
    tassel_compute_s = round(compute, 2), tassel_peak_gb = round(tass_rss / 2^30, 2),
    tassel_cells_per_s = round(cells / compute),
    port_s = round(port_sec, 2), port_threads = eff_threads,
    port_par_s = round(port_par_sec, 2),
    port_cells_per_s = round(cells / best_port),
    speedup_vs_total = round(tass_sec / best_port, 1),
    speedup_vs_compute = round(compute / best_port, 1))
}

res <- rbindlist(rows, fill = TRUE)
cat("\n================ Tier-0 benchmark: FSFHap port vs TASSEL ================\n")
print(res)
csv <- file.path(OUTDIR, "benchmark_fsfhap_results.csv")
fwrite(res, csv)
cat(sprintf("\nwrote %s\n", csv))
cat("NOTE: speedup columns are NA until the port caller exists; startup is the\n",
    "fixed JVM+classload+IO cost TASSEL pays per family and the port avoids.\n")
