# io.R adapters: TASSEL HapMap + PLINK (.bed/.fam) + pedigree -> canonical tables.
# Round-trip each reader from a hand-built fixture with a known genotype matrix.

test_that("read_hapmap maps single-char IUPAC to canonical g", {
  f <- tempfile(fileext = ".hmp.txt")
  writeLines(c(
    paste(c("rs#","alleles","chrom","pos","strand","assembly#","center",
            "protLSID","assayLSID","panelLSID","QCcode","L1","L2","L3"), collapse = "\t"),
    paste(c("m1","A/C",1,100,"+",NA,NA,NA,NA,NA,NA, "A","C","M"), collapse = "\t"),   # A-hom, C-hom, het
    paste(c("m2","A/C",1,200,"+",NA,NA,NA,NA,NA,NA, "N","A","C"), collapse = "\t")),  # missing, A, C
    f)
  d <- read_hapmap(f)
  expect_named(d, c("name", "chr", "pos", "g"))
  expect_equal(d$g[d$name == "L1"], c(0L, 3L))    # A, N
  expect_equal(d$g[d$name == "L2"], c(2L, 0L))    # C, A
  expect_equal(d$g[d$name == "L3"], c(1L, 2L))    # M(het), C
  expect_equal(read_hapmap(f, samples = "L2")$name, rep("L2", 2))   # subset
})

test_that("read_hapmap handles two-char diploid cells", {
  f <- tempfile(fileext = ".hmp.txt")
  writeLines(c(
    paste(c("rs#","alleles","chrom","pos","strand","assembly#","center",
            "protLSID","assayLSID","panelLSID","QCcode","L1"), collapse = "\t"),
    paste(c("m1","A/C",1,100,"+",NA,NA,NA,NA,NA,NA, "AA"), collapse = "\t"),
    paste(c("m2","A/C",1,200,"+",NA,NA,NA,NA,NA,NA, "AC"), collapse = "\t"),
    paste(c("m3","A/C",1,300,"+",NA,NA,NA,NA,NA,NA, "CC"), collapse = "\t"),
    paste(c("m4","A/C",1,400,"+",NA,NA,NA,NA,NA,NA, "NN"), collapse = "\t")), f)
  expect_equal(read_hapmap(f)$g, c(0L, 1L, 2L, 3L))
})

test_that("read_plink decodes a hand-built .bed to A1-dosage g", {
  # known genotype matrix (samples x SNPs), g = A1 dosage {0,1,2,3=missing}
  G <- rbind(c(2L, 1L, 0L, 3L), c(0L, 2L, 1L, 3L)); n <- nrow(G); m <- ncol(G)
  pre <- tempfile()
  # .bim: chr snp cM bp A1 A2 ; .fam: FID IID PID MID sex pheno
  write.table(data.frame(1, paste0("m", 1:m), 0, (1:m) * 100L, "A", "C"),
              paste0(pre, ".bim"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(data.frame(paste0("F", 1:n), paste0("L", 1:n), 0, 0, 0, -9),
              paste0(pre, ".fam"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  # encode .bed (SNP-major): g -> 2-bit code (inverse of the reader's map)
  g2code <- c(`0` = 3L, `1` = 2L, `2` = 0L, `3` = 1L)
  bpS <- ceiling(n / 4); bytes <- as.raw(c(0x6c, 0x1b, 0x01))
  for (j in seq_len(m)) {
    codes <- g2code[as.character(G[, j])]
    for (b in seq_len(bpS)) {
      idx <- (b - 1L) * 4L + 1:4; s <- ifelse(idx <= n, codes[idx], 0L); s[is.na(s)] <- 0L
      bytes <- c(bytes, as.raw(bitwOr(bitwOr(s[1], bitwShiftL(s[2], 2L)),
                                      bitwOr(bitwShiftL(s[3], 4L), bitwShiftL(s[4], 6L)))))
    }
  }
  writeBin(bytes, paste0(pre, ".bed"))
  d <- read_plink(pre)
  expect_named(d, c("name", "chr", "pos", "g"))
  expect_equal(d$g[d$name == "L1"], G[1, ])
  expect_equal(d$g[d$name == "L2"], G[2, ])
  expect_equal(unique(d$pos), (1:m) * 100L)
})

test_that("read_pedigree parses FSFHap TSV and PLINK .fam", {
  p <- tempfile()
  writeLines(c("family\ttaxon\tparent1\tparent2\tp1\tp2\tF",
               "SIM\tL1\tA\tC\t0.75\t0.25\t0.9375",
               "SIM\tL2\tA\tC\t0.75\t0.25\t0.9375"), p)
  ped <- read_pedigree(p)
  expect_equal(ped$taxon, c("L1", "L2"))
  expect_equal(unique(ped$family), "SIM")
  expect_equal(ped$contribution1, c(0.75, 0.75))
  expect_equal(ped$F, c(0.9375, 0.9375))
  # phet derivation from F (matches the caller's .fsfhap_phet)
  expect_equal(.fsfhap_phet(ped$F[1]), 0.03125)

  fam <- tempfile(fileext = ".fam")
  writeLines(c("FAM1 T1 0 0 0 -9", "FAM1 T2 0 0 0 -9", "FAM2 T3 0 0 0 -9"), fam)
  pf <- read_pedigree(fam, format = "fam")
  expect_equal(pf$taxon, c("T1", "T2", "T3"))
  expect_equal(pf$family, c("FAM1", "FAM1", "FAM2"))
  expect_true(all(is.na(pf$F)))
})

test_that("read_hapmap feeds call_ancestry(caller='fsfhap') end to end", {
  # small RIL-ish family written as HapMap -> read -> attach family -> call
  set.seed(1); n <- 30L; nm <- 300L
  g <- matrix(0L, nm, n)
  for (jj in seq_len(n)) { st <- 1L; cur <- sample(c(0L, 2L), 1L)
    for (bp in c(sort(sample(2:(nm - 1L), 3L)), nm)) { g[st:bp, jj] <- cur; cur <- ifelse(cur == 0L, 2L, 0L); st <- bp + 1L } }
  call <- matrix(c("A", "M", "C", "N")[g + 1L], nm, n)
  taxa <- sprintf("L%02d", seq_len(n))
  hd <- c("rs#","alleles","chrom","pos","strand","assembly#","center","protLSID","assayLSID","panelLSID","QCcode")
  rows <- vapply(seq_len(nm), function(i) paste(c(paste0("m", i), "A/C", 1L, i * 100L, "+",
                 rep("NA", 6), call[i, ]), collapse = "\t"), character(1))
  f <- tempfile(fileext = ".hmp.txt"); writeLines(c(paste(c(hd, taxa), collapse = "\t"), rows), f)
  dat <- read_hapmap(f); dat$family <- "F"
  seg <- call_ancestry(dat, caller = "fsfhap", design = "BC2S2")   # RIL-ish -> biparental route
  expect_named(seg, c("source", "donor", "name", "chr", "start_bp", "end_bp", "state"))
  expect_true(nrow(seg) > 0 && all(seg$state %in% c(0L, 1L, 2L)))
})
