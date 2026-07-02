# read_vcf_gt(): VCF GT field -> hard genotype g {0,1,2,3}, fed to the gt
# (categorical) emission via caller="nnil" (the saturated-depth genotype path).

# minimal biallelic VCF: 2 samples, one chr, a donor block for S1
write_test_vcf <- function(path, n = 60L, block = 20:35) {
  con <- file(path, "w")
  writeLines(c("##fileformat=VCFv4.2",
               '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
               paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","S1","S2"),
                     collapse = "\t")), con)
  for (i in seq_len(n)) {
    s1 <- if (i %in% block) (if (i %% 3 == 0) "1/1" else "0/1") else "0/0"  # donor block: het/alt
    s2 <- "0/0"                                                             # pure REF sample
    if (i %% 17 == 0) s1 <- "./."                                          # scattered missing
    writeLines(paste(c(paste0("chr1"), i * 1000000L, ".", "A", "G", ".", "PASS", ".", "GT", s1, s2),
                     collapse = "\t"), con)
  }
  close(con)
}

test_that("read_vcf_gt maps GT to alt-dosage g and returns the long schema", {
  vcf <- tempfile(fileext = ".vcf"); write_test_vcf(vcf)
  d <- read_vcf_gt(vcf)
  expect_named(d, c("name", "chr", "pos", "g"))
  expect_setequal(unique(d$name), c("S1", "S2"))
  expect_true(all(d$g %in% 0:3))
  expect_true(all(d$g[d$name == "S2"] == 0L))                 # pure-REF sample -> all g=0
  expect_true(any(d$g == 3L))                                 # ./. -> missing (3)
  expect_true(any(d$g[d$name == "S1"] == 1L) && any(d$g[d$name == "S1"] == 2L))  # het + alt-hom present
})

test_that("sample subset works", {
  vcf <- tempfile(fileext = ".vcf"); write_test_vcf(vcf)
  expect_setequal(unique(read_vcf_gt(vcf, samples = "S1")$name), "S1")
})

test_that("call_ancestry runs the gt genotype path on a read_vcf_gt table", {
  vcf <- tempfile(fileext = ".vcf"); write_test_vcf(vcf)
  d <- read_vcf_gt(vcf)
  calls <- call_ancestry(d, caller = "nnil", design = "BC2S3")   # g-only -> auto gt emission
  expect_named(calls, c("source", "donor", "name", "chr", "start_bp", "end_bp", "state"))
  expect_true(all(calls$state %in% 0:2))
  # S1's donor block (20-35 Mb) is detected as non-REF; S2 is all REF
  s1nr <- calls[calls$name == "S1" & calls$state > 0L, ]
  expect_true(any(s1nr$start_bp <= 35e6 & s1nr$end_bp >= 20e6))
  expect_true(all(calls$state[calls$name == "S2"] == 0L))
})

test_that("germ passthrough smooths isolated miscalls (higher germ -> fewer non-REF)", {
  pos <- as.integer((1:200) * 1e6)                            # 200 markers on chr1
  g <- rep(0L, 200); g[seq(10, 200, by = 13)] <- 1L           # scattered isolated het miscalls in REF
  d <- data.frame(name = "S", chr = 1L, pos = pos, g = g)
  lo <- call_ancestry(d, caller = "nnil", design = "BC2S3", germ = 0.001)   # trust calls -> fragments
  hi <- call_ancestry(d, caller = "nnil", design = "BC2S3", germ = 0.30)    # absorb as error -> smooth
  nonref <- function(x) sum(x$state > 0L)
  expect_gte(nonref(lo), nonref(hi))                          # raising germ never increases donor calls
  expect_true(all(hi$state == 0L))                            # high germ absorbs the isolated hets -> all REF
})

test_that("a g-only genotype input rejects callers that need read counts", {
  d <- data.frame(name = "S1", chr = 1L, pos = as.integer((1:50) * 1e6), g = 0L)
  expect_error(call_ancestry(d, caller = "rtiger"), "read counts")
  expect_error(call_ancestry(d, caller = "binhmm", design = "BC2S3"), "n_ref|read counts|columns")
})
