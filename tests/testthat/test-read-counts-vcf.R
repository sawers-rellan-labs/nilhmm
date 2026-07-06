# read_counts(format = "vcf_ad"): per-sample allelic depths from a biallelic VCF's
# AD FORMAT field into the common (name, chr, pos, n_ref, n_alt) table.

vcf_lines <- function(records, fmt = "GT:AD:DP", samples = c("NIL1", "NIL2")) {
  c("##fileformat=VCFv4.2",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples),
          collapse = "\t"),
    records)
}
write_vcf <- function(lines) { f <- tempfile(fileext = ".vcf"); writeLines(lines, f); f }

test_that("AD is parsed into n_ref/n_alt for every sample", {
  f <- write_vcf(vcf_lines(c(
    "1\t100000\t.\tA\tT\t.\t.\t.\tGT:AD:DP\t0/0:8,0:8\t0/1:4,5:9",
    "1\t200000\t.\tC\tG\t.\t.\t.\tGT:AD:DP\t1/1:0,7:7\t0/0:6,0:6")))
  d <- read_counts(f, format = "vcf_ad")
  expect_named(d, c("name", "chr", "pos", "n_ref", "n_alt"))
  expect_equal(sort(unique(d$name)), c("NIL1", "NIL2"))
  n1 <- d[d$name == "NIL1", ]; n1 <- n1[order(n1$pos), ]
  expect_equal(n1$n_ref, c(8L, 0L)); expect_equal(n1$n_alt, c(0L, 7L))
  n2 <- d[d$name == "NIL2", ]; n2 <- n2[order(n2$pos), ]
  expect_equal(n2$n_ref, c(4L, 6L)); expect_equal(n2$n_alt, c(5L, 0L))
  expect_equal(sort(unique(d$pos)), c(100000L, 200000L))
})

test_that("missing / '.' AD fields become 0,0", {
  f <- write_vcf(vcf_lines(c(
    "1\t100000\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:5,4\t./.:.",   # NIL2 fully missing
    "1\t200000\t.\tA\tT\t.\t.\t.\tGT:AD\t./.\t1/1:0,9")))    # NIL1 GT-only cell
  d <- read_counts(f, format = "vcf_ad")
  n2_1 <- d[d$name == "NIL2" & d$pos == 100000L, ]
  expect_equal(n2_1$n_ref, 0L); expect_equal(n2_1$n_alt, 0L)
  n1_2 <- d[d$name == "NIL1" & d$pos == 200000L, ]
  expect_equal(n1_2$n_ref, 0L); expect_equal(n1_2$n_alt, 0L)
})

test_that("AD position is found within FORMAT (not assumed first) and varies per record", {
  # AD in different colon positions per record -> per-row index path.
  f <- write_vcf(vcf_lines(c(
    "1\t100000\t.\tA\tT\t.\t.\t.\tGT:AD:DP\t0/0:8,1:9\t0/1:3,3:6",
    "1\t200000\t.\tA\tT\t.\t.\t.\tGT:DP:AD\t0/0:9:7,2\t1/1:8:1,8")))
  d <- read_counts(f, format = "vcf_ad")
  r2 <- d[d$name == "NIL1" & d$pos == 200000L, ]      # AD is the 3rd field here
  expect_equal(r2$n_ref, 7L); expect_equal(r2$n_alt, 2L)
  r2b <- d[d$name == "NIL2" & d$pos == 200000L, ]
  expect_equal(r2b$n_ref, 1L); expect_equal(r2b$n_alt, 8L)
})

test_that("multiallelic AD keeps only the first ALT (biallelic assumption)", {
  f <- write_vcf(vcf_lines(c(
    "1\t100000\t.\tA\tT,C\t.\t.\t.\tGT:AD\t0/1:5,3,2\t1/2:0,4,6")))
  d <- read_counts(f, format = "vcf_ad")
  n1 <- d[d$name == "NIL1", ]; n2 <- d[d$name == "NIL2", ]
  expect_equal(n1$n_ref, 5L); expect_equal(n1$n_alt, 3L)
  expect_equal(n2$n_ref, 0L); expect_equal(n2$n_alt, 4L)
})

test_that("errors: no AD in FORMAT, no #CHROM header", {
  f <- write_vcf(vcf_lines("1\t100000\t.\tA\tT\t.\t.\t.\tGT:DP\t0/0:8\t0/1:9"))
  expect_error(read_counts(f, format = "vcf_ad"), "AD")
  g <- tempfile(fileext = ".vcf"); writeLines(c("##fileformat=VCFv4.2", "no header"), g)
  expect_error(read_counts(g, format = "vcf_ad"), "#CHROM")
})

test_that("vcf_ad output feeds call_ancestry(caller = 'lbimpute')", {
  recs <- sprintf("1\t%d\t.\tA\tT\t.\t.\t.\tGT:AD\t%s\t%s",
                  seq_len(10) * 1e5L,
                  c(rep("0/0:8,0", 5), rep("1/1:0,8", 5)),   # NIL1: REF then ALT
                  rep("0/0:8,0", 10))                          # NIL2: all REF
  d <- read_counts(write_vcf(vcf_lines(recs)), format = "vcf_ad")
  seg <- call_ancestry(d, caller = "lbimpute")
  expect_equal(seg$state[seg$name == "NIL2"], 0L)
  expect_setequal(seg$state[seg$name == "NIL1"], c(0L, 2L))
})
