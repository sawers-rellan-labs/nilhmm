# Helpers for locating bundled test fixtures. Works under devtools::test()
# (source tree) and R CMD check (installed pkg) via testthat::test_path().

golden_slice_path <- function(...) {
  testthat::test_path("fixtures", "golden_slice", ...)
}

# Read a nilHMM count TSV (chr pos ref n_ref alt n_alt; no header) for one chr1
# golden-slice sample.
read_golden_counts <- function(path) {
  df <- utils::read.table(path, sep = "\t", header = FALSE,
                          stringsAsFactors = FALSE,
                          col.names = c("chr", "pos", "ref", "n_ref", "alt", "n_alt"))
  df
}

read_golden_expected <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE)
}
