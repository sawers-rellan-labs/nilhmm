# Recipe for the shipped data/maize_map_v5.rda (see ?maize_map_v5).
#
# SOURCE (cite): the cleaned B73 v5 maize consensus genetic map -- marker -> chr,
# cM, bp; 19,486 loci, 10 chromosomes, ~1783 cM total -- assembled in the
# companion zealhmm pipeline. cM-space is assembly-robust; the bp column is tied
# to B73 v5.
#
# The raw source is NOT vendored in this repo; data/maize_map_v5.rda is the
# authoritative shipped value (this script documents its provenance + recipe).
# To regenerate, point NILHMM_MAP_SRC at the cleaned source .rds (columns:
# locus, chr, cm, bp) and run from the package root:
#   NILHMM_MAP_SRC=/path/to/maize_map_v5_clean.rds Rscript data-raw/make_maize_map.R

SRC <- Sys.getenv("NILHMM_MAP_SRC")
if (!nzchar(SRC))
  stop("set NILHMM_MAP_SRC to the cleaned consensus-map .rds (columns: locus, chr, cm, bp)")
if (!file.exists(SRC)) stop("map source not found: ", SRC)

src <- as.data.frame(readRDS(SRC), stringsAsFactors = FALSE)
stopifnot(all(c("locus", "chr", "cm", "bp") %in% names(src)))

maize_map_v5 <- data.frame(
  locus = as.character(src$locus),
  chr   = as.integer(src$chr),
  cm    = as.numeric(src$cm),
  bp    = as.numeric(src$bp),
  stringsAsFactors = FALSE
)
maize_map_v5 <- maize_map_v5[order(maize_map_v5$chr, maize_map_v5$cm), ]
rownames(maize_map_v5) <- NULL
attr(maize_map_v5, "assembly") <- "B73v5"

usethis::use_data(maize_map_v5, overwrite = TRUE, compress = "xz")
message(sprintf("maize_map_v5: %d loci, %d chr, %.0f cM total",
                nrow(maize_map_v5), length(unique(maize_map_v5$chr)),
                sum(tapply(maize_map_v5$cm, maize_map_v5$chr, max))))
