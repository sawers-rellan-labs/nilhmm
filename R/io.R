# I/O (S8). Readers return marker-level observation tables; the writer emits the
# common segment schema. Data-agnostic: no hardcoded paths. Pipeline scripts
# pass the file lists in.

#' Read allelic counts into the engine's observation table
#'
#' Reads CollectAllelicCounts-style read counts (the GATK-table readcount
#' standard, memory `gatk-table-readcount-standard`) into a per-sample,
#' per-marker `(n_ref, n_alt, chr, pos)` table for the `count` emission. A
#' TSV->counts adapter covers the wideseq-thinned BRB inputs (cf. the Python
#' `agent/brb_nilhmm_counts.py`).
#'
#' @param path A count TSV file, or a directory of them.
#' @param format One of `"tsv"` (the `chr pos ref n_ref alt n_alt` headerless
#'   layout used by the skim/BRB counts), `"gatk_table"`, `"vcf_ad"`.
#' @param name Optional sample name for a single file; defaults to the file's
#'   basename with extensions stripped.
#' @return A long observation table: `name, chr, pos, n_ref, n_alt`.
#' @examples
#' # Headerless "chr pos ref n_ref alt n_alt" TSV, as produced upstream.
#' f <- tempfile(fileext = ".tsv")
#' write.table(
#'   data.frame(chr = "chr1", pos = c(1e5, 2e5, 3e5),
#'              ref = "A", n_ref = c(8, 5, 0), alt = "T", n_alt = c(0, 3, 7)),
#'   f, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
#' read_counts(f, name = "NIL1")
#' @export
read_counts <- function(path, format = c("tsv", "gatk_table", "vcf_ad"), name = NULL) {
  format <- match.arg(format)
  if (format != "tsv") stop("read_counts(): only format='tsv' is implemented (Task 4)")

  if (dir.exists(path)) {
    files <- list.files(path, pattern = "\\.tsv(\\.gz)?$", full.names = TRUE)
    parts <- lapply(files, read_counts, format = "tsv")
    # rbindlist avoids the O(n^2) copying of do.call(rbind, .) over many files.
    if (requireNamespace("data.table", quietly = TRUE)) {
      return(as.data.frame(data.table::rbindlist(parts)))
    }
    return(do.call(rbind, parts))
  }

  # fread is ~10x faster than read.table here (and reads .gz directly); the base
  # reader is the fallback when data.table is not installed.
  cn <- c("chr", "pos", "ref", "n_ref", "alt", "n_alt")
  if (requireNamespace("data.table", quietly = TRUE)) {
    d <- data.table::fread(path, sep = "\t", header = FALSE, col.names = cn,
                           colClasses = list(character = 1L, integer = c(2L, 4L, 6L)),
                           showProgress = FALSE)
  } else {
    d <- utils::read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                           col.names = cn)
  }
  nm <- if (!is.null(name)) name else sub("\\..*$", "", basename(path))
  data.frame(
    name  = nm,
    chr   = as.integer(sub("^chr", "", d$chr)),
    pos   = as.integer(d$pos),
    n_ref = as.integer(d$n_ref),
    n_alt = as.integer(d$n_alt),
    stringsAsFactors = FALSE
  )
}

#' Read hard genotype calls (VCF GT) into the engine's observation table
#'
#' Extracts the per-sample `GT` field of a **biallelic diploid** VCF into a
#' `(name, chr, pos, g)` table for the categorical `gt` emission -- Holland's
#' nNIL genotype path, for the saturated-depth regime (e.g. MolBreeding target sequencing at
#' ~20x or more) where the caller's *called* genotype is trustworthy. `g` is the
#' alt-allele dosage: `0` REF-hom, `1` het, `2` ALT-hom, `3` missing. Unlike
#' [read_counts()] (which feeds the count/BetaBinomial emission from `AD` read
#' depths), this reads the called genotype directly, so it never touches `AD`.
#'
#' Run with `call_ancestry(read_vcf_gt(path), caller = "nnil", emission = "gt",
#' design = ...)` (a `g`-only input auto-selects the gt emission).
#'
#' @param path A `.vcf` or `.vcf.gz` file. `GT` must be the first FORMAT field
#'   (the VCF spec requirement when GT is present).
#' @param samples Optional character vector to restrict to a subset of samples.
#' @return A long observation table: `name, chr, pos, g`.
#' @examples
#' # A minimal biallelic-diploid VCF with a GT field.
#' f <- tempfile(fileext = ".vcf")
#' writeLines(c(
#'   "##fileformat=VCFv4.2",
#'   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNIL1",
#'   "1\t100000\t.\tA\tT\t.\t.\t.\tGT\t0/0",
#'   "1\t200000\t.\tA\tT\t.\t.\t.\tGT\t0/1",
#'   "1\t300000\t.\tA\tT\t.\t.\t.\tGT\t1/1"), f)
#' read_vcf_gt(f)   # g in {0 REF-hom, 1 het, 2 ALT-hom, 3 missing}
#' @export
read_vcf_gt <- function(path, samples = NULL) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  repeat {                                          # skip ## meta lines to the #CHROM header
    ln <- readLines(con, n = 1L)
    if (!length(ln)) stop("read_vcf_gt(): no #CHROM header found in ", path)
    if (startsWith(ln, "#CHROM")) break
  }
  hdr <- strsplit(sub("^#", "", ln), "\t", fixed = TRUE)[[1]]
  if (length(hdr) < 10L) stop("read_vcf_gt(): VCF has no sample (FORMAT) columns")
  snames <- hdr[10:length(hdr)]
  keep <- if (is.null(samples)) seq_along(snames) else which(snames %in% samples)
  if (!length(keep)) stop("read_vcf_gt(): none of `samples` present in the VCF")
  body <- utils::read.table(con, sep = "\t", header = FALSE, quote = "",
                            comment.char = "", stringsAsFactors = FALSE, colClasses = "character")
  chr <- as.integer(sub("^chr", "", body[[1]])); pos <- as.integer(body[[2]])
  # GT is the first ':'-subfield (VCF spec); alt-allele dosage over "a/b"|"a|b".
  to_g <- function(col) {
    gt <- sub(":.*", "", col)
    a1 <- sub("[/|].*", "", gt); a2 <- sub(".*[/|]", "", gt)
    miss <- a1 == "." | a2 == "."
    g <- (a1 != "0" & a1 != ".") + (a2 != "0" & a2 != ".")   # 0 / 1 / 2 alt alleles
    g <- as.integer(g); g[miss] <- 3L; g
  }
  do.call(rbind, lapply(keep, function(j)
    data.frame(name = snames[j], chr = chr, pos = pos, g = to_g(body[[9 + j]]),
               stringsAsFactors = FALSE)))
}

#' Write segment calls in the common schema
#'
#' Columns: `source, donor, name, chr, start_bp, end_bp, state` (0/1/2). `name`
#' is the NIL id (pedigree string) per the doc terminology convention.
#'
#' @param calls A data.frame of segment calls.
#' @param path Output CSV path.
#' @return `path`, invisibly.
#' @examples
#' calls <- data.frame(source = "nilHMM", donor = NA, name = "NIL1", chr = 1L,
#'                     start_bp = 1L, end_bp = 5e6L, state = 2L)
#' write_common_schema(calls, tempfile(fileext = ".csv"))
#' @export
write_common_schema <- function(calls, path) {
  cols <- c("source", "donor", "name", "chr", "start_bp", "end_bp", "state")
  if (!all(cols %in% names(calls)))
    stop("write_common_schema(): calls must have columns ", paste(cols, collapse = ", "))
  utils::write.csv(calls[, cols], path, row.names = FALSE)
  invisible(path)
}