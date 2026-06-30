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
#' @export
read_counts <- function(path, format = c("tsv", "gatk_table", "vcf_ad"), name = NULL) {
  format <- match.arg(format)
  if (format != "tsv") stop("read_counts(): only format='tsv' is implemented (Task 4)")

  if (dir.exists(path)) {
    files <- list.files(path, pattern = "\\.tsv$", full.names = TRUE)
    return(do.call(rbind, lapply(files, read_counts, format = "tsv")))
  }

  d <- utils::read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                         col.names = c("chr", "pos", "ref", "n_ref", "alt", "n_alt"))
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

#' Write segment calls in the common schema
#'
#' Columns: `source, donor, name, chr, start_bp, end_bp, state` (0/1/2). `name`
#' is the NIL id (pedigree string) per the doc terminology convention.
#'
#' @param calls A data.frame of segment calls.
#' @param path Output CSV path.
#' @return `path`, invisibly.
#' @export
write_common_schema <- function(calls, path) {
  cols <- c("source", "donor", "name", "chr", "start_bp", "end_bp", "state")
  if (!all(cols %in% names(calls)))
    stop("write_common_schema(): calls must have columns ", paste(cols, collapse = ", "))
  utils::write.csv(calls[, cols], path, row.names = FALSE)
  invisible(path)
}