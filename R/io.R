# I/O (§8). Readers return marker-level observation tables; the writer emits the
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
#' @param path File or directory of counts.
#' @param format One of `"gatk_table"`, `"tsv"`, `"vcf_ad"`.
#' @return A counts observation table.
#' @export
read_counts <- function(path, format = c("gatk_table", "tsv", "vcf_ad")) {
  format <- match.arg(format)
  stop("nilHMM::read_counts() not yet implemented (Task 4)")
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
  stop("nilHMM::write_common_schema() not yet implemented (Task 4)")
}