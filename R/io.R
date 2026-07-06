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
#' @param path A count TSV file (or directory of them) for `format = "tsv"`; a
#'   single `.vcf`/`.vcf.gz` for `format = "vcf_ad"`.
#' @param format One of `"tsv"` (the `chr pos ref n_ref alt n_alt` headerless
#'   layout used by the skim/BRB counts), `"vcf_ad"` (per-sample allelic depths
#'   from a biallelic VCF's `AD` FORMAT field -- e.g. the LB-Impute example data
#'   and GATK/bcftools output), or `"gatk_table"` (not yet implemented).
#' @param name Optional sample name for a single `tsv` file; defaults to the
#'   file's basename with extensions stripped. Ignored for `vcf_ad` (sample names
#'   come from the VCF `#CHROM` header).
#' @return A long observation table: `name, chr, pos, n_ref, n_alt`.
#' @examples
#' # Headerless "chr pos ref n_ref alt n_alt" TSV, as produced upstream.
#' f <- tempfile(fileext = ".tsv")
#' write.table(
#'   data.frame(chr = "chr1", pos = c(1e5, 2e5, 3e5),
#'              ref = "A", n_ref = c(8, 5, 0), alt = "T", n_alt = c(0, 3, 7)),
#'   f, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
#' read_counts(f, name = "NIL1")
#'
#' # Biallelic VCF with a GT:AD FORMAT: AD = "n_ref,n_alt" per sample.
#' v <- tempfile(fileext = ".vcf")
#' writeLines(c(
#'   "##fileformat=VCFv4.2",
#'   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNIL1\tNIL2",
#'   "1\t100000\t.\tA\tT\t.\t.\t.\tGT:AD\t0/0:8,0\t0/1:4,5",
#'   "1\t200000\t.\tC\tG\t.\t.\t.\tGT:AD\t1/1:0,7\t./.:."), v)
#' read_counts(v, format = "vcf_ad")
#' @export
read_counts <- function(path, format = c("tsv", "gatk_table", "vcf_ad"), name = NULL) {
  format <- match.arg(format)
  if (format == "vcf_ad") return(.read_counts_vcf_ad(path))
  if (format == "gatk_table")
    stop("read_counts(): format = 'gatk_table' is not yet implemented; ",
         "use 'tsv' or 'vcf_ad'.")

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

# Read per-sample allelic depths from a biallelic VCF's AD FORMAT field into the
# common (name, chr, pos, n_ref, n_alt) table. AD is "n_ref,n_alt" (GATK/bcftools);
# a missing/`.` field or absent alt becomes 0. The AD position within FORMAT is
# located per record (usually constant -> a vectorized fast path); multiallelic AD
# keeps only the first ALT (biallelic assumption). Sample names come from #CHROM.
.read_counts_vcf_ad <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  repeat {                                          # skip ## meta to the #CHROM header
    ln <- readLines(con, n = 1L)
    if (!length(ln)) stop("read_counts(vcf_ad): no #CHROM header found in ", path)
    if (startsWith(ln, "#CHROM")) break
  }
  hdr <- strsplit(sub("^#", "", ln), "\t", fixed = TRUE)[[1]]
  if (length(hdr) < 10L) stop("read_counts(vcf_ad): VCF has no sample (FORMAT) columns")
  snames <- hdr[10:length(hdr)]
  body <- utils::read.table(con, sep = "\t", header = FALSE, quote = "",
                            comment.char = "", stringsAsFactors = FALSE, colClasses = "character")
  if (!nrow(body)) stop("read_counts(vcf_ad): VCF has no records")
  chr <- as.integer(sub("^chr", "", body[[1]])); pos <- as.integer(body[[2]])

  # AD index within the colon-delimited FORMAT (col 9), per record.
  ai <- vapply(strsplit(body[[9]], ":", fixed = TRUE),
               function(x) { m <- match("AD", x); if (is.na(m)) NA_integer_ else m }, integer(1))
  if (all(is.na(ai))) stop("read_counts(vcf_ad): no `AD` field found in the FORMAT column")
  if (anyNA(ai))
    warning("read_counts(vcf_ad): ", sum(is.na(ai)), " record(s) have no AD field in ",
            "FORMAT; their genotypes are read as 0 counts.")
  k_const <- length(unique(ai)) == 1L && !anyNA(ai)   # constant FORMAT -> fast vectorized path

  # Extract the k-th ':'-field of each sample cell (k scalar when FORMAT is
  # constant). Both paths guard truncated cells (fewer than k fields) -> NA, so a
  # cell missing AD never falls through to an earlier (possibly comma-bearing, e.g.
  # PL "0,255,255") field being misread as AD.
  kth_field <- function(col, k) {
    if (length(k) == 1L) {
      field  <- sub(sprintf("^(?:[^:]*:){%d}([^:]*).*", k - 1L), "\\1", col, perl = TRUE)
      ncolon <- nchar(col) - nchar(gsub(":", "", col, fixed = TRUE))
      field[ncolon < (k - 1L)] <- NA_character_       # k-th field absent (truncated cell)
      field
    } else {
      vapply(seq_along(col), function(r) {
        p <- strsplit(col[r], ":", fixed = TRUE)[[1]]
        if (is.na(k[r]) || k[r] > length(p)) NA_character_ else p[k[r]]
      }, character(1))
    }
  }
  ad_to_counts <- function(ad) {
    r <- sub(",.*", "", ad)                          # before first comma
    a <- sub(",.*", "", sub("^[^,]*,?", "", ad))     # first ALT depth (drops extra alt cols)
    nref <- suppressWarnings(as.integer(r)); nalt <- suppressWarnings(as.integer(a))
    nref[is.na(nref)] <- 0L; nalt[is.na(nalt)] <- 0L  # "." / missing / GT-only -> 0
    list(n_ref = nref, n_alt = nalt)
  }
  k_use <- if (k_const) ai[1] else ai
  parts <- lapply(seq_along(snames), function(j) {
    cnt <- ad_to_counts(kth_field(body[[9L + j]], k_use))
    data.frame(name = snames[j], chr = chr, pos = pos,
               n_ref = cnt$n_ref, n_alt = cnt$n_alt, stringsAsFactors = FALSE)
  })
  if (requireNamespace("data.table", quietly = TRUE))
    as.data.frame(data.table::rbindlist(parts)) else do.call(rbind, parts)
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

#' Write imputed per-marker genotypes as a VCF (LB-Impute-style deliverable)
#'
#' Optional, decoupled from the engine: turns a per-marker state table (from
#' [call_states()], typically `caller = "lbimpute"`) into an imputed biallelic
#' VCF -- LB-Impute's native output form (missing filled, false homozygotes
#' corrected), as opposed to nilHMM's default segment schema. State -> genotype is
#' `0 -> 0/0` (REF-hom), `1 -> 0/1` (het), `2 -> 1/1` (ALT-hom); any marker x
#' sample absent from `states` is written `./.`. This never re-reads the input
#' VCF; supply per-marker `REF`/`ALT` alleles via `markers` if you want the real
#' alleles rather than the `N` placeholder.
#'
#' @param states A per-marker state table with columns `name, chr, pos, state`
#'   (0/1/2), e.g. from `call_states(..., caller = "lbimpute")`.
#' @param path Output `.vcf` path.
#' @param markers Optional data.frame keyed by `chr, pos` carrying `ref`, `alt`
#'   (and optional `id`) alleles to emit; when `NULL`, REF/ALT default to `N`.
#' @param ref,alt Fallback single-character REF/ALT alleles when `markers` is
#'   `NULL` (default `"N"`).
#' @return `path`, invisibly.
#' @examples
#' st <- data.frame(name = c("NIL1", "NIL1", "NIL2", "NIL2"),
#'                  chr = 1L, pos = c(1e5, 2e5, 1e5, 2e5),
#'                  state = c(0L, 2L, 1L, 0L))
#' write_vcf_impute(st, tempfile(fileext = ".vcf"))
#' @export
write_vcf_impute <- function(states, path, markers = NULL, ref = "N", alt = "N") {
  need <- c("name", "chr", "pos", "state")
  if (!all(need %in% names(states)))
    stop("write_vcf_impute(): `states` needs columns ", paste(need, collapse = ", "))
  st <- as.data.frame(states, stringsAsFactors = FALSE)
  if (anyDuplicated(paste(st$name, st$chr, st$pos, sep = "\r")))
    stop("write_vcf_impute(): duplicate (name, chr, pos) rows in `states`")
  samples <- sort(unique(st$name))
  mk <- unique(st[, c("chr", "pos")])
  mk <- mk[order(mk$chr, mk$pos), , drop = FALSE]
  M <- nrow(mk); S <- length(samples)

  # marker x sample state grid (NA where a sample has no call at a marker)
  mkey <- paste(mk$chr, mk$pos, sep = "\r")
  grid <- matrix(NA_integer_, M, S, dimnames = list(NULL, samples))
  ri <- match(paste(st$chr, st$pos, sep = "\r"), mkey)
  ci <- match(st$name, samples)
  grid[cbind(ri, ci)] <- as.integer(st$state)
  gt <- c("0/0", "0/1", "1/1")[grid + 1L]
  gt[is.na(gt)] <- "./."
  dim(gt) <- c(M, S)

  # per-marker REF/ALT (from `markers` if supplied, else the scalar fallback)
  ref_v <- rep(ref, M); alt_v <- rep(alt, M); id_v <- rep(".", M)
  if (!is.null(markers)) {
    if (!all(c("chr", "pos", "ref", "alt") %in% names(markers)))
      stop("write_vcf_impute(): `markers` needs columns chr, pos, ref, alt")
    j <- match(mkey, paste(markers$chr, markers$pos, sep = "\r"))
    ref_v <- ifelse(is.na(j), ref, markers$ref[j])
    alt_v <- ifelse(is.na(j), alt, markers$alt[j])
    if ("id" %in% names(markers)) id_v <- ifelse(is.na(j), ".", markers$id[j])
  }

  con <- file(path, "wt"); on.exit(close(con))
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##source=nilHMM::write_vcf_impute",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Imputed genotype\">",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            "FORMAT", samples), collapse = "\t")), con)
  body <- paste(mk$chr, mk$pos, id_v, ref_v, alt_v, ".", ".", ".", "GT",
                apply(gt, 1L, paste, collapse = "\t"), sep = "\t")
  writeLines(body, con)
  invisible(path)
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