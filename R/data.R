#' B73 v5 maize consensus genetic map
#'
#' A consensus genetic map for maize: each marker with its chromosome, genetic
#' (cM) and physical (bp, B73 v5 assembly) position. cM-space is assembly-robust;
#' the `bp` column is tied to B73 v5 (see [load_map()], which warns on an assembly
#' override). Used as the default map for [simulate_nil()] and for cM<->Mb
#' position conversion.
#'
#' @format A data.frame with 19,486 rows (loci) and 4 columns, plus an
#'   `"assembly"` attribute (`"B73v5"`):
#' \describe{
#'   \item{locus}{marker id (character)}
#'   \item{chr}{chromosome, 1-10 (integer)}
#'   \item{cm}{genetic position within the chromosome (cM)}
#'   \item{bp}{physical position, B73 v5 (bp)}
#' }
#' @source Cleaned B73 v5 consensus map assembled in the companion zealhmm
#'   pipeline; regenerate with `data-raw/make_maize_map.R`.
#' @seealso [load_map()], [simulate_nil()]
"maize_map_v5"
