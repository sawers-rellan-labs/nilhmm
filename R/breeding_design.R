# Breeding-design axis (S6). Orthogonal to the emission axis: a dataset carries
# BOTH a depth regime (-> emission) and a breeding design (-> generation-derived
# quantities). e.g. skim = count x BC2S2, BRB = count x BC2S3.
#
# The breeding design determines the single-locus genotype-freq priors (f_1, f_2)
# and the expected fragment-size law (the fitted Gamma). The Gamma is CALIBRATION/NULL
# ONLY -- NEVER an engine prior (S7 circularity trap): baking it into the
# transition would make the caller reproduce it by construction and mask the
# open "callers run longer than the model" signal.

#' Parse a `"BCnSm"` breeding-design string into its generation counts
#'
#' @param design Design key of the form `"BC<n>S<m>"` (e.g. `"BC2S2"`, `"BC1S4"`).
#' @return `list(n_bc, n_self)` -- backcross and selfing generation counts.
#' @examples
#' parse_design("BC2S2")   # list(n_bc = 2, n_self = 2)
#' @export
parse_design <- function(design) {
  mm <- regmatches(design, regexec("^BC(\\d+)S(\\d+)$", design))[[1]]
  if (!length(mm))
    stop("parse_design(): design must be 'BC<n>S<m>' (e.g. BC1S4, BC2S2): ", design)
  list(n_bc = as.integer(mm[2]), n_self = as.integer(mm[3]))
}

#' Single-locus (Mendelian) genotype expectation for a breeding design
#'
#' Propagates the single-locus genotype-frequency vector through the backcross and
#' selfing transition matrices for the design named by `design` (`"BC<n>S<m>"`):
#' `n` backcrosses to the recurrent parent, then `m` generations of selfing. The
#' genotype basis is `(AA, Aa, aa)` with `AA` the recurrent (backcross-target)
#' homozygote and `aa` the donor homozygote, so the recurrent allele is `A`
#' (frequency `p_A`) and the donor allele is `a`. The donor allele frequency after
#' the backcrosses is `p_a = (aa + Aa/2) * 0.5^n` in terms of the starting F1
#' composition `f1 = (AA, Aa, aa)` -- `0.5^(n + 1)` for the default fully-het F1 --
#' and selfing leaves it invariant. The backcross matrix `B`
#' mates each generation to `AA` (`aa -> Aa`, `Aa -> half AA / half Aa`) and the
#' selfing matrix `S` splits a quarter of the remaining `Aa` into equal parts
#' `AA`/`aa` each generation.
#'
#' The starting composition `f1` is a free argument, so the expectation is exact
#' for a non-inbred F1, not only the fully-heterozygous `c(0, 1, 0)` case that
#' inbred parents give. At the inbred default, BC2S2 lands at REF 0.844 /
#' HET 0.0625 / ALT 0.0938.
#'
#' @param design Design key of the form `"BC<n>S<m>"` (e.g. `"BC2S2"`, `"BC1S4"`).
#' @param f1 Starting genotype-frequency vector in `(AA, Aa, aa)` order (recurrent
#'   hom, het, donor hom). Defaults to `c(0, 1, 0)`, the fully-heterozygous F1 of
#'   inbred parents.
#' @return Named numeric vector `c(REF, HET, ALT)` summing to 1, mapping the
#'   `(AA, Aa, aa)` genotypes to the REF/HET/ALT ancestry states.
#' @examples
#' single_locus_expectation("BC2S2")                       # inbred parents
#' single_locus_expectation("BC2S2", f1 = c(0.5, 0.5, 0))  # non-inbred F1
#' @seealso [parse_design()], [breeding_prior()]
#' @export
single_locus_expectation <- function(design, f1 = c(0, 1, 0)) {
  d <- parse_design(design)
  B <- rbind(c(1, 0, 0), c(0.5, 0.5, 0), c(0, 1, 0))     # backcross to AA (recurrent)
  S <- rbind(c(1, 0, 0), c(0.25, 0.5, 0.25), c(0, 0, 1)) # selfing
  v <- as.numeric(f1)                                    # (AA, Aa, aa)
  for (i in seq_len(d$n_bc))   v <- as.numeric(v %*% B)
  for (i in seq_len(d$n_self)) v <- as.numeric(v %*% S)
  c(REF = v[1], HET = v[2], ALT = v[3])                  # AA=recurrent=REF, aa=donor=ALT
}

#' Genotype-frequency prior implied by a breeding design
#'
#' The public design-prior contract: the `c(REF, HET, ALT)` vector of single-locus
#' genotype frequencies the breeding scheme implies, consumed directly by
#' [call_gt()] as its `prior` and by the engine as its state frequencies. A thin
#' wrapper over [single_locus_expectation()] -- callers should depend on this stable
#' name rather than on the genetics primitive behind it.
#'
#' @param design Design key of the form `"BC<n>S<m>"` (e.g. `"BC2S2"`, `"BC2S3"`).
#' @param f1 Starting F1 genotype-frequency vector, passed through to
#'   [single_locus_expectation()]; defaults to inbred parents `c(0, 1, 0)`.
#' @return Named numeric length-3 vector `c(REF, HET, ALT)` summing to 1.
#' @examples
#' breeding_prior("BC2S3")                          # c(REF = .8594, HET = .0312, ALT = .1094)
#' call_gt(0, 1, prior = breeding_prior("BC2S3"))   # design prior resists the het flip -> 2
#' @seealso [single_locus_expectation()], [call_gt()]
#' @export
breeding_prior <- function(design, f1 = c(0, 1, 0)) {
  single_locus_expectation(design, f1 = f1)
}

# Resolve the engine's state frequencies (Holland's f_1 = HET, f_2 = donor-hom)
# from either a design (via breeding_prior) or explicit f_1/f_2, else error. One
# place for the pattern the engine call sites share.
.state_freqs <- function(design, f_1 = NULL, f_2 = NULL, ctx = "call_states") {
  if (!is.null(design)) {
    p <- breeding_prior(design)
    list(f_1 = unname(p[["HET"]]), f_2 = unname(p[["ALT"]]))
  } else if (!is.null(f_1) && !is.null(f_2)) {
    if (f_1 + f_2 >= 1) stop(ctx, "(): `f_1` + `f_2` must be < 1")
    list(f_1 = f_1, f_2 = f_2)
  } else {
    stop(ctx, "(): supply `design` or both `f_1` and `f_2`")
  }
}

#' Expected fragment-size distribution (the Null / KS target)
#'
#' Returns the fitted Gamma `(k, lambda)` for a design plus density/CDF
#' closures, for plotting the grey Null and as the KS-vs-sim calibration target.
#' Never feeds the engine transition.
#'
#' @param design Design key.
#' @return `list(k, lambda, density, cdf)`.
#' @examples
#' \dontrun{
#' # Planned (Task 4): the grey Null / KS target for a design.
#' nd <- expected_fragment_dist("BC2S3")
#' }
#' @export
expected_fragment_dist <- function(design) {
  stop("nilHMM::expected_fragment_dist() not yet implemented (Task 4)")
}

#' Fit a design Gamma from simulated segments
#'
#' Extensible entry point used by `data-raw/make_breeding_designs.R` to derive
#' the bundled `(k, lambda)` from simcross output.
#'
#' @param sim_segments data.frame of simulated segment sizes (cM).
#' @return `list(k, lambda)`.
#' @examples
#' \dontrun{
#' # Planned (Task 4): fit the design Gamma from simcross segment sizes.
#' fit_design_gamma(sim_segments)
#' }
#' @export
fit_design_gamma <- function(sim_segments) {
  stop("nilHMM::fit_design_gamma() not yet implemented (Task 4)")
}
