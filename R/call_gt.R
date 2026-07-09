# Closed-form per-site genotype caller. A data-agnostic sibling of
# interpolate_genotype(): counts in -> hard 0/1/2 genotype calls out, with a
# swappable single-locus prior. NOT an HMM -- there is no linkage/duration; each
# (marker, sample) is decided independently from its own (n_ref, n_alt). The
# genotype-likelihood (GL) stage it shares with GATK/ANGSD/bcftools is the
# emission, not the purpose: this calls a *genotype*, hence `call_gt`.
#
# The HWE prior makes this the deliberately het-EXCESS control for the
# coverage-degradation study: under HWE a single ALT read is called HET (the
# 2p(1-p) prior term beats the p^2 hom-ALT term at low depth), so low coverage
# inflates het calls -- the failure mode we want to contrast against the
# duration-aware HMM callers.

#' Per-site genotype caller with a swappable prior
#'
#' A closed-form GATK-style diploid genotype call from allelic read counts, decided
#' independently per (marker, sample) -- no linkage model. ALT = donor (teosinte) =
#' dosage 2. Given per-read error `error` and counts `(n_ref, n_alt)`, the three
#' genotype log-likelihoods (GL) are
#' \deqn{L(0) = (1-\epsilon)^{n_{ref}}\,\epsilon^{n_{alt}}}
#' \deqn{L(2) = \epsilon^{n_{ref}}\,(1-\epsilon)^{n_{alt}}}
#' \deqn{L(1) = 0.5^{\,n_{ref}+n_{alt}}}
#' and the call is `argmax_g L(g) * pi(g)` over the prior `pi` (the MAP / argmax-GP
#' estimate). All arithmetic is in log space and vectorized over the whole matrix
#' (no per-cell loop). Markers with zero total depth carry no signal and are
#' returned as `NA`.
#'
#' @details
#' The prior *is* the argument. `prior` is polymorphic:
#' \describe{
#'   \item{`"flat"`}{Uniform prior -> pure argmax-GL (the ML call). Het-blind at
#'     depth 1 (a single ALT read is called hom-ALT), the pessimistic naive caller.}
#'   \item{`"hwe"`}{Hardy-Weinberg from the per-marker ALT (donor/teosinte) allele
#'     frequency `p`: `pi = {(1-p)^2, 2p(1-p), p^2}`. This is the **het-excess**
#'     control -- a single ALT read is called HET. `af` supplies `p` per marker; if
#'     `NULL`, `p` is self-estimated per marker from the reads
#'     (`rowSums(n_alt) / rowSums(n_ref + n_alt)`) -- discouraged, as the estimate
#'     is circular at low depth.}
#'   \item{numeric `c(f_REF, f_HET, f_ALT)`}{A **fixed** genome-wide prior,
#'     renormalized to sum 1. Covers both the **design** prior (a vector derived
#'     from the cross -- see [design_prior()]) and an arbitrary **custom** prior;
#'     the code path is identical.}
#' }
#' A length-3 numeric vector is genome-wide fixed -- it cannot express a per-marker
#' custom prior. (An `M x 3` matrix could be accepted for that later; out of scope
#' now.)
#'
#' @param n_ref,n_alt Reference / alternate read counts. Integer (or numeric)
#'   matrices `markers x samples`, or plain vectors (treated as one sample =
#'   `markers x 1`). Same shape required.
#' @param prior The single-locus prior: `"hwe"` (default), `"flat"`, or a length-3
#'   numeric vector `c(f_REF, f_HET, f_ALT)` (fixed genome-wide, renormalized).
#' @param af For `prior = "hwe"`: per-marker ALT allele frequency `p`, length
#'   `nrow`. `NULL` (default) self-estimates `p` from the reads (discouraged).
#' @param error Per-read error rate `epsilon` in `(0, 0.5)` (default `0.01`).
#' @param return `"call"` (default): integer 0/1/2 hard call, `NA` at zero depth.
#'   `"dosage"`: the same hard call cast to numeric (convenience for
#'   [interpolate_genotype()]) -- this is the hardcall as a double, **not**
#'   `E[G | reads]`; for the true posterior-mean dosage take `return = "post"` and
#'   compute `gp[,,2]*1 + gp[,,3]*2`. `"post"`: the normalized genotype-posterior
#'   (GP) array `markers x samples x 3` (`NA` rows at zero depth). "GP" is the VCF
#'   FORMAT field for `P(G | reads)`; the hard `"call"` is its MAP (argmax-GP).
#' @return For `"call"`/`"dosage"`: a matrix `markers x samples` (a plain vector if
#'   the inputs were vectors). For `"post"`: the genotype-posterior (GP) array,
#'   `markers x samples x 3` (slices ordered REF/HET/ALT = dosage 0/1/2).
#' @examples
#' # A single ALT read decided under four priors:
#' call_gt(0, 1, prior = "flat")                 # 2 (hom-ALT, het-blind: argmax-GL)
#' call_gt(0, 1, prior = "hwe", af = 0.30)        # 1 (HET, the het-excess control)
#' call_gt(0, 1, prior = design_prior("BC2S3"))   # 2 (design prior resists the het flip)
#' call_gt(0, 1, prior = c(.98, .01, .01))        # custom fixed prior
#' # High depth: the data dominates the prior.
#' call_gt(0, 10, prior = "hwe", af = 0.05)       # 2
#' @seealso [design_prior()], [interpolate_genotype()], [call_ancestry()]
#' @export
call_gt <- function(n_ref, n_alt,
                    prior = "hwe",
                    af = NULL, error = 0.01,
                    return = c("call", "dosage", "post")) {
  return <- match.arg(return)

  if (!is.numeric(error) || length(error) != 1L || error <= 0 || error >= 0.5)
    stop("call_gt(): `error` must be a single value in (0, 0.5).")

  was_vec <- is.null(dim(n_ref)) && is.null(dim(n_alt))
  R <- as.matrix(n_ref); A <- as.matrix(n_alt)
  storage.mode(R) <- "double"; storage.mode(A) <- "double"
  if (!identical(dim(R), dim(A)))
    stop("call_gt(): `n_ref` and `n_alt` must have the same shape.")
  if (anyNA(R) || anyNA(A))
    stop("call_gt(): `n_ref`/`n_alt` must not contain NA (use 0 for no reads).")
  M <- nrow(R); N <- ncol(R)

  depth <- R + A
  le <- log(error); lc <- log1p(-error)          # log(error), log(1 - error)

  # genotype log-likelihoods (markers x samples)
  logL0 <- R * lc + A * le                        # REF hom
  logL2 <- R * le + A * lc                        # ALT hom
  logL1 <- depth * log(0.5)                        # HET

  # per-marker log priors (length M; recycled column-wise across the M x N matrix)
  lp <- .gt_log_prior(prior, M, R, A, af)

  P0 <- logL0 + lp$l0
  P1 <- logL1 + lp$l1
  P2 <- logL2 + lp$l2

  if (return == "post") {
    # Genotype posterior (GP): normalized P(G | reads) per (marker, sample).
    mx <- pmax(P0, P1, P2)
    e0 <- exp(P0 - mx); e1 <- exp(P1 - mx); e2 <- exp(P2 - mx)
    s  <- e0 + e1 + e2
    gp <- array(NA_real_, dim = c(M, N, 3),
                dimnames = list(rownames(R), colnames(R), c("0", "1", "2")))
    gp[, , 1] <- e0 / s; gp[, , 2] <- e1 / s; gp[, , 3] <- e2 / s
    gp[depth == 0] <- NA_real_               # recycles the M x N mask over all 3 slices
    return(gp)
  }

  # MAP call: argmax-GP over the three log-posteriors; ties resolve to the lower
  # state index (which.max semantics: >= keeps the first maximizer).
  call <- ifelse(P0 >= P1 & P0 >= P2, 0L, ifelse(P1 >= P2, 1L, 2L))
  call[depth == 0] <- NA_integer_
  dim(call) <- c(M, N); dimnames(call) <- dimnames(R)

  if (return == "dosage") storage.mode(call) <- "double"
  if (was_vec) call <- call[, 1L]
  call
}

# Per-marker log-prior vectors (each length M) for the three genotypes. `prior` is
# polymorphic: "flat", "hwe", or a length-3 numeric vector (fixed genome-wide).
# R/A are the count matrices, used only to self-estimate `af` when prior == "hwe"
# and `af` is NULL.
.gt_log_prior <- function(prior, M, R, A, af) {
  eps <- 1e-12

  # numeric length-3 vector: a fixed genome-wide prior (design or custom).
  if (is.numeric(prior)) {
    if (length(prior) != 3L || any(!is.finite(prior)) || any(prior < 0) ||
        all(prior == 0))
      stop("call_gt(): a numeric `prior` must be length 3, finite, non-negative, ",
           "and not all zero: c(f_REF, f_HET, f_ALT).")
    f <- prior / sum(prior)
    f <- pmin(pmax(f, eps), 1)
    return(list(l0 = rep(log(f[1]), M), l1 = rep(log(f[2]), M), l2 = rep(log(f[3]), M)))
  }

  prior <- match.arg(prior, c("hwe", "flat"))
  if (prior == "flat") {
    v <- rep(log(1 / 3), M)
    return(list(l0 = v, l1 = v, l2 = v))
  }
  # hwe
  if (is.null(af)) {                             # estimate p per marker from the reads
    tot <- rowSums(R + A)
    p <- ifelse(tot > 0, rowSums(A) / tot, 0.5)
  } else {
    p <- af
  }
  if (length(p) != M)
    stop("call_gt(): `af` length (", length(p),
         ") must equal the number of markers (", M, ").")
  if (anyNA(p) || any(p < 0 | p > 1))
    stop("call_gt(): `af` must be non-NA and in [0, 1].")
  p <- pmin(pmax(p, eps), 1 - eps)               # clamp to keep log finite
  list(l0 = 2 * log1p(-p),
       l1 = log(2) + log(p) + log1p(-p),
       l2 = 2 * log(p))
}