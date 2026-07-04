# Closed-form per-site GATK-style genotype-likelihood (GL) caller. A data-agnostic
# sibling of interpolate_genotype(): counts in -> hard 0/1/2 calls out, with a
# swappable single-locus prior. NOT an HMM -- there is no linkage/duration; each
# (marker, sample) is decided independently from its own (n_ref, n_alt).
#
# This is the deliberately het-EXCESS control for the coverage-degradation study:
# under the HWE prior a single ALT read is called HET (the 2p(1-p) prior term beats
# the p^2 hom-ALT term at low depth), so low coverage inflates het calls -- the
# failure mode we want to contrast against the duration-aware HMM callers.

#' Per-site genotype-likelihood (GL) genotype caller with a swappable prior
#'
#' A closed-form GATK-style diploid genotype call from allelic read counts, decided
#' independently per (marker, sample) -- no linkage model. ALT = donor (teosinte) =
#' dosage 2. Given per-read error `error` and counts `(n_ref, n_alt)`, the three
#' genotype log-likelihoods are
#' \deqn{L(0) = (1-\epsilon)^{n_{ref}}\,\epsilon^{n_{alt}}}
#' \deqn{L(2) = \epsilon^{n_{ref}}\,(1-\epsilon)^{n_{alt}}}
#' \deqn{L(1) = 0.5^{\,n_{ref}+n_{alt}}}
#' and the call is `argmax_g L(g) * pi(g)` over the prior `pi`. All arithmetic is in
#' log space and vectorized over the whole matrix (no per-cell loop). Markers with
#' zero total depth carry no signal and are returned as `NA`.
#'
#' @details
#' Priors (`prior`):
#' \describe{
#'   \item{`"hwe"`}{Hardy-Weinberg from the per-marker ALT (donor/teosinte) allele
#'     frequency `p`: `pi = {(1-p)^2, 2p(1-p), p^2}`. This is the **het-excess**
#'     control -- a single ALT read is called HET. `af` supplies `p` per marker; if
#'     `NULL`, `p` is estimated per marker from the reads
#'     (`rowSums(n_alt) / rowSums(n_ref + n_alt)`).}
#'   \item{`"flat"`}{Uniform prior -> pure argmax-GL. Het-blind at depth 1 (a single
#'     ALT read is called hom-ALT), the pessimistic naive caller.}
#'   \item{`"breeding"`}{Fixed design frequencies `f = c(f_REF, f_HET, f_ALT)`
#'     (e.g. BC1S4 ~ `c(0.84, 0.06, 0.10)`); renormalized to sum 1.}
#' }
#'
#' @param n_ref,n_alt Reference / alternate read counts. Integer (or numeric)
#'   matrices `markers x samples`, or plain vectors (treated as one sample =
#'   `markers x 1`). Same shape required.
#' @param prior One of `"hwe"` (default), `"flat"`, `"breeding"`.
#' @param af For `prior = "hwe"`: per-marker ALT allele frequency `p`, length
#'   `nrow`. `NULL` (default) estimates `p` from the reads.
#' @param error Per-read error rate `epsilon` in `(0, 0.5)` (default `0.01`).
#' @param f For `prior = "breeding"`: length-3 genotype frequencies
#'   `c(f_REF, f_HET, f_ALT)` (renormalized).
#' @param return `"call"` (default): integer 0/1/2 hard call, `NA` at zero depth.
#'   `"dosage"`: the same hard call as numeric (convenience for
#'   [interpolate_genotype()]). `"post"`: the normalized posterior array
#'   `markers x samples x 3` (`NA` rows at zero depth).
#' @return For `"call"`/`"dosage"`: a matrix `markers x samples` (a plain vector if
#'   the inputs were vectors). For `"post"`: a `markers x samples x 3` array.
#' @examples
#' # A single ALT read: het-blind under flat, HET under HWE (het-excess control).
#' call_gl(0, 1, prior = "flat")               # 2 (hom-ALT)
#' call_gl(0, 1, prior = "hwe", af = 0.3)       # 1 (HET)
#' # High depth: the data dominates the prior.
#' call_gl(0, 10, prior = "hwe", af = 0.05)     # 2
#' @seealso [interpolate_genotype()], [call_ancestry()]
#' @export
call_gl <- function(n_ref, n_alt,
                    prior = c("hwe", "flat", "breeding"),
                    af = NULL, error = 0.01, f = NULL,
                    return = c("call", "dosage", "post")) {
  prior  <- match.arg(prior)
  return <- match.arg(return)
  if (!is.numeric(error) || length(error) != 1L || error <= 0 || error >= 0.5)
    stop("call_gl(): `error` must be a single value in (0, 0.5).")

  was_vec <- is.null(dim(n_ref)) && is.null(dim(n_alt))
  R <- as.matrix(n_ref); A <- as.matrix(n_alt)
  storage.mode(R) <- "double"; storage.mode(A) <- "double"
  if (!identical(dim(R), dim(A)))
    stop("call_gl(): `n_ref` and `n_alt` must have the same shape.")
  if (anyNA(R) || anyNA(A))
    stop("call_gl(): `n_ref`/`n_alt` must not contain NA (use 0 for no reads).")
  M <- nrow(R); N <- ncol(R)

  depth <- R + A
  le <- log(error); lc <- log1p(-error)          # log(error), log(1 - error)

  # genotype log-likelihoods (markers x samples)
  logL0 <- R * lc + A * le                        # REF hom
  logL2 <- R * le + A * lc                        # ALT hom
  logL1 <- depth * log(0.5)                        # HET

  # per-marker log priors (length M; recycled column-wise across the M x N matrix)
  lp <- .gl_log_prior(prior, M, R, A, af, f)

  P0 <- logL0 + lp$l0
  P1 <- logL1 + lp$l1
  P2 <- logL2 + lp$l2

  if (return == "post") {
    mx <- pmax(P0, P1, P2)
    e0 <- exp(P0 - mx); e1 <- exp(P1 - mx); e2 <- exp(P2 - mx)
    s  <- e0 + e1 + e2
    post <- array(NA_real_, dim = c(M, N, 3),
                  dimnames = list(rownames(R), colnames(R), c("0", "1", "2")))
    post[, , 1] <- e0 / s; post[, , 2] <- e1 / s; post[, , 3] <- e2 / s
    post[depth == 0] <- NA_real_               # recycles the M x N mask over all 3 slices
    return(post)
  }

  # argmax over the three log-posteriors; ties resolve to the lower state index
  # (which.max semantics: >= keeps the first maximizer).
  call <- ifelse(P0 >= P1 & P0 >= P2, 0L, ifelse(P1 >= P2, 1L, 2L))
  call[depth == 0] <- NA_integer_
  dim(call) <- c(M, N); dimnames(call) <- dimnames(R)

  if (return == "dosage") storage.mode(call) <- "double"
  if (was_vec) call <- call[, 1L]
  call
}

# Per-marker log-prior vectors (each length M) for the three genotypes. R/A are the
# count matrices, used only to estimate `af` from the data when prior == "hwe" and
# `af` is NULL.
.gl_log_prior <- function(prior, M, R, A, af, f) {
  eps <- 1e-12
  if (prior == "flat") {
    v <- rep(log(1 / 3), M)
    return(list(l0 = v, l1 = v, l2 = v))
  }
  if (prior == "hwe") {
    if (is.null(af)) {                           # estimate p per marker from the reads
      tot <- rowSums(R + A)
      p <- ifelse(tot > 0, rowSums(A) / tot, 0.5)
    } else {
      p <- af
    }
    if (length(p) != M)
      stop("call_gl(): `af` length (", length(p),
           ") must equal the number of markers (", M, ").")
    if (anyNA(p) || any(p < 0 | p > 1))
      stop("call_gl(): `af` must be non-NA and in [0, 1].")
    p <- pmin(pmax(p, eps), 1 - eps)             # clamp to keep log finite
    return(list(l0 = 2 * log1p(-p),
                l1 = log(2) + log(p) + log1p(-p),
                l2 = 2 * log(p)))
  }
  # breeding
  if (is.null(f) || length(f) != 3 || any(!is.finite(f)) || any(f < 0))
    stop("call_gl(): prior = 'breeding' needs `f` = c(f_REF, f_HET, f_ALT), ",
         "non-negative and length 3.")
  f <- f / sum(f)
  f <- pmin(pmax(f, eps), 1)
  list(l0 = rep(log(f[1]), M), l1 = rep(log(f[2]), M), l2 = rep(log(f[3]), M))
}