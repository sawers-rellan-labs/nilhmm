# Pedigree-aware ancestry refinement (design of record: design/PEDIGREE_HMM.md).
#
# refine_ancestry() takes a per-marker hard-state mosaic (from call_states()) plus
# a pedigree, couples relatives through the pedigree via loopy belief propagation
# (pedigree_bp_cpp), and returns refined per-marker states. Depth-blind: emission
# is emission_gt() over the hard calls. R owns all IO; the C++ kernel is
# data-agnostic. The generation prior enters ONLY at each family founder (rho/pi =
# pi_0) and propagates downward through the selfing transmission kernel Tsib
# (applying Tsib to pi_0 yields pi_1, ..., so the pi_t marginals emerge without
# double-counting); non-root chains are marginal-neutral (uniform stationary).

# donor genome fraction q and founder prior pi_0 from a design, derived (not typed):
# design_priors() gives terminal f_1 (HET) + f_2 (donor-hom); q = f_2 + f_1/2.
.founder_prior <- function(design) {
  d <- design_priors(design)
  q <- d$f_2 + d$f_1 / 2                       # donor genome fraction (selfing-invariant)
  c(1 - 2 * q, 2 * q, 0)                       # BC founder: all donor content is het
}

# backcross count from "BC{n}S{m}"
.n_bc <- function(design) {
  mm <- regmatches(design, regexec("^BC(\\d+)S(\\d+)$", design))[[1]]
  if (!length(mm)) stop("refine_ancestry(): design must be 'BC<n>S<m>', got '", design, "'")
  as.integer(mm[2])
}

# per-interval recombination fraction from genetic (cm, Haldane) or physical (rrate*bp) gaps
.interval_r <- function(pos, cm = NULL, rrate = 0.01) {
  if (!is.null(cm) && !all(is.na(cm))) {
    d <- pmax(0, diff(cm)) / 100                     # Morgans
    return(pmin(0.5, 0.5 * (1 - exp(-2 * d))))       # Haldane
  }
  pmin(0.5, rrate * pmax(0, diff(as.numeric(pos))))
}

#' Refine per-individual ancestry calls over a pedigree
#'
#' Couples relatives through the pedigree to correct per-individual `call_states()`
#' calls: structured loopy belief propagation over the pedigree x genome grid
#' (see [pedigree_bp_cpp()], design/PEDIGREE_HMM.md). Depth-blind refinement --
#' emission is [emission_gt()] over the hard `state` calls. Families are processed
#' independently; latent ungenotyped ancestors (taxa named as parents but absent
#' from `mosaic`) impose chromosome continuity across siblings.
#'
#' @param mosaic Per-marker state table from [call_states()]: needs `name, chr,
#'   pos, state` (+ optional `source, donor, cm`). `state` in `{0 REF, 1 HET,
#'   2 ALT, 3 missing}`.
#' @param pedigree A pedigree: a path (read via [read_pedigree()]) or a
#'   [read_pedigree()]-shaped data.frame (`taxon, family, parent1, parent2`).
#'   `taxon` joins to `mosaic$name`; a `taxon` used as a parent but absent from
#'   `mosaic` is a latent ancestor.
#' @param design Breeding design `"BC{n}S{m}"` -> founder prior `pi_0` and
#'   per-node `meioses` (via [design_priors()]).
#' @param err,gert [emission_gt()] genotyping-error rates (call-level, not raw
#'   sequencing error -- see the depth caveat in design/PEDIGREE_HMM.md).
#' @param rrate Per-bp recombination fraction applied to `pos` gaps when `mosaic`
#'   has no `cm` column (else Haldane on `cm`).
#' @param maxiter,tol,lambda BP sweeps, convergence tolerance, damping.
#' @param ped_format Passed to [read_pedigree()] when `pedigree` is a path.
#' @return `mosaic` with `state` replaced by the refined per-marker calls (same
#'   columns and row order; genotyped leaves only). Feed to [to_segments()].
#' @seealso [call_states()], [simulate_family()], [pedigree_bp_cpp()], [to_segments()].
#' @examples
#' \dontrun{
#' sim   <- simulate_family("BC2S3", families = 10, sibs = 10, seed = 1)
#' obs   <- simulate_counts(sim$truth, depth = 0.5, seed = 1)
#' calls <- call_states(obs, caller = "rtiger", design = "BC2S3")
#' ref   <- refine_ancestry(calls, sim$pedigree, design = "BC2S3")
#' }
#' @export
refine_ancestry <- function(mosaic, pedigree, design = "BC2S3",
                            err = 0.05, gert = 0.10, rrate = 0.01,
                            maxiter = 30L, tol = 1e-4, lambda = 0.5,
                            ped_format = c("fam", "fsfhap")) {
  mo <- as.data.frame(mosaic, stringsAsFactors = FALSE)
  if (!all(c("name", "chr", "pos", "state") %in% names(mo)))
    stop("refine_ancestry(): `mosaic` needs columns name, chr, pos, state")
  ped <- if (is.character(pedigree)) read_pedigree(pedigree, match.arg(ped_format))
         else as.data.frame(pedigree, stringsAsFactors = FALSE)
  if (!all(c("taxon", "family", "parent1") %in% names(ped)))
    stop("refine_ancestry(): `pedigree` needs columns taxon, family, parent1")

  pi0    <- .founder_prior(design)
  uni    <- c(1, 1, 1) / 3
  fmfounder <- .n_bc(design) + 1L                    # founder meiosis count
  emimat <- .gt_emimat(emission_gt(germ = err, gert = gert))   # 3 states x 4 obs

  fams <- unique(ped$family)
  out_parts <- list()
  for (fam in fams) {
    pf   <- ped[ped$family == fam, , drop = FALSE]
    taxa <- pf$taxon
    par  <- pf$parent1[match(taxa, pf$taxon)]        # parent taxon (selfing: p1 == p2)
    # forest wiring: 0-based parent index within this family; root has parent -1.
    pidx <- match(par, taxa) - 1L
    pidx[is.na(pidx)] <- -1L
    root <- which(pidx < 0L)
    if (length(root) != 1L)
      stop("refine_ancestry(): family '", fam, "' must have exactly one root (got ",
           length(root), ")")
    # meioses = founder count + depth from root along parent links
    depth <- integer(length(taxa))
    repeat {
      upd <- FALSE
      for (v in seq_along(taxa)) if (pidx[v] >= 0L) {
        nd <- depth[pidx[v] + 1L] + 1L
        if (nd != depth[v]) { depth[v] <- nd; upd <- TRUE }
      }
      if (!upd) break
    }
    meioses <- fmfounder + depth
    hasData <- taxa %in% mo$name

    mof <- mo[mo$name %in% taxa, , drop = FALSE]
    for (ch in sort(unique(mof$chr))) {
      moc <- mof[mof$chr == ch, , drop = FALSE]
      pos <- sort(unique(moc$pos)); M <- length(pos)
      if (M < 1L) next
      cmv <- if ("cm" %in% names(moc)) moc$cm[match(pos, moc$pos)] else NULL
      rvec <- if (M > 1L) .interval_r(pos, cmv, rrate) else numeric(0)

      # per-node priors and emission
      rho <- pim <- matrix(rep(uni, each = length(taxa)), ncol = 3)
      rho[root, ] <- pim[root, ] <- pi0
      emit <- vector("list", length(taxa))
      for (v in seq_along(taxa)) {
        if (!hasData[v]) { emit[[v]] <- matrix(1, M, 3); next }
        sv <- moc[moc$name == taxa[v], , drop = FALSE]
        g  <- sv$state[match(pos, sv$pos)]           # 0/1/2, NA -> missing
        g[is.na(g)] <- 3L
        emit[[v]] <- t(emimat[, g + 1L, drop = FALSE])   # M x 3 = P(obs | state)
      }

      bel <- pedigree_bp_cpp(M, as.integer(pidx), as.integer(meioses),
                             as.logical(hasData), emit, rho, pim, rvec,
                             root = as.integer(root - 1L),
                             maxIters = as.integer(maxiter), tol = tol, lambda = lambda)
      # decode genotyped leaves -> refined state at each pos
      for (v in seq_along(taxa)) {
        if (!hasData[v]) next
        state <- max.col(bel[[v]], ties.method = "first") - 1L    # argmax dosage
        sv <- moc[moc$name == taxa[v], , drop = FALSE]
        sv$state <- state[match(sv$pos, pos)]
        out_parts[[length(out_parts) + 1L]] <- sv
      }
    }
  }
  refined <- do.call(rbind, out_parts)
  # restore original row order of the genotyped rows
  key_in  <- paste(mo$name, mo$chr, mo$pos)
  key_out <- paste(refined$name, refined$chr, refined$pos)
  refined[match(key_in[key_in %in% key_out], key_out), , drop = FALSE]
}
