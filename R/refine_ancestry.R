# Pedigree-aware ancestry calling / refinement (design of record: design/PEDIGREE_HMM.md).
#
# One BP engine, two front doors sharing the core `.pedigree_states()`:
#   * call_ancestry(caller = "pedigree", pedigree = ...)  -- a count-based caller
#     in the call_ancestry() family (de novo from read counts; the primary use).
#   * refine_ancestry(mosaic, pedigree, ...)              -- a thin wrapper for the
#     hard-call refinement use (couple relatives to correct an existing mosaic).
#
# `.pedigree_states()` couples relatives through the pedigree via loopy belief
# propagation (pedigree_bp_cpp) over the (pedigree x genome) grid. The emission is
# the swappable observation channel (Text S2 of the methods supplement):
#   "count" -> emission_count() BetaBinomial over read depths; missing = zero depth
#              -> flat, so the pedigree fills it in weighted by relative confidence.
#   "gt"    -> emission_gt() categorical genotype-confusion over hard `state` calls;
#              missing = the explicit g=3 column -> flat.
# The generation prior enters ONLY at each family founder (rho/pi = pi_0) and
# propagates downward through the selfing transmission kernel Tsib; non-root chains
# are marginal-neutral (uniform stationary). R owns all IO; the C++ kernel is
# data-agnostic. Families are processed independently; latent ungenotyped ancestors
# (taxa named as parents but absent from the data) impose chromosome continuity.

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

# per-interval recombination fraction: Haldane on genetic distance where both
# flanking cm are finite, else the physical rrate*bp fallback for that gap.
.interval_r <- function(pos, cm = NULL, rrate = 0.01) {
  phys <- pmin(0.5, rrate * pmax(0, diff(as.numeric(pos))))
  if (is.null(cm) || all(is.na(cm))) return(phys)
  d   <- pmax(0, diff(cm)) / 100                     # Morgans
  gen <- pmin(0.5, 0.5 * (1 - exp(-2 * d)))          # Haldane
  ok  <- is.finite(gen)
  phys[ok] <- gen[ok]
  phys
}

# --- shared BP core (data-agnostic to caller vs refine) ----------------------
# `data`: per-marker table keyed by name/chr/pos, plus `state` (gt) or n_ref/n_alt
#   (count), optionally `cm` and `donor`. `ped`: a resolved read_pedigree() frame.
# Returns one row per genotyped-leaf input marker: name, chr, pos, donor, state
# (marginal-MAP decode of the posterior belief). No IO.
.pedigree_states <- function(data, ped, design = "BC2S3",
                             emission = c("count", "gt"),
                             err = NULL, gert = 0.10, conc = 20, rrate = 0.01,
                             maxiter = 30L, tol = 1e-4, lambda = 0.5) {
  emission <- match.arg(emission)
  if (is.null(err)) err <- if (emission == "gt") 0.05 else 0.01
  mo <- as.data.frame(data, stringsAsFactors = FALSE)
  need <- if (emission == "gt") c("name", "chr", "pos", "state")
          else                  c("name", "chr", "pos", "n_ref", "n_alt")
  if (!all(need %in% names(mo)))
    stop(".pedigree_states(): emission='", emission, "' needs columns ",
         paste(need, collapse = ", "))
  if (!all(c("taxon", "family", "parent1") %in% names(ped)))
    stop(".pedigree_states(): `pedigree` needs columns taxon, family, parent1")
  has_donor <- "donor" %in% names(mo)

  pi0    <- .founder_prior(design)
  uni    <- c(1, 1, 1) / 3
  fmfounder <- .n_bc(design) + 1L                    # founder meiosis count
  emimat <- if (emission == "gt") .gt_emimat(emission_gt(germ = err, gert = gert)) else NULL
  ctheta <- if (emission == "count") .emission_theta(emission_count(err = err, conc = conc)) else NULL

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
      stop(".pedigree_states(): family '", fam, "' must have exactly one root (got ",
           length(root), ")")
    # meioses = founder count + depth from root along parent links
    depth <- integer(length(taxa))
    for (iter in seq_len(length(taxa) + 1L)) {       # bounded: a valid tree converges in <= depth passes
      upd <- FALSE
      for (v in seq_along(taxa)) if (pidx[v] >= 0L) {
        nd <- depth[pidx[v] + 1L] + 1L
        if (nd != depth[v]) { depth[v] <- nd; upd <- TRUE }
      }
      if (!upd) break
      if (iter == length(taxa) + 1L)
        stop(".pedigree_states(): family '", fam, "' pedigree contains a cycle")
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
        if (emission == "gt") {
          g  <- sv$state[match(pos, sv$pos)]           # 0/1/2, NA -> missing
          g[is.na(g)] <- 3L
          emit[[v]] <- t(emimat[, g + 1L, drop = FALSE])   # M x 3 = P(obs | state)
        } else {                                        # count: BetaBinomial over depths
          nref <- sv$n_ref[match(pos, sv$pos)]; nref[is.na(nref)] <- 0L
          nalt <- sv$n_alt[match(pos, sv$pos)]; nalt[is.na(nalt)] <- 0L
          n <- as.integer(nref + nalt); a <- as.integer(nalt)   # depth, alt reads
          emit[[v]] <- exp(count_emission_loglik_cpp(n, a, ctheta, conc))  # M x 3; n=0 -> flat
        }
      }

      bel <- pedigree_bp_cpp(M, as.integer(pidx), as.integer(meioses),
                             as.logical(hasData), emit, rho, pim, rvec,
                             root = as.integer(root - 1L),
                             maxIters = as.integer(maxiter), tol = tol, lambda = lambda)
      # decode genotyped leaves -> refined state at each pos (marginal-MAP)
      for (v in seq_along(taxa)) {
        if (!hasData[v]) next
        state <- max.col(bel[[v]], ties.method = "first") - 1L    # argmax dosage
        sv <- moc[moc$name == taxa[v], , drop = FALSE]
        out_parts[[length(out_parts) + 1L]] <- data.frame(
          name  = sv$name, chr = as.integer(sv$chr), pos = as.integer(sv$pos),
          donor = if (has_donor) sv$donor else NA_character_,
          state = as.integer(state[match(sv$pos, pos)]),
          stringsAsFactors = FALSE)
      }
    }
  }
  refined <- do.call(rbind, out_parts)
  if (is.null(refined))
    return(data.frame(name = character(), chr = integer(), pos = integer(),
                      donor = character(), state = integer(), stringsAsFactors = FALSE))
  refined
}

#' Refine per-individual ancestry calls over a pedigree
#'
#' Thin wrapper on the shared pedigree belief-propagation kernel
#' ([call_states()] `caller = "pedigree"`) for the **hard-call refinement** use:
#' couple a family's relatives to correct a per-marker mosaic produced by another
#' caller. Same engine as `caller = "pedigree"`; the difference is ergonomic ---
#' `refine_ancestry()` takes and returns a mosaic in place (so a
#' `call_states()` output goes in and a same-shape corrected mosaic comes back),
#' whereas the caller takes the raw observation table. See design/PEDIGREE_HMM.md.
#'
#' Two emission modes (the swappable observation channel): `"gt"` (default;
#' depth-blind, [emission_gt()] over the hard `state` calls) or `"count"`
#' (depth-aware, [emission_count()] BetaBinomial over read depths). The count mode
#' is identical to `call_ancestry(caller = "pedigree")` on the same counts; prefer
#' the caller for de novo calling and reserve `refine_ancestry()` for correcting an
#' existing hard-call mosaic.
#'
#' @param mosaic Per-marker table keyed by `name, chr, pos` (+ optional `source,
#'   donor, cm`). For `emission = "gt"` needs `state` in `{0,1,2,3=missing}` (from
#'   [call_states()]); for `emission = "count"` needs `n_ref, n_alt` read depths.
#' @param pedigree A pedigree: a path (read via [read_pedigree()]) or a
#'   [read_pedigree()]-shaped data.frame (`taxon, family, parent1, parent2`).
#'   `taxon` joins to `mosaic$name`; a `taxon` used as a parent but absent from
#'   `mosaic` is a latent ancestor.
#' @param design Breeding design `"BC{n}S{m}"` -> founder prior `pi_0` and
#'   per-node `meioses` (via [design_priors()]).
#' @param emission `"gt"` (depth-blind, over hard states) or `"count"`
#'   (depth-aware BetaBinomial over `n_ref`/`n_alt`).
#' @param err Genotyping/read error: [emission_gt()] `germ` when `emission="gt"`
#'   (default 0.05), per-read error when `emission="count"` (default 0.01).
#' @param gert [emission_gt()] genotyping-error-rate-in-transmission (`"gt"` only).
#' @param conc [emission_count()] BetaBinomial concentration (`"count"` only).
#' @param rrate Per-bp recombination fraction applied to `pos` gaps when `mosaic`
#'   has no `cm` column (else Haldane on `cm`).
#' @param maxiter,tol,lambda BP sweeps, convergence tolerance, damping.
#' @param ped_format Passed to [read_pedigree()] when `pedigree` is a path.
#' @return `mosaic` with `state` replaced by the refined per-marker calls (same
#'   columns and row order; genotyped leaves only). Feed to [to_segments()].
#' @seealso [call_states()] (`caller = "pedigree"`), [simulate_family()],
#'   [pedigree_bp_cpp()], [to_segments()].
#' @examples
#' \dontrun{
#' sim   <- simulate_family("BC2S3", families = 10, sibs = 10, seed = 1)
#' obs   <- simulate_counts(sim$truth, depth = 0.5, seed = 1)
#' calls <- call_states(obs, caller = "rtiger", design = "BC2S3")
#' ref   <- refine_ancestry(calls, sim$pedigree, design = "BC2S3")
#' }
#' @export
refine_ancestry <- function(mosaic, pedigree, design = "BC2S3",
                            emission = c("gt", "count"),
                            err = NULL, gert = 0.10, conc = 20, rrate = 0.01,
                            maxiter = 30L, tol = 1e-4, lambda = 0.5,
                            ped_format = c("fam", "fsfhap")) {
  emission <- match.arg(emission)
  ped <- if (is.character(pedigree)) read_pedigree(pedigree, match.arg(ped_format))
         else as.data.frame(pedigree, stringsAsFactors = FALSE)
  mo  <- as.data.frame(mosaic, stringsAsFactors = FALSE)
  st  <- .pedigree_states(mo, ped, design = design, emission = emission, err = err,
                          gert = gert, conc = conc, rrate = rrate,
                          maxiter = maxiter, tol = tol, lambda = lambda)
  if (!nrow(st)) return(mo[0, , drop = FALSE])
  # restore the input mosaic's rows/columns/order; replace state on genotyped rows
  key_mo <- paste(mo$name, mo$chr, mo$pos)
  key_st <- paste(st$name, st$chr, st$pos)
  keep   <- key_mo %in% key_st
  out    <- mo[keep, , drop = FALSE]
  out$state <- st$state[match(key_mo[keep], key_st)]
  out
}
