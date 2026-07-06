# Flanking-marker genotype interpolation (Tian 2011 / Chen 2019). Densifies a
# COMPLETE dense genotype block onto a target marker grid by linear interpolation
# in genetic distance (cM). Data-agnostic: the caller supplies cM (from a map on
# its side) and, for block-sparse per-family data, owns the family loop. This is
# deterministic hard-call interpolation, NOT ancestry projection.

#' Interpolate genotypes onto a target marker grid (Tian 2011 / Chen 2019)
#'
#' Densifies a complete genotype block onto a common target marker grid by
#' flanking-marker interpolation in genetic distance (cM). For a target between
#' called flanking markers L and R, the interpolation weight is
#' `w = (cM_x - cM_L) / (cM_R - cM_L)` and the continuous dosage is
#' `v_L + w * (v_R - v_L)` -- the Tian 2011 rule (the weighted average
#' interpreted as the probability the SNP comes from the non-recurrent parent).
#' Ends are clamped to the terminal observed value (equivalent to
#' [stats::approx()] with `rule = 2`).
#'
#' @details
#' Genotypes are the alt/teosinte-allele dosage: `0` = REF/recurrent hom,
#' `1` = HET, `2` = ALT/donor hom (`continuous` returns real values in `[0, 2]`).
#' Three modes share one interpolation core:
#' \describe{
#'   \item{`continuous`}{dosage ramp (Tian 2011); returns fractional dosage.}
#'   \item{`step`}{nearest flanking value (`w < 0.5 ? v_L : v_R`, hard step at the
#'     cM midpoint, tie `w == 0.5` -> `v_R`); the faithful Chen/TeoNAM discrete
#'     form -- fabricates no het across a `0<->2` gap.}
#'   \item{`round`}{`round(continuous)` to 0/1/2; fabricates a het band across a
#'     `0<->2` gap (width proportional to the cM gap) -- provided to demonstrate
#'     that artifact.}
#'   \item{`chen2019`}{the composite genetic-MAP construction rule of Chen et al.
#'     2019 (TeoNAM): *"the missing genotypes were imputed according to the
#'     flanking markers. If the flanking markers had same genotypes, the missing
#'     genotype was imputed as the same with flanking markers, or otherwise left
#'     as missing."* So concordant flanks fill; **discordant flanks -> NA**; and
#'     targets past the terminal observed marker (chromosome ends) -> **NA** (no
#'     rule=2 clamping). Order-based: the coordinate only locates the flanks, no
#'     distance weight enters the decision. Use this for MAP construction; use
#'     `continuous`/`step`/`round` for the Tian 2011 JLM/GWAS densification, which
#'     imputes by the genetic distance of the flanking markers. Unlike the other
#'     three modes, `chen2019` **can emit `NA`**.}
#' }
#'
#' The function densifies **one dense block** to a target grid. For block-sparse
#' per-family data (each family genotyped on its own marker subset), loop over
#' families -- interpolating each family's dense block onto the shared grid -- and
#' `cbind` the results; that loop lives in the consumer pipeline, not here.
#'
#' Tied target positions: `obs$cm` must be strictly increasing (source flanks
#' must be distinct), but **`target$cm` may repeat** -- the target grid does NOT
#' need unique cM. Every target row gets its own output row, so pass the full
#' marker set to densify all of it; do not pre-collapse to unique cM (that
#' silently drops markers). Two target markers at the same cM are, by the map, at
#' one genetic location (zero modelled recombination between them, i.e. perfect
#' LD), so they receive an **identical genotype vector** -- deterministic, and
#' correct for the interpolation model, not a bug. This is common in
#' recombination-cold centromeric regions (a whole Mb mapping to one cM). If a
#' downstream step needs non-duplicate columns, thin the twins there
#' ([select_independent()] on an `r2`/[position_distance()] matrix); interpolation
#' cannot manufacture sub-cM resolution the map lacks.
#'
#' `geno` columns are assumed COMPLETE (no `NA`) -- true for imputed truth and for
#' HMM-caller output. There is no `NA`-aware flanking search.
#'
#' @param geno Numeric matrix, `nrow(geno) == nrow(obs)`, column = sample;
#'   COMPLETE (no `NA`). Values are the alt-allele dosage in `[0, 2]`.
#' @param obs `data.frame` with columns `chr` and the coordinate column named by
#'   `coord` (default `cm`), aligned row-for-row to `geno`, sorted by
#'   `(chr, coord)` with `coord` strictly increasing within each chromosome.
#' @param target `data.frame` with columns `chr` and the coordinate column named
#'   by `coord`, sorted by `(chr, coord)`. The coordinate may repeat (tied
#'   positions allowed); each row yields one output row, and markers sharing a
#'   coordinate get identical genotypes (see Details).
#' @param mode One of `"continuous"` (Tian 2011 dosage ramp), `"step"`, `"round"`
#'   (all three = distance-based JLM/GWAS densification), or `"chen2019"` (the
#'   composite genetic-map rule: concordant flanks fill, discordant flanks or
#'   chromosome ends -> `NA`; see Details).
#' @param coord Name of the coordinate column to interpolate along; default
#'   `"cm"` (genetic distance). Use e.g. `"bp"` for physical-distance
#'   interpolation -- appropriate when interpolating in cM would be circular, such
#'   as building a native genetic map from bp positions. The coordinate is just
#'   the axis to interpolate along; the arithmetic is unit-agnostic.
#' @return Numeric matrix `nrow(target)` x `ncol(geno)`; rownames from `target`
#'   (if any), colnames from `geno`. Modes `continuous`/`step`/`round` never
#'   introduce `NA`; mode `chen2019` returns `NA` at discordant-flank targets and
#'   chromosome ends.
#' @examples
#' obs <- data.frame(chr = 1L, cm = c(0, 1))
#' geno <- matrix(c(0, 2), nrow = 2, dimnames = list(NULL, "S1"))
#' target <- data.frame(chr = 1L, cm = c(0, 0.5, 1))
#' interpolate_genotype(geno, obs, target, "continuous")  # 0, 1, 2
#'
#' # Interpolate along physical distance (bp) instead of cM.
#' obs_bp    <- data.frame(chr = 1L, bp = c(1e6, 3e6))
#' target_bp <- data.frame(chr = 1L, bp = c(1e6, 2e6, 3e6))
#' interpolate_genotype(geno, obs_bp, target_bp, "continuous", coord = "bp")  # 0, 1, 2
#'
#' # Chen 2019 composite-map rule: concordant flanks fill, discordant/ends -> NA.
#' obs2    <- data.frame(chr = 1L, cm = c(0, 2))
#' geno2   <- matrix(c(0, 0,  0, 2), nrow = 2,
#'                   dimnames = list(NULL, c("concord", "discord")))
#' target2 <- data.frame(chr = 1L, cm = c(0, 1, 2, 3))
#' interpolate_genotype(geno2, obs2, target2, "chen2019")
#' #            concord discord
#' # cm0 (obs)        0       0
#' # cm1 (mid)        0      NA   # discordant flanks -> NA
#' # cm2 (obs)        0       2
#' # cm3 (end)       NA      NA   # chromosome end -> NA
#' @export
interpolate_genotype <- function(geno, obs, target,
                                 mode = c("continuous", "step", "round", "chen2019"),
                                 coord = "cm") {
  mode <- match.arg(mode)
  mode_int <- match(mode, c("continuous", "step", "round", "chen2019")) - 1L  # 0/1/2/3
  if (!is.character(coord) || length(coord) != 1L)
    stop("interpolate_genotype(): `coord` must be a single column name.")

  geno <- as.matrix(geno)
  if (!is.numeric(geno)) storage.mode(geno) <- "double"
  if (anyNA(geno))
    stop("interpolate_genotype(): `geno` must be complete (no NA).")
  if (!all(c("chr", coord) %in% names(obs)))
    stop("interpolate_genotype(): `obs` must have columns `chr` and `", coord, "`.")
  if (!all(c("chr", coord) %in% names(target)))
    stop("interpolate_genotype(): `target` must have columns `chr` and `", coord, "`.")
  if (nrow(geno) != nrow(obs))
    stop("interpolate_genotype(): nrow(geno) (", nrow(geno),
         ") must equal nrow(obs) (", nrow(obs), ").")

  # obs must be sorted by (chr, coord) and strictly increasing in coord within a chr.
  ord <- order(obs$chr, obs[[coord]])
  if (!identical(ord, seq_len(nrow(obs))))
    stop("interpolate_genotype(): `obs` must be sorted by (chr, ", coord, ").")
  for (ch in unique(obs$chr)) {
    c_ch <- obs[[coord]][obs$chr == ch]
    if (any(diff(c_ch) <= 0))
      stop("interpolate_genotype(): `obs$", coord, "` must be strictly increasing ",
           "within chromosome ", ch, " (collapse tied positions upstream).")
  }
  if (!identical(order(target$chr, target[[coord]]), seq_len(nrow(target))))
    stop("interpolate_genotype(): `target` must be sorted by (chr, ", coord, ").")

  # Empty target grid -> a 0-row matrix with geno's columns (rbind of no blocks
  # would be NULL and break the colnames<- assignment below).
  if (nrow(target) == 0L)
    return(matrix(numeric(0), nrow = 0L, ncol = ncol(geno),
                  dimnames = list(NULL, colnames(geno))))

  # Split by chr, interpolate each chr independently (no cross-chr bleed), rbind.
  chrs <- unique(target$chr)
  blocks <- vector("list", length(chrs))
  for (i in seq_along(chrs)) {
    ch <- chrs[i]
    orows <- which(obs$chr == ch)
    trows <- which(target$chr == ch)
    if (!length(orows))
      stop("interpolate_genotype(): target chromosome ", ch,
           " has no observed markers in `obs`.")
    blocks[[i]] <- interp_geno_cpp(
      obs_cm    = as.numeric(obs[[coord]][orows]),
      G         = geno[orows, , drop = FALSE],
      target_cm = as.numeric(target[[coord]][trows]),
      mode      = mode_int
    )
  }
  out <- do.call(rbind, blocks)

  colnames(out) <- colnames(geno)
  # Carry explicit target rownames only; a data.frame's automatic row names
  # (integer 1:n) are not meaningful marker ids and are dropped.
  rn <- attr(target, "row.names")
  if (is.character(rn)) rownames(out) <- rn
  out
}
