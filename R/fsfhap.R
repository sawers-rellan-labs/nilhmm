# FSFHap caller (design/FSFHAP_PORT.md) — native port of TASSEL's FSFHap, WIRED
# into call_ancestry(caller = "fsfhap"). BOTH parent-calling routes are ported and
# TASSEL-validated end-to-end (per-cell imputed_parents concordance vs TASSEL):
#   - stage 1a   .fsfhap_segregating        (whichSitesSegregateCorrectly)
#   - stage 1b   .fsfhap_call_parents_bc    BC route: same-tag filter, A/C recode,
#                coverage. BC1 (contribution 0.75). e2e parity 1.00000.
#   - stage 2b   .fsfhap_biparental_call    NON-BC1 route (BiparentalHaplotypeFinder):
#                preFilterSites -> reconstruct 2 parental haplotypes -> recode.
#                e2e parity 0.99997 (F2). Uses fsfhap_prefilter_sites_cpp,
#                fsfhap_biparental_alleles_cpp, fsfhap_cluster_window_cpp.
#   - stage 3    .fsfhap_impute             (imputeUsingViterbiFiveState + fillGaps)
#   - .fsfhap_phet    design-derived phet = (1 - F)/2
#   - .fsfhap_states  per family x chromosome driver; dispatches route by design
#                (BC1 -> "bc", else -> "biparental") -> common REF/HET/ALT schema.
# NOT ported (separate -ImputeProgenyStatesPlugin tool, off every FSFHap path): the
# 4-state ImputeCrossProgeny / CrossProgenyEmissionMatrix lineage.
# See design/FSFHAP_PORT.md for the route dispatch and remaining checklist items.

#' FSFHap stage 1a: segregating-site test for one family x chromosome
#'
#' Thin wrapper over [fsfhap_segregating_sites_cpp()] — a faithful port of
#' TASSEL's `whichSitesSegregateCorrectly`. Given one family's genotype matrix on
#' one chromosome, flags the sites whose minor-allele count fits the family's
#' segregation model (backcross vs F2) better than a monomorphic/error model.
#'
#' `ratio` selects the branch, matching TASSEL: `0.25`/`0.75` → backcross (keep
#' iff `pquarter > phalf && pquarter > pmono`), anything else → F2 (`phalf /
#' (pmono + pquarter) > 2`). Note TASSEL's backcross branch always tests the
#' **0.25** model regardless of backcross depth, so deeper backcrosses (e.g. a
#' BC2 donor fraction ~0.125) are tested against 0.25 too; the port reproduces
#' this deliberately for parity.
#'
#' @param G Integer matrix, taxa x sites, values 0/1/2/3 = REF-hom / het /
#'   ALT-hom / missing (engine canonical `g`); sites sorted by position.
#' @param max_missing Max missing-genotype proportion for a site to be tested
#'   (TASSEL default 0.9; the TeoNAM run passes `1.0`).
#' @param ratio Expected minor-allele frequency selecting the segregation model:
#'   `0.25`/`0.75` (backcross) or `0.5` (F2). **Design-derived, no default** — it
#'   must be supplied by the caller from the breeding design (the routing
#'   dispatcher, mirroring TASSEL's `contribution1` → route → ratio), never
#'   guessed or defaulted here.
#' @return A list: `seg` (logical, kept sites), `Mj`/`Mn` (major/minor allele
#'   counts), `p_missing`, and model probabilities `pmono`/`pquarter`/`phalf`.
#' @keywords internal
.fsfhap_segregating <- function(G, max_missing = 0.9, ratio) {
  if (!is.matrix(G)) stop(".fsfhap_segregating(): G must be a taxa x sites matrix")
  if (!is.integer(G)) storage.mode(G) <- "integer"
  fsfhap_segregating_sites_cpp(G, max_missing, ratio)
}

#' FSFHap stage 1b: backcross parent-allele calling for one family x chromosome
#'
#' Faithful orchestration of TASSEL's `callParentAllelesByWindowForBackcrosses`:
#' keep sites that both segregate as backcross (stage 1a) and survive the same-tag
#' filter ([fsfhap_same_tag_keep_cpp()]); assign A = major / C = minor per kept
#' site; drop low-coverage taxa (`<= min_gametes` non-missing gametes across kept
#' sites); recode genotypes into the parent-origin A/C frame
#' (`0` = A-hom / `1` = het / `2` = C-hom / `3` = missing), where A is the site's
#' major allele. The result is the observation matrix the 5-state smoother
#' consumes next.
#'
#' @param G Integer matrix, taxa x sites, canonical `g` in `{0,1,2,3}`; one family,
#'   one chromosome, sites sorted by position.
#' @param pos Integer marker positions (bp), length = ncol(G).
#' @param max_missing Max missing proportion for the segregating-site test
#'   (TASSEL default 0.9; TeoNAM run 1.0).
#' @param min_rsq Same-tag R^2 threshold (TASSEL 0.8).
#' @param min_gametes Coverage floor: taxa with `> min_gametes` non-missing
#'   gametes (2 x non-missing sites) are kept (TASSEL 200).
#' @param min_r LD-filter threshold; `> 0` enables TASSEL's `ldfilter`, **not yet
#'   ported** (warns and ignores). The TeoNAM run uses `0`.
#' @return List: `G` (recoded, kept taxa x kept sites), `keep_sites`/`keep_taxa`
#'   (indices into the input), `pos` (kept-site positions), `major_is_ref`
#'   (per kept site), and counts `n_seg`/`n_sametag`/`n_kept_sites`.
#' @keywords internal
.fsfhap_call_parents_bc <- function(G, pos, max_missing = 1.0, min_rsq = 0.8,
                                    min_gametes = 200L, min_r = 0.0) {
  if (!is.matrix(G)) stop(".fsfhap_call_parents_bc(): G must be a taxa x sites matrix")
  if (!is.integer(G)) storage.mode(G) <- "integer"
  if (length(pos) != ncol(G)) stop(".fsfhap_call_parents_bc(): length(pos) != ncol(G)")
  if (is.unsorted(pos)) stop(".fsfhap_call_parents_bc(): pos must be sorted ascending ",
                             "(one chromosome; the same-tag/window kernels assume it).")
  if (min_r > 0) warning(".fsfhap_call_parents_bc(): LD filter (min_r > 0) not yet ",
                         "ported; ignoring (matches the TeoNAM min_r = 0 run).")

  # per-site allele counts -> major allele (ties -> REF, deterministic)
  nRef <- colSums(G == 0L); nHet <- colSums(G == 1L); nAlt <- colSums(G == 2L)
  major_is_ref <- (2 * nRef + nHet) >= (2 * nAlt + nHet)

  # 0.25 is the BC route's design-FIXED segregation ratio: this function *is* the
  # backcross route (TASSEL enters it only for contribution1 in {0.75,0.25}), so
  # the ratio is fixed by the design that routed here, not guessed.
  seg     <- fsfhap_segregating_sites_cpp(G, max_missing, 0.25)$seg   # stage 1a
  sametag <- fsfhap_same_tag_keep_cpp(G, as.integer(pos), major_is_ref, min_rsq)
  keep_sites <- which(seg & sametag)                                 # filteredBits & polybits

  Gk  <- G[, keep_sites, drop = FALSE]
  mir <- major_is_ref[keep_sites]
  # recode to A/C frame: where major is ALT, flip 0<->2 (het/missing unchanged)
  flip <- which(!mir)
  if (length(flip)) {
    sub <- Gk[, flip, drop = FALSE]
    sub[Gk[, flip, drop = FALSE] == 0L] <- 2L
    sub[Gk[, flip, drop = FALSE] == 2L] <- 0L
    Gk[, flip] <- sub
  }
  # coverage filter (on kept sites): keep taxa with > min_gametes non-missing gametes
  keep_taxa <- which(2L * rowSums(Gk != 3L) > min_gametes)

  list(G = Gk[keep_taxa, , drop = FALSE], keep_sites = keep_sites,
       keep_taxa = keep_taxa, pos = pos[keep_sites], major_is_ref = mir,
       n_seg = sum(seg), n_sametag = sum(sametag), n_kept_sites = length(keep_sites))
}

#' FSFHap stage 3: 5-state EM imputation for one family x chromosome
#'
#' Thin wrapper over [fsfhap_impute_five_state_cpp()] — the real FSFHap imputation
#' (via TASSEL's `ViterbiAlgorithmPlugin` → `imputeUsingViterbiFiveState`), then
#' [fsfhap_fill_gaps_cpp()] (`fillGapsInAlignment`) if `fill_gaps`. Runs on the
#' parent-called A/het/C frame from [.fsfhap_call_parents_bc()].
#'
#' @param G Integer matrix, taxa x sites, parent-called `g` in `{0,1,2,3}`.
#' @param pos Integer marker positions (bp), length = ncol(G).
#' @param phet Expected heterozygosity. **Design-derived, no default** — the
#'   caller supplies `(1 - F)/2` from the pedigree inbreeding coefficient
#'   (BC1S4 F=0.9375 → 0.03125), mirroring TASSEL's `ViterbiAlgorithmPlugin`.
#' @param max_iter EM iteration cap (TASSEL 50).
#' @param fill_gaps Forward-fill missing runs between equal flanking calls.
#' @return The imputed matrix (taxa x sites, `0`/`1`/`2`/`3`), with attributes
#'   `iters` and `emission` (final 5x3) from the EM.
#' @keywords internal
.fsfhap_impute <- function(G, pos, phet, max_iter = 50L, fill_gaps = TRUE) {
  if (!is.matrix(G)) stop(".fsfhap_impute(): G must be a taxa x sites matrix")
  if (!is.integer(G)) storage.mode(G) <- "integer"
  if (length(pos) != ncol(G)) stop(".fsfhap_impute(): length(pos) != ncol(G)")
  if (is.unsorted(pos)) stop(".fsfhap_impute(): pos must be sorted ascending.")
  res <- fsfhap_impute_five_state_cpp(G, as.integer(pos), phet, as.integer(max_iter))
  out <- if (isTRUE(fill_gaps)) fsfhap_fill_gaps_cpp(res$imputed) else res$imputed
  attr(out, "iters") <- res$iters
  attr(out, "emission") <- res$emission
  out
}

#' Design-derived expected heterozygosity from the pedigree inbreeding coefficient
#'
#' `phet = (1 - F)/2` when `F` in `[0, 1]` (TASSEL's `ViterbiAlgorithmPlugin`),
#' else the `default` fallback (TASSEL `probHeterozygous` = 0.07). Keeps `phet`
#' design-derived, not a magic constant (no hardcoded design parameters).
#'
#' @param F Inbreeding coefficient (pedigree column 7 / `1 - 0.5^nself`).
#' @param default Fallback when `F` is outside `[0, 1]`.
#' @return The expected-heterozygosity scalar for [.fsfhap_impute()].
#' @keywords internal
.fsfhap_phet <- function(F, default = 0.07) {
  if (length(F) == 1 && !is.na(F) && F >= 0 && F <= 1) (1 - F) / 2 else default
}

# seed window used by fsfhap_biparental_alleles_cpp (C++ BHF_WINDOW); the finder
# cannot seed on fewer preFiltered sites than one window.
.FSFHAP_BHF_WINDOW <- 100L

#' FSFHap stage 2b: BiparentalHaplotypeFinder parent-calling (non-BC1 route)
#'
#' The non-backcross route: [fsfhap_prefilter_sites_cpp()] → reconstruct two
#' parental haplotypes ([fsfhap_biparental_alleles_cpp()]) → recode each genotype
#' into the parent-origin A/het/C frame using the per-site `alleleA`/`alleleC` →
#' coverage filter. Returns the same shape as [.fsfhap_call_parents_bc()] so the
#' shared stage-3 imputation can consume it.
#'
#' @param G Integer matrix, taxa x sites, canonical `g` in `{0,1,2,3}`; one family,
#'   one chromosome, sites sorted.
#' @param pos Integer marker positions (bp), length = ncol(G).
#' @param min_maf,min_coverage,max_het_deviation,min_r2 preFilterSites knobs
#'   (TASSEL 0.05 / 0.2 / 5 / 0.2).
#' @param min_gametes Coverage floor (TASSEL 200).
#' @return List `G` (recoded, kept taxa x kept sites), `keep_sites`, `keep_taxa`,
#'   `pos`; empty if preFilterSites leaves too few sites or no seed is found.
#' @keywords internal
.fsfhap_biparental_call <- function(G, pos, min_maf = 0.05, min_coverage = 0.2,
                                    max_het_deviation = 5, min_r2 = 0.2, min_gametes = 200L) {
  if (!is.matrix(G)) stop(".fsfhap_biparental_call(): G must be a taxa x sites matrix")
  if (!is.integer(G)) storage.mode(G) <- "integer"
  if (length(pos) != ncol(G)) stop(".fsfhap_biparental_call(): length(pos) != ncol(G)")
  if (is.unsorted(pos)) stop(".fsfhap_biparental_call(): pos must be sorted ascending.")
  empty <- list(G = G[, integer(0), drop = FALSE], keep_sites = integer(0),
                keep_taxa = integer(0), pos = integer(0))
  keepPre <- fsfhap_prefilter_sites_cpp(G, as.integer(pos), min_maf, min_coverage,
                                        max_het_deviation, min_r2)
  fsites <- which(keepPre)
  if (length(fsites) < .FSFHAP_BHF_WINDOW) return(empty)       # fewer than one seed window
  Gf <- G[, fsites, drop = FALSE]
  al <- fsfhap_biparental_alleles_cpp(Gf)
  if (!isTRUE(al$seeded)) return(empty)
  aA <- al$alleleA; aC <- al$alleleC
  # usable sites: both parent alleles known and distinct (0 vs 2)
  keep2 <- which(aA != 3L & aC != 3L & aA != aC)
  if (!length(keep2)) return(empty)
  Gk <- Gf[, keep2, drop = FALSE]; aAk <- aA[keep2]; aCk <- aC[keep2]
  rec <- matrix(3L, nrow(Gk), ncol(Gk))                        # recode to A(0)/het(1)/C(2)/N(3)
  for (j in seq_along(keep2)) {
    col <- Gk[, j]
    rec[col == aAk[j], j] <- 0L
    rec[col == aCk[j], j] <- 2L
    rec[col == 1L,     j] <- 1L
  }
  storage.mode(rec) <- "integer"
  keep_taxa <- which(2L * rowSums(rec != 3L) > min_gametes)
  list(G = rec[keep_taxa, , drop = FALSE], keep_sites = fsites[keep2],
       keep_taxa = keep_taxa, pos = pos[fsites[keep2]])
}

#' FSFHap caller driver: per-family x chromosome parent-calling + imputation
#'
#' The `caller = "fsfhap"` path for [call_states()]. Unlike the per-sample callers,
#' FSFHap **pools each family**: for every `family x chr` it builds the taxa x sites
#' genotype matrix, runs BC parent-calling ([.fsfhap_call_parents_bc()], stage 1)
#' then the 5-state EM imputation + gap-fill ([.fsfhap_impute()], stage 3), and
#' emits per-marker states in the common REF/HET/ALT frame (major allele = the
#' recurrent parent in a backcross → `0` REF, `1` HET, `2` ALT/donor).
#'
#' @param data Long table with `name, chr, pos, g` (`g` in `{0,1,2,3}`) and a
#'   `family` grouping column; optional `donor`.
#' @param phet Design-derived expected heterozygosity ([.fsfhap_phet()]).
#' @param source,donor Output `source`/`donor` labels (donor per-`name` if
#'   `has_donor`).
#' @param has_donor Whether `data` carries a per-row `donor` column.
#' @param route `"bc"` (backcross, BC1) or `"biparental"` (BiparentalHaplotypeFinder,
#'   non-BC1) — the design-selected parent-calling route.
#' @param min_gametes,max_missing,max_iter Stage-1/3 knobs (TASSEL 200 / 1.0 / 50).
#' @param threads Fan-out width for the family x chromosome units. Each unit is
#'   fully independent (chromosome-separate processing, family-pooled), so it maps
#'   cleanly onto [parallel::mclapply()] (unix only; serial [lapply()] elsewhere).
#'   This is FSFHap's coarse-grained parallelism — the single-threaded per-unit
#'   stages 1+3 are untouched, so results are identical to `threads = 1L`.
#' @return `data.frame(source, donor, name, chr, pos, state)`, `state` in `{0,1,2}`;
#'   uncalled markers (dropped by the filters, or still missing after gap-fill) are
#'   simply absent. Feed to [to_segments()].
#' @keywords internal
.fsfhap_states <- function(data, phet, source = "nilHMM", donor = NA_character_,
                           has_donor = FALSE, route = "bc", min_gametes = 200L,
                           max_missing = 1.0, max_iter = 50L, threads = 1L) {
  if (!"family" %in% names(data)) stop(".fsfhap_states(): caller='fsfhap' needs a `family` column")
  if (!"g" %in% names(data)) stop(".fsfhap_states(): caller='fsfhap' needs a `g` genotype column")
  empty <- data.frame(source = character(), donor = character(), name = character(),
                      chr = integer(), pos = integer(), state = integer(), stringsAsFactors = FALSE)

  # Independent work units: one per (family, chromosome). Chromosome-separate
  # processing + family pooling means no unit shares state with another, so the
  # units fan out with no coordination (see `threads`).
  units <- list()
  for (fam in unique(data$family)) {
    fi <- data[data$family == fam, , drop = FALSE]
    for (ch in unique(fi$chr))
      units[[length(units) + 1L]] <- fi[fi$chr == ch, , drop = FALSE]
  }
  if (!length(units)) return(empty)

  # One (family, chromosome) unit: build the taxa x sites matrix, run stage 1
  # (route-dispatched parent-calling) then stage 3 (5-state EM + gap-fill), and
  # emit the per-marker states. Returns NULL when the unit has nothing to phase.
  one_unit <- function(ci) {
    ch   <- ci$chr[1L]
    taxa <- unique(ci$name); sites <- sort(unique(ci$pos))
    if (length(sites) < 2L || length(taxa) < 2L) return(NULL)  # nothing to phase
    G <- matrix(3L, length(taxa), length(sites))               # taxa x sites, default missing
    G[cbind(match(ci$name, taxa), match(ci$pos, sites))] <- as.integer(ci$g)
    storage.mode(G) <- "integer"
    s1 <- if (route == "biparental")
      .fsfhap_biparental_call(G, sites, min_gametes = min_gametes)
    else
      .fsfhap_call_parents_bc(G, sites, max_missing = max_missing,
                              min_gametes = min_gametes, min_r = 0.0)
    if (!length(s1$keep_sites) || !length(s1$keep_taxa)) return(NULL)
    imp <- .fsfhap_impute(s1$G, s1$pos, phet = phet, max_iter = max_iter, fill_gaps = TRUE)
    kt <- taxa[s1$keep_taxa]; kpos <- s1$pos
    nr <- length(kt); nc <- length(kpos)
    dvec <- if (has_donor) ci$donor[match(kt, ci$name)] else rep(donor, nr)
    st <- as.integer(imp)                                      # column-major: taxa fastest
    keep <- st != 3L
    data.frame(
      source = source, donor = rep(dvec, times = nc)[keep], name = rep(kt, times = nc)[keep],
      chr = as.integer(ch), pos = rep(kpos, each = nr)[keep], state = st[keep],
      stringsAsFactors = FALSE)
  }

  out <- if (threads > 1L && .Platform$OS.type == "unix")
           parallel::mclapply(units, one_unit, mc.cores = threads) else lapply(units, one_unit)
  if (any(vapply(out, function(x) inherits(x, "try-error"), logical(1))))
    stop(".fsfhap_states(): a family x chromosome unit failed to decode")
  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) return(empty)
  res <- do.call(rbind, out)
  # Radix on INTEGER factor codes, not a character multi-key order(): the latter
  # profiled at ~50% of full-genome wall-clock (a serial tail that capped the
  # per-unit parallel speedup). Factor levels sort in the session locale, so the
  # ordering matches the previous character sort. Same fix as `.lbimpute_prep`.
  o <- order(as.integer(factor(res$donor)), as.integer(factor(res$name)),
             as.integer(res$chr), as.integer(res$pos), method = "radix")
  res[o, , drop = FALSE]
}
