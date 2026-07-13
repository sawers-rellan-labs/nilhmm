# Design-driven cross simulator (simcross truth generator + generic degradation).
#
# The data-agnostic core: given a breeding design ("BC{n}S{m}") and a genetic map,
# simulate the true donor mosaic of one or more lines via kbroman/simcross (real
# meiosis / recombination), then optionally degrade the truth to observed allelic
# read counts under generic (depth, error, missingness) parameters. No file paths,
# no hardcoded map, no experiment-specific coverage regimes -- pipelines
# (e.g. zealhmm's simulate_source()) layer a real map + named source regimes on top.
# Uses simcross (Suggests).

# --- pedigree ---------------------------------------------------------------
# Founders 1 = recurrent (REF), 2 = donor (ALT). F1 = 1 x 2; then `n_bc` backcrosses
# to the recurrent (1); then `n_self` generations of selfing. Builds the pedigree
# for ANY BC_n S_m design. simcross selfing is sex-agnostic (mom == dad). Returns a
# base data.frame (simcross::check_pedigree/sim_from_pedigree index it as ped$col)
# plus the id of the terminal (fully-inbred) individual to sample.
.bcsft_pedigree <- function(n_bc, n_self) {
  id <- c(1L, 2L); mom <- c(0L, 0L); dad <- c(0L, 0L); gen <- c(0L, 0L)
  add <- function(mm, dd, g) {                       # append one progeny, return its id
    i <- length(id) + 1L
    id[i] <<- i; mom[i] <<- mm; dad[i] <<- dd; gen[i] <<- g
    i
  }
  cur <- add(1L, 2L, 1L)                             # F1
  for (k in seq_len(n_bc)) cur <- add(1L, cur, 1L + k)   # BCk = recurrent x prev
  g <- 1L + n_bc
  for (j in seq_len(n_self)) { g <- g + 1L; cur <- add(cur, cur, g) }  # self x n_self
  # sex is cosmetic: sim_from_pedigree takes the maternal gamete from `mom` and the
  # paternal from `dad` regardless of sex, and backcrosses put the recurrent as mom.
  sex <- ifelse(id %in% mom[mom > 0L], 0L, 1L)
  list(ped = data.frame(id = id, mom = mom, dad = dad, sex = sex, gen = gen),
       nil_id = as.character(cur))
}

# "BC1S4" -> list(n_bc = 1, n_self = 4)
.parse_bcsft <- function(design) {
  mm <- regmatches(design, regexec("^BC(\\d+)S(\\d+)$", design))[[1]]
  if (!length(mm))
    stop(".parse_bcsft(): `design` must be 'BC<n>S<m>' (e.g. 'BC1S4', 'BC2S2'), got '", design, "'")
  list(n_bc = as.integer(mm[2]), n_self = as.integer(mm[3]))
}

#' Build a physical marker grid from a genetic map
#'
#' Spread `n_markers` markers across the chromosomes of `map` in proportion to
#' physical span, interpolating each marker's cM from the map (monotone Hyman
#' spline, bp -> cM). This is how [simulate_nil()] samples the bundled consensus
#' map ([load_map()]) at a chosen density.
#'
#' @param map A map with columns `chr`, `cm`, `bp` (e.g. from [load_map()]).
#' @param n_markers Approximate total marker count across all chromosomes.
#' @return data.frame `chr`, `pos` (bp), `cm`, sorted by `chr, pos`.
#' @seealso [load_map()], [simulate_nil()]
#' @examples
#' g <- build_marker_grid(load_map(), n_markers = 500)
#' table(g$chr)
#' @export
build_marker_grid <- function(map, n_markers = 2000L) {
  map <- as.data.frame(map, stringsAsFactors = FALSE)
  if (!all(c("chr", "cm", "bp") %in% names(map)))
    stop("build_marker_grid(): `map` needs columns chr, cm, bp")
  chrs <- sort(unique(map$chr))
  span <- vapply(chrs, function(ch) { b <- map$bp[map$chr == ch]; max(b) - min(b) }, numeric(1))
  quota <- pmax(2L, as.integer(round(n_markers * span / sum(span))))
  do.call(rbind, lapply(seq_along(chrs), function(i) {
    sub <- map[map$chr == chrs[i], , drop = FALSE]
    cm_by_bp <- tapply(sub$cm, sub$bp, mean)                # collapse duplicate bp
    bp_u <- as.numeric(names(cm_by_bp)); o <- order(bp_u)
    f <- stats::splinefun(bp_u[o], as.numeric(cm_by_bp)[o], method = "hyman")  # monotone bp -> cM
    gbp <- round(seq(min(bp_u), max(bp_u), length.out = quota[i]))
    data.frame(chr = as.integer(chrs[i]), pos = as.integer(gbp),
               cm = as.numeric(f(gbp)), stringsAsFactors = FALSE)
  }))
}

# allele carried by a simcross haplotype (locations + alleles) at each cM position
.hap_allele_at <- function(hap, cm) {
  idx <- findInterval(cm, hap$locations, left.open = TRUE) + 1L
  idx[idx > length(hap$alleles)] <- length(hap$alleles)
  idx[idx < 1L] <- 1L
  hap$alleles[idx]
}

#' Simulate the true donor mosaic of a cross via simcross
#'
#' Design-driven truth generator: build the pedigree for a `"BC{n}S{m}"` design
#' (F1, then `n` backcrosses to the recurrent parent, then `m` selfings), simulate
#' it with \pkg{simcross} (real meiosis / crossover interference on a genetic map),
#' and read off each line's **true** donor dosage at the markers. The result is a
#' per-marker table with `state` in `{0 REF, 1 HET, 2 ALT}` (ALT = donor) that feeds
#' [to_segments()] (truth segments), [paint_calls()] (a truth track), or
#' [simulate_counts()] (degrade to observed read counts). Recurrent parent = REF = A;
#' donor = ALT.
#'
#' Data-agnostic: defaults to the bundled B73 v5 consensus map ([load_map()]),
#' sampled to `n_markers` by [build_marker_grid()]; pass your own `map` to override.
#' Hardcodes no design parameters (the design token drives the pedigree). For
#' designs beyond backcross-self (RIL, AIL, MAGIC, ...), pass a ready simcross
#' `pedigree` + `nil_id` and `design` is used only to label output.
#'
#' @param design Breeding design `"BC{n}S{m}"` (e.g. `"BC2S2"`, `"BC1S4"`). Ignored
#'   when `pedigree` is supplied.
#' @param n Number of lines to simulate (independent meioses).
#' @param map Reference genetic map with `chr`, `cm`, `bp`. `NULL` (default) uses
#'   the bundled [load_map()] (`maize_map_v5`).
#' @param n_markers Approximate marker count to sample from `map` ([build_marker_grid()]).
#' @param chr Optional integer vector restricting to a subset of chromosomes.
#' @param m,p simcross Stahl interference (`m = 10`, `p = 0` ~ maize-like).
#' @param seed Optional RNG seed for reproducibility.
#' @param donor Donor/ALT label for the output `donor` column.
#' @param names Optional length-`n` sample names (default `sim0001` ...).
#' @param pedigree,nil_id Escape hatch: a ready simcross pedigree
#'   (`id, mom, dad, sex, gen`) and the id to sample; bypasses `design`.
#' @return data.frame `source, donor, name, chr, pos, cm, state` (`state` 0/1/2),
#'   one row per (line, marker), sorted by `name, chr, pos`.
#' @seealso [load_map()], [build_marker_grid()], [simulate_counts()], [to_segments()].
#' @examples
#' if (requireNamespace("simcross", quietly = TRUE)) {
#'   truth <- simulate_nil("BC2S2", n = 2, chr = 1:2, n_markers = 200, seed = 1)
#'   head(truth)
#'   to_segments(truth)          # true donor segments
#' }
#' @export
simulate_nil <- function(design = "BC2S2", n = 1L, map = NULL, n_markers = 2000L,
                         chr = NULL, m = 10L, p = 0, seed = NULL, donor = "B",
                         names = NULL, pedigree = NULL, nil_id = NULL) {
  if (!requireNamespace("simcross", quietly = TRUE))
    stop("simulate_nil() needs the 'simcross' package (kbroman/simcross)")
  if (!is.null(seed)) set.seed(seed)
  if (is.null(map)) map <- load_map()
  map <- as.data.frame(map, stringsAsFactors = FALSE)
  if (!all(c("chr", "cm", "bp") %in% names(map)))
    stop("simulate_nil(): `map` needs columns chr, cm, bp")
  if (!is.null(chr)) map <- map[map$chr %in% chr, , drop = FALSE]
  if (!nrow(map)) stop("simulate_nil(): no markers left after the `chr` subset")
  grid <- build_marker_grid(map, n_markers)                # sampled marker positions
  chrs <- sort(unique(grid$chr))
  L <- vapply(chrs, function(ch) max(map$cm[map$chr == ch]), numeric(1))  # true genetic length

  if (!is.null(pedigree)) {
    if (is.null(nil_id)) stop("simulate_nil(): supply `nil_id` with a custom `pedigree`")
    ped <- as.data.frame(pedigree); nid <- as.character(nil_id)
  } else {
    pd <- .parse_bcsft(design); bp <- .bcsft_pedigree(pd$n_bc, pd$n_self)
    ped <- bp$ped; nid <- bp$nil_id
  }
  nms <- if (!is.null(names)) as.character(names) else sprintf("sim%04d", seq_len(n))
  if (length(nms) != n) stop("simulate_nil(): `names` must have length `n`")

  one <- function(nm) {
    sim <- simcross::sim_from_pedigree(ped, L = L, m = m, p = p)
    # simcross nests per-individual for a single chromosome (sim[[id]]) but
    # per-chromosome for several (sim[[chr]][[id]]) -- handle both.
    hap_of <- if (length(L) == 1L) function(ci) sim[[nid]] else function(ci) sim[[ci]][[nid]]
    dose <- integer(nrow(grid))
    for (ci in seq_along(chrs)) {
      rows <- which(grid$chr == chrs[ci]); if (!length(rows)) next
      hap <- hap_of(ci); cm <- grid$cm[rows]
      dose[rows] <- (.hap_allele_at(hap$mat, cm) == 2L) + (.hap_allele_at(hap$pat, cm) == 2L)
    }
    data.frame(source = "sim", donor = donor, name = nm,
               chr = as.integer(grid$chr), pos = as.integer(grid$pos), cm = grid$cm,
               state = as.integer(dose), stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, lapply(nms, one))
  out[order(out$name, out$chr, out$pos), , drop = FALSE]
}

#' Degrade simulated truth to observed allelic read counts
#'
#' Turn a [simulate_nil()] truth table into observed data under a generic
#' sequencing model: per marker draw depth `~ Poisson(depth)`, mask a fraction
#' `p_missing` (and zero-depth markers) to missing, and draw ALT reads
#' `~ Binomial(depth, p)` where `p` is the state's donor fraction (`0` / `0.5` / `1`)
#' flipped by the per-read `error`. Returns the `read_counts()`-style columns plus a
#' hard genotype `g` (masked to `3` where missing), so the result drives both the
#' count callers (`n_ref`/`n_alt`) and the genotype callers (`g`). The named
#' experiment-specific coverage regimes (skim / BRB / target) live in the pipeline,
#' not here.
#'
#' @param truth A [simulate_nil()] table (needs `name, donor, chr, pos, state`).
#' @param depth Mean per-marker read depth (Poisson mean).
#' @param error Per-read error rate (homozygote miscall probability).
#' @param p_missing Extra fraction of markers set missing (on top of zero-depth).
#' @param seed Optional RNG seed.
#' @return data.frame `name, donor, chr, pos, n_ref, n_alt, g` (`g` in `{0,1,2,3}`).
#' @seealso [simulate_nil()], [read_counts()].
#' @examples
#' if (requireNamespace("simcross", quietly = TRUE)) {
#'   truth <- simulate_nil("BC2S2", n = 2, chr = 1:2, n_markers = 100, seed = 1)
#'   obs <- simulate_counts(truth, depth = 6, seed = 1)
#'   head(obs)
#' }
#' @export
simulate_counts <- function(truth, depth = 1, error = 0.01, p_missing = 0, seed = NULL) {
  need <- c("name", "chr", "pos", "state")
  if (!all(need %in% names(truth)))
    stop("simulate_counts(): `truth` needs columns ", paste(need, collapse = ", "),
         " (from simulate_nil())")
  if (!is.null(seed)) set.seed(seed)
  N <- nrow(truth)
  d <- stats::rpois(N, depth)
  d[stats::runif(N) < p_missing] <- 0L                 # extra missingness
  present <- d > 0L
  p_alt <- c(0, 0.5, 1)[truth$state + 1L]
  p_eff <- p_alt * (1 - error) + (1 - p_alt) * error
  alt <- stats::rbinom(N, d, p_eff)
  g <- truth$state; g[!present] <- 3L
  data.frame(
    name = truth$name, donor = if ("donor" %in% names(truth)) truth$donor else NA_character_,
    chr = as.integer(truth$chr), pos = as.integer(truth$pos),
    n_ref = as.integer(ifelse(present, d - alt, 0L)),
    n_alt = as.integer(ifelse(present, alt, 0L)),
    g = as.integer(g), stringsAsFactors = FALSE)
}

#' Simulate a tracked forest of BC{n}S{m} sibling families (shared ancestry)
#'
#' Unlike [simulate_nil()] (independent lines with no relatives), this draws ONE
#' \pkg{simcross} meiosis per family and reads off `sibs` terminal siblings from
#' that single draw, so the sibs genuinely share the IBD blocks inherited through
#' the (latent, ungenotyped) founder -> selfing chain. This is the ground-truth
#' generator for validating a pedigree-aware refinement ([refine_ancestry()]):
#' the coupling refine exploits only exists when relatives share ancestry.
#'
#' Returns both the per-marker `truth` for the genotyped sibs (a [simulate_nil()]
#' table, so it feeds [simulate_counts()] / [to_segments()]) and a
#' [read_pedigree()]-shaped `pedigree`: one row per taxon, with the founder and
#' the intermediate selfing nodes present as **parent-only latent rows** (no
#' genotype), the sibs as leaves. Family `f` is a star: founder `g0` -> `g1` ->
#' ... -> `g{m-1}` -> `sibs` leaves; families are independent (independent founders).
#'
#' @param design "BC{n}S{m}" (default "BC2S3"); needs >= 1 selfing generation.
#' @param families Number of independent families (independent founders).
#' @param sibs Genotyped siblings per family (they share the founder + chain).
#' @param map,n_markers,chr,m,p,seed,donor As in [simulate_nil()].
#' @param prefix Family-id prefix for taxon/family names (default "fam").
#' @return `list(truth = <per-marker table: source, donor, name, family, chr,
#'   pos, cm, state>, pedigree = <data.frame taxon, family, parent1, parent2>)`.
#' @seealso [simulate_nil()], [simulate_counts()], [read_pedigree()], [refine_ancestry()].
#' @examples
#' if (requireNamespace("simcross", quietly = TRUE)) {
#'   fam <- simulate_family("BC2S3", families = 2, sibs = 3, chr = 1,
#'                          n_markers = 100, seed = 1)
#'   head(fam$pedigree)
#'   table(fam$truth$name)
#' }
#' @export
simulate_family <- function(design = "BC2S3", families = 10L, sibs = 10L,
                            map = NULL, n_markers = 1000L, chr = NULL,
                            m = 10L, p = 0, seed = NULL, donor = "B",
                            prefix = "fam") {
  if (!requireNamespace("simcross", quietly = TRUE))
    stop("simulate_family() needs the 'simcross' package (kbroman/simcross)")
  posint <- function(x) length(x) == 1L && is.finite(x) && x >= 1 && x == as.integer(x)
  if (!posint(families)) stop("simulate_family(): `families` must be a positive integer")
  if (!posint(sibs))     stop("simulate_family(): `sibs` must be a positive integer")
  if (!is.null(seed)) set.seed(seed)
  if (is.null(map)) map <- load_map()
  map <- as.data.frame(map, stringsAsFactors = FALSE)
  if (!all(c("chr", "cm", "bp") %in% names(map)))
    stop("simulate_family(): `map` needs columns chr, cm, bp")
  if (!is.null(chr)) map <- map[map$chr %in% chr, , drop = FALSE]
  if (!nrow(map)) stop("simulate_family(): no markers left after the `chr` subset")
  pd <- .parse_bcsft(design)
  if (pd$n_self < 1L)
    stop("simulate_family(): `design` needs >= 1 selfing generation (got ", design, ")")
  grid <- build_marker_grid(map, n_markers)
  chrs <- sort(unique(grid$chr))
  L <- vapply(chrs, function(ch) max(map$cm[map$chr == ch]), numeric(1))

  # simcross pedigree for ONE family: 1 = RP, 2 = donor, F1, then n_bc backcrosses
  # to the recurrent -> founder (selfing depth 0); then (n_self - 1) SINGLE selfings
  # (the latent chain g1..g{n_self-1}); then `sibs` terminal sibs off the last node.
  id <- c(1L, 2L); mom <- c(0L, 0L); dad <- c(0L, 0L)
  addn <- function(mm, dd) { i <- length(id) + 1L; id[i] <<- i; mom[i] <<- mm; dad[i] <<- dd; i }
  cur <- addn(1L, 2L)                                  # F1
  for (k in seq_len(pd$n_bc)) cur <- addn(1L, cur)     # BCk = recurrent x prev
  for (j in seq_len(pd$n_self - 1L)) cur <- addn(cur, cur)  # single selfings g1..g{m-1}
  sib_ids <- vapply(seq_len(sibs), function(s) addn(cur, cur), integer(1))  # terminal sibs
  ped_sc <- data.frame(id = id, mom = mom, dad = dad,
                       sex = ifelse(id %in% mom[mom > 0L], 0L, 1L), gen = 0L)

  truth_rows <- vector("list", families * sibs); ti <- 0L
  ped_rows <- vector("list", families)
  for (f in seq_len(families)) {
    fam    <- sprintf("%s%02d", prefix, f)
    g_lab  <- sprintf("%s_g%d", fam, seq_len(pd$n_self) - 1L)  # depth 0 (founder) .. n_self-1
    sib_lab<- sprintf("%s_L%02d", fam, seq_len(sibs))
    n_lat  <- length(g_lab)
    ped_rows[[f]] <- data.frame(
      taxon   = c(g_lab, sib_lab), family = fam,
      parent1 = c(NA_character_, g_lab[-n_lat], rep(g_lab[n_lat], sibs)),
      parent2 = c(NA_character_, g_lab[-n_lat], rep(g_lab[n_lat], sibs)),
      stringsAsFactors = FALSE)

    sim <- simcross::sim_from_pedigree(ped_sc, L = L, m = m, p = p)  # one draw / family
    hap_of <- if (length(L) == 1L) function(ci, iid) sim[[as.character(iid)]]
              else                 function(ci, iid) sim[[ci]][[as.character(iid)]]
    for (s in seq_len(sibs)) {
      dose <- integer(nrow(grid))
      for (ci in seq_along(chrs)) {
        rows <- which(grid$chr == chrs[ci]); if (!length(rows)) next
        hap <- hap_of(ci, sib_ids[s]); cm <- grid$cm[rows]
        dose[rows] <- (.hap_allele_at(hap$mat, cm) == 2L) + (.hap_allele_at(hap$pat, cm) == 2L)
      }
      ti <- ti + 1L
      truth_rows[[ti]] <- data.frame(
        source = "sim", donor = donor, name = sib_lab[s], family = fam,
        chr = as.integer(grid$chr), pos = as.integer(grid$pos), cm = grid$cm,
        state = as.integer(dose), stringsAsFactors = FALSE)
    }
  }
  truth <- do.call(rbind, truth_rows)
  list(truth    = truth[order(truth$name, truth$chr, truth$pos), , drop = FALSE],
       pedigree = do.call(rbind, ped_rows))
}
