# binhmm: per-bin ancestry from binned alt-freq. bin the genome (default 1 Mb;
# per bin ALT_FREQ = alt reads / total reads), then genotype each bin by one of
# two backends (cluster_method):
#   "gauss" (DEFAULT): a 3-state Gaussian-emission HMM with REF anchored at the
#     baseline floor, HET/ALT seeded + hard-EM, collapse guards, de-speckle
#     (.binhmm_gauss_states). Fixes the HET over-call and high-coverage
#     fragmentation of the forced K=3 clustering.
#   "gmm"/"kmeans"/"rebmix": the original rpubs route
#     (rpubs.com/faustovrz, "...Ancestry Analysis by bins") -- K=3 1-D cluster the
#     non-zero bins, relabel REF/HET/ALT by ascending mean, then a per-chromosome
#     fixed-confusion sticky Viterbi. Supports joint_clust / obs_weights.
# Dependency-light: base-R GMM-EM + stats::kmeans (rebmix optional/Suggests, HMM
# reuses the package's viterbi_log_cpp).

# --- vectorized binning (the "genome walker", done as a base-R group-by) -------
.binhmm_bin <- function(data, bin_size) {
  bin <- data$pos %/% bin_size
  key <- paste(data$name, data$chr, bin, sep = "\r")
  tot <- data$n_ref + data$n_alt
  sums <- rowsum(cbind(alt = data$n_alt, tot = tot, ninf = as.integer(tot > 0)), key)
  sb <- tapply(data$pos, key, min); eb <- tapply(data$pos, key, max)
  k  <- rownames(sums); p <- do.call(rbind, strsplit(k, "\r", fixed = TRUE))
  data.frame(name = p[, 1], chr = as.integer(p[, 2]), bin = as.integer(p[, 3]),
             start_bp = as.integer(sb[k]), end_bp = as.integer(eb[k]),
             alt = sums[, "alt"], tot = sums[, "tot"], ninf = sums[, "ninf"],
             alt_freq = ifelse(sums[, "tot"] > 0, sums[, "alt"] / sums[, "tot"], 0),
             stringsAsFactors = FALSE)
}

# --- 1-D 3-component Gaussian mixture via EM (replaces rebmix) ------------------
# `weights` (optional): per-point observation weights, folded into the E-step
# responsibilities so high-weight points pull the component means/vars harder.
# It weights the FIT, not the x values (the alt-freq axis is untouched), and the
# final assignment is the ordinary argmax posterior. weights = NULL reproduces
# the unweighted EM exactly.
.gmm1d <- function(x, k = 3L, iter = 200L, tol = 1e-6, weights = NULL) {
  n <- length(x); ow <- if (is.null(weights)) rep(1, n) else weights; sw <- sum(ow)
  mu <- as.numeric(stats::quantile(x, seq(0, 1, length.out = k + 2L)[2:(k + 1L)]))
  v  <- rep(stats::var(x) + 1e-8, k); w <- rep(1 / k, k)
  dens_mat <- function() vapply(1:k, function(j) w[j] * stats::dnorm(x, mu[j], sqrt(v[j])), numeric(n))
  for (it in 1:iter) {
    d <- dens_mat(); tot <- rowSums(d); tot[tot == 0] <- 1e-300
    resp <- (d / tot) * ow; Nk <- colSums(resp); Nk[Nk == 0] <- 1e-300
    mu_new <- colSums(resp * x) / Nk
    v_new  <- vapply(1:k, function(j) sum(resp[, j] * (x - mu_new[j])^2) / Nk[j], numeric(1))
    v_new[v_new < 1e-8] <- 1e-8; w <- Nk / sw
    if (max(abs(mu_new - mu)) < tol) { mu <- mu_new; v <- v_new; break }
    mu <- mu_new; v <- v_new
  }
  max.col(dens_mat(), ties.method = "first")   # cluster per point
}

# --- optional rebmix GMM backend (Suggests) -----------------------------------
# The rpubs Kgmm used rebmix REBMIX/RCLRMIX. It is NOT a hard dependency (the
# base-R .gmm1d reproduces the rpubs Kgmm_HMM calls to ~99% after smoothing, and
# a 1-D 3-component mixture does not warrant rebmix's weight). This backend is
# offered only for users who want bit-exact rpubs reproduction; it is gated on
# rebmix being installed and falls back to .gmm1d when <3 distinct values make
# REBMIX's forced 3 components infeasible.
.binhmm_rebmix <- function(x, k) {
  if (!requireNamespace("rebmix", quietly = TRUE))
    stop("cluster_method='rebmix' needs the 'rebmix' package (Suggests); ",
         "install.packages('rebmix') or use 'gmm'/'kmeans'.")
  if (k < 3L) return(.gmm1d(x, k))                # REBMIX forces cmin=cmax=3
  quiet <- function(e) utils::capture.output(suppressMessages(suppressWarnings(v <- e)))
  quiet(est <- rebmix::REBMIX(
    Dataset = list(data.frame(Value = x)), Preprocessing = "histogram",
    cmin = 3, cmax = 3, Criterion = "BIC", pdf = "normal"))
  quiet(cl <- rebmix::RCLRMIX(x = est))
  as.integer(cl@Zp)
}

# --- cluster a sample's bins into 0/1/2 (REF/HET/ALT) --------------------------
# zero-ALT -> REF; non-zero -> K=3 cluster on alt_freq, relabel by ascending mean.
# Backend: "gmm" (base-R .gmm1d, default), "kmeans" (stats::kmeans), or "rebmix"
# (the rpubs GMM, optional Suggests). Degenerate clusterings are handled
# explicitly: K is capped at the number of DISTINCT non-zero values (so K-means
# never asks for more centres than distinct points), and if fewer than 3 clusters
# actually emerge (e.g. a clean bimodal sample where a GMM component wins
# nothing), the survivors are spread across the extremes -- 2 clusters ->
# {REF, ALT}, 1 -> {REF} -- so a clearly-ALT group is never silently collapsed
# into REF.
.binhmm_cluster <- function(alt_freq, method = "gmm", weights = NULL) {
  st <- integer(length(alt_freq))                 # all REF (0) by default
  nz <- which(alt_freq > 0)
  x  <- alt_freq[nz]
  if (length(nz) >= 3L && length(unique(x)) >= 2L) {
    k  <- min(3L, length(unique(x)))
    wz <- if (is.null(weights)) NULL else weights[nz]
    cl <- switch(method,
                 kmeans = { if (!is.null(weights)) stop("cluster_method='kmeans' has no observation-weight support; use 'gmm'")
                            stats::kmeans(x, centers = k, nstart = 10L)$cluster },
                 rebmix = { if (!is.null(weights)) stop("cluster_method='rebmix' has no observation-weight support; use 'gmm'")
                            .binhmm_rebmix(x, k) },
                 .gmm1d(x, k, weights = wz))       # "gmm" default
    m       <- tapply(x, cl, mean)                # emergent cluster means
    present <- as.integer(names(m))
    lab_by  <- switch(as.character(length(present)),
                      "1" = 0L, "2" = c(0L, 2L), c(0L, 1L, 2L))  # REF / REF+ALT / REF+HET+ALT
    labels  <- integer(length(present)); labels[order(m)] <- lab_by  # ascending mean
    map     <- integer(max(present)); map[present] <- labels
    st[nz]  <- map[cl]
  }
  st
}

# --- HMM smoothing over a chromosome's bins (reuses viterbi_log_cpp) -----------
# obs = cluster labels 0/1/2; start = design REF/HET/ALT freqs; fixed confusion
# emission + sticky transition, exactly as the rpubs smooth_ancestry_with_hmm.
.BINHMM_EMISS <- matrix(c(0.90, 0.08, 0.02,
                          0.10, 0.80, 0.10,
                          0.02, 0.08, 0.90), nrow = 3, byrow = TRUE)  # rows=state, cols=obs
.binhmm_smooth <- function(obs, start, stay = 0.995) {
  sw <- (1 - stay) / 2
  trans <- matrix(sw, 3, 3); diag(trans) <- stay
  log_emit <- t(log(.BINHMM_EMISS)[, obs + 1L, drop = FALSE])       # T x 3
  viterbi_log_cpp(log(start), log(trans), log_emit)                # 0-indexed 0/1/2
}

# --- corrected binhmm stage: 3-state GAUSSIAN-emission HMM on per-bin alt_freq --
# The default backend (cluster_method = "gauss"). Replaces the K=3 GMM, whose
# forced third component carves a spurious HET cluster out of REF-baseline noise
# (HET over-call), and which fragments at high coverage. Instead:
#   * REF is ANCHORED at the empirical baseline floor (the bulk-low mean), NOT 0
#     -- the diluted read alt_freq REF floor is a small positive value;
#   * the REF Gaussian width is FLOORED to the signal scale so a hyper-tight
#     baseline (high coverage) can't reject mildly-elevated single bins into donor;
#   * HET/ALT means are SEEDED from the data (ALT = upper quantile of clearly
#     elevated bins, HET = midway) and refined by Viterbi hard-EM with REF fixed;
#   * sticky transitions come from the design start priors;
#   * collapse guards fold a state not separated from the floor back to REF, and
#     donor runs shorter than `min_run` bins are dropped (de-speckle).
# Ports the nilhifi chromosome_painting.R anchoring to the read-alt_freq scale.
# `alt_freq`/`chr` must be ordered by (chr, bin). Returns final 0/1/2 states.
.binhmm_gauss_states <- function(alt_freq, chr, start, stay = 0.995, min_run = 2L) {
  n <- length(alt_freq); af <- alt_freq
  base <- af[af <= stats::quantile(af, 0.6)]
  ref_mean <- mean(base); sd_base <- max(1e-3, stats::sd(base))
  hi <- af[af > ref_mean + 4 * sd_base]
  if (length(hi) < 4L) return(integer(n))                # collapse guard: no donor -> all REF
  alt_mean <- as.numeric(stats::quantile(hi, 0.75)); het_mean <- ref_mean + (alt_mean - ref_mean) / 2
  sd_ref <- max(sd_base, (alt_mean - ref_mean) / 6)      # floor REF width to the signal scale
  sd_hi  <- max(sd_ref, stats::sd(hi))
  means <- c(ref_mean, het_mean, alt_mean); sds <- c(sd_ref, sd_hi, sd_hi)
  sw <- (1 - stay) / 2; ltrans <- log({ tr <- matrix(sw, 3, 3); diag(tr) <- stay; tr }); lstart <- log(start)
  chrs <- unique(chr)
  emit <- function() vapply(1:3, function(k) stats::dnorm(af, means[k], max(sds[k], 1e-3), log = TRUE), numeric(n))
  path <- integer(n)
  for (it in 1:8) {                                      # hard-EM: Viterbi -> refit HET/ALT means (REF fixed)
    le <- emit(); np <- integer(n)
    for (cc in chrs) { i <- which(chr == cc); np[i] <- viterbi_log_cpp(lstart, ltrans, le[i, , drop = FALSE]) }
    if (identical(np, path)) break
    path <- np
    for (k in 2:3) if (sum(path == k - 1L) > 0L) {
      means[k] <- mean(af[path == k - 1L])
      v <- stats::sd(af[path == k - 1L]); sds[k] <- max(sd_ref, if (is.na(v)) sd_ref else v)
    }
  }
  if (means[2] < ref_mean + 4 * sd_ref) path[path == 1L] <- 0L   # HET indistinct from floor -> REF
  if (means[3] < ref_mean + 4 * sd_ref) path[path == 2L] <- 0L   # ALT indistinct from floor -> REF
  for (cc in chrs) { i <- which(chr == cc); rr <- rle(path[i])   # drop < min_run donor specks
    rr$values[rr$values > 0L & rr$lengths < min_run] <- 0L; path[i] <- inverse.rle(rr) }
  path
}

# --- RLE a chromosome's smoothed bin states into common-schema segments --------
.binhmm_segments <- function(states, start_bp, end_bp) {
  n <- length(states); brk <- which(states[-1L] != states[-n])
  s <- c(1L, brk + 1L); e <- c(brk, n)
  data.frame(start_bp = as.integer(start_bp[s]), end_bp = as.integer(end_bp[e]),
             state = as.integer(states[s]))
}

# --- caller backend: bin -> genotype bins -> common schema ---------------------
# Two per-bin genotyping backends, selected by `cluster_method`:
#   "gauss" (DEFAULT): the anchored 3-state Gaussian-emission HMM
#     (.binhmm_gauss_states) -- REF pinned at the baseline floor, HET/ALT seeded +
#     hard-EM, collapse guards, de-speckle. Fixes the K=3 HET over-call and the
#     high-coverage fragmentation. Per-sample; `joint_clust`/`obs_weights` do not
#     apply.
#   "gmm"/"kmeans"/"rebmix": the older cluster-then-smooth route (the rpubs
#     pipeline) -- K=3 cluster the per-bin alt-freq, then a fixed-confusion sticky
#     Viterbi (.binhmm_smooth). Supports `joint_clust` (pool all samples, one
#     shared fit on RAW alt-freq, a la get_joint_ancestry_calls.R) and, for "gmm",
#     `obs_weights` (informative-count weights in the fit -- influence, not value).
# HMM smoothing/decoding is always per (sample, chromosome).
.binhmm_states <- function(data, bin_size, cluster_method, priors,
                           source, donor, has_donor, stay = 0.995,
                           joint_clust = FALSE, obs_weights = FALSE) {
  # Accept EITHER raw per-SNP counts (binned here) OR a table that is already
  # binned upstream, carrying `alt_freq` + `start_bp`/`end_bp` (e.g. the bzeaseq
  # `*_bin_genotypes.tsv`). Bin summarization is a data-prep concern, so both
  # paths are supported at the input boundary; see call_states().
  bins <- if ("alt_freq" %in% names(data)) {
    if (!all(c("start_bp", "end_bp") %in% names(data))) {
      stop("binhmm: pre-binned input (`alt_freq`) also needs `start_bp` and `end_bp`")
    }
    data.frame(
      name = data$name, chr = as.integer(data$chr), bin = as.integer(data$pos),
      start_bp = as.integer(data$start_bp), end_bp = as.integer(data$end_bp),
      alt_freq = as.numeric(data$alt_freq),
      ninf = if ("ninf" %in% names(data)) as.integer(data$ninf) else 1L,
      stringsAsFactors = FALSE
    )
  } else {
    .binhmm_bin(data, bin_size)
  }
  start <- c(1 - priors$f_1 - priors$f_2, priors$f_1, priors$f_2)   # REF/HET/ALT (Mendelian)
  donor_of <- if (has_donor) tapply(as.character(data$donor), data$name, `[`, 1L) else NULL
  gauss <- identical(cluster_method, "gauss")

  # cluster backends: set per-bin RAW cluster labels bins$state now (the gauss
  # backend computes its FINAL states per-sample in the loop below).
  if (!gauss) {
    if (joint_clust) {
      ow <- if (obs_weights) as.numeric(bins$ninf) else NULL   # influence weights, never a value rescale
      bins$state <- .binhmm_cluster(bins$alt_freq, cluster_method, weights = ow)  # ONE pooled fit
    } else {
      if (obs_weights) warning("binhmm: obs_weights applies only to joint_clust = TRUE; ignored")
      st <- integer(nrow(bins))
      for (nm in unique(bins$name)) {                          # each sample clustered independently
        j <- which(bins$name == nm)
        w <- bins$alt_freq[j] * (bins$ninf[j] / mean(bins$ninf[j])); w <- pmin(pmax(w, 0), 1)
        st[j] <- .binhmm_cluster(w, cluster_method)
      }
      bins$state <- st
    }
  } else if (joint_clust || obs_weights) {
    warning("binhmm: joint_clust/obs_weights apply only to the cluster backends (gmm/kmeans/rebmix); ignored for cluster_method = 'gauss'")
  }

  out <- list()
  for (nm in unique(bins$name)) {
    b <- bins[bins$name == nm, , drop = FALSE]
    b <- b[order(b$chr, b$bin), , drop = FALSE]              # (chr, bin) order for per-chr decode
    dn <- if (has_donor) donor_of[[nm]] else donor
    fstate <- if (gauss) {
      .binhmm_gauss_states(b$alt_freq, b$chr, start, stay)   # final (Gaussian-HMM) states
    } else {
      fs <- integer(nrow(b))                                 # smooth the cluster labels per chr
      for (cc in unique(b$chr)) { i <- which(b$chr == cc); fs[i] <- .binhmm_smooth(b$state[i], start, stay) }
      fs
    }
    for (cc in unique(b$chr)) {
      i <- which(b$chr == cc)                                # b is (chr, bin)-ordered -> i is bin-ordered
      # coordinate-free per-bin states; bins are interval units so carry their
      # own start_bp/end_bp for segments() to reinstate.
      out[[length(out) + 1L]] <- data.frame(
        source = source, donor = dn, name = nm, chr = as.integer(cc),
        pos = as.integer(b$bin[i]), start_bp = as.integer(b$start_bp[i]),
        end_bp = as.integer(b$end_bp[i]), state = as.integer(fstate[i]),
        stringsAsFactors = FALSE)
    }
  }
  states <- do.call(rbind, out)
  states[order(states$donor, states$name, states$chr, states$pos), , drop = FALSE]
}

