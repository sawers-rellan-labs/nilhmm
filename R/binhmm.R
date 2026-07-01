# binhmm: per-bin ancestry by clustering + HMM smoothing. Faithful to the rpubs
# pipeline (rpubs.com/faustovrz/1306822, "...Ancestry Analysis by bins"):
#   1. bin the genome (default 1 Mb); per bin, ALT_FREQ = alt reads / total reads,
#      weighted by informative-variant count.
#   2. zero-ALT bins -> REF; the rest get K=3 1-D clustering (GMM or K-means),
#      relabeled REF/HET/ALT by ascending mean ALT_FREQ.
#   3. per-chromosome HMM smoothing (Viterbi) with breeding-design (Mendelian)
#      start probabilities, a fixed confusion emission, and a sticky transition.
# Dependency-light: a small 1-D GMM-EM + base stats::kmeans (no rebmix / HMM pkg);
# the HMM reuses the package's viterbi_log_cpp.

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
.gmm1d <- function(x, k = 3L, iter = 200L, tol = 1e-6) {
  mu <- as.numeric(stats::quantile(x, seq(0, 1, length.out = k + 2L)[2:(k + 1L)]))
  v  <- rep(stats::var(x) + 1e-8, k); w <- rep(1 / k, k); n <- length(x)
  dens_mat <- function() vapply(1:k, function(j) w[j] * stats::dnorm(x, mu[j], sqrt(v[j])), numeric(n))
  for (it in 1:iter) {
    d <- dens_mat(); tot <- rowSums(d); tot[tot == 0] <- 1e-300
    resp <- d / tot; Nk <- colSums(resp); Nk[Nk == 0] <- 1e-300
    mu_new <- colSums(resp * x) / Nk
    v_new  <- vapply(1:k, function(j) sum(resp[, j] * (x - mu_new[j])^2) / Nk[j], numeric(1))
    v_new[v_new < 1e-8] <- 1e-8; w <- Nk / n
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
.binhmm_cluster <- function(alt_freq, method = "gmm") {
  st <- integer(length(alt_freq))                 # all REF (0) by default
  nz <- which(alt_freq > 0)
  x  <- alt_freq[nz]
  if (length(nz) >= 3L && length(unique(x)) >= 2L) {
    k  <- min(3L, length(unique(x)))
    cl <- switch(method,
                 kmeans = stats::kmeans(x, centers = k, nstart = 10L)$cluster,
                 rebmix = .binhmm_rebmix(x, k),
                 .gmm1d(x, k))                     # "gmm" default
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

# --- RLE a chromosome's smoothed bin states into common-schema segments --------
.binhmm_segments <- function(states, start_bp, end_bp) {
  n <- length(states); brk <- which(states[-1L] != states[-n])
  s <- c(1L, brk + 1L); e <- c(brk, n)
  data.frame(start_bp = as.integer(start_bp[s]), end_bp = as.integer(end_bp[e]),
             state = as.integer(states[s]))
}

# --- caller backend: bin -> cluster -> smooth -> common schema -----------------
.call_ancestry_binhmm <- function(data, bin_size, cluster_method, priors,
                                  source, donor, has_donor, stay = 0.995) {
  bins <- .binhmm_bin(data, bin_size)
  start <- c(1 - priors$f_1 - priors$f_2, priors$f_1, priors$f_2)   # REF/HET/ALT (Mendelian)
  donor_of <- if (has_donor) tapply(as.character(data$donor), data$name, `[`, 1L) else NULL
  out <- list()
  for (nm in unique(bins$name)) {
    b <- bins[bins$name == nm, , drop = FALSE]
    w <- b$alt_freq * (b$ninf / mean(b$ninf)); w <- pmin(pmax(w, 0), 1)  # informative-count weighting
    st <- .binhmm_cluster(w, cluster_method)                            # per-sample K=3 cluster
    dn <- if (has_donor) donor_of[[nm]] else donor
    for (cc in unique(b$chr)) {
      i <- which(b$chr == cc); i <- i[order(b$bin[i])]
      sm  <- .binhmm_smooth(st[i], start, stay)                        # HMM-smoothed bins
      seg <- .binhmm_segments(sm, b$start_bp[i], b$end_bp[i])
      out[[length(out) + 1L]] <- data.frame(source = source, donor = dn, name = nm,
                                            chr = cc, seg, stringsAsFactors = FALSE)
    }
  }
  calls <- do.call(rbind, out)
  calls[order(calls$donor, calls$name, calls$chr, calls$start_bp), , drop = FALSE]
}
