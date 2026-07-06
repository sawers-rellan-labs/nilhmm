# Engine: the duration-aware 3-state HMM (REF/HET/ALT). Layer 1 of the
# architecture (REFACTOR_R_PACKAGE.md S4). Emission and duration are pluggable
# interfaces (see emissions.R, duration.R). Hot loops live in src/ (Rcpp:
# count_emission_loglik_cpp, viterbi_log_cpp).

# --- internal: transition + init from r and state-frequency priors ----------
# Identical parameterization to the Python _build_transition: self-transition
# 1-r, off-diagonal recombination weighted by expected state frequencies. The
# generation enters ONLY through f_1/f_2 here (a prior), never a baked-in
# duration (S7).
.build_transition <- function(r, f_1, f_2) {
  f_0 <- 1 - f_1 - f_2
  if (f_0 <= 0) stop(".build_transition: f_1 + f_2 must be < 1")
  tmat <- matrix(c(
    1 - r,                       r * f_1 / (f_1 + f_2),       r * f_2 / (f_1 + f_2),
    r * f_0 / (f_0 + f_2),       1 - r,                       r * f_2 / (f_0 + f_2),
    r * f_0 / (f_0 + f_1),       r * f_1 / (f_0 + f_1),       1 - r
  ), nrow = 3, byrow = TRUE)
  list(log_start = log(c(f_0, f_1, f_2)), log_trans = log(tmat), n_sub = 1L)
}

# --- internal: rigidity (RTIGER) transition via state expansion (S4, S7) ------
# Each macro-state (REF/HET/ALT) becomes a chain of `r` sub-states: positions
# 1..r-1 advance deterministically (forced minimum run length r), position r is
# a free state that either self-loops (stay in the macro-state) or exits to
# another macro-state's sub-state 1, with probabilities from the macro
# transition `A` (built from `p_switch` + priors). r = 1 collapses to `A`
# itself. Expanded state index = m * r + k (m, k 0-based). All sub-states of a
# macro-state share its emission (handled in decode()).
.build_transition_rigidity <- function(r, p_switch, f_1, f_2) {
  macro <- .build_transition(p_switch, f_1, f_2)
  A     <- exp(macro$log_trans)              # 3 x 3 macro transition
  start <- exp(macro$log_start)              # c(f_0, f_1, f_2)
  n <- as.integer(r); K <- 3L; S <- K * n
  idx <- function(m, k) m * n + k + 1L       # m,k 0-based -> 1-based
  Tm <- matrix(0, S, S)
  for (m in 0:2) for (k in 0:(n - 1L)) {
    if (k < n - 1L) {
      Tm[idx(m, k), idx(m, k + 1L)] <- 1                       # forced advance
    } else {
      Tm[idx(m, k), idx(m, n - 1L)] <- A[m + 1L, m + 1L]       # free state: stay
      for (mp in 0:2) if (mp != m) Tm[idx(m, k), idx(mp, 0L)] <- A[m + 1L, mp + 1L]  # exit
    }
  }
  s0 <- numeric(S); for (m in 0:2) s0[idx(m, 0L)] <- start[m + 1L]
  list(log_start = log(s0), log_trans = log(Tm), n_sub = n)
}

# Dispatch: (duration, priors) -> list(log_start, log_trans, n_sub).
.duration_transition <- function(duration, priors) {
  if (inherits(duration, "nilHMM_duration_geometric"))
    return(.build_transition(duration$r, priors$f_1, priors$f_2))
  if (inherits(duration, "nilHMM_duration_rigidity"))
    return(.build_transition_rigidity(duration$r, duration$p_switch, priors$f_1, priors$f_2))
  stop(".duration_transition(): unsupported duration '", duration$type,
       "' (hsmm is reserved, S4)")
}

# --- internal: run-length-encode a state path into common-schema segments ----
# Boundaries where the state changes; start_bp/end_bp = first/last marker
# position of each contiguous run (same convention as the Python rle_segments
# and the RTIGER segments).
.rle_segments <- function(path, pos, chr, name, source, donor) {
  if (length(path) == 0) return(NULL)
  r <- rle_segments_cpp(as.integer(path), as.integer(pos))   # boundaries in C++
  data.frame(
    source = source, donor = donor, name = name, chr = as.integer(chr),
    start_bp = r$start_bp, end_bp = r$end_bp, state = r$state,
    stringsAsFactors = FALSE
  )
}

# Batched fixed-means count calls: per chromosome, pivot the (rectangular)
# cohort to T x S count matrices, memoize the emission over distinct (n,a)
# pairs, and run one C++ viterbi_batch over all samples. Returns NULL if the
# data is not a clean rectangular grid (every sample x every position once),
# so call_ancestry can fall back to the safe per-sample loop.
.batched_states <- function(data, emission, td, theta, source, donor, has_donor,
                            parallel = FALSE) {
  viterbi_fun <- if (parallel) viterbi_batch_par_cpp else viterbi_batch_cpp
  by_chr <- split(seq_len(nrow(data)), data$chr)   # row indices per chr, O(n)
  out <- vector("list", length(by_chr))
  for (i in seq_along(by_chr)) {
    ri    <- by_chr[[i]]
    nm_v  <- data$name[ri]; pos_v <- data$pos[ri]
    samples <- unique(nm_v); pos_lv <- sort(unique(pos_v))
    Tn <- length(pos_lv); S <- length(samples)
    si <- match(nm_v, samples); pii <- match(pos_v, pos_lv)
    lin <- (si - 1L) * Tn + pii                    # column-major cell id
    if (length(ri) != Tn * S || anyDuplicated(lin)) return(NULL)   # not a rectangular grid
    n_mat <- integer(Tn * S); n_mat[lin] <- as.integer(data$n_ref[ri] + data$n_alt[ri]); dim(n_mat) <- c(Tn, S)
    a_mat <- integer(Tn * S); a_mat[lin] <- as.integer(data$n_alt[ri]);                   dim(a_mat) <- c(Tn, S)
    base  <- max(n_mat) + 1
    key   <- as.double(n_mat) * base + as.double(a_mat)
    u     <- unique(key)
    em_u  <- count_emission_loglik_cpp(as.integer(u %/% base), as.integer(u %% base),
                                       theta, emission$conc)
    inv   <- matrix(match(key, u) - 1L, Tn, S)
    paths <- viterbi_fun(td$log_start, td$log_trans, em_u, inv)         # Tn x S (positions x samples)
    # coordinate-free: one state per original observation row (align via pos x sample)
    dn_v <- if (has_donor) data$donor[ri] else rep(donor, length(ri))
    out[[i]] <- data.frame(
      source = source, donor = dn_v, name = nm_v,
      chr = as.integer(names(by_chr)[i]), pos = as.integer(pos_v),
      state = as.integer(paths[cbind(pii, si)]),
      stringsAsFactors = FALSE)
  }
  states <- do.call(rbind, out)
  states[order(states$donor, states$name, states$chr, states$pos), , drop = FALSE]
}

#' Fit HMM emission/transition parameters
#'
#' For a fixed-mean emission (the default, reproducing the Python count caller)
#' this resolves `theta` and the transition and is otherwise a no-op. For a
#' fittable-mean emission (`fit_means = TRUE`; required for RNA / reference-biased
#' data per S10) it Baum-Welch EM-fits the emission means. Viterbi/EM hot loops
#' are in Rcpp.
#'
#' @param obs A per-(sample, chromosome) observation table from an emission's
#'   reader (for `count`: columns `n`, `a`).
#' @param emission Emission spec ([emission_count()] / [emission_gt()]).
#' @param duration Duration spec ([duration_geometric()] /
#'   [duration_rigidity()] / [duration_hsmm()]).
#' @param priors Single-locus genotype-frequency priors `list(f_1, f_2)`.
#' @param control EM control (`max_iter`, `tol`).
#' @return A fitted model: `list(theta, log_start, log_trans, emission, duration)`.
#' @examples
#' # Per-(sample, chromosome) count obs: n = depth, a = alt count.
#' obs <- data.frame(n = c(10, 9, 11, 8, 12), a = c(0, 0, 1, 4, 6))
#' model <- fit(obs, emission_count(), duration_geometric(1e-4),
#'              priors = design_priors("BC2S2"))
#' @export
fit <- function(obs, emission, duration, priors, control = list()) {
  td <- .duration_transition(duration, priors)
  theta <- .emission_theta(emission)                 # initial / fixed means
  if (isTRUE(emission$fit_means)) {
    td_fit <- if (td$n_sub == 1L) td                 # fit means on collapsed 3-state model
              else .build_transition(duration$p_switch, priors$f_1, priors$f_2)
    theta <- .em_fit_means(list(obs), emission, td_fit, theta, control)
  }
  structure(list(theta = theta, log_start = td$log_start, log_trans = td$log_trans,
                 n_sub = td$n_sub, emission = emission, duration = duration),
            class = "nilHMM_model")
}

# Replicate a T x 3 macro emission across each macro-state's n_sub sub-states
# (expanded column index = m * n_sub + k). n_sub = 1 is a no-op.
.expand_emission <- function(em3, n_sub) {
  if (n_sub == 1L) return(em3)
  out <- matrix(0, nrow(em3), 3L * n_sub)
  for (m in 0:2) for (k in 0:(n_sub - 1L)) out[, m * n_sub + k + 1L] <- em3[, m + 1L]
  out
}

#' Decode the most-likely state path (Viterbi)
#'
#' @param model A fitted model from [fit()].
#' @param obs A per-(sample, chromosome) observation table.
#' @return Integer macro-state path over markers (0 = REF, 1 = HET, 2 = ALT);
#'   for rigidity the expanded sub-states are mapped back to macro-states.
#' @examples
#' obs <- data.frame(n = c(10, 9, 11, 8, 12), a = c(0, 0, 1, 4, 6))
#' model <- fit(obs, emission_count(), duration_geometric(1e-4),
#'              priors = design_priors("BC2S2"))
#' decode(model, obs)
#' @export
decode <- function(model, obs) {
  em <- .expand_emission(.emission_loglik(model$emission, obs, model$theta), model$n_sub)
  path <- viterbi_log_cpp(model$log_start, model$log_trans, em)
  if (model$n_sub > 1L) path %/% model$n_sub else path
}

#' Coordinate-free ancestry decode (per-observation states)
#'
#' Resolves a named caller (and optional design preset) into emission + duration
#' specs and runs [fit()] / [decode()] per (sample, chromosome), returning one
#' ancestry **state per observation unit** --- coordinate-free: no genomic
#' intervals. Feed the result to [to_segments()] to reinstate coordinates and
#' collapse into common-schema segments; [call_ancestry()] chains the two.
#' Data-agnostic: pass the observations in; this function never reads paths.
#'
#' @param data Long observation table with columns `name, chr, pos` and either
#'   `n_ref, n_alt` read counts (from [read_counts()]; for the count/rtiger/binhmm/
#'   atlas paths) or a pre-called hard-genotype column `g` in `{0,1,2,3}` (from
#'   [read_vcf_gt()]; auto-selects `caller = "nnil"`, `emission = "gt"`).
#'   Optionally `donor`.
#' @param caller One of `"nnil"`, `"rtiger"`, `"binhmm"`, `"atlas"`, `"lbimpute"`.
#' @param design Breeding-design key for priors (e.g. `"BC2S2"`, `"BC2S3"`).
#'   Required unless `f_1`/`f_2` are supplied.
#' @param rrate Count/geometric callers (`nnil`, `atlas`): expected per-marker
#'   **recombination rate** for the geometric duration (self-stay = `1 - rrate`).
#'   A resolution hyperparameter, not an MLE. Holland's nNIL sets it to
#'   `2 * total_cM / (100 * n_markers)` (~`30 / n_markers` for a 1500 cM maize
#'   map; the factor 2 is the RIL-like doubling for the backcross + self meioses,
#'   Haldane & Waddington), and it is optimizable from data by KS-vs-sim
#'   ([calibrate_r()]).
#' @param rigidity `rtiger` caller only: integer minimum run length (e.g. `5`).
#' @param xrate Exit rate of nilHMM's **rigidity duration** ([duration_rigidity()]):
#'   the per-marker switch probability at the free (post-minimum-run) state — the
#'   geometric tail beyond the enforced run. A nilHMM construct, **not** a RTIGER
#'   parameter; the faithful `caller = "rtiger"` port uses only `rigidity` (plus
#'   `threads`/`seed`/`postprocess`) and ignores `xrate`.
#' @param err,conc,fit_means Count-emission parameters forwarded to [caller_spec()].
#' @param f_1,f_2 Single-locus priors, used when `design` is `NULL`.
#' @param source Value for the output `source` column.
#' @param donor Donor/taxon label when `data` has no `donor` column.
#' @param parallel If `TRUE`, decode samples across cores via RcppParallel
#'   (only the fixed-means batched path; control thread count with
#'   [RcppParallel::setThreadOptions()]). Identical results to serial.
#' @param threads,seed RTIGER caller only: E-step threads and the seed for its
#'   randomized init.
#' @param postprocess RTIGER caller only: apply the border re-placement (default TRUE).
#' @param min_cov Drop no-coverage units before decoding (default `1L`; `0L` keeps
#'   everything, the old behaviour). "No coverage" is caller-specific but the
#'   intent is uniform — never make a confident call from no data, and keep all
#'   callers on the same support for comparability:
#'   * `nnil` (count) and `rtiger`: markers with `n_ref + n_alt < min_cov`;
#'   * `nnil` gt path: missing genotypes (`g == 3`);
#'   * `binhmm`: bins with fewer than `min_cov` informative markers (`ninf`).
#'   Zero-coverage units carry no emission signal — they only slow decoding,
#'   marginally inflate `nnil` fragmentation, and dilute `rtiger`'s rigidity run.
#'   (The `atlas` caller has its own `atlas_min_reads` gate and is unaffected.)
#' @param emission Optional emission override (`"count"`, `"gt"`) for the `nnil`
#'   caller; `NULL` uses the caller's default.
#' @param bin_size,cluster_method `binhmm` caller only: genomic bin width in bp
#'   (default 1 Mb) and the per-bin genotyping backend --- `"gauss"` (default: the
#'   anchored 3-state Gaussian-emission HMM; fixes the HET over-call and
#'   high-coverage fragmentation of the K=3 clustering), or the original
#'   cluster-then-smooth route `"gmm"` / `"kmeans"` / `"rebmix"` (the last needs
#'   the suggested \pkg{rebmix} package, for bit-exact rpubs reproduction).
#'   `joint_clust`/`obs_weights` apply only to the cluster backends.
#' @param joint_clust `binhmm` caller only: if `TRUE`, pool all samples' bins and
#'   learn one shared set of REF/HET/ALT clusters on raw alt-freq (borrows
#'   strength across the cohort; cf. `get_joint_ancestry_calls.R`); if `FALSE`
#'   (default), cluster each sample independently on its informative-count-weighted
#'   alt-freq. This is joint *clustering* only; the HMM stays per-sample.
#' @param obs_weights `binhmm` caller only, `joint_clust = TRUE` only: if `TRUE`,
#'   weight the (gmm) clustering fit by each bin's informative-variant count
#'   (weights the influence, not the alt-freq value).
#' @param genotypeerr,recombdist,drp,unit `lbimpute` caller only (Fragoso et al.
#'   2014): `genotypeerr` is the coverage-independent genotyping-error rate
#'   (LB-Impute `genotypeerr`, default 0.05) that bounds every emission to
#'   `[genotypeerr, 1 - genotypeerr]`; `err` doubles as LB-Impute's per-read
#'   `readerr` (its default is 0.05, higher than the shared `err = 0.01`).
#'   `unit` chooses the coordinate the transition decays over: `"bp"` (default,
#'   the faithful LB-Impute model -- a uniform genome-wide recombination rate
#'   over physical distance) or `"cm"` (map-aware -- uses a `cm` column of genetic
#'   positions per marker, so the local recombination rate, e.g. maize
#'   centromeric suppression, is captured). Output coordinates are always bp
#'   (`pos`); cM only feeds the transition. `recombdist` is the coordinate
#'   distance over which the recombination probability equalizes (LB-Impute
#'   `recombdist`); it shares units with the coordinate, so its default is
#'   unit-aware (`1e7` bp = ~50 cM in maize, or `50` cM) and larger values mean
#'   stiffer paths -- a validation layer warns on a bp/cM `recombdist` mismatch.
#'   `drp` (LB-Impute `-dr`): when `TRUE`, a homozygous->homozygous switch is
#'   priced as a single recombination rather than a double event (use for inbred
#'   / RIL populations). For `lbimpute`, `design`/`f_1,f_2` only seed the start
#'   distribution (flat if absent) and zero-coverage markers are kept so the
#'   transition sees true marker spacing.
#' @param germ,gert,p,mr,nir Genotype-error rates for the `gt` (categorical)
#'   emission (Holland's nNIL model): `germ` error on true homozygotes, `gert` on
#'   true heterozygotes, `p` fraction of hom errors called het, `mr` missing rate,
#'   `nir` non-informative-marker rate. Raising `germ`/`gert` toward the platform's
#'   real error rate is the native cure for over-fragmentation (isolated
#'   miscalled markers are absorbed as errors rather than opening 1-marker
#'   segments); calibrate to a clean control, don't crank blindly. Used by the gt
#'   emission (`emission = "gt"`, a `g` genotype input, or the `atlas` caller).
#' @param atlas_thresh,atlas_het,atlas_min_reads `atlas` caller only: GOOGA
#'   genotype-call thresholds on the donor read fraction --- homozygous call when
#'   a parent's fraction >= `atlas_thresh` (0.95), HET when both parents >=
#'   `atlas_het` (0.25), and a minimum of `atlas_min_reads` (5) informative reads
#'   per gene (else missing). For `atlas`, `n_ref`/`n_alt` are the recurrent/donor
#'   competitive-alignment read counts (ambiguous excluded upstream).
#' @return data.frame of per-unit states
#'   (`source, donor, name, chr, pos, state`; the `binhmm` caller additionally
#'   carries per-bin `start_bp, end_bp`). Pass to [to_segments()] for intervals.
#' @examples
#' # Toy single-NIL count table: a REF stretch then a donor (ALT) block.
#' set.seed(1)
#' toy <- data.frame(
#'   name  = "NIL1", chr = 1L, pos = seq_len(40L) * 1e5L,
#'   n_ref = c(rpois(20, 8), rpois(20, 4)),
#'   n_alt = c(rpois(20, 0), rpois(20, 4)))
#' st <- call_states(toy, caller = "nnil", design = "BC2S2", rrate = 1e-4, err = 0.01)
#' to_segments(st)                       # or, in one step:
#' call_ancestry(toy, caller = "nnil", design = "BC2S2", rrate = 1e-4, err = 0.01)
#' @export
call_states <- function(data, caller = c("nnil", "rtiger", "binhmm", "atlas", "lbimpute"),
                        design = NULL, rrate = 0.01, rigidity = NULL,
                        err = 0.01, conc = 20,
                        fit_means = FALSE, xrate = 0.01,
                        f_1 = NULL, f_2 = NULL,
                        source = "nilHMM", donor = NA_character_,
                        parallel = FALSE, threads = 1L, seed = 1L,
                        postprocess = TRUE, min_cov = 1L, emission = NULL,
                        bin_size = 1e6, cluster_method = c("gauss", "gmm", "kmeans", "rebmix"),
                        joint_clust = FALSE, obs_weights = FALSE,
                        atlas_thresh = 0.95, atlas_het = 0.25, atlas_min_reads = 5L,
                        genotypeerr = 0.05, recombdist = NULL, drp = FALSE,
                        unit = c("bp", "cm"),
                        germ = 0.05, gert = 0.10, p = 0.5, mr = 0.10, nir = 0.01) {
  caller <- match.arg(caller)
  has_counts <- all(c("n_ref", "n_alt") %in% names(data))
  has_gt     <- "g" %in% names(data)          # pre-called hard genotype (0/1/2/3), e.g. from read_vcf_gt()
  has_af     <- "alt_freq" %in% names(data)   # pre-binned binhmm input (bin summary done upstream)
  if (!all(c("name", "chr", "pos") %in% names(data)))
    stop("call_states(): data needs columns name, chr, pos")
  if (caller == "binhmm") {
    if (!(has_counts || has_af))
      stop("call_states(): caller = 'binhmm' needs (n_ref, n_alt) read counts to bin, ",
           "or a pre-binned `alt_freq` column (with start_bp, end_bp)")
  } else if (!(has_counts || has_gt)) {
    stop("call_states(): caller = '", caller, "' needs (n_ref, n_alt) read counts ",
         "or a `g` genotype column")
  }
  has_donor <- "donor" %in% names(data)

  # A hard-genotype (`g`-only, no counts) input is the categorical GT path
  # (Holland's nNIL genotype model on called genotypes; the saturated-depth /
  # MolBreeding regime). It only makes sense for caller = "nnil" + the gt emission
  # (the count/rtiger/binhmm/atlas paths all need read counts).
  if (has_gt && !has_counts) {
    if (caller != "nnil")
      stop("call_states(): a `g` genotype input (no read counts) requires caller = 'nnil'; ",
           "'", caller, "' needs n_ref/n_alt read counts")
    emission <- "gt"
  }

  # Required priors are a caller-argument problem, not a data problem, so validate
  # them BEFORE any coverage filtering below. Otherwise an all-zero-coverage input
  # would trip the min_cov "no covered markers" error and mask a genuinely missing
  # `design`/`f_1`/`f_2`. Every caller but rtiger (which fits its own start freqs)
  # needs priors; each path re-resolves them below, this just fails fast up front.
  if (!caller %in% c("rtiger", "lbimpute") &&
      is.null(design) && !(!is.null(f_1) && !is.null(f_2)))
    stop("call_states(): supply `design` or both `f_1` and `f_2`")

  # gt path: "no coverage" is a missing genotype (g == 3). Drop missing calls so
  # the gt (nnil hard-genotype) caller ignores uncovered positions too, mirroring
  # the count/rtiger covered-marker filter below.
  if (identical(emission, "gt") && "g" %in% names(data) &&
      !is.null(min_cov) && min_cov > 0L) {
    data <- data[data$g != 3L, , drop = FALSE]
    if (!nrow(data)) stop("call_states(): no non-missing genotypes after the min_cov filter")
  }

  # Covered-marker filter for the count-emission callers (`nnil` count path and
  # `rtiger`): drop markers with fewer than `min_cov` total reads. A zero-coverage
  # panel position carries no emission signal — it only slows decoding, marginally
  # inflates fragmentation for `nnil`, and dilutes the rigidity run for `rtiger`
  # (which requires covered markers anyway). This makes the two callers run on the
  # SAME marker support, so their calls are directly comparable. `binhmm` bins its
  # own way and the `gt`/`alt_freq` inputs have no read counts to threshold, so
  # they are exempt. `min_cov = 0L` restores decoding every input marker.
  if (has_counts && !identical(emission, "gt") && caller %in% c("nnil", "rtiger") &&
      !is.null(min_cov) && min_cov > 0L) {
    data <- data[data$n_ref + data$n_alt >= min_cov, , drop = FALSE]
    if (!nrow(data))
      stop("call_states(): no markers with coverage >= min_cov (", min_cov, ")")
  }

  # RTIGER caller: its own EM/Viterbi (src/rtiger.cpp, R/rtiger.R) — a faithful
  # port of the RTIGER fork, not the count engine. `r` is the integer rigidity;
  # `postprocess` applies the border re-placement (on by default, as RTIGER does).
  if (caller == "rtiger") {
    rig <- if (is.null(rigidity)) 5L else as.integer(rigidity)   # minimum run length
    if (rig < 1L) stop("rtiger: `rigidity` must be an integer >= 1")
    return(.rtiger_states(data, rig, source, donor, has_donor, threads, seed, postprocess))
  }

  # binhmm caller: the bin -> K=3 cluster -> HMM-smooth pipeline (R/binhmm.R),
  # architecturally distinct from the count engine (it clusters per-bin ALT
  # fraction rather than modelling reads). Only priors (start freqs) feed it.
  if (caller == "binhmm") {
    cluster_method <- match.arg(cluster_method)
    priors <- if (!is.null(design)) design_priors(design)
              else if (!is.null(f_1) && !is.null(f_2)) list(f_1 = f_1, f_2 = f_2)
              else stop("call_states(): supply `design` or both `f_1` and `f_2`")
    return(.binhmm_states(data, bin_size, cluster_method, priors,
                          source, donor, has_donor,
                          joint_clust = joint_clust, obs_weights = obs_weights,
                          min_cov = min_cov))
  }

  # LB-Impute caller (R/lbimpute.R, src/lbimpute.cpp): a native port of
  # Fragoso et al. 2014. Coverage-aware emission + distance-dependent transition
  # (double-recomb penalty), decoded with a full-chromosome Viterbi. Routed
  # separately because its transition is per-marker (distance-based), not the
  # single time-homogeneous matrix the count/gt callers share. `design`/`f_1,f_2`
  # only seed the start distribution (LB-Impute has no state-frequency prior in
  # its transition); absent, the start is flat. Zero-coverage markers are kept
  # (flat emission) so the distance transition sees true marker spacing.
  #
  # The transition decays over a coordinate: physical bp (`unit = "bp"`, the
  # faithful LB-Impute model, a uniform genome-wide recombination rate) or
  # genetic cM (`unit = "cm"`, map-aware, so local recombination rate -- e.g.
  # maize centromeric suppression -- is captured). Output coordinates are always
  # bp (`pos`); cM only feeds the transition. `recombdist` and the coordinate
  # must share units, so its default is unit-aware and a validation layer warns
  # on the dangerous unit/`recombdist` mismatches (see the design discussion).
  if (caller == "lbimpute") {
    unit <- match.arg(unit)
    if (!has_counts)
      stop("call_states(): caller = 'lbimpute' needs (n_ref, n_alt) read counts")
    # Transition coordinate (bp `pos` / cM `cm`) + value-based validation, then the
    # unit-aware `recombdist` default and mismatch warnings, then the start seed --
    # all shared with caller_sweep via helpers in R/lbimpute.R.
    tcol <- .lbimpute_tcol(data, unit)
    if (is.null(recombdist)) recombdist <- if (unit == "cm") 50 else 1e7
    if (recombdist <= 0) stop("call_states(): `recombdist` must be > 0")
    .lbimpute_warn_recombdist(recombdist, unit, data[[tcol]])
    log_init <- .lbimpute_log_init(design, f_1, f_2)
    return(.lbimpute_states(data, err, genotypeerr, recombdist, drp, log_init,
                            source, donor, has_donor, tcol, threads))
  }

  # ATLAS caller (R/atlas.R): GOOGA-style per-gene ancestry from competitive-
  # alignment read counts (n_ref = recurrent, n_alt = donor; ambiguous excluded
  # upstream by the cassini pipeline). Same engine as nnil+gt, but the per-unit
  # genotype call uses GOOGA's hard fraction thresholds + a min-read gate rather
  # than the 1/3-2/3 cutoffs. It is RNA/read-competition data, so the categorical
  # gt (confusion) emission is used, NOT the count/BetaBinomial mean model.
  googa <- caller == "atlas"

  spec <- caller_spec(caller, rrate = rrate, err = err, conc = conc, fit_means = fit_means,
                      xrate = xrate, germ = germ, gert = gert, p = p, mr = mr, nir = nir)
  if (!is.null(emission)) spec$emission <- switch(match.arg(emission, c("count","gt")),
    count  = emission_count(err, conc, fit_means),
    gt     = emission_gt(germ, gert, p, mr, nir))

  priors <- if (!is.null(design)) design_priors(design)
            else if (!is.null(f_1) && !is.null(f_2)) list(f_1 = f_1, f_2 = f_2)
            else stop("call_states(): supply `design` or both `f_1` and `f_2`")

  td     <- .duration_transition(spec$duration, priors)   # same for all samples
  theta0 <- .emission_theta(spec$emission)

  # Fast path: fixed-means count caller on a rectangular cohort -> one batched
  # Viterbi per chromosome across all samples in C++ (mirrors the numpy-batched
  # Python). Falls back to the per-sample loop when fit_means / rigidity / a
  # non-rectangular grid make batching unsafe.
  if (!isTRUE(spec$emission$fit_means) && td$n_sub == 1L &&
      inherits(spec$emission, "nilHMM_emission_count")) {
    batched <- .batched_states(data, spec$emission, td, theta0, source, donor, has_donor, parallel)
    if (!is.null(batched)) return(batched)
  }

  # Split by sample ONCE (O(n)) rather than re-scanning the frame per sample
  # (which is O(samples * rows) and dominates wall-clock on large cohorts).
  by_name <- split(data, data$name, drop = TRUE)

  # Per-sample work is independent (each sample builds its own chains, optionally
  # fits its own emission means, and decodes), so parallelize across samples when
  # threads > 1 on unix. This is the win for the fit_means path (whose per-sample
  # EM dominates) and for any non-rectangular cohort that misses the batched fast
  # path. Each call returns this sample's stacked per-chromosome state rows.
  one_sample <- function(nm) {
    dn <- by_name[[nm]]
    donor_nm <- if (has_donor) dn$donor[1] else donor
    # per-chromosome observation sequences (HMM decodes each chr independently;
    # emission means are fit pooled across a sample's chromosomes).
    obs_list <- lapply(split(dn, dn$chr, drop = TRUE), function(dc) {
      dc <- dc[order(dc$pos), , drop = FALSE]
      if (has_gt && !has_counts) {                 # pre-called genotypes (read_vcf_gt): use g directly
        n <- rep(NA_integer_, nrow(dc)); a <- n; g <- as.integer(dc$g)
      } else {
        n <- dc$n_ref + dc$n_alt; a <- dc$n_alt; f <- ifelse(n == 0, NA_real_, a / n)
        # derived hard genotype call for the gt emission: ATLAS uses GOOGA's
        # fraction thresholds + min-read gate; otherwise the 1/3-2/3 cutoffs.
        g <- if (googa) .googa_gt_call(a, n, atlas_thresh, atlas_het, atlas_min_reads)
             else ifelse(is.na(f), 3L, ifelse(f < 1/3, 0L, ifelse(f > 2/3, 2L, 1L)))
      }
      list(chr = dc$chr[1], pos = dc$pos, n = n, a = a, g = g)
    })
    # Emission means are fit on the COLLAPSED 3-state model (means are
    # ~duration-independent; the expanded rigidity FB is K^2 = (3r)^2 and far
    # slower). Decoding still uses the full duration transition `td`.
    theta <- theta0
    if (isTRUE(spec$emission$fit_means)) {
      td_fit <- if (td$n_sub == 1L) td
                else .build_transition(spec$duration$p_switch, priors$f_1, priors$f_2)
      theta <- .em_fit_means(obs_list, spec$emission, td_fit, theta0)
    }
    model <- structure(list(theta = theta, log_start = td$log_start,
                            log_trans = td$log_trans, n_sub = td$n_sub,
                            emission = spec$emission, duration = spec$duration),
                       class = "nilHMM_model")
    do.call(rbind, lapply(obs_list, function(o) data.frame(
      source = source, donor = donor_nm, name = nm, chr = as.integer(o$chr),
      pos = as.integer(o$pos), state = as.integer(decode(model, o)),
      stringsAsFactors = FALSE)))
  }

  nms <- names(by_name)
  per <- if (threads > 1L && .Platform$OS.type == "unix")
           parallel::mclapply(nms, one_sample, mc.cores = threads)
         else lapply(nms, one_sample)
  bad <- vapply(per, function(x) inherits(x, "try-error") || is.null(x), logical(1))
  if (any(bad)) stop("call_states(): decoding failed for sample(s): ",
                     paste(nms[bad], collapse = ", "))
  # rbindlist avoids the O(n^2)-ish copying of do.call(rbind, .) over many
  # per-sample frames — it's the serial tail that caps the mclapply speedup.
  states <- if (requireNamespace("data.table", quietly = TRUE))
              as.data.frame(data.table::rbindlist(per)) else do.call(rbind, per)
  states[order(states$donor, states$name, states$chr, states$pos), , drop = FALSE]
}

#' Collapse per-unit ancestry states into genomic segments
#'
#' The coordinate-reinstatement step: run-length-collapses the equal-state runs
#' within each `(name, chr)` from [call_states()] into segments. **Point** units
#' (markers) use `pos` for both boundaries; **interval** units (bins) carry their
#' own per-unit `start_bp`/`end_bp` (as the `binhmm` caller emits). This is the
#' only place genomic coordinates re-enter — the decode itself is coordinate-free.
#'
#' @param states A per-unit state table from [call_states()]: `name, chr, pos,
#'   state` (+ optional `source, donor`; + optional `start_bp, end_bp` for
#'   interval/bin units).
#' @return data.frame in the common schema
#'   (`source, donor, name, chr, start_bp, end_bp, state`).
#' @examples
#' st <- data.frame(name = "NIL1", chr = 1L, pos = (1:6) * 1e5L,
#'                  state = c(0L, 0L, 2L, 2L, 0L, 0L))
#' to_segments(st)
#' @export
to_segments <- function(states) {
  if (is.null(states) || !nrow(states)) {
    return(states)
  }
  st <- as.data.frame(states, stringsAsFactors = FALSE)
  if (!all(c("name", "chr", "pos", "state") %in% names(st))) {
    stop("to_segments(): `states` needs columns name, chr, pos, state")
  }
  if (!"source" %in% names(st)) st$source <- "nilHMM"
  if (!"donor" %in% names(st)) st$donor <- NA_character_
  interval <- all(c("start_bp", "end_bp") %in% names(st))
  # Single vectorized RLE across all (name, chr) groups: order by (group, pos),
  # then a run starts wherever the group or the state changes. Avoids the per-group
  # loop + do.call(rbind) that dominated wall-clock on large cohorts (the RLE and
  # segment build are then O(n) over the whole table).
  gkey <- paste(st$name, st$chr, sep = "\r")
  o <- order(gkey, st$pos)
  st <- st[o, , drop = FALSE]
  gkey <- gkey[o]
  state <- as.integer(st$state)
  n <- length(state)
  s <- which(c(TRUE, state[-1L] != state[-n] | gkey[-1L] != gkey[-n]))  # run-start rows
  e <- c(s[-1L] - 1L, n)                                                # run-end rows
  calls <- data.frame(
    source   = st$source[s], donor = st$donor[s], name = st$name[s],
    chr      = as.integer(st$chr[s]),
    start_bp = as.integer(if (interval) st$start_bp[s] else st$pos[s]),
    end_bp   = as.integer(if (interval) st$end_bp[e] else st$pos[e]),
    state    = state[s], stringsAsFactors = FALSE
  )
  calls[order(calls$donor, calls$name, calls$chr, calls$start_bp), , drop = FALSE]
}

#' Top-level ancestry-calling API (decode + segment)
#'
#' Convenience wrapper that chains the two stages: coordinate-free decode
#' ([call_states()]) followed by coordinate reinstatement ([to_segments()]).
#' `call_ancestry(...)` is exactly `to_segments(call_states(...))`. See
#' [call_states()] for the full argument list.
#'
#' @param data,... Passed to [call_states()].
#' @return data.frame in the common schema
#'   (`source, donor, name, chr, start_bp, end_bp, state`).
#' @examples
#' toy <- data.frame(name = "NIL1", chr = 1L, pos = (1:6) * 1e5L,
#'                   n_ref = c(9, 8, 4, 5, 9, 8), n_alt = c(0, 0, 4, 5, 0, 0))
#' call_ancestry(toy, caller = "nnil", design = "BC2S2", rrate = 1e-4, err = 0.01)
#' @export
call_ancestry <- function(data, ...) {
  to_segments(call_states(data, ...))
}
