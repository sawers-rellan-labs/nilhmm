# googa / atlas callers -- per-gene ancestry from COMPETITIVE ALIGNMENT read
# counts to the two parents (recurrent vs donor transcriptomes; the cassini
# pipeline produces these). NOT a SNP/dosage caller and NOT the count/BetaBinomial
# mean model -- RNA read counts are expression-driven and the allele fraction is
# confounded by allele-specific expression, so a hard threshold + categorical (gt)
# confusion HMM is used.
#
# Two callers share this GOOGA thresholding, differing only in the duration:
#   `googa` -- gt + GEOMETRIC. The faithful reproduction of GOOGA (Flagel et al.
#     2019, pcbi.1006949) and its ASE application (Veltsos et al. 2024,
#     pgen.1011072), whose HMM is a recombination-fraction (geometric) F2 model
#     with NO minimum-run/rigidity duration (confirmed in the GOOGA source and both
#     papers' methods).
#   `atlas` -- gt + RIGIDITY. This work's transcript caller: the same GOOGA
#     thresholding decoded with the engine's rigidity (minimum-run-length) duration.
#
# Input to call_ancestry(caller = "atlas"): the common count schema
# `name, chr, pos, n_ref, n_alt`, where per gene:
#   n_ref = reads assigned to the RECURRENT parent (e.g. B73),
#   n_alt = reads assigned to the DONOR parent.
# AMBIGUOUS/tie reads (map equally to both parents = conserved, non-informative)
# are EXCLUDED upstream (not a state, not HET). Genes are positioned by their B73
# gene coordinate (a consumer-side staging concern; the package stays data-agnostic).
#
# GOOGA genotype call (Program_Set_2/s2.py), from donor fraction f = n_alt / ntot
# with ntot = n_ref + n_alt (thresh = 0.95, het = 0.25, min_reads = 5):
#   ntot < min_reads              -> NN / missing (3)
#   recurrent frac (1-f) >= thresh -> REF (0)
#   donor frac (f)      >= thresh -> ALT (2)
#   both fracs          >= het    -> HET (1)
#   otherwise (ambiguous)         -> NN / missing (3)
# The resulting per-gene calls feed the engine's gt (categorical confusion)
# emission HMM: with the geometric duration for `googa` (faithful GOOGA) or the
# rigidity duration for `atlas` (this work). See caller_spec()/call_ancestry.

# vectorized GOOGA hard call -> {0 REF, 1 HET, 2 ALT, 3 missing}
.googa_gt_call <- function(a, n, thresh = 0.95, het = 0.25, min_reads = 5L) {
  f <- ifelse(n == 0, NA_real_, a / n)
  g <- rep(3L, length(a))                                 # NN / missing default
  ok <- n >= min_reads & !is.na(f)
  g[ok & f <= 1 - thresh]            <- 0L                # recurrent >= thresh -> REF
  g[ok & f >= thresh]                <- 2L                # donor >= thresh -> ALT
  g[ok & f >= het & f <= 1 - het]    <- 1L                # both parents >= het -> HET
  g                                                       # ambiguous ok-bins stay 3 (missing)
}
