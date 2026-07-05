// LB-Impute caller (Fragoso et al. 2014, G3; dellaporta-laboratory/LB-Impute).
// A native nilHMM port: the coverage-aware emission and the distance-dependent
// transition (with the double-recombination penalty) are LB-Impute's model,
// verbatim from the reference Java (ImputeOffspring.getprobabilities2 and
// FindPath2). What DIFFERS from the original is the decoder: LB-Impute uses an
// iterating fixed-window best-path plus a forward/reverse consensus -- a
// computational compromise the paper itself flags ("it would be ideal for this
// chain to stretch the entire length of the chromosome"). nilHMM does that
// ideal: one full-chromosome log-space Viterbi. So the emission/transition are
// faithful; the decode is the optimal path the window approximates.
//
// State order is the common schema REF/HET/ALT = 0/1/2. LB-Impute internally
// orders [parent1, parent2, het] with het LAST; the two homozygous "parent"
// states map to REF (ref allele) and ALT (alt allele) here.

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' LB-Impute coverage-aware emission (log-space, REF/HET/ALT)
//'
//' Per-marker emission from allelic read depths, the model of
//' `ImputeOffspring.getprobabilities2`. Raw state likelihoods are
//' \eqn{E_{REF}\propto(1-err)^{n_{ref}}err^{n_{alt}}},
//' \eqn{E_{ALT}\propto(1-err)^{n_{alt}}err^{n_{ref}}},
//' \eqn{E_{HET}\propto 0.5^{n_{ref}+n_{alt}}}; each is divided by the per-marker
//' maximum, scaled by \eqn{1-2\,err_g} and offset by \eqn{err_g}, bounding every
//' emission to \eqn{[err_g,\,1-err_g]} so a single artifactual marker cannot
//' dominate the path. A zero-coverage marker emits flat (all states equal).
//'
//' @param nref,nalt Integer per-marker reference / alternate read counts.
//' @param err Per-read sequencing-error probability (LB-Impute `readerr`).
//' @param errg Coverage-independent genotyping-error probability (LB-Impute
//'   `genotypeerr`); sets the emission floor/ceiling.
//' @return A T x 3 matrix of log emission probabilities (columns REF/HET/ALT).
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix lb_emission_loglik_cpp(IntegerVector nref, IntegerVector nalt,
                                     double err, double errg) {
  const int T = nref.size();
  if (nalt.size() != T) stop("lb_emission_loglik_cpp: nref/nalt length mismatch");
  if (err <= 0.0 || err >= 0.5) stop("lb_emission_loglik_cpp: err must be in (0, 0.5)");
  if (errg < 0.0 || errg >= 0.5) stop("lb_emission_loglik_cpp: errg must be in [0, 0.5)");
  NumericMatrix out(T, 3);
  const double log_ok  = std::log(1.0 - err);
  const double log_bad = std::log(err);
  const double log_half = std::log(0.5);
  const double maxprob = 1.0 - 2.0 * errg;   // emission range scale
  const double minprob = errg;               // emission floor
  for (int t = 0; t < T; ++t) {
    const double r = nref[t], a = nalt[t];
    // raw log-likelihoods (computed in log space so high depth cannot underflow)
    const double lr = r * log_ok  + a * log_bad;   // REF hom
    const double la = a * log_ok  + r * log_bad;   // ALT hom
    const double lh = (r + a) * log_half;          // HET
    double m = lr; if (la > m) m = la; if (lh > m) m = lh;
    // normvalue = (raw/max) * maxprob + minprob  in [errg, 1-errg]  -> log
    out(t, 0) = std::log(std::exp(lr - m) * maxprob + minprob);
    out(t, 1) = std::log(std::exp(lh - m) * maxprob + minprob);
    out(t, 2) = std::log(std::exp(la - m) * maxprob + minprob);
  }
  return out;
}

//' LB-Impute distance-aware full-chromosome Viterbi
//'
//' Decodes the most-likely REF/HET/ALT path under LB-Impute's distance-dependent
//' transition (`FindPath2`). Between markers a physical distance `d` bp apart,
//' the stay probability is \eqn{p_s = 0.5(1 + e^{-d/recombdist})} and the single
//' recombination probability is \eqn{p_r = 0.5(1 - e^{-d/recombdist})}. Homozygous
//' <-> heterozygous transitions cost one recombination (\eqn{p_r}); homozygous
//' -> the OTHER homozygous state costs two (\eqn{p_r^2}) unless `drp = TRUE`
//' (LB-Impute `-dr`), which prices it as a single event (for inbred / RIL
//' populations). Transition weights are LB-Impute's exact (un-normalized) model;
//' the Viterbi argmax over full-chromosome paths supersedes the original's
//' windowed best-path + forward/reverse consensus.
//'
//' @param log_init Length-3 vector of log initial-state probabilities (REF/HET/ALT).
//' @param log_emit T x 3 matrix of log emissions (from [lb_emission_loglik_cpp()]).
//' @param pos Integer length-T marker positions in bp (sorted ascending).
//' @param recombdist Distance in bp over which recombination probability
//'   equalizes (LB-Impute `recombdist`, default 1e7).
//' @param drp If `TRUE`, a homozygous->homozygous switch is priced as a single
//'   recombination rather than a double event.
//' @return Integer length-T most-likely state path (0 = REF, 1 = HET, 2 = ALT).
//' @keywords internal
// [[Rcpp::export]]
IntegerVector lb_viterbi_cpp(NumericVector log_init, NumericMatrix log_emit,
                             IntegerVector pos, double recombdist, bool drp) {
  const int T = log_emit.nrow();
  if (log_emit.ncol() != 3) stop("lb_viterbi_cpp: log_emit must be T x 3");
  if (log_init.size() != 3) stop("lb_viterbi_cpp: log_init must have length 3");
  if (pos.size() != T) stop("lb_viterbi_cpp: length(pos) must equal nrow(log_emit)");
  if (recombdist <= 0.0) stop("lb_viterbi_cpp: recombdist must be > 0");

  IntegerVector path(T);
  if (T == 0) return path;

  NumericMatrix delta(T, 3);
  IntegerMatrix psi(T, 3);
  for (int k = 0; k < 3; ++k) { delta(0, k) = log_init[k] + log_emit(0, k); psi(0, k) = 0; }

  // per-edge log transition, rebuilt from the marker gap (0 = REF, 1 = HET, 2 = ALT).
  double lt[3][3];
  for (int t = 1; t < T; ++t) {
    double d = (double)pos[t] - (double)pos[t - 1];
    if (d < 0) stop("lb_viterbi_cpp: `pos` must be sorted ascending");
    const double e  = std::exp(-d / recombdist);
    const double ps = 0.5 * (1.0 + e);   // stay
    const double pr = 0.5 * (1.0 - e);   // single recombination
    const double lps = std::log(ps);
    const double lpr = std::log(pr);           // -Inf when d == 0 (coincident markers): forbids recomb
    const double lhh = drp ? lpr : (2.0 * lpr);// hom -> other hom (single vs double event)
    lt[0][0] = lps; lt[0][1] = lpr; lt[0][2] = lhh;   // from REF
    lt[1][0] = lpr; lt[1][1] = lps; lt[1][2] = lpr;   // from HET
    lt[2][0] = lhh; lt[2][1] = lpr; lt[2][2] = lps;   // from ALT
    for (int k = 0; k < 3; ++k) {
      double best = R_NegInf; int arg = 0;
      for (int j = 0; j < 3; ++j) {
        double v = delta(t - 1, j) + lt[j][k];
        if (v > best) { best = v; arg = j; }
      }
      delta(t, k) = best + log_emit(t, k);
      psi(t, k) = arg;
    }
  }

  double best = R_NegInf; int last = 0;
  for (int k = 0; k < 3; ++k) if (delta(T - 1, k) > best) { best = delta(T - 1, k); last = k; }
  path[T - 1] = last;
  for (int t = T - 1; t > 0; --t) path[t - 1] = psi(t, path[t]);
  return path;
}