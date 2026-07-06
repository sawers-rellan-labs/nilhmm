// Flanking-marker genotype interpolation (Tian 2011 / Chen 2019), one chromosome.
// Densifies a COMPLETE k x n genotype block onto a target cM grid by linear
// interpolation in genetic distance. The speed win: the flanking index jj[m] and
// weight w[m] depend only on the two cM vectors, so they are precomputed once via
// std::upper_bound and reused across all n samples (O(M log k + M*n) total).
// This is deterministic hard-call interpolation, NOT ancestry inference.

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

//' Interpolate a complete genotype block onto a target cM grid (single chr)
//'
//' Linear flanking-marker interpolation in genetic distance for one chromosome.
//' For a target at cM `t` between called flanking markers L (at `obs_cm[jj]`) and
//' R (at `obs_cm[jj+1]`), the weight is `w = (t - cM_L) / (cM_R - cM_L)` and the
//' continuous dosage is `vL + w*(vR - vL)`. Ends are clamped to the terminal
//' observed value (`stats::approx` rule = 2); tied flanking positions
//' (`denom == 0`) collapse to `w = 0`.
//'
//' @param obs_cm Numeric cM of the observed markers, ascending, length k.
//' @param G Numeric k x n genotype matrix (row = observed marker, col = sample),
//'   COMPLETE (no NA); values are the alt/teosinte-allele dosage in 0 to 2.
//' @param target_cm Numeric cM of the target grid, ascending, length M.
//' @param mode Interpolation mode: 0 = continuous dosage ramp (Tian 2011),
//'   1 = step / nearest flanking value (`w < 0.5 ? vL : vR`, tie `w == 0.5` -> vR;
//'   Chen/TeoNAM densification), 2 = round(continuous) to 0/1/2, 3 = Chen 2019
//'   composite-map rule (concordant flanks fill, discordant flanks or chromosome
//'   ends -> NA; distance-independent).
//' @return Numeric M x n matrix of interpolated genotypes. Mode 3 may contain NA
//'   (discordant flanks / chromosome ends); modes 0-2 never introduce NA.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix interp_geno_cpp(NumericVector obs_cm,
                              NumericMatrix G,
                              NumericVector target_cm,
                              int mode) {
  const int k = obs_cm.size();
  const int n = G.ncol();
  const int M = target_cm.size();

  NumericMatrix out(M, n);
  if (k == 0 || M == 0) return out;   // nothing to interpolate onto / from

  // mode 3 = Chen 2019 composite-map rule: concordant flanks -> fill, discordant
  // OR chromosome ends -> NA. Order-based; the coordinate only locates flanks (no
  // distance weight). Unlike modes 0/1/2 this can EMIT NA (ends are not clamped).
  if (mode == 3) {
    std::vector<int> lf(M), rf(M);                 // flank indices; lf==rf -> exact obs; lf<0 -> end
    for (int m = 0; m < M; ++m) {
      const double t = target_cm[m];
      if (t < obs_cm[0] || t > obs_cm[k - 1]) { lf[m] = -1; rf[m] = -1; continue; }  // end -> NA
      NumericVector::iterator up = std::upper_bound(obs_cm.begin(), obs_cm.end(), t);
      int r = (int)(up - obs_cm.begin());          // first obs strictly greater than t
      if (r > 0 && obs_cm[r - 1] == t) { lf[m] = r - 1; rf[m] = r - 1; }  // exact observed marker
      else { lf[m] = r - 1; rf[m] = r; }                                   // strictly between l, r
    }
    for (int i = 0; i < n; ++i) {
      for (int m = 0; m < M; ++m) {
        if (lf[m] < 0) { out(m, i) = NA_REAL; continue; }             // chromosome end
        if (lf[m] == rf[m]) { out(m, i) = G(lf[m], i); continue; }    // exact observed value
        const double vL = G(lf[m], i), vR = G(rf[m], i);
        out(m, i) = (vL == vR) ? vL : NA_REAL;                        // concordant fill, else NA
      }
    }
    return out;
  }

  // Precompute (jj, w) once: shared across every sample.
  std::vector<int> jj(M);
  std::vector<double> w(M);
  for (int m = 0; m < M; ++m) {
    const double t = target_cm[m];
    if (t <= obs_cm[0]) {                       // clamp left end (rule = 2)
      jj[m] = 0; w[m] = 0.0;
    } else if (t >= obs_cm[k - 1]) {            // clamp right end (rule = 2)
      jj[m] = k - 1; w[m] = 0.0;
    } else {
      // first obs_cm strictly greater than t -> right flank; left flank is one before.
      NumericVector::iterator up = std::upper_bound(obs_cm.begin(), obs_cm.end(), t);
      int r = (int)(up - obs_cm.begin());       // 1 <= r <= k-1 here
      int l = r - 1;
      jj[m] = l;
      const double denom = obs_cm[r] - obs_cm[l];
      w[m] = (denom == 0.0) ? 0.0 : (t - obs_cm[l]) / denom;   // tie -> w = 0
    }
  }

  for (int i = 0; i < n; ++i) {
    for (int m = 0; m < M; ++m) {
      const int l = jj[m];
      const double vL = G(l, i);
      const double vR = (l + 1 < k) ? G(l + 1, i) : vL;   // right end: no right flank
      const double ww = w[m];
      double val;
      if (mode == 1) {                    // step: nearest flanking value (tie -> vR)
        val = (ww < 0.5) ? vL : vR;
      } else {
        val = vL + ww * (vR - vL);        // continuous dosage
        if (mode == 2) val = std::round(val);   // round: rasterize to 0/1/2
      }
      out(m, i) = val;
    }
  }
  return out;
}