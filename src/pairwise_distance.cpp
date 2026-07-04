// Pluggable pairwise-relatedness matrix builder for marker thinning. Given a
// markers x samples genotype matrix (dosages, typically {0,1,2}), compute a
// symmetric markers x markers matrix under a chosen measure, to be fed to
// `fast_indep_cpp` for LD-based independent-set selection.
//
// Measures (all computed over the samples where BOTH markers are non-missing):
//   * r2 (method 0, a SIMILARITY): squared Pearson correlation of the dosage
//     rows -- the PLINK --indep-pairwise convention.
//   * mi (method 1, a SIMILARITY): plug-in mutual information from the joint
//     genotype contingency table with the Miller-Madow bias correction.
//   * vi (method 2, a DISTANCE, a true metric): variation of information
//     VI(X,Y) = H(X) + H(Y) - 2*I(X,Y), on MM-corrected entropies/MI.
//
// The measure is measure-agnostic downstream: `fast_indep_cpp` only needs the
// matrix and the edge sense (>= for similarity, <= for distance). No r2<->MI
// conversion layer -- each measure is thresholded in its own units.
//
// mi/vi are reported in nats (`base = 0`) or bits (`base = 1`). Entropies are
// estimated in nats (the MM correction (m-1)/(2N) is a nats-unit bias term),
// then the final MI/VI is divided by ln(2) for bits. Missing genotypes are any
// non-finite entry (NA/NaN); the R wrapper maps a caller's sentinel to NA first.

#include <Rcpp.h>
#include <vector>
#include <map>
#include <cmath>
using namespace Rcpp;

namespace {

// Squared Pearson correlation of rows i, j of `g` over complete sample pairs.
// Zero-variance (constant) marker -> 0 (no linkage information).
double r2_pair(const NumericMatrix& g, int i, int j) {
  const int S = g.ncol();
  double n = 0, sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
  for (int s = 0; s < S; ++s) {
    double x = g(i, s), y = g(j, s);
    if (!R_finite(x) || !R_finite(y)) continue;
    n += 1; sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
  }
  if (n < 2) return 0.0;
  double varx = n * sxx - sx * sx, vary = n * syy - sy * sy, cov = n * sxy - sx * sy;
  if (varx <= 0.0 || vary <= 0.0) return 0.0;
  return (cov * cov) / (varx * vary);
}

// Miller-Madow corrected MI and VI (both in nats) for rows i, j over complete
// pairs. Fills `mi_out` and `vi_out`. Builds the joint genotype contingency
// table with std::map so any integer dosage states are handled (not just 0/1/2).
void mi_vi_pair(const NumericMatrix& g, int i, int j, double& mi_out, double& vi_out) {
  const int S = g.ncol();
  std::map<int, double> cx, cy;
  std::map<std::pair<int, int>, double> cxy;
  double N = 0;
  for (int s = 0; s < S; ++s) {
    double xv = g(i, s), yv = g(j, s);
    if (!R_finite(xv) || !R_finite(yv)) continue;
    int x = (int)std::floor(xv + 0.5), y = (int)std::floor(yv + 0.5);   // categorize dosage
    cx[x] += 1; cy[y] += 1; cxy[std::make_pair(x, y)] += 1; N += 1;
  }
  if (N < 1) { mi_out = 0.0; vi_out = 0.0; return; }

  double Hx = 0, Hy = 0, Hxy = 0;
  for (std::map<int, double>::const_iterator it = cx.begin(); it != cx.end(); ++it) {
    double p = it->second / N; Hx -= p * std::log(p);
  }
  for (std::map<int, double>::const_iterator it = cy.begin(); it != cy.end(); ++it) {
    double p = it->second / N; Hy -= p * std::log(p);
  }
  for (std::map<std::pair<int, int>, double>::const_iterator it = cxy.begin();
       it != cxy.end(); ++it) {
    double p = it->second / N; Hxy -= p * std::log(p);
  }
  // Miller-Madow: add (m - 1)/(2N) per entropy, m = number of non-empty cells.
  double mX = (double)cx.size(), mY = (double)cy.size(), mXY = (double)cxy.size();
  double Hx_mm  = Hx  + (mX  - 1.0) / (2.0 * N);
  double Hy_mm  = Hy  + (mY  - 1.0) / (2.0 * N);
  double Hxy_mm = Hxy + (mXY - 1.0) / (2.0 * N);
  double I_mm   = Hx_mm + Hy_mm - Hxy_mm;                 // = I_plugin + (mX+mY-mXY-1)/(2N)
  mi_out = I_mm;
  vi_out = Hx_mm + Hy_mm - 2.0 * I_mm;                    // = Hxy_mm - I_mm  (>= 0, metric)
}

}  // namespace

//' Pairwise marker relatedness matrix (r2 / MI / VI)
//'
//' Build a symmetric markers x markers relatedness matrix from a
//' markers x samples dosage matrix, under a pluggable measure, for feeding
//' `fast_indep_cpp`. Computed over the samples where both markers are
//' non-missing (any non-finite entry is missing).
//'
//' @param geno Numeric matrix, rows = markers, cols = samples; dosages
//'   (typically 0/1/2). Missing entries are `NA`/`NaN`.
//' @param method 0 = r2 (squared Pearson correlation, a similarity),
//'   1 = mi (Miller-Madow mutual information, a similarity),
//'   2 = vi (variation of information, a distance / metric).
//' @param base For mi/vi only: 0 = nats, 1 = bits. Ignored for r2.
//' @return Symmetric numeric matrix (markers x markers). Diagonal: r2 = 1,
//'   mi = H(X), vi = 0.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix pairwise_distance_cpp(NumericMatrix geno, int method, int base) {
  const int m = geno.nrow();
  NumericMatrix out(m, m);
  const double logbase = (base == 1) ? std::log(2.0) : 1.0;   // nats -> bits divisor

  for (int i = 0; i < m; ++i) {
    for (int j = i; j < m; ++j) {
      double v;
      if (method == 0) {
        v = (i == j) ? 1.0 : r2_pair(geno, i, j);   // diagonal is 1 even for a constant marker
      } else {
        double mi, vi;
        mi_vi_pair(geno, i, j, mi, vi);
        v = (method == 2 ? vi : mi) / logbase;
      }
      out(i, j) = v;
      out(j, i) = v;
    }
  }
  return out;
}
