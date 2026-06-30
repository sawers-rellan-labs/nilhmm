// Faithful port of RTIGER's rHMM_methods.jl (faustovrz/RTIGER fork, commit
// 649cbf6, branch optimize-julia-core). Each kernel reproduces the SERIAL path's
// arithmetic (the fork marks it bit-identical to upstream). Verified against the
// installed Julia via agent/rtref.jl.
//
// Data convention (from generateObject): k = P1/ref-allele count, n = total.
// States ordered [pat (ref~0.95), het (0.5), mat/donor (ref~0.05)].
//
// This file is built incrementally; kernels are added + Julia-verified one at a
// time. Done so far: getlogpsi, productpsi.

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <unordered_map>
using namespace Rcpp;

//' RTIGER emission log-probabilities (getlogpsi)
//'
//' logpsi[i,t] = logpdf(BetaBinomial(n_t, a_i, b_i), k_t). Memoized over distinct
//' (k,n) pairs (same as the fork's getlogpsi cache; bit-identical values).
//' @param k Integer vector of ref-allele counts (length T).
//' @param n Integer vector of totals (length T).
//' @param a,b Per-state BetaBinomial shape vectors (length s).
//' @return s x T matrix of log emission probabilities.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix rtiger_getlogpsi_cpp(IntegerVector k, IntegerVector n,
                                   NumericVector a, NumericVector b) {
  const int T = k.size(), s = a.size();
  std::vector<double> lbeta_ab(s);
  for (int i = 0; i < s; ++i) lbeta_ab[i] = R::lbeta(a[i], b[i]);
  NumericMatrix logpsi(s, T);
  std::unordered_map<long long, int> cache;          // (k,n) -> column already filled
  for (int t = 0; t < T; ++t) {
    const int kt = k[t], nt = n[t];
    long long key = (long long)nt * 2147483647LL + kt;
    auto it = cache.find(key);
    if (it != cache.end()) {
      const int src = it->second;
      for (int i = 0; i < s; ++i) logpsi(i, t) = logpsi(i, src);
    } else {
      const double lchoose_t = R::lchoose((double)nt, (double)kt);
      for (int i = 0; i < s; ++i)
        logpsi(i, t) = lchoose_t + R::lbeta(kt + a[i], nt - kt + b[i]) - lbeta_ab[i];
      cache[key] = t;
    }
  }
  return logpsi;
}

//' RTIGER windowed emission product (productpsi)
//'
//' PSI[i,t] = sum of psi[i, t-r+1 .. t] (sliding window of r), as a cumulative
//' sum; PSI is s x (T+r) with the tail columns T+1..T+r-1 carrying the trailing
//' partial sums and column T+r left 0 (exactly as the fork's productpsi).
//' @param logpsi s x T matrix from rtiger_getlogpsi_cpp.
//' @param r Rigidity.
//' @return s x (T+r) matrix.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix rtiger_productpsi_cpp(NumericMatrix logpsi, int r) {
  const int s = logpsi.nrow(), T = logpsi.ncol();
  NumericMatrix PSI(s, T + r);                        // zero-initialized
  for (int i = 0; i < s; ++i) {
    PSI(i, 0) = logpsi(i, 0);
    for (int t = 1; t < r; ++t)        PSI(i, t) = PSI(i, t - 1) + logpsi(i, t);
    for (int t = r; t < T; ++t)        PSI(i, t) = PSI(i, t - 1) + logpsi(i, t) - logpsi(i, t - r);
    for (int t = T; t < T + r - 1; ++t) PSI(i, t) = PSI(i, t - 1) - logpsi(i, t - r);
  }
  return PSI;
}

// log-add-exp of two scalars, matching the fork's logaddexp (handles -Inf).
static inline double logaddexp2(double x, double y) {
  double m = x > y ? x : y, mn = x > y ? y : x;
  return m == R_NegInf ? R_NegInf : m + std::log1p(std::exp(mn - m));
}

//' RTIGER forward pass (rigidity-aware Baum-Welch forward)
//'
//' Literal port of the fork's `forward` (rHMM_methods.jl): only a diagonal
//' "stay" each step, or an "enter" from another state r positions back. Indices
//' are 1-based to mirror the Julia exactly.
//' @param logPI length-s log start probabilities.
//' @param logPSI s x (T+r) windowed-emission matrix (rtiger_productpsi_cpp).
//' @param logA s x s log transition matrix.
//' @param logpsi s x T log emission matrix (rtiger_getlogpsi_cpp).
//' @param r Rigidity.
//' @return s x T matrix of log forward probabilities (cols 1 and >T-r+1 stay -Inf).
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix rtiger_forward_cpp(NumericVector logPI, NumericMatrix logPSI,
                                 NumericMatrix logA, NumericMatrix logpsi, int r) {
  const int s = logpsi.nrow(), T = logpsi.ncol();
  NumericMatrix alpha(s, T);
  std::fill(alpha.begin(), alpha.end(), R_NegInf);
  for (int k = 1; k <= s; ++k) alpha(k - 1, r - 1) = logPI[k - 1] + logPSI(k - 1, r - 1);
  for (int t = r + 1; t <= 2 * r - 1; ++t)
    for (int k = 1; k <= s; ++k)
      alpha(k - 1, t - 1) = logpsi(k - 1, t - 1) + logA(k - 1, k - 1) + alpha(k - 1, t - 2);
  for (int t = 2 * r; t <= T - r + 1; ++t) {
    for (int k = 1; k <= s; ++k) {
      double stay = logpsi(k - 1, t - 1) + logA(k - 1, k - 1) + alpha(k - 1, t - 2);
      double maxv = R_NegInf; int amax = 0;
      for (int i = 1; i <= s; ++i) if (i != k) {
        double x = logA(i - 1, k - 1) + alpha(i - 1, t - r - 1);
        if (x > maxv) { maxv = x; amax = i; }
      }
      double enter = R_NegInf;
      if (maxv != R_NegInf) {
        double acc = 0.0;
        for (int i = 1; i <= s; ++i) if (i != k && i != amax)
          acc += std::exp(logA(i - 1, k - 1) + alpha(i - 1, t - r - 1) - maxv);
        enter = logPSI(k - 1, t - 1) + maxv + std::log1p(acc);
      }
      alpha(k - 1, t - 1) = logaddexp2(stay, enter);
    }
  }
  return alpha;
}

//' RTIGER backward pass (rigidity-aware Baum-Welch backward)
//'
//' Literal port of the fork's `backward`. @inheritParams rtiger_forward_cpp
//' @return s x T matrix of log backward probabilities.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix rtiger_backward_cpp(NumericMatrix logPSI, NumericMatrix logA,
                                  NumericMatrix logpsi, int r) {
  const int s = logpsi.nrow(), T = logpsi.ncol();
  NumericMatrix beta(s, T);
  std::fill(beta.begin(), beta.end(), R_NegInf);
  for (int c = T - r + 1; c <= T; ++c)
    for (int j = 1; j <= s; ++j) beta(j - 1, c - 1) = logPSI(j - 1, c + r - 1);
  for (int i = 0; i <= T - 2 * r; ++i) {
    const int t = T - r - i;
    for (int j = 1; j <= s; ++j) {
      double stay = logpsi(j - 1, t) + logA(j - 1, j - 1) + beta(j - 1, t);  // (t+1)-1
      double maxv = R_NegInf; int amax = 0;
      for (int k = 1; k <= s; ++k) if (k != j) {
        double x = logA(j - 1, k - 1) + logPSI(k - 1, t + r - 1) + beta(k - 1, t + r - 1);
        if (x > maxv) { maxv = x; amax = k; }
      }
      double leave = R_NegInf;
      if (maxv != R_NegInf) {
        double acc = 0.0;
        for (int k = 1; k <= s; ++k) if (k != j && k != amax)
          acc += std::exp(logA(j - 1, k - 1) + logPSI(k - 1, t + r - 1) + beta(k - 1, t + r - 1) - maxv);
        leave = maxv + std::log1p(acc);
      }
      beta(j - 1, t - 1) = logaddexp2(stay, leave);
    }
  }
  return beta;
}
