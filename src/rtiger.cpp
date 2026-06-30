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
