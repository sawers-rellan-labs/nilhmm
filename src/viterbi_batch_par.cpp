// Threaded batched Viterbi (RcppParallel/TBB). The per-sample decode is
// embarrassingly parallel (samples share the transition + memoized emission and
// write disjoint output columns), so this splits the sample axis across cores.
// On arm64 this is the higher-ROI lever than NEON SIMD (NEON = 2 doubles/reg;
// here = up to #cores). Prototype: opt-in (call_ancestry(parallel = TRUE));
// the serial viterbi_batch_cpp remains the default validated path.

#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
using namespace Rcpp;
using namespace RcppParallel;

struct ViterbiWorker : public Worker {
  const RVector<double> log_init;
  const RMatrix<double> log_trans;
  const RMatrix<double> em_uniq;
  const RMatrix<int>    inv;      // T x S, 0-based emission row per (t, sample)
  RMatrix<int>          paths;    // T x S output
  const int T, S, K;

  ViterbiWorker(const NumericVector& li, const NumericMatrix& lt,
                const NumericMatrix& eu, const IntegerMatrix& iv, IntegerMatrix& pth)
    : log_init(li), log_trans(lt), em_uniq(eu), inv(iv), paths(pth),
      T(iv.nrow()), S(iv.ncol()), K(eu.ncol()) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<double> prev(K), cur(K);
    std::vector<int> psi(static_cast<size_t>(T) * K);
    for (std::size_t s = begin; s < end; ++s) {
      for (int k = 0; k < K; ++k) prev[k] = log_init[k] + em_uniq(inv(0, s), k);
      for (int t = 1; t < T; ++t) {
        const int row = inv(t, s);
        for (int k = 0; k < K; ++k) {
          double best = R_NegInf; int arg = 0;
          for (int j = 0; j < K; ++j) {
            double v = prev[j] + log_trans(j, k);
            if (v > best) { best = v; arg = j; }
          }
          cur[k] = best + em_uniq(row, k);
          psi[static_cast<size_t>(t) * K + k] = arg;
        }
        std::swap(prev, cur);
      }
      double best = R_NegInf; int last = 0;
      for (int k = 0; k < K; ++k) if (prev[k] > best) { best = prev[k]; last = k; }
      paths(T - 1, s) = last;
      for (int t = T - 1; t > 0; --t)
        paths(t - 1, s) = psi[static_cast<size_t>(t) * K + paths(t, s)];
    }
  }
};

//' Threaded batched log-space Viterbi (RcppParallel)
//'
//' Identical result to [viterbi_batch_cpp]; splits the sample axis across cores.
//' @inheritParams viterbi_batch_cpp
//' @return T x S integer matrix of most-likely state paths (0-based).
//' @keywords internal
// [[Rcpp::export]]
IntegerMatrix viterbi_batch_par_cpp(NumericVector log_init, NumericMatrix log_trans,
                                    NumericMatrix em_uniq, IntegerMatrix inv) {
  const int T = inv.nrow(), S = inv.ncol();
  IntegerMatrix paths(T, S);
  if (T == 0 || S == 0) return paths;
  if (log_init.size() != em_uniq.ncol() || log_trans.nrow() != em_uniq.ncol())
    stop("viterbi_batch_par_cpp: dimension mismatch (K)");
  ViterbiWorker w(log_init, log_trans, em_uniq, inv, paths);
  parallelFor(0, S, w);
  return paths;
}
