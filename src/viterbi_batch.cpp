// Batched Viterbi over many samples sharing one transition and one memoized
// emission table (REFACTOR S4; the perf path mirroring the Python numpy-batched
// caller). Looping samples in C++ removes the per-(sample,chromosome) R
// orchestration (unique/match/alloc/split/GC) that dominates call_ancestry
// wall-clock; emission is read through a precomputed index so the BetaBinomial
// is evaluated once per distinct (n,a) pair, not per cell.

#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

//' Batched log-space Viterbi (shared transition, memoized emission)
//'
//' @param log_init Length-K log initial-state probabilities.
//' @param log_trans K x K log transition matrix (row = from, col = to).
//' @param em_uniq U x K log emission table, one row per distinct observation.
//' @param inv T x S integer matrix; inv(t,s) is the 0-based row of `em_uniq`
//'   giving the emission for sample s at marker t.
//' @return T x S integer matrix of most-likely state paths (0-based states).
//' @keywords internal
// [[Rcpp::export]]
IntegerMatrix viterbi_batch_cpp(NumericVector log_init, NumericMatrix log_trans,
                                NumericMatrix em_uniq, IntegerMatrix inv) {
  const int T = inv.nrow();
  const int S = inv.ncol();
  const int K = em_uniq.ncol();
  IntegerMatrix paths(T, S);
  if (T == 0 || S == 0) return paths;
  if (log_init.size() != K || log_trans.nrow() != K || log_trans.ncol() != K)
    stop("viterbi_batch_cpp: dimension mismatch (K)");

  std::vector<double> prev(K), cur(K);
  std::vector<int> psi(static_cast<size_t>(T) * K);

  for (int s = 0; s < S; ++s) {
    int r0 = inv(0, s);
    for (int k = 0; k < K; ++k) prev[k] = log_init[k] + em_uniq(r0, k);
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
  return paths;
}
