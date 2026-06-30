// Hot loops for the nilHMM engine (REFACTOR_R_PACKAGE.md §4, §8).
// Rcpp recovers Julia-class speed within a CRAN/Bioconductor-friendly package.
//
// This file currently provides the generic log-space Viterbi workhorse. The
// state-expansion (rigidity) decoder, BetaBinomial emission evaluation, and the
// EM hot loop land in Task 4; they build their (log_init, log_trans, log_emit)
// arguments in R and call viterbi_log_cpp here.

#include <Rcpp.h>
using namespace Rcpp;

//' Generic log-space Viterbi decode (time-homogeneous transitions)
//'
//' @param log_init Length-K vector of log initial-state probabilities.
//' @param log_trans K x K matrix of log transition probabilities
//'   (row = from-state, column = to-state).
//' @param log_emit T x K matrix of log emission probabilities per marker/state.
//' @return Integer length-T most-likely state path, 0-indexed.
//' @keywords internal
// [[Rcpp::export]]
IntegerVector viterbi_log_cpp(NumericVector log_init,
                              NumericMatrix log_trans,
                              NumericMatrix log_emit) {
  const int T = log_emit.nrow();
  const int K = log_emit.ncol();
  if (log_init.size() != K)
    stop("viterbi_log_cpp: length(log_init) must equal ncol(log_emit)");
  if (log_trans.nrow() != K || log_trans.ncol() != K)
    stop("viterbi_log_cpp: log_trans must be K x K");

  IntegerVector path(T);
  if (T == 0) return path;

  NumericMatrix delta(T, K);   // best log-prob ending in state k at t
  IntegerMatrix psi(T, K);     // argmax backpointer

  for (int k = 0; k < K; ++k) {
    delta(0, k) = log_init[k] + log_emit(0, k);
    psi(0, k) = 0;
  }

  for (int t = 1; t < T; ++t) {
    for (int k = 0; k < K; ++k) {
      double best = R_NegInf;
      int arg = 0;
      for (int j = 0; j < K; ++j) {
        double v = delta(t - 1, j) + log_trans(j, k);
        if (v > best) { best = v; arg = j; }
      }
      delta(t, k) = best + log_emit(t, k);
      psi(t, k) = arg;
    }
  }

  // terminate + backtrace
  double best = R_NegInf;
  int last = 0;
  for (int k = 0; k < K; ++k) {
    if (delta(T - 1, k) > best) { best = delta(T - 1, k); last = k; }
  }
  path[T - 1] = last;
  for (int t = T - 1; t > 0; --t) path[t - 1] = psi(t, path[t]);

  return path;
}