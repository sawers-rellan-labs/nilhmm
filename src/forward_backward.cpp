// Forward-backward posteriors for the Baum-Welch E-step (REFACTOR S4, S10).
// Log-space alpha/beta; returns gamma_t(k) = P(state_t = k | obs) in normal
// space. Used by the EM that fits the count-emission means (.em_fit_means),
// the S10 fix for reference-biased data where the fixed theta_ALT collapses
// ALT into HET. -Inf entries (forbidden transitions in the expanded rigidity
// model) are handled in the log-sum-exp.

#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

static inline double logsumexp(const std::vector<double>& v) {
  double m = R_NegInf;
  for (double x : v) if (x > m) m = x;
  if (!R_FINITE(m)) return R_NegInf;        // all -Inf
  double s = 0.0;
  for (double x : v) s += std::exp(x - m);
  return m + std::log(s);
}

//' Forward-backward state posteriors (log-space; time-homogeneous transitions)
//'
//' @param log_init Length-K log initial-state probabilities.
//' @param log_trans K x K log transition matrix (row = from, col = to).
//' @param log_emit T x K log emission matrix.
//' @return T x K matrix of posterior state probabilities (rows sum to 1).
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix forward_backward_cpp(NumericVector log_init,
                                   NumericMatrix log_trans,
                                   NumericMatrix log_emit) {
  const int T = log_emit.nrow();
  const int K = log_emit.ncol();
  NumericMatrix g(T, K);
  if (T == 0) return g;

  NumericMatrix la(T, K), lb(T, K);

  for (int k = 0; k < K; ++k) la(0, k) = log_init[k] + log_emit(0, k);
  for (int t = 1; t < T; ++t) {
    for (int k = 0; k < K; ++k) {
      std::vector<double> terms(K);
      for (int j = 0; j < K; ++j) terms[j] = la(t - 1, j) + log_trans(j, k);
      la(t, k) = logsumexp(terms) + log_emit(t, k);
    }
  }

  for (int k = 0; k < K; ++k) lb(T - 1, k) = 0.0;
  for (int t = T - 2; t >= 0; --t) {
    for (int k = 0; k < K; ++k) {
      std::vector<double> terms(K);
      for (int j = 0; j < K; ++j) terms[j] = log_trans(k, j) + log_emit(t + 1, j) + lb(t + 1, j);
      lb(t, k) = logsumexp(terms);
    }
  }

  std::vector<double> last(K);
  for (int k = 0; k < K; ++k) last[k] = la(T - 1, k);
  const double ll = logsumexp(last);        // total log-likelihood

  for (int t = 0; t < T; ++t)
    for (int k = 0; k < K; ++k)
      g(t, k) = std::exp(la(t, k) + lb(t, k) - ll);

  return g;
}
