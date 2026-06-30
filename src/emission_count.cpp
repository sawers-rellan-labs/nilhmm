// BetaBinomial count emission (REFACTOR_R_PACKAGE.md §4, §5). Reproduces the
// Python introgression_hmm_counts emission: for marker t and state s with
// expected alt-fraction theta_s, emission = BetaBinomial(a_t | n_t,
// theta_s*conc, (1-theta_s)*conc). Depth-0 markers (n=0,a=0) emit log 0 (flat /
// uninformative), since lchoose(0,0)=0 and lbeta(alpha,beta) cancels.
//
// logpmf(a|n,A,B) = lchoose(n,a) + lbeta(a+A, n-a+B) - lbeta(A,B)

#include <Rcpp.h>
using namespace Rcpp;

//' Log BetaBinomial emission matrix over REF/HET/ALT states
//'
//' @param n Integer vector of total depths (n_ref + n_alt) per marker.
//' @param a Integer vector of alt counts (n_alt) per marker.
//' @param theta Length-K vector of expected alt-fractions per state
//'   (e.g. c(err, 0.5, 1 - err)).
//' @param conc BetaBinomial concentration (alpha + beta); larger -> closer to
//'   Binomial.
//' @return T x K matrix of log emission probabilities.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix count_emission_loglik_cpp(IntegerVector n, IntegerVector a,
                                        NumericVector theta, double conc) {
  const int T = n.size();
  const int K = theta.size();
  if (a.size() != T) stop("count_emission_loglik_cpp: length(a) must equal length(n)");
  if (conc <= 0.0)   stop("count_emission_loglik_cpp: conc must be > 0");

  std::vector<double> alpha(K), beta(K), lbeta_ab(K);
  for (int s = 0; s < K; ++s) {
    if (theta[s] <= 0.0 || theta[s] >= 1.0)
      stop("count_emission_loglik_cpp: theta must be strictly in (0, 1)");
    alpha[s] = theta[s] * conc;
    beta[s]  = (1.0 - theta[s]) * conc;
    lbeta_ab[s] = R::lbeta(alpha[s], beta[s]);
  }

  NumericMatrix out(T, K);
  for (int t = 0; t < T; ++t) {
    const double nt = n[t], at = a[t];
    if (nt < 0 || at < 0 || at > nt)
      stop("count_emission_loglik_cpp: require 0 <= a <= n");
    const double lchoose_t = R::lchoose(nt, at);
    for (int s = 0; s < K; ++s) {
      out(t, s) = lchoose_t + R::lbeta(at + alpha[s], nt - at + beta[s]) - lbeta_ab[s];
    }
  }
  return out;
}
