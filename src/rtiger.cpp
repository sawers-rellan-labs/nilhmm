// Faithful port of RTIGER's rHMM_methods.jl (faustovrz/RTIGER fork, commit
// 649cbf6, branch optimize-julia-core). Each kernel reproduces the SERIAL path's
// arithmetic (the fork marks it bit-identical to upstream) and mirrors the
// Julia's structure line-for-line so it can be reviewed against the source.
// Verified against the installed Julia via agent/rtref.jl (max|Δ| ~5e-9).
//
// Data convention (from generateObject): k = P1/ref-allele count, n = total.
// States ordered [pat (ref~0.95), het (0.5), mat/donor (ref~0.05)].
//
// Indices follow the Julia as 1-based loop variables; matrix/array accesses
// subtract 1. (s, T) = (#states, #markers).

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <unordered_map>
using namespace Rcpp;

// log-add-exp of two scalars, matching the fork's logaddexp (handles -Inf).
static inline double logaddexp2(double x, double y) {
  double hi = x > y ? x : y;
  double lo = x > y ? y : x;
  if (hi == R_NegInf) return R_NegInf;
  return hi + std::log1p(std::exp(lo - hi));
}


//' RTIGER emission log-probabilities (getlogpsi)
//'
//' logpsi[i,t] = logpdf(BetaBinomial(n_t, a_i, b_i), k_t). Memoized over distinct
//' (k,n) pairs (as in the fork's getlogpsi cache; bit-identical values).
//' @param k Integer vector of ref-allele counts (length T).
//' @param n Integer vector of totals (length T).
//' @param a,b Per-state BetaBinomial shape vectors (length s).
//' @return s x T matrix of log emission probabilities.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix rtiger_getlogpsi_cpp(IntegerVector k, IntegerVector n,
                                   NumericVector a, NumericVector b) {
  const int T = k.size();
  const int s = a.size();

  std::vector<double> lbeta_ab(s);
  for (int i = 0; i < s; ++i) {
    lbeta_ab[i] = R::lbeta(a[i], b[i]);
  }

  NumericMatrix logpsi(s, T);
  std::unordered_map<long long, int> seen;             // (k,n) -> a column already computed

  for (int t = 0; t < T; ++t) {
    const int kt = k[t];
    const int nt = n[t];
    const long long key = (long long)nt * 2147483647LL + kt;

    auto it = seen.find(key);
    if (it != seen.end()) {                            // reuse a previously computed column
      const int src = it->second;
      for (int i = 0; i < s; ++i) logpsi(i, t) = logpsi(i, src);
      continue;
    }

    const double lchoose_t = R::lchoose((double)nt, (double)kt);
    for (int i = 0; i < s; ++i) {
      logpsi(i, t) = lchoose_t + R::lbeta(kt + a[i], nt - kt + b[i]) - lbeta_ab[i];
    }
    seen[key] = t;
  }
  return logpsi;
}


//' RTIGER windowed emission product (productpsi)
//'
//' PSI[i,t] = sum of psi[i, t-r+1 .. t] (sliding window of r), as a running
//' cumulative sum. PSI is s x (T+r): the tail columns T+1..T+r-1 carry the
//' trailing partial sums and column T+r is left 0 (exactly as the fork).
//' @param logpsi s x T matrix from rtiger_getlogpsi_cpp.
//' @param r Rigidity.
//' @return s x (T+r) matrix.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix rtiger_productpsi_cpp(NumericMatrix logpsi, int r) {
  const int s = logpsi.nrow();
  const int T = logpsi.ncol();
  NumericMatrix PSI(s, T + r);                          // zero-initialized

  for (int i = 0; i < s; ++i) {
    PSI(i, 0) = logpsi(i, 0);

    // grow the window to size r
    for (int t = 1; t < r; ++t) {
      PSI(i, t) = PSI(i, t - 1) + logpsi(i, t);
    }
    // slide the window: add the entering marker, drop the leaving one
    for (int t = r; t < T; ++t) {
      PSI(i, t) = PSI(i, t - 1) + logpsi(i, t) - logpsi(i, t - r);
    }
    // drain the window past the last marker
    for (int t = T; t < T + r - 1; ++t) {
      PSI(i, t) = PSI(i, t - 1) - logpsi(i, t - r);
    }
  }
  return PSI;
}


//' RTIGER forward pass (rigidity-aware Baum-Welch forward)
//'
//' Literal port of the fork's `forward`: at each step a state either "stays"
//' (diagonal transition) or was "entered" from another state r positions back.
//' @param logPI length-s log start probabilities.
//' @param logPSI s x (T+r) windowed-emission matrix (rtiger_productpsi_cpp).
//' @param logA s x s log transition matrix.
//' @param logpsi s x T log emission matrix (rtiger_getlogpsi_cpp).
//' @param r Rigidity.
//' @return s x T log forward matrix (cols 1 and >T-r+1 stay -Inf).
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix rtiger_forward_cpp(NumericVector logPI, NumericMatrix logPSI,
                                 NumericMatrix logA, NumericMatrix logpsi, int r) {
  const int s = logpsi.nrow();
  const int T = logpsi.ncol();

  NumericMatrix alpha(s, T);
  std::fill(alpha.begin(), alpha.end(), R_NegInf);

  // alpha[:, r] = logPI + logPSI[:, r]
  for (int k = 1; k <= s; ++k) {
    alpha(k - 1, r - 1) = logPI[k - 1] + logPSI(k - 1, r - 1);
  }

  // first window: only the diagonal "stay" is possible
  for (int t = r + 1; t <= 2 * r - 1; ++t) {
    for (int k = 1; k <= s; ++k) {
      alpha(k - 1, t - 1) = logpsi(k - 1, t - 1) + logA(k - 1, k - 1) + alpha(k - 1, t - 2);
    }
  }

  for (int t = 2 * r; t <= T - r + 1; ++t) {
    for (int k = 1; k <= s; ++k) {
      // stay in state k (diagonal transition)
      const double stay = logpsi(k - 1, t - 1) + logA(k - 1, k - 1) + alpha(k - 1, t - 2);

      // enter state k from any i != k, r positions back; log-sum-exp over i,
      // factoring out the max (amax) for numerical stability (matches the fork)
      double maxv = R_NegInf;
      int amax = 0;
      for (int i = 1; i <= s; ++i) {
        if (i == k) continue;
        const double x = logA(i - 1, k - 1) + alpha(i - 1, t - r - 1);
        if (x > maxv) { maxv = x; amax = i; }
      }

      double enter = R_NegInf;
      if (maxv != R_NegInf) {
        double acc = 0.0;
        for (int i = 1; i <= s; ++i) {
          if (i == k || i == amax) continue;
          acc += std::exp(logA(i - 1, k - 1) + alpha(i - 1, t - r - 1) - maxv);
        }
        enter = logPSI(k - 1, t - 1) + maxv + std::log1p(acc);
      }

      alpha(k - 1, t - 1) = logaddexp2(stay, enter);
    }
  }
  return alpha;
}


//' RTIGER backward pass (rigidity-aware Baum-Welch backward)
//'
//' Literal port of the fork's `backward` (the time-reversed mirror of forward).
//' @inheritParams rtiger_forward_cpp
//' @return s x T log backward matrix.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix rtiger_backward_cpp(NumericMatrix logPSI, NumericMatrix logA,
                                  NumericMatrix logpsi, int r) {
  const int s = logpsi.nrow();
  const int T = logpsi.ncol();

  NumericMatrix beta(s, T);
  std::fill(beta.begin(), beta.end(), R_NegInf);

  // beta[:, (T-r+1):T] = logPSI[:, (T+1):(T+r)]   (col c <- logPSI col c+r)
  for (int c = T - r + 1; c <= T; ++c) {
    for (int j = 1; j <= s; ++j) {
      beta(j - 1, c - 1) = logPSI(j - 1, c + r - 1);
    }
  }

  for (int i = 0; i <= T - 2 * r; ++i) {
    const int t = T - r - i;
    for (int j = 1; j <= s; ++j) {
      // stay in state j (diagonal transition)
      const double stay = logpsi(j - 1, t) + logA(j - 1, j - 1) + beta(j - 1, t);  // index (t+1)-1

      // leave state j to any k != j, r positions ahead
      double maxv = R_NegInf;
      int amax = 0;
      for (int k = 1; k <= s; ++k) {
        if (k == j) continue;
        const double x = logA(j - 1, k - 1) + logPSI(k - 1, t + r - 1) + beta(k - 1, t + r - 1);
        if (x > maxv) { maxv = x; amax = k; }
      }

      double leave = R_NegInf;
      if (maxv != R_NegInf) {
        double acc = 0.0;
        for (int k = 1; k <= s; ++k) {
          if (k == j || k == amax) continue;
          acc += std::exp(logA(j - 1, k - 1) + logPSI(k - 1, t + r - 1) + beta(k - 1, t + r - 1) - maxv);
        }
        leave = maxv + std::log1p(acc);
      }

      beta(j - 1, t - 1) = logaddexp2(stay, leave);
    }
  }
  return beta;
}


//' RTIGER zeta (pairwise posteriors over the rigidity window)
//'
//' Literal port of the fork's `zeta`: build the log values, then normalize by
//' the global max and exponentiate (so exp(-Inf - PO) = 0). Returned as a
//' T x s x s array (column-major), non-log, as in the Julia.
//' @return T x s x s numeric array.
//' @keywords internal
// [[Rcpp::export]]
NumericVector rtiger_zeta_cpp(NumericMatrix logalpha, NumericMatrix logbeta,
                              NumericMatrix logA, NumericMatrix logPSI,
                              NumericMatrix logpsi, int r) {
  const int s = logpsi.nrow();
  const int T = logpsi.ncol();

  NumericVector zeta(T * s * s, R_NegInf);
  // flat index into a T x s x s column-major array (0-based t, j, k)
  auto at = [&](int t, int j, int k) -> double& { return zeta[t + j * T + k * (T * s)]; };

  for (int k = 1; k <= s; ++k) {
    for (int j = 1; j <= s; ++j) {
      if (j != k) {
        for (int t = r + 1; t <= T - r + 1; ++t) {
          at(t - 1, j - 1, k - 1) = logalpha(j - 1, t - 2) + logA(j - 1, k - 1)
                                  + logbeta(k - 1, t + r - 2) + logPSI(k - 1, t + r - 2);
        }
      } else {
        for (int t = r + 1; t <= T - r + 1; ++t) {
          at(t - 1, k - 1, k - 1) = logalpha(k - 1, t - 2) + logA(k - 1, k - 1)
                                  + logpsi(k - 1, t - 1) + logbeta(k - 1, t - 1);
        }
      }
    }
  }

  double PO = R_NegInf;
  for (R_xlen_t i = 0; i < zeta.size(); ++i) {
    if (zeta[i] > PO) PO = zeta[i];
  }
  for (R_xlen_t i = 0; i < zeta.size(); ++i) {
    zeta[i] = std::exp(zeta[i] - PO);
  }

  zeta.attr("dim") = IntegerVector::create(T, s, s);
  return zeta;
}


//' RTIGER gamma (state posteriors from zeta)
//'
//' Literal port of the fork's `gamma`. Works in linear space (zeta is already
//' exponentiated): the first r columns are the normalized alpha*beta at the
//' window edge; interior columns add the windowed cumulative sum of the
//' off-diagonal zeta to the diagonal term; the tail repeats; columns are then
//' normalized to sum to 1.
//' @param zeta T x s x s array from rtiger_zeta_cpp.
//' @return s x T matrix of state posteriors.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix rtiger_gamma_cpp(NumericVector zeta, NumericMatrix logalpha,
                               NumericMatrix logbeta, int r) {
  const int s = logalpha.nrow();
  const int T = logalpha.ncol();
  auto at = [&](int t, int j, int k) -> double { return zeta[t + j * T + k * (T * s)]; };

  NumericMatrix gam(s, T);

  // gammar = normalize(exp(logalpha[:,r] + logbeta[:,r]))
  std::vector<double> gammar(s);
  double mx = R_NegInf;
  for (int k = 0; k < s; ++k) {
    gammar[k] = logalpha(k, r - 1) + logbeta(k, r - 1);
    if (gammar[k] > mx) mx = gammar[k];
  }
  double gsum = 0.0;
  for (int k = 0; k < s; ++k) {
    gammar[k] = std::exp(gammar[k] - mx);
    gsum += gammar[k];
  }
  for (int k = 0; k < s; ++k) {
    gammar[k] /= gsum;
  }

  // gam[:, 1:r] = repeat(gammar)
  for (int t = 1; t <= r; ++t) {
    for (int k = 0; k < s; ++k) gam(k, t - 1) = gammar[k];
  }

  std::vector<double> L(T);   // L[t] = sum_{j != k} zeta[t, j, k]
  std::vector<double> M(T);   // M    = cumsum(L)
  for (int k = 1; k <= s; ++k) {
    for (int t = 1; t <= T; ++t) {
      double acc = 0.0;
      for (int j = 1; j <= s; ++j) {
        if (j != k) acc += at(t - 1, j - 1, k - 1);
      }
      L[t - 1] = acc;
    }

    double run = 0.0;
    for (int t = 1; t <= T; ++t) {
      run += L[t - 1];
      M[t - 1] = run;
    }

    for (int t = r + 1; t <= T - r + 1; ++t) {
      const double Mm = (t <= 2 * r) ? M[r - 1] : M[t - r - 1];   // window-back sum
      gam(k - 1, t - 1) = at(t - 1, k - 1, k - 1) + M[t - 1] - Mm;
    }
  }

  // gam[:, T-r+2:T] = repeat(gam[:, T-r+1])
  for (int t = T - r + 2; t <= T; ++t) {
    for (int k = 0; k < s; ++k) gam(k, t - 1) = gam(k, T - r);    // col T-r+1 is index T-r
  }

  // normalize each column to sum to 1
  for (int t = 0; t < T; ++t) {
    double cs = 0.0;
    for (int k = 0; k < s; ++k) cs += gam(k, t);
    for (int k = 0; k < s; ++k) gam(k, t) /= cs;
  }
  return gam;
}
