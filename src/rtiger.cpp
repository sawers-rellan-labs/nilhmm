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
#include <RcppParallel.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include <functional>
#include <unordered_map>
using namespace Rcpp;
using namespace RcppParallel;

// forward declarations (kernels defined below; used by the EM driver)
NumericMatrix rtiger_getlogpsi_cpp(IntegerVector, IntegerVector, NumericVector, NumericVector);
NumericMatrix rtiger_productpsi_cpp(NumericMatrix, int);
NumericMatrix rtiger_forward_cpp(NumericVector, NumericMatrix, NumericMatrix, NumericMatrix, int);
NumericMatrix rtiger_backward_cpp(NumericMatrix, NumericMatrix, NumericMatrix, int);
NumericVector rtiger_zeta_cpp(NumericMatrix, NumericMatrix, NumericMatrix, NumericMatrix, NumericMatrix, int);
NumericMatrix rtiger_gamma_cpp(NumericVector, NumericMatrix, NumericMatrix, int);

// log-add-exp of two scalars, matching the fork's logaddexp (handles -Inf).
static inline double logaddexp2(double x, double y) {
  double hi = x > y ? x : y;
  double lo = x > y ? y : x;
  if (hi == R_NegInf) return R_NegInf;
  return hi + std::log1p(std::exp(lo - hi));
}


//' RTIGER emission log-probabilities (getlogpsi)
//'
//' logpsi(i,t) = logpdf(BetaBinomial(n_t, a_i, b_i), k_t). Memoized over distinct
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
//' PSI(i,t) = sum of psi over the window (t-r+1 .. t) (sliding window of r), as a running
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


//' RTIGER Viterbi (rigidity max-product + rigid backtrace)
//'
//' Literal port of the fork's `viterbi`. Each step is a max over: "stay" in the
//' target state (diagonal), or "enter" it from another state r positions back
//' (using the windowed PSI). The backtrace fills a whole r-block with the
//' current state on a switch, giving the r-rigid path. Tie-break is first state
//' (seed at state 1, replace only on strict >), matching the fork's findmax.
//' @param PI length-s log start probabilities.
//' @param PSI s x (T+r) windowed-emission matrix.
//' @param psi s x T log emission matrix.
//' @param A s x s log transition matrix.
//' @param r Rigidity.
//' @return length-T integer state path (1-based states).
//' @keywords internal
// [[Rcpp::export]]
IntegerVector rtiger_viterbi_cpp(NumericVector PI, NumericMatrix PSI,
                                 NumericMatrix psi, NumericMatrix A, int r) {
  const int s = psi.nrow();
  const int T = psi.ncol();

  NumericMatrix phi(s, T);
  std::fill(phi.begin(), phi.end(), R_NegInf);
  IntegerMatrix b(s, T);

  // phi[:, r] = PSI[:, r] + PI ; b[:, 1:r] = 1:s
  for (int k = 1; k <= s; ++k) phi(k - 1, r - 1) = PSI(k - 1, r - 1) + PI[k - 1];
  for (int t = 1; t <= r; ++t)
    for (int k = 1; k <= s; ++k) b(k - 1, t - 1) = k;

  for (int t = r + 1; t <= T; ++t) {
    for (int j = 1; j <= s; ++j) {
      // candidate score of reaching state j at t coming from row i:
      //   i == j : stay (diagonal);  i != j : enter from i, r positions back
      auto cand = [&](int i) -> double {
        if (i == j) return psi(j - 1, t - 1) + A(j - 1, j - 1) + phi(j - 1, t - 2);
        else        return A(i - 1, j - 1) + PSI(j - 1, t - 1) + phi(i - 1, t - r - 1);
      };
      double best = cand(1);
      int bestidx = 1;
      for (int i = 2; i <= s; ++i) {
        const double val = cand(i);
        if (val > best) { best = val; bestidx = i; }
      }
      phi(j - 1, t - 1) = best;
      b(j - 1, t - 1) = bestidx;
    }
  }

  IntegerVector v(T);
  // v[T] = argmax_k phi[k, T]  (first max)
  {
    double mx = R_NegInf; int arg = 1;
    for (int k = 1; k <= s; ++k) if (phi(k - 1, T - 1) > mx) { mx = phi(k - 1, T - 1); arg = k; }
    v[T - 1] = arg;
  }

  int t = T;
  while (t > 1) {
    const int pointer = b(v[t - 1] - 1, t - 1);
    if (pointer != v[t - 1] && r > 1) {
      // on a switch, fill the trailing r-block with the current state
      const int lo = std::max(t - r + 1, 1);
      for (int p = lo; p <= t - 1; ++p) v[p - 1] = v[t - 1];
      t = t - r + 1;
      if (t < 2) break;
    }
    v[t - 2] = pointer;   // v[t-1] (1-based) <- pointer
    t = t - 1;
    if (t < 2) break;
  }
  return v;
}


//' RTIGER one-EM-iteration sufficient statistics (E-step + fold, all in C++)
//'
//' Runs the per-chain E-step (getlogpsi..gamma) and accumulates the pooled
//' sufficient statistics, replacing the slow R-level fold (apply()/tapply()).
//' Mirrors the fork's EM accumulation: transition band-sum of zeta, start =
//' Sum of gamma column 1, and per-state emission weights grouped by distinct (k,n).
//' @param ks_list,ns_list Lists of per-chain integer (k, n) vectors.
//' @param logPI,logA log start / log transition.
//' @param alpha,beta current BetaBinomial shape vectors.
//' @param r,nstates Rigidity and number of states.
//' @return list(sumZeta, startAcc, nOb, kvals, nvals, wmat, sumk, sumn).
//' @keywords internal
// [[Rcpp::export]]
List rtiger_em_suffstats_cpp(List ks_list, List ns_list, NumericVector logPI,
                             NumericMatrix logA, NumericVector alpha,
                             NumericVector beta, int r, int nstates) {
  const int s = nstates;
  NumericMatrix sumZeta(s, s);
  NumericVector startAcc(s);
  int nOb = 0;
  std::map<long long, std::vector<double> > W;          // (k,n) key -> per-state Σγ
  std::vector<double> sumk(s, 0.0), sumn(s, 0.0);

  const int nchains = ks_list.size();
  for (int c = 0; c < nchains; ++c) {
    IntegerVector O_k = ks_list[c];
    IntegerVector O_n = ns_list[c];
    const int Tc = O_k.size();

    NumericMatrix lp  = rtiger_getlogpsi_cpp(O_k, O_n, alpha, beta);
    NumericMatrix LP  = rtiger_productpsi_cpp(lp, r);
    NumericMatrix al  = rtiger_forward_cpp(logPI, LP, logA, lp, r);
    NumericMatrix be  = rtiger_backward_cpp(LP, logA, lp, r);
    NumericVector z   = rtiger_zeta_cpp(al, be, logA, LP, lp, r);   // T x s x s flat
    NumericMatrix gam = rtiger_gamma_cpp(z, al, be, r);            // s x T

    // transition: Σ over band (r+1):(Tc-r+1) of zeta[t, i, j]
    for (int i = 0; i < s; ++i) {
      for (int j = 0; j < s; ++j) {
        double acc = 0.0;
        for (int t = r + 1; t <= Tc - r + 1; ++t) acc += z[(t - 1) + i * Tc + j * (Tc * s)];
        sumZeta(i, j) += acc;
      }
    }
    // start: Σ gamma[,1]
    for (int st = 0; st < s; ++st) startAcc[st] += gam(st, 0);
    nOb += 1;
    // emission: per marker, accumulate per-state weights by distinct (k,n)
    for (int t = 0; t < Tc; ++t) {
      const int kk = O_k[t], nn = O_n[t];
      const long long key = (long long)nn * 2147483647LL + kk;
      std::vector<double>& w = W[key];
      if (w.empty()) w.assign(s, 0.0);
      for (int st = 0; st < s; ++st) {
        const double g = gam(st, t);
        w[st]   += g;
        sumk[st] += kk * g;
        sumn[st] += nn * g;
      }
    }
  }

  const int np = (int)W.size();
  IntegerVector kvals(np), nvals(np);
  NumericMatrix wmat(np, s);
  int p = 0;
  for (std::map<long long, std::vector<double> >::iterator it = W.begin(); it != W.end(); ++it) {
    nvals[p] = (int)(it->first / 2147483647LL);
    kvals[p] = (int)(it->first % 2147483647LL);
    for (int st = 0; st < s; ++st) wmat(p, st) = it->second[st];
    ++p;
  }
  return List::create(_["sumZeta"] = sumZeta, _["startAcc"] = startAcc, _["nOb"] = nOb,
                      _["kvals"] = kvals, _["nvals"] = nvals, _["wmat"] = wmat,
                      _["sumk"] = NumericVector(sumk.begin(), sumk.end()),
                      _["sumn"] = NumericVector(sumn.begin(), sumn.end()));
}


// Brent's method (1-D minimization on [ax,bx]) — the algorithm R's optimize and
// Julia's Optim.Brent both implement, so the emission M-step is now in C++ and
// matches the fork's Brent (not the R optimize approximation). Minimizes f.
static double brent_min(double ax, double bx, std::function<double(double)> f,
                        double tol = 1e-8, int max_it = 200) {
  const double gold = 0.3819660112501051;   // (3 - sqrt(5)) / 2
  double a = ax, b = bx;
  double x = a + gold * (b - a), w = x, v = x;
  double fx = f(x), fw = fx, fv = fx;
  double d = 0.0, e = 0.0;
  for (int it = 0; it < max_it; ++it) {
    double xm = 0.5 * (a + b);
    double tol1 = tol * std::fabs(x) + 1e-12, tol2 = 2.0 * tol1;
    if (std::fabs(x - xm) <= tol2 - 0.5 * (b - a)) break;
    bool use_golden = true;
    if (std::fabs(e) > tol1) {                 // try parabolic interpolation
      double r = (x - w) * (fx - fv);
      double q = (x - v) * (fx - fw);
      double p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0.0) p = -p;
      q = std::fabs(q);
      double etmp = e; e = d;
      if (std::fabs(p) < std::fabs(0.5 * q * etmp) && p > q * (a - x) && p < q * (b - x)) {
        d = p / q; double u = x + d;
        if (u - a < tol2 || b - u < tol2) d = (xm - x >= 0 ? tol1 : -tol1);
        use_golden = false;
      }
    }
    if (use_golden) { e = (x >= xm ? a - x : b - x); d = gold * e; }
    double u = (std::fabs(d) >= tol1) ? x + d : x + (d >= 0 ? tol1 : -tol1);
    double fu = f(u);
    if (fu <= fx) {
      if (u >= x) a = x; else b = x;
      v = w; w = x; x = u; fv = fw; fw = fx; fx = fu;
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w == x) { v = w; w = u; fv = fw; fw = fu; }
      else if (fu <= fv || v == x || v == w) { v = u; fv = fu; }
    }
  }
  return x;
}

// Per-state emission M-step (port of emissionUpdateState), Brent in C++.
static void emission_update_state(int i0, const std::vector<int>& ks,
                                  const std::vector<int>& ns, const std::vector<double>& ws,
                                  double sumk, double sumn,
                                  const NumericVector& alpha_old, const NumericVector& beta_old,
                                  double& a_out, double& b_out) {
  double mi = sumk / sumn;
  if (!R_FINITE(mi)) {                                   // sumn == 0: state unsupported
    a_out = alpha_old[i0]; b_out = beta_old[i0]; return;
  }
  if (mi < 0.01) mi = 0.01;
  if (mi > 0.99) mi = 0.99;
  double tau_i = alpha_old[i0] / mi + beta_old[i0] / (1.0 - mi);
  if (tau_i > 100.0) tau_i = 100.0;
  auto negQ = [&](double t) -> double {
    double a = std::max(t * mi, 1e-6), b = std::max(t * (1.0 - mi), 1e-6);
    double lbeta_ab = R::lbeta(a, b), acc = 0.0;
    for (size_t p = 0; p < ws.size(); ++p)
      acc += ws[p] * (R::lchoose((double)ns[p], (double)ks[p]) + R::lbeta(ks[p] + a, ns[p] - ks[p] + b) - lbeta_ab);
    return -acc;
  };
  double tau = brent_min(std::max(1e-6, tau_i - 100.0), std::max(tau_i + 1.0, 100.0), negQ);
  a_out = tau * mi; b_out = tau * (1.0 - mi);
}

// ---------------------------------------------------------------------------
// In-place buffered E-step kernels (fork optimization #4/#5: ~5x fwd/bwd,
// ~12x viterbi constant factor). Each writes into caller-provided buffers
// reused across all chains/iterations, eliminating per-call allocation. Same
// arithmetic as the verified rtiger_*_cpp kernels. Column-major strides: an
// s x T matrix indexes (i,t)=i+t*s; the zeta buffer is Tmax x s x s and indexes
// (t,j,k)=t + j*Tmax + k*Tmax*s (Tmax stride, so one buffer fits every chain).
// ---------------------------------------------------------------------------
static void getlogpsi_buf(std::vector<double>& lp, const int* O_k,
                          const int* O_n, const NumericVector& a,
                          const NumericVector& b, int T, int s) {
  std::vector<double> lbeta_ab(s);
  for (int i = 0; i < s; ++i) lbeta_ab[i] = R::lbeta(a[i], b[i]);
  std::unordered_map<long long, int> seen;
  for (int t = 0; t < T; ++t) {
    const int kt = O_k[t], nt = O_n[t];
    const long long key = (long long)nt * 2147483647LL + kt;
    std::unordered_map<long long, int>::iterator it = seen.find(key);
    if (it != seen.end()) {
      const int src = it->second;
      for (int i = 0; i < s; ++i) lp[i + t * s] = lp[i + src * s];
      continue;
    }
    const double lc = R::lchoose((double)nt, (double)kt);
    for (int i = 0; i < s; ++i) lp[i + t * s] = lc + R::lbeta(kt + a[i], nt - kt + b[i]) - lbeta_ab[i];
    seen[key] = t;
  }
}

static void productpsi_buf(std::vector<double>& PSI, const std::vector<double>& lp, int T, int r, int s) {
  for (int i = 0; i < s; ++i) {
    for (int t = 0; t < T + r; ++t) PSI[i + t * s] = 0.0;
    PSI[i] = lp[i];
    for (int t = 1; t < r; ++t)         PSI[i + t * s] = PSI[i + (t - 1) * s] + lp[i + t * s];
    for (int t = r; t < T; ++t)         PSI[i + t * s] = PSI[i + (t - 1) * s] + lp[i + t * s] - lp[i + (t - r) * s];
    for (int t = T; t < T + r - 1; ++t) PSI[i + t * s] = PSI[i + (t - 1) * s] - lp[i + (t - r) * s];
  }
}

static void forward_buf(std::vector<double>& al, int T, int r, int s, const NumericVector& logPI,
                        const std::vector<double>& PSI, const NumericMatrix& logA, const std::vector<double>& lp) {
  for (int t = 0; t < T; ++t) for (int k = 0; k < s; ++k) al[k + t * s] = R_NegInf;
  for (int k = 1; k <= s; ++k) al[(k - 1) + (r - 1) * s] = logPI[k - 1] + PSI[(k - 1) + (r - 1) * s];
  for (int t = r + 1; t <= 2 * r - 1; ++t) for (int k = 1; k <= s; ++k)
    al[(k - 1) + (t - 1) * s] = lp[(k - 1) + (t - 1) * s] + logA(k - 1, k - 1) + al[(k - 1) + (t - 2) * s];
  for (int t = 2 * r; t <= T - r + 1; ++t) for (int k = 1; k <= s; ++k) {
    const double stay = lp[(k - 1) + (t - 1) * s] + logA(k - 1, k - 1) + al[(k - 1) + (t - 2) * s];
    double maxv = R_NegInf; int amax = 0;
    for (int i = 1; i <= s; ++i) if (i != k) { double x = logA(i - 1, k - 1) + al[(i - 1) + (t - r - 1) * s]; if (x > maxv) { maxv = x; amax = i; } }
    double enter = R_NegInf;
    if (maxv != R_NegInf) {
      double acc = 0.0;
      for (int i = 1; i <= s; ++i) if (i != k && i != amax) acc += std::exp(logA(i - 1, k - 1) + al[(i - 1) + (t - r - 1) * s] - maxv);
      enter = PSI[(k - 1) + (t - 1) * s] + maxv + std::log1p(acc);
    }
    al[(k - 1) + (t - 1) * s] = logaddexp2(stay, enter);
  }
}

static void backward_buf(std::vector<double>& be, int T, int r, int s, const std::vector<double>& PSI,
                         const NumericMatrix& logA, const std::vector<double>& lp) {
  for (int t = 0; t < T; ++t) for (int j = 0; j < s; ++j) be[j + t * s] = R_NegInf;
  for (int c = T - r + 1; c <= T; ++c) for (int j = 1; j <= s; ++j) be[(j - 1) + (c - 1) * s] = PSI[(j - 1) + (c + r - 1) * s];
  for (int i = 0; i <= T - 2 * r; ++i) {
    const int t = T - r - i;
    for (int j = 1; j <= s; ++j) {
      const double stay = lp[(j - 1) + t * s] + logA(j - 1, j - 1) + be[(j - 1) + t * s];
      double maxv = R_NegInf; int amax = 0;
      for (int k = 1; k <= s; ++k) if (k != j) { double x = logA(j - 1, k - 1) + PSI[(k - 1) + (t + r - 1) * s] + be[(k - 1) + (t + r - 1) * s]; if (x > maxv) { maxv = x; amax = k; } }
      double leave = R_NegInf;
      if (maxv != R_NegInf) {
        double acc = 0.0;
        for (int k = 1; k <= s; ++k) if (k != j && k != amax) acc += std::exp(logA(j - 1, k - 1) + PSI[(k - 1) + (t + r - 1) * s] + be[(k - 1) + (t + r - 1) * s] - maxv);
        leave = maxv + std::log1p(acc);
      }
      be[(j - 1) + (t - 1) * s] = logaddexp2(stay, leave);
    }
  }
}

static void zeta_buf(std::vector<double>& z, const std::vector<double>& al, const std::vector<double>& be,
                     const NumericMatrix& logA, const std::vector<double>& PSI, const std::vector<double>& lp,
                     int r, int T, int s, int Tmax) {
  const long long sT = (long long)Tmax * s;
  for (int k = 0; k < s; ++k) for (int j = 0; j < s; ++j) for (int t = 0; t < T; ++t) z[t + j * Tmax + k * sT] = R_NegInf;
  for (int k = 1; k <= s; ++k) for (int j = 1; j <= s; ++j) {
    if (j != k) {
      for (int t = r + 1; t <= T - r + 1; ++t)
        z[(t - 1) + (j - 1) * Tmax + (k - 1) * sT] = al[(j - 1) + (t - 2) * s] + logA(j - 1, k - 1) + be[(k - 1) + (t + r - 2) * s] + PSI[(k - 1) + (t + r - 2) * s];
    } else {
      for (int t = r + 1; t <= T - r + 1; ++t)
        z[(t - 1) + (k - 1) * Tmax + (k - 1) * sT] = al[(k - 1) + (t - 2) * s] + logA(k - 1, k - 1) + lp[(k - 1) + (t - 1) * s] + be[(k - 1) + (t - 1) * s];
    }
  }
  double PO = R_NegInf;
  for (int k = 0; k < s; ++k) for (int j = 0; j < s; ++j) for (int t = 0; t < T; ++t) { double v = z[t + j * Tmax + k * sT]; if (v > PO) PO = v; }
  for (int k = 0; k < s; ++k) for (int j = 0; j < s; ++j) for (int t = 0; t < T; ++t) { long long ix = t + j * Tmax + k * sT; z[ix] = std::exp(z[ix] - PO); }
}

static void gamma_buf(std::vector<double>& gam, const std::vector<double>& z, const std::vector<double>& al,
                      const std::vector<double>& be, int r, int T, int s, int Tmax,
                      std::vector<double>& L, std::vector<double>& M) {
  const long long sT = (long long)Tmax * s;
  std::vector<double> gammar(s); double mx = R_NegInf;
  for (int k = 0; k < s; ++k) { gammar[k] = al[k + (r - 1) * s] + be[k + (r - 1) * s]; if (gammar[k] > mx) mx = gammar[k]; }
  double gsum = 0.0; for (int k = 0; k < s; ++k) { gammar[k] = std::exp(gammar[k] - mx); gsum += gammar[k]; }
  for (int k = 0; k < s; ++k) gammar[k] /= gsum;
  for (int t = 1; t <= r; ++t) for (int k = 0; k < s; ++k) gam[k + (t - 1) * s] = gammar[k];
  for (int k = 1; k <= s; ++k) {
    for (int t = 1; t <= T; ++t) { double acc = 0.0; for (int j = 1; j <= s; ++j) if (j != k) acc += z[(t - 1) + (j - 1) * Tmax + (k - 1) * sT]; L[t - 1] = acc; }
    double run = 0.0; for (int t = 1; t <= T; ++t) { run += L[t - 1]; M[t - 1] = run; }
    for (int t = r + 1; t <= T - r + 1; ++t) { double Mm = (t <= 2 * r) ? M[r - 1] : M[t - r - 1]; gam[(k - 1) + (t - 1) * s] = z[(t - 1) + (k - 1) * Tmax + (k - 1) * sT] + M[t - 1] - Mm; }
  }
  for (int t = T - r + 2; t <= T; ++t) for (int k = 0; k < s; ++k) gam[k + (t - 1) * s] = gam[k + (T - r) * s];
  for (int t = 0; t < T; ++t) { double cs = 0.0; for (int k = 0; k < s; ++k) cs += gam[k + t * s]; for (int k = 0; k < s; ++k) gam[k + t * s] /= cs; }
}

// Parallel E-step worker: each chunk processes a contiguous range of chains
// into its OWN partial accumulators (parts[ch]) with its OWN buffer set, so
// there are no races. The serial ordered reduce of parts[0..K-1] then makes the
// result deterministic for a fixed thread count (Viterbi-identical to serial;
// parameters differ only by float-summation order). Read-only Rcpp params
// (logPI/logA/alpha/beta) are shared — concurrent reads are safe.
struct EStepChunk : public Worker {
  const std::vector<std::vector<int> >& KS;
  const std::vector<std::vector<int> >& NS;
  const NumericVector& logPI; const NumericMatrix& logA;
  const NumericVector& alpha; const NumericVector& beta;
  const std::vector<int>& edges;          // chunk boundaries, size nchunks+1
  int r, s, Tmax;
  std::vector<std::vector<double> >& sumZeta_p;   // [ch] -> s*s (col-major)
  std::vector<std::vector<double> >& startAcc_p;  // [ch] -> s
  std::vector<std::vector<double> >& sumk_p;      // [ch] -> s
  std::vector<std::vector<double> >& sumn_p;      // [ch] -> s
  std::vector<std::map<long long, std::vector<double> > >& W_p;

  EStepChunk(const std::vector<std::vector<int> >& KS_, const std::vector<std::vector<int> >& NS_,
             const NumericVector& logPI_, const NumericMatrix& logA_,
             const NumericVector& alpha_, const NumericVector& beta_,
             const std::vector<int>& edges_, int r_, int s_, int Tmax_,
             std::vector<std::vector<double> >& sZ, std::vector<std::vector<double> >& sA,
             std::vector<std::vector<double> >& sk, std::vector<std::vector<double> >& sn,
             std::vector<std::map<long long, std::vector<double> > >& Wp)
    : KS(KS_), NS(NS_), logPI(logPI_), logA(logA_), alpha(alpha_), beta(beta_),
      edges(edges_), r(r_), s(s_), Tmax(Tmax_),
      sumZeta_p(sZ), startAcc_p(sA), sumk_p(sk), sumn_p(sn), W_p(Wp) {}

  void operator()(std::size_t begin, std::size_t end) {
    const size_t sT = (size_t)Tmax * s;
    std::vector<double> lp(s * Tmax), PSI(s * (Tmax + r)), al(s * Tmax),
                        be(s * Tmax), gam(s * Tmax), L(Tmax), M(Tmax), z((size_t)Tmax * s * s);
    for (std::size_t ch = begin; ch < end; ++ch) {
      std::vector<double>& sZ = sumZeta_p[ch]; std::vector<double>& sA = startAcc_p[ch];
      std::vector<double>& sk = sumk_p[ch];    std::vector<double>& sn = sumn_p[ch];
      std::map<long long, std::vector<double> >& W = W_p[ch];
      for (int c = edges[ch]; c < edges[ch + 1]; ++c) {
        const int* O_k = KS[c].data(); const int* O_n = NS[c].data();
        const int Tc = (int)KS[c].size();
        getlogpsi_buf(lp, O_k, O_n, alpha, beta, Tc, s);
        productpsi_buf(PSI, lp, Tc, r, s);
        forward_buf(al, Tc, r, s, logPI, PSI, logA, lp);
        backward_buf(be, Tc, r, s, PSI, logA, lp);
        zeta_buf(z, al, be, logA, PSI, lp, r, Tc, s, Tmax);
        gamma_buf(gam, z, al, be, r, Tc, s, Tmax, L, M);
        for (int i = 0; i < s; ++i) for (int j = 0; j < s; ++j) {
          double acc = 0.0;
          for (int t = r + 1; t <= Tc - r + 1; ++t) acc += z[(t - 1) + i * Tmax + j * sT];
          sZ[i + j * s] += acc;
        }
        for (int st = 0; st < s; ++st) sA[st] += gam[st];
        for (int t = 0; t < Tc; ++t) {
          const int kk = O_k[t], nn = O_n[t];
          std::vector<double>& w = W[(long long)nn * 2147483647LL + kk];
          if (w.empty()) w.assign(s, 0.0);
          for (int st = 0; st < s; ++st) { double g = gam[st + t * s]; w[st]+=g; sk[st]+=kk*g; sn[st]+=nn*g; }
        }
      }
    }
  }
};

//' RTIGER full EM fit (E-step + M-step + convergence loop, all in C++)
//'
//' The entire fit (port of the fork's `fit`/`EM`): per-chain rigidity E-step,
//' pooled M-steps (transition, start, emission via C++ Brent), iterated until
//' max(|Δα|,|Δβ|) <= eps or max.iter. Deterministic init (generate_params forms,
//' randomize off). Returns the fitted parameters and iteration count.
//' @param ks_list,ns_list Lists of per-chain integer (k, n) vectors.
//' @param r,nstates Rigidity and number of states.
//' @param eps,max_iter Convergence tolerance and iteration cap.
//' @return list(A, pi, alpha, beta, iterations).
//' @keywords internal
// [[Rcpp::export]]
List rtiger_fit_cpp(List ks_list, List ns_list, int r, int nstates,
                    double eps, int max_iter, int threads,
                    NumericVector init_alpha, NumericVector init_beta) {
  const int s = nstates;
  const int nchains = ks_list.size();

  // deterministic init (generate_params forms, randomize off)
  NumericMatrix A(s, s);
  for (int i = 0; i < s; ++i) { for (int j = 0; j < s; ++j) A(i, j) = 0.1; A(i, i) += 10.0; }
  for (int i = 0; i < s; ++i) { double rs = 0; for (int j = 0; j < s; ++j) rs += A(i, j); for (int j = 0; j < s; ++j) A(i, j) /= rs; }
  NumericVector PI(s, 1.0 / s);
  NumericVector alpha(s), beta(s);
  if (init_alpha.size() == s && init_beta.size() == s) {              // caller-supplied init
    for (int i=0;i<s;++i){ alpha[i]=init_alpha[i]; beta[i]=init_beta[i]; }
  } else if (s == 3) { alpha[0]=20; alpha[1]=20; alpha[2]=1; beta[0]=1; beta[1]=20; beta[2]=20; }
  else { for (int i=0;i<s;++i){ alpha[i]=20; beta[i]=20; } }

  // Extract chains to plain C++ ONCE (thread-safe + avoids per-iter List re-wrap).
  std::vector<std::vector<int> > KS(nchains), NS(nchains);
  int Tmax = 0;
  for (int c = 0; c < nchains; ++c) {
    KS[c] = as<std::vector<int> >(ks_list[c]);
    NS[c] = as<std::vector<int> >(ns_list[c]);
    if ((int)KS[c].size() > Tmax) Tmax = (int)KS[c].size();
  }
  // Hard stop on under-covered chains, as the RTIGER (Julia) fit does: a chain
  // shorter than 2*rigidity has an empty zeta window and its E-step degenerates
  // to NaN. RTIGER aborts rather than fit such data; do the same, up front on the
  // main thread (throwing from inside the parallel E-step would be unsafe). The
  // caller must drop low-coverage (sample, chromosome) chains first.
  for (int c = 0; c < nchains; ++c)
    if ((int)KS[c].size() < 2 * r)
      Rcpp::stop("rtiger: chain %d has %d covered markers, below the 2*rigidity = %d "
                 "floor (RTIGER requires >= 2*rigidity per chromosome). Drop "
                 "low-coverage samples/chromosomes before fitting.",
                 c, (int)KS[c].size(), 2 * r);
  // Contiguous chunking for the parallel E-step: K = min(threads, nchains)
  // chunks, each folded into its own partial, reduced in chunk order (so the
  // result is deterministic for a fixed thread count). threads=1 -> 1 chunk =
  // serial order = bit-identical to the single-threaded path.
  const int nthr = std::max(1, threads);
  const int nchunks = std::max(1, std::min(nthr, nchains));
  std::vector<int> edges(nchunks + 1);
  for (int ch = 0; ch <= nchunks; ++ch) edges[ch] = (int)((long long)ch * nchains / nchunks);

  int iter = 0;
  for (;;) {
    NumericVector logPI(s); for (int k=0;k<s;++k) logPI[k] = std::log(PI[k]);
    NumericMatrix logA(s, s); for (int i=0;i<s;++i) for (int j=0;j<s;++j) logA(i,j)=std::log(A(i,j));

    // ---- E-step: parallel per-chunk fold, then ordered serial reduce ----
    std::vector<std::vector<double> > sumZeta_p(nchunks, std::vector<double>(s * s, 0.0));
    std::vector<std::vector<double> > startAcc_p(nchunks, std::vector<double>(s, 0.0));
    std::vector<std::vector<double> > sumk_p(nchunks, std::vector<double>(s, 0.0));
    std::vector<std::vector<double> > sumn_p(nchunks, std::vector<double>(s, 0.0));
    std::vector<std::map<long long, std::vector<double> > > W_p(nchunks);
    EStepChunk worker(KS, NS, logPI, logA, alpha, beta, edges, r, s, Tmax,
                      sumZeta_p, startAcc_p, sumk_p, sumn_p, W_p);
    // nchunks == 1 (threads <= 1) runs the worker DIRECTLY, without TBB. Going
    // through parallelFor even for a single chunk was nondeterministic across
    // repeated calls: it did not reliably invoke the worker over the whole range,
    // leaving the pooled sufficient statistics empty -> NaN emission M-step ->
    // silent return of the init as a "converged" fit. A direct call is
    // deterministic. threads > 1 keeps the per-chunk parallel fold.
    if (nchunks <= 1) worker(0, (std::size_t)nchunks);
    else parallelFor(0, nchunks, worker, 1);         // grain 1 -> each chunk independent

    // reduce partials in fixed chunk order (deterministic for fixed thread count)
    NumericMatrix sumZeta(s, s);
    std::vector<double> startAcc(s, 0.0), sumk(s, 0.0), sumn(s, 0.0);
    std::map<long long, std::vector<double> > W;
    int nOb = nchains;
    for (int ch = 0; ch < nchunks; ++ch) {
      for (int i = 0; i < s; ++i) for (int j = 0; j < s; ++j) sumZeta(i, j) += sumZeta_p[ch][i + j * s];
      for (int st = 0; st < s; ++st) { startAcc[st] += startAcc_p[ch][st]; sumk[st] += sumk_p[ch][st]; sumn[st] += sumn_p[ch][st]; }
      for (std::map<long long, std::vector<double> >::iterator it = W_p[ch].begin(); it != W_p[ch].end(); ++it) {
        std::vector<double>& w = W[it->first];
        if (w.empty()) w.assign(s, 0.0);
        for (int st = 0; st < s; ++st) w[st] += it->second[st];
      }
    }

    // Guard: a non-empty E-step must assign some emission weight. If the pooled
    // total is zero the E-step produced nothing (e.g. it never ran, or the input
    // has no covered markers) — every state's mi = sumk/sumn would be NaN and the
    // emission M-step would silently keep the init, reporting a bogus 1-iteration
    // "convergence". Fail loudly instead.
    double totw = 0.0; for (int st = 0; st < s; ++st) totw += sumn[st];
    if (!(totw > 0.0))
      Rcpp::stop("rtiger_fit_cpp: E-step produced zero total emission weight over %d "
                 "chains (degenerate/empty fit); check that the input has covered markers.",
                 nchains);

    // ---- M-step ----
    NumericMatrix Anew(s, s);
    for (int i = 0; i < s; ++i) {
      double rs = 0; for (int j = 0; j < s; ++j) rs += sumZeta(i, j);
      if (rs == 0) rs = 1;
      for (int j = 0; j < s; ++j) Anew(i, j) = sumZeta(i, j) / rs;
    }
    NumericVector PInew(s); for (int st = 0; st < s; ++st) PInew[st] = startAcc[st] / nOb;
    // flatten W to per-state (k,n,w) and update emission
    std::vector<int> ks, ns; ks.reserve(W.size()); ns.reserve(W.size());
    for (std::map<long long, std::vector<double> >::iterator it = W.begin(); it != W.end(); ++it) {
      ns.push_back((int)(it->first / 2147483647LL)); ks.push_back((int)(it->first % 2147483647LL));
    }
    NumericVector anew(s), bnew(s);
    for (int st = 0; st < s; ++st) {
      std::vector<double> ws; ws.reserve(W.size());
      for (std::map<long long, std::vector<double> >::iterator it = W.begin(); it != W.end(); ++it) ws.push_back(it->second[st]);
      double aa, bb; emission_update_state(st, ks, ns, ws, sumk[st], sumn[st], alpha, beta, aa, bb);
      anew[st] = aa; bnew[st] = bb;
    }

    double er = 0.0;
    for (int st = 0; st < s; ++st) { er = std::max(er, std::fabs(alpha[st]-anew[st])); er = std::max(er, std::fabs(beta[st]-bnew[st])); }
    A = Anew; PI = PInew; alpha = anew; beta = bnew;
    iter += 1;
    if (std::floor(er * 1e6 + 0.5) / 1e6 <= eps || iter >= max_iter) break;
  }
  return List::create(_["A"]=A, _["pi"]=PI, _["alpha"]=alpha, _["beta"]=beta, _["iterations"]=iter);
}
