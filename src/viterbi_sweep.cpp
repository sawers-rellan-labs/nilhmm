// Unified batched Viterbi sweep over a transition grid with a fixed emission
// (Scope A of the sweep-kernel consolidation; design/TODO_viterbi_kernel_consolidation.md).
//
// One entry point, viterbi_sweep_cpp(), with two transition modes:
//   "const"    -- time-homogeneous: one S x S log-transition matrix per grid value
//                 (geometric / switch-rate callers: nnil, bbnil, ...). General S.
//   "distance" -- LB-Impute's per-gap distance transition (+ double-recombination
//                 penalty drp), 3 states. Byte-identical to lb_viterbi_sweep_cpp
//                 (that kernel is retired in favour of this one; parity-gated).
//
// The emission axis of a full parameter grid is driven at the R level
// (caller_grid(): build log_emit once per emission combo, one sweep call per
// combo), so this kernel only ever sweeps the TRANSITION. All grid values reuse
// the delta/psi buffers; only the transition is rebuilt per value (const) or per
// gap (distance).

#include <Rcpp.h>
#include <cmath>
#include <string>
#include <vector>
using namespace Rcpp;

//' Batched Viterbi over a transition grid, shared (fixed) emission
//'
//' @param log_init Length-S log initial-state probabilities.
//' @param log_emit T x S log emissions (fixed across the transition grid).
//' @param mode `"const"` (time-homogeneous; one S x S matrix per grid value via
//'   `trans_list`) or `"distance"` (LB-Impute per-gap transition; 3 states, via
//'   `tpos` + `recombdists` + `drp`).
//' @param trans_list `"const"` mode: a length-V list of S x S log-transition
//'   matrices (row = from, col = to), one per grid value.
//' @param tpos,recombdists,drp `"distance"` mode: non-decreasing length-T
//'   coordinate, the length-V recombdist grid (same units), and the
//'   double-recombination penalty flag. Ignored in `"const"` mode.
//' @param tie_break Transition-backpointer tie policy, `"const"` mode only (see
//'   [viterbi_log_cpp()]): `0` (default) keeps the first (lowest-index)
//'   predecessor; `1` ("incumbent") keeps the last, needed to match the GT/nnil
//'   decode where emission-degenerate positions (missing / het) tie. `"distance"`
//'   mode is always strict-first (byte-identical to `lb_viterbi_sweep_cpp`).
//' @return A T x V integer matrix of 0-indexed state paths; column k is the decode
//'   at grid value k. Terminal argmax always keeps the first max (as
//'   `numpy.argmax`).
//' @keywords internal
// [[Rcpp::export]]
IntegerMatrix viterbi_sweep_cpp(NumericVector log_init, NumericMatrix log_emit,
                                std::string mode,
                                Nullable<List> trans_list = R_NilValue,
                                Nullable<NumericVector> tpos = R_NilValue,
                                Nullable<NumericVector> recombdists = R_NilValue,
                                bool drp = false,
                                int tie_break = 0) {
  const int T = log_emit.nrow();
  const int S = log_emit.ncol();
  if (log_init.size() != S)
    stop("viterbi_sweep_cpp: length(log_init) must equal ncol(log_emit)");

  // ---- distance mode: byte-identical to lb_viterbi_sweep_cpp (3-state) -------
  if (mode == "distance") {
    if (S != 3)
      stop("viterbi_sweep_cpp: mode='distance' requires a 3-state (REF/HET/ALT) emission");
    if (tpos.isNull() || recombdists.isNull())
      stop("viterbi_sweep_cpp: mode='distance' needs `tpos` and `recombdists`");
    NumericVector tp(tpos.get());
    NumericVector rd(recombdists.get());
    const int V = rd.size();
    if (tp.size() != T)
      stop("viterbi_sweep_cpp: length(tpos) must equal nrow(log_emit)");
    IntegerMatrix paths(T, V);
    if (T == 0 || V == 0) return paths;
    NumericMatrix delta(T, 3);
    IntegerMatrix psi(T, 3);
    double lt[3][3];
    for (int vi = 0; vi < V; ++vi) {
      const double recombdist = rd[vi];
      if (recombdist <= 0.0) stop("viterbi_sweep_cpp: recombdist must be > 0");
      for (int k = 0; k < 3; ++k) { delta(0, k) = log_init[k] + log_emit(0, k); psi(0, k) = 0; }
      for (int t = 1; t < T; ++t) {
        double d = tp[t] - tp[t - 1];
        if (d < 0) stop("viterbi_sweep_cpp: `tpos` must be non-decreasing");
        const double e   = std::exp(-d / recombdist);
        const double lps = std::log(0.5 * (1.0 + e));
        const double lpr = std::log(0.5 * (1.0 - e));
        const double lhh = drp ? lpr : (2.0 * lpr);
        lt[0][0] = lps; lt[0][1] = lpr; lt[0][2] = lhh;   // from REF
        lt[1][0] = lpr; lt[1][1] = lps; lt[1][2] = lpr;   // from HET
        lt[2][0] = lhh; lt[2][1] = lpr; lt[2][2] = lps;   // from ALT
        for (int k = 0; k < 3; ++k) {
          double best = R_NegInf; int arg = 0;
          for (int j = 0; j < 3; ++j) {
            double v = delta(t - 1, j) + lt[j][k];
            if (v > best) { best = v; arg = j; }
          }
          delta(t, k) = best + log_emit(t, k);
          psi(t, k) = arg;
        }
      }
      double best = R_NegInf; int last = 0;
      for (int k = 0; k < 3; ++k) if (delta(T - 1, k) > best) { best = delta(T - 1, k); last = k; }
      paths(T - 1, vi) = last;
      for (int t = T - 1; t > 0; --t) paths(t - 1, vi) = psi(t, paths(t, vi));
    }
    return paths;
  }

  // ---- const mode: general S, one time-homogeneous matrix per grid value -----
  if (mode == "const") {
    if (trans_list.isNull())
      stop("viterbi_sweep_cpp: mode='const' needs `trans_list`");
    List tl(trans_list.get());
    const int V = tl.size();
    IntegerMatrix paths(T, V);
    if (T == 0 || V == 0) return paths;
    std::vector<double> dprev(S), dcur(S);
    IntegerMatrix psi(T, S);
    const bool incumbent = (tie_break == 1); // keep last predecessor on ties (GT/nnil)
    for (int vi = 0; vi < V; ++vi) {
      NumericMatrix lt = tl[vi];
      if (lt.nrow() != S || lt.ncol() != S)
        stop("viterbi_sweep_cpp: each trans_list matrix must be S x S (S = ncol(log_emit))");
      for (int k = 0; k < S; ++k) { dprev[k] = log_init[k] + log_emit(0, k); psi(0, k) = 0; }
      for (int t = 1; t < T; ++t) {
        for (int k = 0; k < S; ++k) {
          double best = R_NegInf; int arg = 0;
          for (int j = 0; j < S; ++j) {
            double v = dprev[j] + lt(j, k);
            if (v > best || (incumbent && v == best)) { best = v; arg = j; }
          }
          dcur[k] = best + log_emit(t, k);
          psi(t, k) = arg;
        }
        dprev.swap(dcur);
      }
      double best = R_NegInf; int last = 0;
      for (int k = 0; k < S; ++k) if (dprev[k] > best) { best = dprev[k]; last = k; }
      paths(T - 1, vi) = last;
      for (int t = T - 1; t > 0; --t) paths(t - 1, vi) = psi(t, paths(t, vi));
    }
    return paths;
  }

  stop("viterbi_sweep_cpp: mode must be 'const' or 'distance'");
}
