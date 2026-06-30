// Run-length-encode a decoded state path into segment boundaries (REFACTOR S8).
// Doing this in C++ avoids the R-level diff()/which()/negative-index allocations
// that dominate call_ancestry wall-clock when looping over many (sample, chr)
// sequences. Boundaries where the state changes; start_bp/end_bp = first/last
// marker position of each contiguous run.

#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

//' Run-length-encode a state path into (start_bp, end_bp, state) segments
//'
//' @param path Integer state path (0/1/2), length T.
//' @param pos Integer marker positions, length T (same order as path).
//' @return List of equal-length integer vectors: start_bp, end_bp, state.
//' @keywords internal
// [[Rcpp::export]]
List rle_segments_cpp(IntegerVector path, IntegerVector pos) {
  const int n = path.size();
  std::vector<int> sbp, ebp, st;
  if (n > 0) {
    int run_start = 0;
    for (int t = 1; t < n; ++t) {
      if (path[t] != path[t - 1]) {
        sbp.push_back(pos[run_start]); ebp.push_back(pos[t - 1]); st.push_back(path[run_start]);
        run_start = t;
      }
    }
    sbp.push_back(pos[run_start]); ebp.push_back(pos[n - 1]); st.push_back(path[run_start]);
  }
  return List::create(_["start_bp"] = sbp, _["end_bp"] = ebp, _["state"] = st);
}

//' Run-length-encode a batch of state paths (one per column) into segments
//'
//' @param paths T x S integer matrix of state paths (column = sample).
//' @param pos Integer marker positions, length T (shared across samples).
//' @return List of equal-length vectors: sample (1-based column), start_bp,
//'   end_bp, state.
//' @keywords internal
// [[Rcpp::export]]
List rle_segments_batch_cpp(IntegerMatrix paths, IntegerVector pos) {
  const int T = paths.nrow(), S = paths.ncol();
  std::vector<int> col, sbp, ebp, st;
  for (int s = 0; s < S; ++s) {
    if (T == 0) continue;
    int run = 0;
    for (int t = 1; t < T; ++t) {
      if (paths(t, s) != paths(t - 1, s)) {
        col.push_back(s + 1); sbp.push_back(pos[run]); ebp.push_back(pos[t - 1]); st.push_back(paths(run, s));
        run = t;
      }
    }
    col.push_back(s + 1); sbp.push_back(pos[run]); ebp.push_back(pos[T - 1]); st.push_back(paths(run, s));
  }
  return List::create(_["sample"] = col, _["start_bp"] = sbp, _["end_bp"] = ebp, _["state"] = st);
}
