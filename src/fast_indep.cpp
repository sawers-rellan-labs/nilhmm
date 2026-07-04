// Maximal independent set selection -- a faithful port of FastIndep
// (Joseph Abraham 2013; CLI options added by F. Rodriguez 2017,
// https://github.com/faustovrz/FastIndep). Given a symmetric matrix and a
// threshold, find a large independent set: a set of nodes in which no two are
// "related" (adjacent). Used for LD-based marker thinning before joint-linkage
// mapping -- prune a per-chromosome LD matrix to a set with all pairwise
// r^2 (or MI / VI) below/above a cutoff.
//
// Edge convention (matches FastIndep):
//   * similarity input (r2, mi; higher = more related): edge if value >= threshold;
//   * distance   input (vi;      lower  = more related): edge if value <= threshold.
// An independent set has NO edges among its members. The diagonal is ignored.
//
// Following FastIndep, each node i carries `seps[i]` = the ids of nodes it is
// NOT related to (the "separated" set = non-neighbours). The greedy heuristic
// repeatedly extends the set by the candidate of highest remaining degree of
// separation; the stochastic runs pick candidates with probability proportional
// to that degree. `singletons` (separated from everyone) join every set;
// `all_connectors` (related to everyone) can never join one.
//
// RNG FIDELITY: the greedy heuristic is deterministic and does NOT touch the
// RNG, so `fast_indep_cpp`'s greedy set is bit-identical to the FastIndep CLI.
// The stochastic runs are algorithmically faithful but use a self-contained
// splitmix64 generator (seeded by `seed`), NOT FastIndep's Mersenne Twister, so
// individual random sets are reproducible here but not bit-identical to the CLI.
// (Reproducing the CLI's exact random stream would require vendoring its
// LGPL MersenneTwister.h, whose seeding relies on platform `unsigned long`
// width; that license does not mix cleanly with this MIT package.)

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <map>
#include <cstdint>
using namespace Rcpp;

namespace {

// splitmix64: a tiny, self-contained, reproducible PRNG for the stochastic runs.
// Independent of R's global RNG state, so `seed` alone fixes the result.
struct RNG {
  uint64_t s;
  explicit RNG(uint64_t seed) : s(seed) {}
  uint64_t next_u64() {
    uint64_t z = (s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
  }
  double unif() {                         // uniform in [0, 1)
    return (next_u64() >> 11) * (1.0 / 9007199254740992.0);
  }
};

// Greedy start: the eligible node (not a singleton, not a connector) with the
// largest separation set; ties broken by lowest id (strict `>`, ascending scan),
// matching FastIndep::FindstartSNP.
int find_start(const std::vector<std::vector<int> >& seps,
               const std::vector<char>& is_singleton,
               const std::vector<char>& is_connector) {
  int best = 0;
  std::size_t largest = 0;
  for (std::size_t i = 0; i < seps.size(); ++i) {
    if (is_connector[i] || is_singleton[i]) continue;
    if (seps[i].size() > largest) { largest = seps[i].size(); best = (int)i + 1; }
  }
  return best;
}

// Greedy next: the candidate with the largest separation set; ties -> lowest id
// (candVec is kept ascending). Mirrors FastIndep::FindnextSNP.
int find_next(const std::vector<int>& cand,
              const std::vector<std::vector<int> >& seps) {
  int best = 0;
  std::size_t largest = 0;
  for (std::size_t i = 0; i < cand.size(); ++i) {
    std::size_t d = seps[cand[i] - 1].size();
    if (d > largest) { largest = d; best = cand[i]; }
  }
  return best;
}

// Stochastic pick over `pool` (node ids), probability proportional to each
// node's separation-set size -- FastIndep's FindstartSNP_rand / FindnextSNP_rand
// biased coin toss, with our splitmix64 draw in place of its Mersenne Twister.
int weighted_pick(const std::vector<int>& pool,
                  const std::vector<std::vector<int> >& seps,
                  RNG& rng) {
  double total = 0.0;
  for (std::size_t i = 0; i < pool.size(); ++i) total += (double)seps[pool[i] - 1].size();
  if (total <= 0.0) return 0;
  double u = rng.unif();
  double cum = 0.0;
  for (std::size_t i = 0; i < pool.size(); ++i) {
    double p = (double)seps[pool[i] - 1].size() / total;
    if (u < cum + p) return pool[i];
    cum += p;
  }
  return pool.back();                     // guard against fp round-off
}

// Sorted set-intersection of two ascending id vectors.
std::vector<int> intersect_sorted(const std::vector<int>& a, const std::vector<int>& b) {
  std::vector<int> out;
  out.reserve(std::min(a.size(), b.size()));
  std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(out));
  return out;
}

// One independent-set search (FastIndep::FindSeparatedVec). `rand_mode` selects
// the greedy vs stochastic candidate rule. Returns 1-based node ids, ascending.
std::vector<int> find_set(const std::vector<std::vector<int> >& seps,
                          const std::vector<int>& singletons,
                          const std::vector<char>& is_singleton,
                          const std::vector<char>& is_connector,
                          bool rand_mode, RNG& rng) {
  const int n = (int)seps.size();
  int start = rand_mode
    ? [&]() {                             // weighted start over eligible nodes
        std::vector<int> pool;
        for (int i = 0; i < n; ++i)
          if (!is_connector[i] && !is_singleton[i]) pool.push_back(i + 1);
        return pool.empty() ? 0 : weighted_pick(pool, seps, rng);
      }()
    : find_start(seps, is_singleton, is_connector);

  std::vector<int> chosen;

  // No eligible start: the graph is all singletons and/or one big clique. Return
  // the singletons if any (each is independent of everything); otherwise a single
  // node -- a complete graph's maximum independent set has size 1.
  if (start == 0) {
    if (!singletons.empty()) return singletons;             // already ascending
    chosen.push_back(1);
    return chosen;
  }

  chosen.push_back(start);
  std::vector<int> cand = seps[start - 1];                  // candidates = start's non-neighbours
  if (!singletons.empty()) {                                // singletons handled separately
    std::vector<int> tmp;
    std::set_difference(cand.begin(), cand.end(), singletons.begin(), singletons.end(),
                        std::back_inserter(tmp));
    cand.swap(tmp);
    chosen.insert(chosen.end(), singletons.begin(), singletons.end());
  }

  while (!cand.empty()) {
    int nxt = rand_mode ? weighted_pick(cand, seps, rng) : find_next(cand, seps);
    if (nxt == 0) break;
    chosen.push_back(nxt);
    cand.erase(std::find(cand.begin(), cand.end(), nxt));   // drop the picked node
    cand = intersect_sorted(cand, seps[nxt - 1]);           // keep only common non-neighbours
  }

  std::sort(chosen.begin(), chosen.end());
  return chosen;
}

}  // namespace

//' Find a large independent set in a thresholded similarity/distance matrix
//'
//' Faithful port of FastIndep (Abraham 2013): the deterministic greedy heuristic
//' plus `n_runs - 1` stochastic runs. Two nodes are "related" (an edge) when the
//' matrix entry crosses `threshold` in the sense given by `distance`; an
//' independent set has no edges among its members. The matrix diagonal is ignored.
//'
//' @param sim Symmetric numeric matrix. A similarity (higher = more related,
//'   `distance = false`) or a distance (lower = more related, `distance = true`).
//' @param threshold Edge cutoff. Similarity: edge if `sim[i,j] >= threshold`.
//'   Distance: edge if `sim[i,j] <= threshold`.
//' @param n_runs Total runs: 1 = greedy only; `> 1` adds `n_runs - 1` stochastic
//'   runs and returns the distinct sets found.
//' @param seed Seed for the stochastic runs (self-contained splitmix64 PRNG;
//'   independent of R's RNG). The greedy run is deterministic and seed-free.
//' @param distance If `true`, `sim` is a distance and the edge sense is inverted
//'   (edge if `<= threshold`). Default `false` (similarity).
//' @return A list: `best` (1-based indices of the largest independent set found,
//'   the greedy set when it is largest), `sets` (list of distinct sets, greedy
//'   first then by decreasing size), `size_dist` (named integer vector, set
//'   size -> number of distinct sets of that size).
//' @keywords internal
// [[Rcpp::export]]
List fast_indep_cpp(NumericMatrix sim, double threshold, int n_runs, int seed,
                    bool distance = false) {
  const int n = sim.nrow();
  if (sim.ncol() != n) stop("fast_indep_cpp(): `sim` must be square");
  if (n_runs < 1) stop("fast_indep_cpp(): `n_runs` must be >= 1");

  // Build separation sets: seps[i] = ascending ids (1-based) of nodes NOT related
  // to i under the edge sense. Flag singletons (separated from all: |seps| = n-1)
  // and connectors (related to all: |seps| = 0), which FastIndep treats specially.
  std::vector<std::vector<int> > seps(n);
  std::vector<char> is_singleton(n, 0), is_connector(n, 0);
  std::vector<int> singletons;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == j) continue;               // diagonal ignored
      double v = sim(i, j);
      bool edge = distance ? (v <= threshold) : (v >= threshold);
      if (!edge) seps[i].push_back(j + 1);
    }
    if (seps[i].empty()) is_connector[i] = 1;
    else if ((int)seps[i].size() == n - 1) { is_singleton[i] = 1; singletons.push_back(i + 1); }
  }
  std::sort(singletons.begin(), singletons.end());

  RNG rng((uint64_t)(unsigned int)seed);

  // Collect distinct sets. Stochastic runs first (mirrors the CLI order), then
  // the deterministic greedy set. Dedup on the sorted membership vector.
  std::map<std::vector<int>, bool> seen;
  std::vector<std::vector<int> > rand_sets;
  for (int r = 0; r < n_runs - 1; ++r) {
    std::vector<int> s = find_set(seps, singletons, is_singleton, is_connector, true, rng);
    if (!s.empty() && seen.find(s) == seen.end()) { seen[s] = true; rand_sets.push_back(s); }
  }
  std::vector<int> greedy =
    find_set(seps, singletons, is_singleton, is_connector, false, rng);
  bool greedy_is_new = !greedy.empty() && seen.find(greedy) == seen.end();

  // size_dist over the distinct sets (random uniques + greedy if new).
  std::map<int, int> sd;
  for (std::size_t i = 0; i < rand_sets.size(); ++i) sd[(int)rand_sets[i].size()]++;
  if (greedy_is_new) sd[(int)greedy.size()]++;

  // `sets`: greedy first, then the remaining distinct sets by decreasing size
  // (stable so equal sizes keep discovery order). When the greedy set was also
  // found by a stochastic run, drop that copy from `others` so it is not listed
  // twice (sets[0] is always the greedy set).
  std::vector<std::vector<int> > others = rand_sets;
  if (!greedy_is_new)
    others.erase(std::remove(others.begin(), others.end(), greedy), others.end());
  std::stable_sort(others.begin(), others.end(),
                   [](const std::vector<int>& a, const std::vector<int>& b) {
                     return a.size() > b.size();
                   });
  List sets(1 + others.size());
  sets[0] = wrap(greedy);
  for (std::size_t i = 0; i < others.size(); ++i) sets[i + 1] = wrap(others[i]);

  // `best`: the largest set; prefer greedy on a size tie (deterministic).
  const std::vector<int>* best = &greedy;
  for (std::size_t i = 0; i < others.size(); ++i)
    if (others[i].size() > best->size()) best = &others[i];

  IntegerVector size_dist(sd.size());
  CharacterVector nm(sd.size());
  int k = 0;
  for (std::map<int, int>::const_iterator it = sd.begin(); it != sd.end(); ++it, ++k) {
    size_dist[k] = it->second;
    nm[k] = std::to_string(it->first);
  }
  size_dist.attr("names") = nm;

  return List::create(_["best"] = wrap(*best),
                      _["sets"] = sets,
                      _["size_dist"] = size_dist);
}
